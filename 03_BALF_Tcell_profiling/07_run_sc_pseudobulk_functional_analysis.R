# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: run_sc_pseudobulk_functional_analysis.R
# Figure/table scope: Fig3; FigS6; Supplementary Table 13
# Data cohort: balf. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(DESeq2)
  library(clusterProfiler)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
})

# ==============================================================================
# 模块一：数据基建与预处理 (Infrastructure & Preprocessing)
# ==============================================================================

# 1. 检查并获取物种数据库
.get_orgdb <- function(organism = c("human", "mouse")) {
  organism <- match.arg(organism)
  
  if (organism == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("请先安装 org.Hs.eg.db: BiocManager::install('org.Hs.eg.db')")
    }
    OrgDb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    kegg_organism <- "hsa"
  } else if (organism == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("请先安装 org.Mm.eg.db: BiocManager::install('org.Mm.eg.db')")
    }
    OrgDb <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
    kegg_organism <- "mmu"
  }
  
  list(OrgDb = OrgDb, kegg_organism = kegg_organism)
}

# 2. 兼容 Seurat v4/v5 的 counts 提取
.get_counts <- function(seu, assay = "RNA", layer = "counts") {
  tryCatch(
    { SeuratObject::GetAssayData(seu, assay = assay, layer = layer) },
    error = function(e) { SeuratObject::GetAssayData(seu, assay = assay, slot = layer) }
  )
}

# 3. 构建 Pseudobulk counts 矩阵与元数据
make_pseudobulk_counts <- function(
    seu,
    group_col = "Group",
    sample_col = "HLA_Splitter",
    celltype_col = NULL,
    default_celltype = NULL,
    assay = "RNA",
    layer = "counts"
) {
  counts <- .get_counts(seu, assay = assay, layer = layer)
  meta <- seu@meta.data[colnames(counts), , drop = FALSE]
  
  # 校验缺失列
  required_cols <- c(group_col, sample_col)
  if (!is.null(celltype_col)) required_cols <- c(required_cols, celltype_col)
  
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("metadata 中缺少这些列: ", paste(missing_cols, collapse = ", "))
  }
  
  # 统一元数据命名与默认值填充
  meta$.group  <- as.character(meta[[group_col]])
  meta$.sample <- as.character(meta[[sample_col]])
  meta$.celltype <- if (is.null(celltype_col)) default_celltype else as.character(meta[[celltype_col]])
  
  # 过滤无效细胞
  keep_cells <- complete.cases(meta[, c(".group", ".sample", ".celltype")])
  meta <- meta[keep_cells, , drop = FALSE]
  counts <- counts[, rownames(meta), drop = FALSE]
  
  # 校验样本分组唯一性
  sample_group_check <- meta %>%
    dplyr::distinct(.sample, .group) %>%
    dplyr::count(.sample) %>%
    dplyr::filter(n > 1)
  
  if (nrow(sample_group_check) > 0) {
    stop("检测到同一个样本对应多个 Group，请检查 metadata。问题样本: ", 
         paste(sample_group_check$.sample, collapse = ", "))
  }
  
  # 矩阵转换计算 Pseudobulk
  meta$.pb_id <- paste(meta$.celltype, meta$.sample, sep = "__")
  pb_factor <- factor(meta$.pb_id)
  
  design_mat <- Matrix::sparse.model.matrix(~ 0 + pb_factor)
  colnames(design_mat) <- levels(pb_factor)
  
  pb_counts <- as(counts %*% design_mat, "dgCMatrix")
  
  pb_meta <- meta %>%
    dplyr::group_by(.pb_id) %>%
    dplyr::summarise(
      celltype = dplyr::first(.celltype),
      sample   = dplyr::first(.sample),
      group    = dplyr::first(.group),
      n_cells  = dplyr::n(),
      .groups  = "drop"
    ) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(".pb_id")
  
  # 对齐矩阵与元数据
  pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]
  
  list(counts = pb_counts, metadata = pb_meta)
}


# ==============================================================================
# 模块二：核心统计与功能注释 (Core Stats & Functional Annotation)
# ==============================================================================

# 4. 基因 ID 转换基础函数
map_genes_to_entrez <- function(genes, OrgDb, from_type = "SYMBOL") {
  genes <- unique(genes[complete.cases(genes) & genes != ""])
  if (length(genes) == 0) return(data.frame())
  
  suppressMessages({
    mp <- clusterProfiler::bitr(
      genes, fromType = from_type, toType = "ENTREZID", OrgDb = OrgDb
    )
  })
  
  if (is.null(mp) || nrow(mp) == 0) return(data.frame())
  unique(mp[, c(from_type, "ENTREZID")])
}

# 5. 单个细胞类型的 DESeq2 分析 (LRT 全局 + Wald 两两比较)
run_deseq_for_one_celltype <- function(
    pb_counts, pb_meta, celltype, comparisons,
    covariates = NULL, min_cells_per_sample = 20,
    min_samples_per_group = 2, min_count = 10, min_samples_expr = 2
) {
  # 样本过滤
  meta_ct <- pb_meta %>% filter(celltype == !!celltype, n_cells >= min_cells_per_sample)
  if (nrow(meta_ct) == 0) {
    message("跳过 ", celltype, ": 样本数不足")
    return(NULL)
  }
  
  group_tab <- table(meta_ct$group)
  if (sum(group_tab >= min_samples_per_group) < 2) {
    message("跳过 ", celltype, ": 有效分组不足两组")
    return(NULL)
  }
  
  # 数据准备与基因过滤
  cts <- as.matrix(pb_counts[, rownames(meta_ct), drop = FALSE])
  storage.mode(cts) <- "integer"
  meta_ct$group <- factor(meta_ct$group)
  
  keep_gene <- rowSums(cts >= min_count) >= min_samples_expr
  cts <- cts[keep_gene, , drop = FALSE]
  if (nrow(cts) < 100) {
    message("跳过 ", celltype, ": 过滤后可用基因数过少")
    return(NULL)
  }
  
  # 构建公式
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(meta_ct))
    if (length(missing_cov) > 0) stop("缺失协变量: ", paste(missing_cov, collapse = ", "))
  }
  
  full_formula <- as.formula(paste("~", paste(c(covariates, "group"), collapse = " + ")))
  reduced_formula <- as.formula(if(length(covariates) > 0) paste("~", paste(covariates, collapse = " + ")) else "~ 1")
  
  # 执行 LRT
  overall_res <- tryCatch({
    dds_lrt <- DESeq2::DESeqDataSetFromMatrix(cts, colData = meta_ct, design = full_formula)
    dds_lrt <- DESeq2::DESeq(dds_lrt, test = "LRT", reduced = reduced_formula, quiet = TRUE)
    as.data.frame(DESeq2::results(dds_lrt)) %>%
      rownames_to_column("gene") %>%
      as_tibble() %>%
      mutate(celltype = celltype, test = "LRT_group_overall") %>%
      relocate(celltype, test, gene)
  }, error = function(e) { message("LRT 失败: ", e$message); tibble() })
  
  # 执行 Wald 测试
  dds_wald <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(cts, colData = meta_ct, design = full_formula)
    DESeq2::DESeq(dds, quiet = TRUE)
  }, error = function(e) { message("Wald DESeq2 失败: ", e$message); NULL })
  
  pairwise_list <- list()
  if (!is.null(dds_wald)) {
    for (comp in comparisons) {
      num <- comp[1]; den <- comp[2]
      contrast_name <- paste0(num, "_vs_", den)
      
      sub_group <- as.character(meta_ct$group[meta_ct$group %in% c(num, den)])
      
      if (!all(c(num, den) %in% sub_group) || any(table(sub_group) < min_samples_per_group)) {
        message("提示：跳过 ", contrast_name, "，因为参与比较的有效样本数不足。")
        next
      }
      
      res <- tryCatch({
        as.data.frame(DESeq2::results(dds_wald, contrast = c("group", num, den))) %>%
          rownames_to_column("gene") %>%
          as_tibble() %>%
          mutate(celltype = celltype, contrast = contrast_name, numerator = num, denominator = den) %>%
          relocate(celltype, contrast, numerator, denominator, gene) %>%
          arrange(padj)
      }, error = function(e) { tibble() })
      
      pairwise_list[[contrast_name]] <- res
    }
  }
  
  list(overall = overall_res, pairwise = pairwise_list, dds = dds_wald, metadata = meta_ct)
}

# 6. 单个 contrast 的 GO/KEGG (ORA + GSEA) 分析
run_functional_for_contrast <- function(
    deg_table, OrgDb, kegg_organism = "hsa", gene_id_type = "SYMBOL",
    padj_cutoff = 0.05, lfc_cutoff = 0.25, go_ont = "BP",
    min_gs_size = 10, max_gs_size = 500
) {
  deg_table <- deg_table %>% filter(!is.na(gene))
  
  # 提取背景基因与差异基因
  universe_entrez <- unique(map_genes_to_entrez(deg_table$gene[!is.na(deg_table$pvalue)], OrgDb, gene_id_type)$ENTREZID)
  
  get_entrez <- function(direction_condition) {
    genes <- deg_table %>% filter(!!rlang::parse_expr(direction_condition)) %>% pull(gene) %>% unique()
    unique(map_genes_to_entrez(genes, OrgDb, gene_id_type)$ENTREZID)
  }
  
  up_entrez   <- get_entrez(paste0("!is.na(padj) & padj < ", padj_cutoff, " & log2FoldChange > ", lfc_cutoff))
  down_entrez <- get_entrez(paste0("!is.na(padj) & padj < ", padj_cutoff, " & log2FoldChange < -", lfc_cutoff))
  
  safe_call <- function(expr) tryCatch(expr, error = function(e) NULL)
  
  # --- ORA 分析 ---
  run_ora <- function(genes, func, ...) {
    if (length(genes) >= 5) safe_call(func(gene = genes, universe = universe_entrez, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, ...)) else NULL
  }
  
  ego_up   <- run_ora(up_entrez, clusterProfiler::enrichGO, OrgDb = OrgDb, keyType = "ENTREZID", ont = go_ont, readable = TRUE)
  ego_down <- run_ora(down_entrez, clusterProfiler::enrichGO, OrgDb = OrgDb, keyType = "ENTREZID", ont = go_ont, readable = TRUE)
  ekegg_up   <- run_ora(up_entrez, clusterProfiler::enrichKEGG, organism = kegg_organism)
  ekegg_down <- run_ora(down_entrez, clusterProfiler::enrichKEGG, organism = kegg_organism)
  
  # --- GSEA 分析 ---
  rank_df <- deg_table %>%
    mutate(score = case_when(
      !is.na(stat) ~ stat,
      is.na(stat) & !is.na(pvalue) & !is.na(log2FoldChange) ~ sign(log2FoldChange) * (-log10(pmax(pvalue, .Machine$double.xmin))),
      TRUE ~ NA_real_
    )) %>%
    filter(!is.na(score)) %>% select(gene, score)
  
  rank_map <- map_genes_to_entrez(rank_df$gene, OrgDb, gene_id_type)
  
  gene_list <- rank_df %>%
    inner_join(rank_map, by = c("gene" = gene_id_type)) %>%
    arrange(desc(abs(score))) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    { setNames(.$score, .$ENTREZID) } %>%
    sort(decreasing = TRUE)
  
  gsea_go <- gsea_kegg <- NULL
  if (length(gene_list) >= 100) {
    gsea_go   <- safe_call(clusterProfiler::gseGO(geneList = gene_list, OrgDb = OrgDb, keyType = "ENTREZID", ont = go_ont, minGSSize = min_gs_size, maxGSSize = max_gs_size, pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE))
    gsea_kegg <- safe_call(clusterProfiler::gseKEGG(geneList = gene_list, organism = kegg_organism, minGSSize = min_gs_size, maxGSSize = max_gs_size, pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE))
  }
  
  list(
    GO_UP = ego_up, GO_DOWN = ego_down, KEGG_UP = ekegg_up, KEGG_DOWN = ekegg_down,
    GSEA_GO = gsea_go, GSEA_KEGG = gsea_kegg
  )
}

# 7. 富集结果对象转长数据框
.enrich_to_table <- function(enrich_obj, method, direction, celltype, contrast) {
  if (is.null(enrich_obj)) return(tibble())
  df <- tryCatch(as_tibble(as.data.frame(enrich_obj)), error = function(e) tibble())
  if (nrow(df) == 0) return(tibble())
  
  df %>%
    mutate(
      direction = if ("NES" %in% colnames(.)) ifelse(NES > 0, "UP_IN_NUMERATOR", "DOWN_IN_NUMERATOR") else direction,
      celltype = celltype, contrast = contrast, method = method
    ) %>%
    relocate(celltype, contrast, method, direction)
}


# ==============================================================================
# 模块三：下游模式挖掘 (Advanced Pattern Mining)
# ==============================================================================

# 8. 通路热图矩阵提取
make_pathway_direction_matrix <- function(pathway_table, method_keep = c("GSEA_GO", "GSEA_KEGG"), padj_cutoff = 0.05) {
  pathway_table %>%
    filter(method %in% method_keep, !is.na(p.adjust), p.adjust < padj_cutoff, !is.na(NES)) %>%
    mutate(
      pathway = paste(method, ID, Description, sep = "__"),
      celltype_contrast = paste(celltype, contrast, sep = "__")
    ) %>%
    select(pathway, celltype_contrast, NES) %>%
    distinct() %>%
    pivot_wider(names_from = celltype_contrast, values_from = NES, values_fill = 0)
}

# 9. 挖掘通路状态逆转模块 (Reversed Pathways)
find_reversed_pathways <- function(pathway_table, contrast_a, contrast_b, method_keep = c("GSEA_GO", "GSEA_KEGG"), padj_cutoff = 0.05) {
  extract_contrast <- function(ctr) {
    pathway_table %>%
      filter(contrast == ctr, method %in% method_keep, !is.na(p.adjust), p.adjust < padj_cutoff, !is.na(NES)) %>%
      select(celltype, method, ID, Description, NES, padj = p.adjust)
  }
  
  inner_join(
    extract_contrast(contrast_a) %>% rename(NES_a = NES, padj_a = padj),
    extract_contrast(contrast_b) %>% rename(NES_b = NES, padj_b = padj),
    by = c("celltype", "method", "ID", "Description")
  ) %>%
    filter(sign(NES_a) != sign(NES_b)) %>%
    mutate(reversal_score = abs(NES_a) + abs(NES_b)) %>%
    arrange(desc(reversal_score))
}

# 10. 挖掘细胞类型特异性通路
find_celltype_specific_pathways <- function(pathway_table, contrast_name, method_keep = c("GSEA_GO", "GSEA_KEGG", "GO_ORA", "KEGG_ORA"), padj_cutoff = 0.05) {
  pathway_table %>%
    filter(contrast == contrast_name, method %in% method_keep, !is.na(p.adjust), p.adjust < padj_cutoff) %>%
    group_by(method, ID, Description) %>%
    mutate(n_celltypes = n_distinct(celltype)) %>%
    ungroup() %>%
    filter(n_celltypes == 1) %>%
    arrange(celltype, method, p.adjust)
}


# ==============================================================================
# 模块四：主控调度引擎 (Main Pipeline Wrapper)
# ==============================================================================

# 11. 主控执行函数
run_sc_pseudobulk_functional_analysis <- function(
    seu, group_col = "Group", sample_col = "HLA_Splitter", celltype_col = NULL, default_celltype = NULL,
    comparisons = NULL, organism = c("human", "mouse"),
    gene_id_type = "SYMBOL", assay = "RNA", layer = "counts", covariates = NULL,
    min_cells_per_sample = 20, min_samples_per_group = 2, min_count = 10, min_samples_expr = 2,
    padj_cutoff = 0.05, lfc_cutoff = 0.25, go_ont = "BP", output_dir = NULL
) {
  # 1. 基础配置构建
  org_info <- .get_orgdb(organism)
  pb <- make_pseudobulk_counts(seu, group_col, sample_col, celltype_col, default_celltype, assay, layer)
  
  celltypes <- unique(pb$metadata$celltype)
  all_overall <- list(); all_pairwise <- list(); functional_objects <- list(); functional_tables <- list()
  
  # 2. 遍历执行细胞层级分析
  for (ct in celltypes) {
    message("========== 分析细胞类型: ", ct, " ==========")
    de_ct <- run_deseq_for_one_celltype(pb$counts, pb$metadata, ct, comparisons, covariates, min_cells_per_sample, min_samples_per_group, min_count, min_samples_expr)
    if (is.null(de_ct)) next
    
    all_overall[[ct]] <- de_ct$overall
    functional_objects[[ct]] <- list()
    
    for (contrast_name in names(de_ct$pairwise)) {
      deg_tbl <- de_ct$pairwise[[contrast_name]]
      all_pairwise[[paste(ct, contrast_name, sep = "__")]] <- deg_tbl
      
      if (nrow(deg_tbl) > 0) {
        message("  -> 功能富集: ", contrast_name)
        fun_res <- run_functional_for_contrast(deg_tbl, org_info$OrgDb, org_info$kegg_organism, gene_id_type, padj_cutoff, lfc_cutoff, go_ont)
        
        functional_objects[[ct]][[contrast_name]] <- fun_res
        functional_tables[[paste(ct, contrast_name, sep = "__")]] <- bind_rows(
          .enrich_to_table(fun_res$GO_UP, "GO_ORA", "UP", ct, contrast_name),
          .enrich_to_table(fun_res$GO_DOWN, "GO_ORA", "DOWN", ct, contrast_name),
          .enrich_to_table(fun_res$KEGG_UP, "KEGG_ORA", "UP", ct, contrast_name),
          .enrich_to_table(fun_res$KEGG_DOWN, "KEGG_ORA", "DOWN", ct, contrast_name),
          .enrich_to_table(fun_res$GSEA_GO, "GSEA_GO", "GSEA", ct, contrast_name),
          .enrich_to_table(fun_res$GSEA_KEGG, "GSEA_KEGG", "GSEA", ct, contrast_name)
        )
      }
    }
  }
  
  # 3. 统计细胞比例
  cell_proportion <- NULL
  if (!is.null(celltype_col)) {
    cell_proportion <- seu@meta.data %>%
      dplyr::filter(!is.na(!!sym(group_col)), !is.na(!!sym(sample_col)), !is.na(!!sym(celltype_col))) %>%
      dplyr::count(sample = !!sym(sample_col), group = !!sym(group_col), celltype = !!sym(celltype_col), name = "n_cells") %>%
      dplyr::group_by(sample, group) %>%
      dplyr::mutate(total_cells = sum(n_cells), proportion = n_cells / total_cells) %>%
      dplyr::ungroup()
  }
  
  # 4. 数据汇编与落盘
  overall_table <- bind_rows(all_overall)
  pairwise_table <- bind_rows(all_pairwise)
  pathway_table <- bind_rows(functional_tables)
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(overall_table, file.path(output_dir, "01_overall_LRT_group_test.csv"))
    readr::write_csv(pairwise_table, file.path(output_dir, "02_pairwise_DEG_all_contrasts.csv"))
    readr::write_csv(pathway_table, file.path(output_dir, "03_GO_KEGG_GSEA_pathway_table.csv"))
    if (!is.null(cell_proportion)) readr::write_csv(cell_proportion, file.path(output_dir, "04_celltype_proportion_by_sample.csv"))
  }
  
  list(
    pseudobulk = pb, overall_LRT = overall_table, pairwise_DEG = pairwise_table,
    functional_objects = functional_objects, pathway_table = pathway_table, cell_proportion = cell_proportion
  )
}

# ==============================================================================
# 模块五：动态变化基因的模式聚类与可视化 (LRT Heatmap & Trend Plot)
# ==============================================================================

# 12. 内部函数：获取特定细胞类型的 Z-score 矩阵
.get_scaled_lrt_matrix <- function(
    lrt_table, pb_counts, pb_meta, target_celltype, 
    filter_by = c("padj", "pvalue"), p_cutoff = 0.01, top_n = 1000
) {
  # 校验并匹配参数输入
  filter_by <- match.arg(filter_by)
  
  # 1. 动态选择用于过滤和排序的列
  sig_genes <- lrt_table %>%
    dplyr::filter(
      celltype == target_celltype, 
      !is.na(.data[[filter_by]]), 
      .data[[filter_by]] < p_cutoff
    ) %>%
    dplyr::arrange(.data[[filter_by]]) %>%
    head(top_n) %>%
    dplyr::pull(gene)
  
  if (length(sig_genes) < 10) {
    stop(
      target_celltype, " 中满足 ", filter_by, " < ", p_cutoff, 
      " 的基因少于 10 个，无法画图。\n提示：请尝试放宽 p_cutoff 阈值（如 0.05），或者将 filter_by 设为 'pvalue'。"
    )
  }
  
  # 2. 提取对应细胞类型的样本
  meta_ct <- pb_meta %>%
    dplyr::filter(celltype == target_celltype) %>%
    dplyr::arrange(group) 
  
  # 3. 提取 counts 矩阵并做 log2(CPM + 1) 标准化
  cts <- pb_counts[sig_genes, rownames(meta_ct), drop = FALSE]
  cpm_mat <- sweep(cts, 2, colSums(cts), FUN = "/") * 1e6
  logcpm_mat <- log2(cpm_mat + 1)
  
  # 4. 行标准化 (Z-score)：为了抹平高低表达基因的绝对差异，只看趋势
  z_mat <- t(scale(t(logcpm_mat)))
  z_mat <- z_mat[complete.cases(z_mat), , drop = FALSE]
  
  list(
    z_matrix = z_mat,
    metadata = meta_ct,
    filter_metric = filter_by
  )
}


# 13. 绘制全局 LRT 热图
plot_lrt_global_heatmap <- function(
    lrt_table, pb_counts, pb_meta, target_celltype, 
    filter_by = "padj", p_cutoff = 0.01, top_n = 1000, 
    k_clusters = 4, show_rownames = FALSE
) {
  # 获取 Z-score 矩阵
  mat_info <- .get_scaled_lrt_matrix(lrt_table, pb_counts, pb_meta, target_celltype, filter_by, p_cutoff, top_n)
  z_mat <- mat_info$z_matrix
  meta_ct <- mat_info$metadata
  metric <- mat_info$filter_metric
  
  # 准备注释条
  annotation_col <- data.frame(
    Group = meta_ct$group,
    row.names = rownames(meta_ct)
  )
  
  # 绘制带 K-means 分块的热图
  p <- pheatmap::pheatmap(
    mat = z_mat,
    cluster_cols = FALSE,     
    cluster_rows = TRUE,      
    kmeans_k = NA,            
    cutree_rows = k_clusters, 
    annotation_col = annotation_col,
    show_colnames = FALSE,
    show_rownames = show_rownames,
    scale = "none",           
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    main = paste0(target_celltype, " - Dynamic Genes (", metric, " < ", p_cutoff, ")\nTop ", nrow(z_mat), " genes")
  )
  
  return(p)
}


# 14. 提取聚类模块，并画出“模式趋势折线图”
plot_lrt_cluster_trends <- function(
    lrt_table, pb_counts, pb_meta, target_celltype, 
    filter_by = "padj", p_cutoff = 0.01, top_n = 1000, 
    k_clusters = 4
) {
  mat_info <- .get_scaled_lrt_matrix(lrt_table, pb_counts, pb_meta, target_celltype, filter_by, p_cutoff, top_n)
  z_mat <- mat_info$z_matrix
  meta_ct <- mat_info$metadata
  metric <- mat_info$filter_metric
  
  # 1. 运行 K-means 对基因进行聚类
  set.seed(42) 
  km <- kmeans(z_mat, centers = k_clusters, iter.max = 50)
  
  # 2. 将 Z-score 矩阵转为长格式，并融合分组信息
  long_df <- as.data.frame(z_mat) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(cluster = paste0("Cluster ", km$cluster[gene])) %>%
    tidyr::pivot_longer(
      cols = -c(gene, cluster), 
      names_to = ".pb_id",    # 修改点 1：转换出来的列名改叫 .pb_id
      values_to = "z_score"
    ) %>%
    dplyr::left_join(
      meta_ct %>% tibble::rownames_to_column(".pb_id") %>% dplyr::select(.pb_id, group), 
      by = ".pb_id"
    )
  
  # 3. 计算每个 Cluster 在不同 Group 的平均表达水平和置信区间
  trend_df <- long_df %>%
    dplyr::group_by(cluster, group) %>%
    dplyr::summarise(
      mean_z = mean(z_score),
      sd_z = sd(z_score),
      n = dplyr::n(),
      se_z = sd_z / sqrt(n),
      .groups = "drop"
    ) %>%
    dplyr::mutate(group = factor(group)) 
  
  # 4. 使用 ggplot2 绘制趋势图
  p <- ggplot(trend_df, aes(x = group, y = mean_z, group = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(color = "#0073C2FF", linewidth = 1.2) +
    geom_point(color = "#0073C2FF", size = 3) +
    geom_errorbar(aes(ymin = mean_z - se_z, ymax = mean_z + se_z), width = 0.2, color = "#0073C2FF") +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(target_celltype, " - Gene Expression Patterns"),
      subtitle = paste0("K-means clustering (k=", k_clusters, ") | Filter: ", metric, " < ", p_cutoff),
      x = "Disease Group",
      y = "Mean Z-score of Expression"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )
  
  gene_clusters <- data.frame(
    gene = names(km$cluster),
    cluster = paste0("Cluster ", km$cluster),
    row.names = NULL
  )
  
  list(
    plot = p,
    gene_clusters = gene_clusters
  )
}

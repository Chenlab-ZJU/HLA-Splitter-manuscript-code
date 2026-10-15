# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: plot_top_pseudobulk.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# 必备 R 包加载
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# ==============================================================================
# 函数 1：自动化输出 Pseudobulk Lollipop 图 + 分组专属 CSV + Bar&Dot 图
# ==============================================================================
#' @param seurat_obj Seurat 对象
#' @param deg_table FindAllMarkers 的初步结果
#' @param group_by 临床/实验分组列 (如 "Group")
#' @param sample_col 样本ID列 (如 "HLA_Splitter")
#' @param top_n_lollipop Lollipop 总图每组展示多少个基因 (默认 10)
#' @param top_n_bardot 每个分组的 Bar+Dot 图展示多少个 Top 基因 (默认 4)
#' @param out_dir 输出文件夹路径
#' @param my_colors 自定义调色板
plot_top_pseudobulk <- function(
    seurat_obj, 
    group_by = "Group", 
    sample_col = "HLA_Splitter", 
    top_n_lollipop = 10, 
    top_n_bardot = 4,
    out_dir = "./Figures/Marker_Analysis",
    my_colors = NULL
) {
  
  # 0. 准备输出目录
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----------------------------------------------------------------------------
  # 1. 核心数据准备与 Pseudobulk 均值测算
  # ----------------------------------------------------------------------------
  message("▶ 步骤 1: 提取候选基因并计算样本级(Pseudobulk)表达量及样本级覆盖度...")
  Idents(seurat_obj) <- group_by
  
  deg_table <- FindAllMarkers(
    seurat_obj, 
    only.pos = TRUE, 
    min.pct = 0.1, 
    logfc.threshold = 0.1,
    return.thresh = 0.05
  )
  
  if (nrow(deg_table) == 0) {
    stop("❌ 错误：FindAllMarkers 没有找到任何符合条件的显著差异基因。")
  }
  
  noise_pattern <- "^MT-|^HB[ABDEGQZ]|^RP[0-9]+-|^AC[0-9]+\\.|^AL[0-9]+\\.|^AP[0-9]+\\.|^CT[ABCD]-|^RPS|^RPL|MALAT1|NEAT1"
  
  candidate_genes <- deg_table %>%
    filter(p_val_adj < 0.05 & avg_log2FC > 0.1) %>%
    filter(!grepl(noise_pattern, gene)) %>%
    group_by(cluster) %>% 
    pull(gene) %>% 
    unique()
  
  expr_data <- FetchData(seurat_obj, vars = unique(c(candidate_genes, group_by, sample_col)))
  colnames(expr_data)[which(colnames(expr_data) == group_by)] <- "Group_Col"
  colnames(expr_data)[which(colnames(expr_data) == sample_col)] <- "Sample_Col"
  
  # 计算每个样本的真实算数均值
  df_sample <- expr_data %>%
    pivot_longer(cols = all_of(candidate_genes), names_to = "gene", values_to = "Expression") %>%
    group_by(Group_Col, Sample_Col, gene) %>%
    summarise(Sample_Mean_Exp = mean(expm1(Expression)), .groups = "drop")
  
  # --- [新增] 计算“样本级”的基因覆盖度 (Coverage = 该组内表达均值 > 0 的样本比例) ---
  df_sample_coverage <- df_sample %>%
    group_by(gene, Group_Col) %>%
    summarise(
      Group_Express_Samples = sum(Sample_Mean_Exp > 0), # 该组内表达量大于0的样本数
      Group_Total_Samples = n(),                        # 该组的总样本数
      .groups = "drop"
    ) %>%
    group_by(gene) %>%
    mutate(
      Global_Express_Samples = sum(Group_Express_Samples),
      Global_Total_Samples = sum(Group_Total_Samples),
      # 目标组样本覆盖度
      Target_Sample_Coverage = Group_Express_Samples / Group_Total_Samples,
      # 其他组样本覆盖度
      Other_Sample_Coverage = ifelse((Global_Total_Samples - Group_Total_Samples) > 0,
                                     (Global_Express_Samples - Group_Express_Samples) / (Global_Total_Samples - Group_Total_Samples),
                                     0)
    ) %>%
    ungroup() %>%
    dplyr::select(gene, Group_Col, Target_Sample_Coverage, Other_Sample_Coverage)
  # ----------------------------------------------------------------------------------
  
  df_fc <- df_sample %>%
    group_by(gene, Group_Col) %>%
    summarise(
      Group_Sum = sum(Sample_Mean_Exp), 
      Group_N = n(),
      .groups = "drop"
    ) %>%
    group_by(gene) %>%
    mutate(
      Global_Sum = sum(Group_Sum), 
      Global_N = sum(Group_N),
      Group_Mean = Group_Sum / Group_N,
      Other_Mean = ifelse(Global_N > Group_N, (Global_Sum - Group_Sum) / (Global_N - Group_N), 0),
      PB_Log2FC = log2(Group_Mean + 1) - log2(Other_Mean + 1)
    ) %>%
    ungroup()
  
  # ----------------------------------------------------------------------------
  # 2. 生成并保存 Lollipop 总图
  # ----------------------------------------------------------------------------
  message("▶ 步骤 2: 绘制并导出差异倍数棒棒糖总图 (Lollipop Plot)...")
  plot_df_lolli <- df_fc %>%
    group_by(Group_Col) %>%
    slice_max(n = top_n_lollipop, order_by = PB_Log2FC, with_ties = FALSE) %>%
    ungroup() %>% 
    arrange(Group_Col, PB_Log2FC) %>%
    mutate(gene_label = factor(paste0(gene, "_", Group_Col), levels = paste0(gene, "_", Group_Col)))
  
  p_lolli <- ggplot(plot_df_lolli, aes(x = PB_Log2FC, y = gene_label, color = Group_Col)) +
    geom_segment(aes(x = 0, xend = PB_Log2FC, y = gene_label, yend = gene_label), linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 4) +
    geom_text(aes(label = round(PB_Log2FC, 2)), hjust = -0.4, size = 3.5, color = "black") +
    facet_wrap(~ Group_Col, scales = "free_y", ncol = 3) +
    scale_y_discrete(labels = function(x) gsub("_.*", "", x)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme_classic(base_size = 14) +
    labs(x = "Sample-level Log2 Fold Change", y = "Top Marker genes") +
    theme(
      axis.text.y = element_text(color = "black", face = "bold.italic"), 
      axis.text.x = element_text(color = "black"),
      strip.text = element_text(size = 14, face = "bold"), 
      strip.background = element_rect(fill = "grey95", color = "black"),
      panel.border = element_rect(color = "black", fill = NA), 
      panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"), 
      legend.position = "none"
    )
  if (!is.null(my_colors)) p_lolli <- p_lolli + scale_color_manual(values = my_colors)
  
  ggsave(file.path(out_dir, "Overall_FoldChange_Lollipop.pdf"), plot = p_lolli, height = 6, width = 12)
  
  # ----------------------------------------------------------------------------
  # 3. 循环遍历每个分组：计算差异检验 + 导出 CSV + 专属 BarDot 图
  # ----------------------------------------------------------------------------
  message("▶ 步骤 3: 为每个分组进行样本级差异检验、生成专属 CSV 文件和 Bar+Dot 验证图...")
  unique_groups <- unique(df_fc$Group_Col)
  
  for (grp in unique_groups) {
    message("  - 正在处理分组: ", grp)
    
    grp_df <- df_fc %>% filter(Group_Col == grp)
    
    # 差异检验
    p_values <- sapply(grp_df$gene, function(g) {
      sub_data <- df_sample %>% filter(gene == g)
      grp_vector <- sub_data$Sample_Mean_Exp[sub_data$Group_Col == grp]
      other_vector <- sub_data$Sample_Mean_Exp[sub_data$Group_Col != grp]
      
      if (length(grp_vector) < 2 || length(other_vector) < 2 || 
          all(grp_vector == grp_vector[1]) && all(other_vector == other_vector[1])) {
        return(NA)
      }
      test_res <- wilcox.test(grp_vector, other_vector, exact = FALSE)
      return(test_res$p.value)
    })
    
    # --- 将“样本级”覆盖度与P值合并，并生成百分比 ---
    grp_coverage <- df_sample_coverage %>% filter(Group_Col == grp) %>% dplyr::select(-Group_Col)
    
    grp_df <- grp_df %>%
      mutate(
        PB_Pval = p_values,
        PB_FDR = p.adjust(PB_Pval, method = "fdr")
      ) %>%
      left_join(grp_coverage, by = "gene") %>%
      mutate(
        Target_Sample_Coverage = paste0(round(Target_Sample_Coverage * 100, 2), "%"),
        Other_Sample_Coverage = paste0(round(Other_Sample_Coverage * 100, 2), "%")
      ) %>%
      arrange(desc(PB_Log2FC)) 
    
    # 导出 CSV 
    csv_name <- file.path(out_dir, paste0("Marker_Stats_", grp, ".csv"))
# Export the analysis source data for figure reproduction.
    write.csv(grp_df, csv_name, row.names = FALSE)
    
    # B. 提取该组 Top N 基因画 Bar+Dot 图
    top_genes_grp <- head(grp_df$gene, top_n_bardot)
    
    plot_df_sample <- df_sample %>% filter(gene %in% top_genes_grp)
    plot_df_sample$gene <- factor(plot_df_sample$gene, levels = top_genes_grp)
    
    plot_df_summary <- plot_df_sample %>%
      group_by(Group_Col, gene) %>%
      summarise(Mean = mean(Sample_Mean_Exp), SD = sd(Sample_Mean_Exp), SE = SD/sqrt(n()), .groups = "drop")
    
    p_bardot <- ggplot() +
      geom_jitter(data = plot_df_sample, aes(x = Group_Col, y = Sample_Mean_Exp, color = Group_Col), 
                  width = 0.2, alpha = 0.8, size = 2) +
      geom_col(data = plot_df_summary, aes(x = Group_Col, y = Mean, fill = Group_Col), 
               color = "black", alpha = 0.5, linewidth = 0.4, width = 0.6) +
      geom_errorbar(data = plot_df_summary, aes(x = Group_Col, ymin = Mean, ymax = Mean + SE), 
                    width = 0.2, color = "black", linewidth = 0.6) +
      facet_wrap(~ gene, scales = "free_y", ncol = 4) +
      theme_classic(base_size = 14) +
      labs(title = paste0("Top Markers for: ", grp), x = "Clinical Group", y = "Average Expression (expm1)") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(size = 14, face = "bold.italic"), 
        strip.background = element_rect(fill = "grey90", color = "black"),
        panel.border = element_rect(color = "black", fill = NA), 
        legend.position = "none"
      )
    if (!is.null(my_colors)) {
      p_bardot <- p_bardot + scale_fill_manual(values = my_colors) + scale_color_manual(values = my_colors)
    }
    
    pdf_name <- file.path(out_dir, paste0("BarDot_", grp, "_Top", top_n_bardot, ".pdf"))
    ggsave(pdf_name, plot = p_bardot, height = top_n_bardot * 0.8, width = 8)
  }
  
  message("✅ 全部流程执行完毕！所有结果已保存至: ", out_dir)
}

# ==============================================================================
# 函数 2：绘制并保存特定亚群的 Top 基因小提琴图
# ==============================================================================
plot_top_violins <- function(
    seu_obj,           # Seurat 对象
    deg_data,          # 差异基因数据框
    target_group,      # 当前要分析的亚群名称 (如 "CAP")
    out_dir,           # 输出文件夹
    comparisons_list,  # 用于 stat_compare_means 的比较列表
    plot_colors = NULL,# 自定义颜色
    group_col = "Group_Col", # deg_data 中代表亚群的列名
    lfc_col = "PB_Log2FC",   # deg_data 中代表 Log2FC 的列名
    gene_col = "gene",       # deg_data 中代表基因名的列名
    vln_group_by = "Group",  # Seurat 对象中用于分组画小提琴的列 (如 "Group")
    lfc_threshold = 0.1,     # 提取基因的最低 Log2FC
    top_n = 8                # 绘制小提琴图的基因数量
) {
  # 1. 过滤并提取 Top 基因
  grp_markers <- deg_data %>% 
    filter(!!sym(group_col) == target_group & !!sym(lfc_col) > lfc_threshold)
  
  if (nrow(grp_markers) == 0) return(NULL)
  
  top_vln_genes <- head(grp_markers[[gene_col]], top_n)
  valid_vln_genes <- intersect(top_vln_genes, rownames(seu_obj))
  
  if (length(valid_vln_genes) == 0) return(NULL)
  
  # 2. 生成原生 VlnPlot 列表
  plots_violins <- VlnPlot(
    seu_obj, features = valid_vln_genes, group.by = vln_group_by, 
    cols = plot_colors, pt.size = 0, ncol = 4, combine = FALSE
  )
  
  # 3. 循环美化并添加显著性检验
  for (i in seq_along(plots_violins)) {
    plots_violins[[i]] <- plots_violins[[i]] + 
      theme_classic() + 
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      theme(
        axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(), 
        legend.position = 'none'
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
      stat_compare_means(aes(label = ..p.format..), comparisons = comparisons_list, step.increase = 0.2, vjust = -0.05)
  }
  
  # 4. 拼接并保存
  combined_vln_plot <- patchwork::wrap_plots(plots_violins, ncol = 4)
  filename <- file.path(out_dir, paste0(target_group, "_VlnPlot.pdf"))
  ggsave(filename = filename, plot = combined_vln_plot, width = 16, height = 8)
  
  message("  - Saved Violin plots to: ", filename)
  return(combined_vln_plot)
}


# ==============================================================================
# 函数 3：执行并保存特定亚群的 GO 富集分析
# ==============================================================================
run_go_enrichment <- function(
    deg_data, 
    target_group, 
    out_dir,
    group_col = "Group_Col",
    lfc_col = "PB_Log2FC",
    gene_col = "gene",
    lfc_threshold = 0.8,    # 富集分析通常需要更严苛的倍数限制
    top_n_pathways = 8      # 气泡图展示的通路数
) {
  # 1. 严格过滤基因
  go_genes_symbol <- deg_data %>% 
    filter(!!sym(group_col) == target_group & !!sym(lfc_col) > lfc_threshold) %>% 
    pull(!!sym(gene_col)) %>% 
    unique()
  
  if (length(go_genes_symbol) < 5) {
    message("  - Not enough highly upregulated genes (LFC > ", lfc_threshold, ") for GO enrichment.")
    return(NULL)
  }
  
  # 2. 基因名转换
  gene_mapped <- tryCatch({
    bitr(go_genes_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)
  
  if (is.null(gene_mapped) || nrow(gene_mapped) == 0) {
    message("  - Gene ID conversion failed.")
    return(NULL)
  }
  
  # 3. 执行 GO 富集
  ego <- enrichGO(
    gene          = gene_mapped$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  # 4. 结果落盘与出图
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    
    # 保存 CSV
    csv_filename <- file.path(out_dir, paste0(target_group, "_enriched_GO.csv"))
    write.csv(as.data.frame(ego), file = csv_filename, row.names = FALSE)
    message("  - Saved GO Enrichment table to: ", csv_filename)
    
    # 绘制气泡图并保存
    p_enrich <- dotplot(ego, showCategory = top_n_pathways) +
      ggtitle(paste0("GO BP Enrichment: ", target_group)) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.y = element_text(size = 10)
      )
    
    plot_filename <- file.path(out_dir, paste0(target_group, "_Enrichment_Plot.pdf"))
    ggsave(filename = plot_filename, plot = p_enrich, width = 6, height = 5)
    message("  - Saved GO Enrichment plot to: ", plot_filename)
    
    return(list(ego_obj = ego, plot = p_enrich))
    
  } else {
    message("  - No significant GO terms enriched.")
    return(NULL)
  }
}

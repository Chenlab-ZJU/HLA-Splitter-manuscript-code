# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Upstream urine processing and sample-level signature analysis.
# Original source: 01_data_integration_and_QC.R
# Figure/table scope: Fig4; FigS7; FigS8
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# File: 01_data_integration_and_QC.R
# Description: 数据合并、HLA矩阵分析、QC、细胞分群与全局特征可视化
# ==============================================================================
source(file.path(getOption("hla.code_root"), "04_Urine_immune_programs/01_setup_and_helpers.R"))

ensure_dir(FIGURE_DIR)
ensure_dir(RESULT_DIR)
ensure_dir(QC_DIR)
ensure_dir(DEMULTIPLEX_RESULT_DIR)

# ======================= A. 读入数据与基础 QC =======================
samplenames <- ANALYSIS_BATCHES
scRNA_list <- list()
counter <- 1

# 读取所有 10X filtered_feature_bc_matrix 目录；兼容 GEM 下直接存放或多批次子目录存放。
for (batch in samplenames) {
  tenx_dirs <- find_10x_dirs(batch)
  if (nrow(tenx_dirs) == 0) {
    warning("No 10X matrix found for batch: ", batch)
    next
  }

  for (row_id in seq_len(nrow(tenx_dirs))) {
    data_path <- tenx_dirs$data_dir[row_id]
    run_index <- tenx_dirs$run_index[row_id]
    message("Reading 10X matrix: ", data_path)

    Counts <- Read10X(data.dir = data_path)
    project_name <- paste0(batch, "_", run_index)
    obj <- CreateSeuratObject(counts = Counts, project = project_name, assay = "RNA")
    new_barcodes <- paste0(str_sub(colnames(obj), 1, 16), "-", batch, "-", run_index)
    obj <- RenameCells(obj, new.names = new_barcodes)
    scRNA_list[[counter]] <- obj
    counter <- counter + 1
  }
}
# 合并所有的 Seurat 对象
if (length(scRNA_list) > 1) {
  scRNA_Intergrate <- merge(scRNA_list[[1]], y = scRNA_list[-1], project = "Urine_cell")
  message("合并成功！共合并了 ", length(scRNA_list), " 个样本。")
} else if (length(scRNA_list) == 1) {
  scRNA_Intergrate <- scRNA_list[[1]]
  message("只读取到一个样本，无需合并。")
} else {
  stop("未读取到任何有效数据，请检查文件夹结构和路径是否正确！")
}

scRNA_Intergrate <- scRNA_QC(scRNA_Intergrate, min_RNA = 200, max_MT = 20, save_path = QC_DIR)

# ======================= B. 去除 Doublets (scrubletR) =======================
options(future.globals.maxSize = 10000 * 1024^2)

scRNA_Intergrate <- NormalizeData(scRNA_Intergrate) %>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData() %>% RunPCA() %>% 
  RunHarmony(group.by.vars = "orig.ident") %>% 
  JoinLayers()

scRNA_Intergrate <- scrublet(scRNA_Intergrate)
scRNA_Intergrate <- subset(scRNA_Intergrate, subset = doublet_scores < 0.2)

# Export the analysis source data for figure reproduction.
saveRDS(scRNA_Intergrate, object_path("Urine_scRNA_Intergrate.rds"))

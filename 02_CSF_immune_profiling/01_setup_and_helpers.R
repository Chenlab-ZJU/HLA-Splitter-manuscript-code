# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Upstream cohort processing and reusable helper definitions.
# Original source: 00_setup_and_helpers.R
# Figure/table scope: Fig2; FigS4; FigS5
# Data cohort: 脑脊液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.

# ==============================================================================
# File: 00_setup_and_helpers.R
# Description: 加载依赖包、设置全局变量、加载外部工具脚本及定义辅助函数
# ==============================================================================

WORK_DIR <- file.path(getOption("hla.data_root"), "\u8111\u810a\u6db2\u6837\u672c")
setwd(WORK_DIR)
set.seed(2026) # 设定随机种子，保证UMAP和聚类结果可重复

# 1. 加载依赖包 ----------------------------------------------------------------
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(harmony)
  library(reshape2)
  library(pheatmap)
  library(sctransform)
  library(RColorBrewer)
  library(grid)
  library(stringr)
  library(gtools)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  library(scrubletR)
  library(paletteer)
  library(ggExtra)
  library(enrichR)
  library(GSVA)
  library(msigdbr)
  library(clusterProfiler)
  library(limma)
  library(ggrepel)
  library(presto)
})

# 2. 统一的颜色盘配置 ----------------------------------------------------------
GLOBAL_COLORS <- c(
  "#E42320","#699ECA","#4DAF4A","#FF8C00","#F898CB","#D65190","#731A73","#FFCB5B",
  "#E89B1E","#0076B9","#3D505A","#0098B2","#8A7355","#CC5B45","#008B8B","#4E5689",
  "#F9B3AD","#EDC66A","#9FDAF7","#C9A1CA","#D9DEE7","#C74546","#EB7E60","#8FB4DC"
)
GROUP_COLORS <- c("Metastases"= "#E42320", "Hemorrhage" = "#699ECA", "Encephalitis" = "#4DAF4A", "Doublets" = "black", "Undefined" = "grey")

comparison <- list(c("Metastases","Hemorrhage"), c("Metastases","Encephalitis"), c("Hemorrhage","Encephalitis"))

# 3. 引入自定义外部脚本 --------------------------------------------------------
source(file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks/01_HLA_demultiplex_v2.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_grouped_umap.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_cell_proportion_barplot.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_HLA_tSNE.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_celltype_ratio.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_cell_ratio_comparison.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_volcano.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/Enrichment_analysis.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_pseudobulk_heatmap.R"))

# 4. 辅助函数 (Helper Functions) -----------------------------------------------
# 提取亚群Top10差异基因并导出
export_top_markers <- function(seurat_obj, ident_col, out_prefix) {
  Idents(seurat_obj) <- ident_col
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Export the analysis source data for figure reproduction.
  write.csv(markers, paste0("./Result/", out_prefix, "_makers.csv"), row.names = FALSE)
  
  df_summary <- seurat_obj@meta.data %>%
    select(all_of(ident_col)) %>% distinct() %>% arrange(across(all_of(ident_col)))
  
  df_top_genes <- markers %>%
    filter(p_val_adj < 0.01) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>%
    summarise(top_genes = paste0(gene, collapse = ", ")) %>%
    rename(!!ident_col := cluster)
  
  final_df <- left_join(df_summary, df_top_genes, by = ident_col)
  final_df$Celltype = c("")
  write.csv(final_df, paste0("./Result/", out_prefix, "_makers_summary.csv"), row.names = FALSE)
  
  message("请手动打开 ./Result/", out_prefix, "_makers_summary.csv 文件，添加 Celltype 列后保存！")
}

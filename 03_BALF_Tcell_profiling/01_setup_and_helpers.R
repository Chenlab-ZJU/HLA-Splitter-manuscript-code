# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: 00_setup_and_helpers.R
# Figure/table scope: Fig3; FigS6; Supplementary Table 13
# Data cohort: 肺泡灌洗液. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# File: 00_setup_and_helpers.R
# Description: 全局依赖加载、颜色配置、自定义辅助函数
# ==============================================================================

WORK_DIR <- file.path(getOption("hla.data_root"), "\u80ba\u6ce1\u704c\u6d17\u6db2")
setwd(WORK_DIR)
set.seed(2026)

# 1. 加载所有依赖包 ------------------------------------------------------------
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(harmony)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
  library(stringr)
  library(gtools)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  library(scrubletR)
  library(paletteer)
  library(ggExtra)
  library(ggrepel)
  library(enrichR)
  library(GSVA)
  #remotes::install_github("BorchLab/scRepertoire@v2.1.0", upgrade = "never")
  library(scRepertoire)
  library(immunarch)
  library(ggseqlogo)
  library(circlize)
  library(tidyr)
  library(VennDiagram)
  library(SCIPAC)
  library(MAST)
  #library(turboGliph)
  library(AUCell)
  library(glmnet)
  library(pROC)
  library(scCustomize)
  library(viridisLite)
})

# 2. 全局调色板 ----------------------------------------------------------------
GLOBAL_COLORS <- c("#E42320","#699ECA","#4DAF4A","#FF8C00","#F898CB","#D65190","#731A73","#FFCB5B",
                   "#E89B1E","#0076B9","#3D505A","#0098B2","#8A7355","#CC5B45","#008B8B","#4E5689",
                   "#F9B3AD","#EDC66A","#9FDAF7","#C9A1CA","#D9DEE7","#C74546","#EB7E60","#8FB4DC")
GROUP_COLORS <- c("LC"="#E42320", "CAP"="#699ECA", "ILD"="#4DAF4A", "Others"="#FF8C00", "Doublets"="black")
CLONE_SIZE_COLORS <- c("#f9f871","#ffc75f","#ff6f91","#845ec2","grey")


comparison <- list(c("CAP","LC"), c("CAP","ILD"), c("ILD","LC"))

# 3. 引入自定义脚本 ------------------------------------------------------------
source(file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks/01_HLA_demultiplex_v2.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_grouped_umap.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_cell_proportion_barplot.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_HLA_tSNE.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_celltype_ratio.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_cell_ratio_comparison.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_volcano.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/Enrichment_analysis.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_pseudobulk_heatmap.R"))
source(file.path(getOption("hla.code_root"), "00_Shared/plot_dotplot.R"))

# 4. 自定义辅助函数 ------------------------------------------------------------
Feature_Plot_Multi <- function(cell_subset, features, ncol = 3, reduction = "umap"){
  p_merge <- list()
  for (i in seq_along(features)) {
    p_merge[[i]] <- FeaturePlot(cell_subset, features = features[i], order = TRUE, pt.size = 0.4, reduction = reduction) +
      scale_color_gradientn(colors = c("#fef9ef","#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),
                            limits = c(0, 6), name = "Expression") +
      theme_void() + labs(title = features[i]) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), legend.position = "none")
  }
  
  legend_plot <- FeaturePlot(cell_subset, features = features[1]) +
    scale_color_gradientn(colors = c("#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"), limits = c(0, 6)) +
    theme(legend.position = "right")
  
  legend <- cowplot::get_legend(legend_plot)
  plot_grid <- wrap_plots(p_merge, ncol = ncol)
  return(cowplot::plot_grid(plot_grid, legend, rel_widths = c(1, 0.2)))
}

export_top_markers <- function(seurat_obj, ident_col, out_prefix) {
  Idents(seurat_obj) <- ident_col
  markers <- SeuratWrappers::RunPrestoAll(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Export the analysis source data for figure reproduction.
  write.csv(markers, paste0("./Rawdata/", out_prefix, "_makers.csv"), row.names = FALSE)
  
  df_summary <- seurat_obj@meta.data %>% select(all_of(ident_col)) %>% distinct() %>% arrange(across(all_of(ident_col)))
  df_top_genes <- markers %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>%
    summarise(top_genes = paste0(gene, collapse = ", ")) %>% rename(!!ident_col := cluster)
  
  final_df <- left_join(df_summary, df_top_genes, by = ident_col)
  write.csv(final_df, paste0("./Rawdata/", out_prefix, "_makers_summary.csv"), row.names = FALSE)
}

# 提取亚群Top10差异基因并导出
export_top_markers <- function(seurat_obj, ident_col, out_prefix) {
  Idents(seurat_obj) <- ident_col
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
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

my_theme <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

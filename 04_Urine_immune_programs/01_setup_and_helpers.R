# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Upstream urine processing and sample-level signature analysis.
# Original source: 00_setup_and_helpers.R
# Figure/table scope: Fig4; FigS7; FigS8
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# File: 00_setup_and_helpers.R
# Description: Shared paths, package loading, plotting style, and helper functions.
# ===============================================================================

PROJECT_DIR <- file.path(getOption("hla.data_root"), "\u5c3f\u6db2\u6837\u672c")
setwd(PROJECT_DIR)
set.seed(2026)

CODE_DIR <- file.path(getOption("hla.code_root"), "04_Urine_immune_programs")
DATA_DIR <- file.path(PROJECT_DIR, "Data")
FIGURE_DIR <- file.path(PROJECT_DIR, "Figures")
RESULT_DIR <- file.path(PROJECT_DIR, "Result")
QC_DIR <- file.path(FIGURE_DIR, "QC")
DEMULTIPLEX_FIGURE_DIR <- file.path(FIGURE_DIR, "Demultiplex")
DEMULTIPLEX_RESULT_DIR <- file.path(RESULT_DIR, "Demultiplex")
HELPER_DIR <- file.path(getOption("hla.code_root"), "00_Shared")
SYN_CODE_DIR <- file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks")

project_path <- function(...) file.path(PROJECT_DIR, ...)
data_path <- function(...) file.path(DATA_DIR, ...)
figure_path <- function(...) file.path(FIGURE_DIR, ...)
result_path <- function(...) file.path(RESULT_DIR, ...)
code_path <- function(...) file.path(CODE_DIR, ...)
object_path <- function(...) result_path("Objects", ...)
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

fig_file <- function(section, number, label, ext = "pdf") {
  clean_label <- gsub("[^A-Za-z0-9_+.-]+", "_", label)
  clean_label <- gsub("_+", "_", clean_label)
  clean_label <- gsub("^_|_$", "", clean_label)
  file.path(ensure_dir(figure_path(section)), sprintf("%02d_%s.%s", number, clean_label, ext))
}

result_file <- function(section, ...) {
  file.path(ensure_dir(result_path(section)), ...)
}

invisible(lapply(c(
  DATA_DIR, FIGURE_DIR, RESULT_DIR, QC_DIR,
  DEMULTIPLEX_FIGURE_DIR, DEMULTIPLEX_RESULT_DIR,
  result_path("Objects"),
  result_path("scRNA", "Markers"),
  result_path("scRNA", "Group_Markers"),
  result_path("Myeloid", "Markers"),
  result_path("Myeloid", "Differential_Expression"),
  result_path("Myeloid", "Other"),
  result_path("Tcell", "Markers"),
  result_path("Tcell", "Differential_Expression"),
  result_path("TCR", "Contigs_Metadata"),
  result_path("TCR", "GLIPH_HLA"),
  result_path("TCR", "Statistics")
), ensure_dir))

source_if_exists <- function(path, required = FALSE) {
  if (file.exists(path)) {
    source(path)
  } else if (required) {
    stop("Required helper script not found: ", path)
  } else {
    warning("Helper script not found, skipped: ", path)
  }
}

load_packages <- function(pkgs) {
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(rlang)
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
  library(scRepertoire)
  library(scales)
  library(immunarch)
  library(ggseqlogo)
  library(circlize)
  library(tidyr)
  library(VennDiagram)
  library(SCIPAC)
  library(MAST)
  library(AUCell)
  library(glmnet)
  library(pROC)
  library(scCustomize)
  library(viridisLite)
})

GLOBAL_COLORS <- c("#E42320", "#699ECA", "#4DAF4A", "#FF8C00", "#F898CB", "#D65190", "#731A73", "#FFCB5B",
                   "#E89B1E", "#0076B9", "#3D505A", "#0098B2", "#8A7355", "#CC5B45", "#008B8B", "#4E5689",
                   "#F9B3AD", "#EDC66A", "#9FDAF7", "#C9A1CA", "#D9DEE7", "#C74546", "#EB7E60", "#8FB4DC",
                   "#2E8B57", "#8B008B", "#FFD700", "#A0522D", "#4682B4", "#DDA0DD", "#556B2F", "#CD5C5C")
GROUP_COLORS <- c("NR" = "#E42320", "TCMR" = "#699ECA", "ABMR" = "#4DAF4A", "Mix" = "#FF8C00", "BK" = "#F898CB", "Undefined" = "grey", "Doublets" = "black")
CLONE_SIZE_COLORS <- c("#f9f871", "#ffc75f", "#ff6f91", "#845ec2", "grey")
comparison <- list(c("NR", "TCMR"), c("NR", "ABMR"), c("NR", "Mix"), c("NR", "BK"))
ANALYSIS_BATCHES <- c("Group1", "Group2", "Group3")

source_if_exists(file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks/01_HLA_demultiplex_v2.R"))
for (helper_file in c(
  "plot_grouped_umap.R", "plot_cell_proportion_barplot.R", "plot_HLA_tSNE.R",
  "plot_celltype_ratio.R", "plot_cell_ratio_comparison.R", "plot_volcano.R",
  "plot_violin_comparison.R", "Enrichment_analysis.R", "plot_pseudobulk_heatmap.R",
  "plot_dotplot.R", "plot_top_pseudobulk.R"
)) {
  source_if_exists(file.path(HELPER_DIR, helper_file))
}

find_10x_dirs <- function(batch) {
  gem_path <- data_path(batch, "GEM")
  if (!dir.exists(gem_path)) return(data.frame())
  candidate_dirs <- unique(c(
    if (dir.exists(file.path(gem_path, "filtered_feature_bc_matrix"))) gem_path else character(),
    dirname(list.files(gem_path, pattern = "^matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE))
  ))
  candidate_dirs <- candidate_dirs[basename(candidate_dirs) == "filtered_feature_bc_matrix"]
  if (!length(candidate_dirs)) return(data.frame())
  data.frame(
    data_dir = normalizePath(candidate_dirs, winslash = "/", mustWork = FALSE),
    run_index = seq_along(candidate_dirs),
    stringsAsFactors = FALSE
  )
}

find_tcr_contig_files <- function(batch, full_length = FALSE) {
  gem_path <- data_path(batch, "GEM")
  if (!dir.exists(gem_path)) return(data.frame())
  pattern <- if (full_length) "filtered_contig_annotations_full_length\\.csv$" else "filtered_contig_annotations\\.csv$"
  files <- list.files(gem_path, pattern = pattern, recursive = TRUE, full.names = TRUE)
  files <- files[grepl("(TCR|vdj_t)", files, ignore.case = TRUE)]
  if (!length(files)) return(data.frame())
  data.frame(file = normalizePath(files, winslash = "/", mustWork = FALSE), run_index = seq_along(files), stringsAsFactors = FALSE)
}

find_hla_demultiplex_files <- function(batch) {
  hla_path <- data_path(batch, "HLA_counts")
  if (!dir.exists(hla_path)) return(data.frame())
  files <- list.files(hla_path, pattern = "HLA_demultiplex\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(data.frame())
  method <- ifelse(grepl("2field|2-field", files, ignore.case = TRUE), "2-field", "1-field")
  data.frame(file = normalizePath(files, winslash = "/", mustWork = FALSE), method = method, run_index = seq_along(files), stringsAsFactors = FALSE)
}

find_existing_file <- function(candidates) {
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) hit[1] else candidates[1]
}

find_object_file <- function(primary, fallback = character()) {
  candidates <- object_path(c(primary, fallback))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) return(hit[1])
  stop("Required object file not found. Checked: ", paste(candidates, collapse = "; "))
}

find_hla_matrix_sets <- function(batch) {
  hla_path <- data_path(batch, "HLA_counts")
  if (!dir.exists(hla_path)) return(data.frame())
  mtx_files <- list.files(hla_path, pattern = "HLA_type\\.mtx$", recursive = TRUE, full.names = TRUE)
  if (!length(mtx_files)) return(data.frame())

  barcode_files <- vapply(mtx_files, function(mtx_file) {
    find_existing_file(c(
      file.path(dirname(mtx_file), "output.barcodes.txt"),
      file.path(dirname(mtx_file), "filter_counts", "output.barcodes.txt"),
      file.path(dirname(dirname(mtx_file)), "output.barcodes.txt"),
      file.path(dirname(dirname(mtx_file)), "filter_counts", "output.barcodes.txt")
    ))
  }, character(1))

  label_files <- vapply(mtx_files, function(mtx_file) {
    find_existing_file(c(
      file.path(dirname(mtx_file), "sorted_genotypes.txt"),
      file.path(dirname(dirname(mtx_file)), "sorted_genotypes.txt"),
      file.path(dirname(dirname(dirname(mtx_file))), "sorted_genotypes.txt")
    ))
  }, character(1))

  data.frame(
    mtx_file = normalizePath(mtx_files, winslash = "/", mustWork = FALSE),
    barcode_file = normalizePath(barcode_files, winslash = "/", mustWork = FALSE),
    label_file = normalizePath(label_files, winslash = "/", mustWork = FALSE),
    run_index = seq_along(mtx_files),
    stringsAsFactors = FALSE
  )
}
# ------------------------------------------------------------------------------
# Custom helper functions
# ------------------------------------------------------------------------------

scRNA_QC <- function(scRNA, min_RNA = 200, max_MT = 20, save_path = figure_path("QC")) {
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
  scRNA <- subset(scRNA, subset = nCount_RNA > min_RNA & percent.mt < max_MT)
  
  # 批量绘制QC小提琴图
  qc_metrics <- list(nFeature = c("nFeature_RNA"), 
                     nCount = c("nCount_RNA"), 
                     MT = c("percent.mt"))
  
  dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  
  for(name in names(qc_metrics)) {
    p <- ggplot(scRNA@meta.data, aes(x = orig.ident, y = .data[[qc_metrics[[name]][1]]])) +
      geom_violin(aes(fill = orig.ident), width = 1, alpha = 0.5) +
      geom_jitter(aes(color = orig.ident), width = 0.2, size = 0.01) +
      geom_boxplot(fill = "white", width = 0.3, outlier.size = 0) +
      theme_bw() +
      theme(legend.position = "none")
    
    ggsave(file.path(ensure_dir(save_path), sprintf("%02d_QC_%s.pdf", match(name, names(qc_metrics)), name)), plot = p, height = 5, width = 8)
  }
  
  return(scRNA)
}


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


# 提取亚群Top10差异基因并导出
export_top_markers <- function(seurat_obj, ident_col, out_prefix) {
  Seurat::Idents(seurat_obj) <- ident_col
  markers <- SeuratWrappers::RunPrestoAll(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  marker_category <- sub("_all$", "", out_prefix)
  marker_dir <- ensure_dir(result_path(marker_category, "Markers"))
# Export the analysis source data for figure reproduction.
  write.csv(markers, file.path(marker_dir, paste0(out_prefix, "_makers.csv")), row.names = FALSE)
  
  gene_col <- dplyr::case_when("gene" %in% colnames(markers) ~ "gene", "feature" %in% colnames(markers) ~ "feature", TRUE ~ NA_character_)
  cluster_col <- dplyr::case_when("cluster" %in% colnames(markers) ~ "cluster", "group" %in% colnames(markers) ~ "group", TRUE ~ NA_character_)
  p_col <- dplyr::case_when("p_val_adj" %in% colnames(markers) ~ "p_val_adj", "padj" %in% colnames(markers) ~ "padj", TRUE ~ NA_character_)
  fc_col <- dplyr::case_when("avg_log2FC" %in% colnames(markers) ~ "avg_log2FC", "logFC" %in% colnames(markers) ~ "logFC", TRUE ~ NA_character_)
  
  if (any(is.na(c(gene_col, cluster_col, p_col, fc_col)))) {
    stop("markers 中缺少必要列，请检查 RunPrestoAll 输出列名：", paste(colnames(markers), collapse = ", "))
  }
  
  df_summary <- seurat_obj@meta.data %>%
    dplyr::select(dplyr::all_of(ident_col)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data[[ident_col]])
  
  df_top_genes <- markers %>%
    dplyr::filter(.data[[p_col]] < 0.01) %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::slice_max(order_by = .data[[fc_col]], n = 10, with_ties = FALSE) %>%
    dplyr::summarise(top_genes = paste0(.data[[gene_col]], collapse = ", "), .groups = "drop") %>%
    dplyr::rename_with(~ ident_col, .cols = dplyr::all_of(cluster_col))
  
  final_df <- dplyr::left_join(df_summary, df_top_genes, by = ident_col)
  write.csv(final_df, file.path(marker_dir, paste0(out_prefix, "_makers_summary.csv")), row.names = FALSE)
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

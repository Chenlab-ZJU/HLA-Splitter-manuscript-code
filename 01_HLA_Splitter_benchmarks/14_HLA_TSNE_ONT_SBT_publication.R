# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Final nanopore typing and reference-concordance figures.
# Original source: 15_revision_HLA_TSNE_ONT_SBT_publication.R
# Figure/table scope: FigS1c; FigS1d; FigS1e; FigS1f
# Data cohort: 30PBMC. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(Rtsne)
  library(patchwork)
})

project_dir <- Sys.getenv("PBMC30_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "30PBMC")
}
project_dir <- gsub("\\\\", "/", project_dir)
data_dir <- file.path(project_dir, "Data")
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

sample_levels <- paste0("Donor", 1:30)
sample_palette <- c(
  "#E42320", "#699ECA", "#4DAF4A", "#FF8C00", "#F898CB", "#D65190",
  "#731A73", "#FFCB5B", "#E89B1E", "#0076B9", "#3D505A", "#0098B2",
  "#7E6148", "#B09C85", "#00A087", "#91D1C2", "#8491B4", "#DC0000",
  "#7CAE00", "#C77CFF", "#00BFC4", "#F8766D", "#619CFF", "#00BA38",
  "#B79F00", "#00B0F6", "#F564E3", "#A3A500", "#39B600", "#9590FF"
)
sample_colors <- setNames(sample_palette, sample_levels)
sample_colors <- c(sample_colors, "Undefined" = "#BFBFBF")

save_plot_safe <- function(plot, path_no_ext, width, height) {
  tmp_pdf <- file.path(tempdir(), paste0(basename(path_no_ext), "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(basename(path_no_ext), "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
  file.copy(tmp_pdf, paste0(path_no_ext, ".pdf"), overwrite = TRUE)
  file.copy(tmp_png, paste0(path_no_ext, ".png"), overwrite = TRUE)
}

read_hla_matrix <- function() {
  hla_dir <- file.path(data_dir, "HLA_counts")
  genotype_file <- file.path(hla_dir, "sorted_genotypes.txt")
  barcode_file <- file.path(hla_dir, "output.barcodes.txt")
  matrix_file <- file.path(hla_dir, "HLA_type.mtx")
  if (!all(file.exists(genotype_file, barcode_file, matrix_file))) {
    stop("Missing HLA count input files in: ", hla_dir)
  }
  genes <- read.delim(genotype_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  barcodes <- read.delim(barcode_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  mat <- Matrix::readMM(matrix_file)
  if (nrow(mat) != length(barcodes) || ncol(mat) != length(genes)) {
    stop("HLA matrix dimensions do not match barcode/gene files.")
  }
  mat <- as.matrix(mat)
  colnames(mat) <- genes
  rownames(mat) <- paste0(barcodes, "-1")
  mat
}

read_demux <- function(method) {
  f <- file.path(result_dir, paste0("HLA_Splitter_", method), "HLA_demultiplex.csv")
  if (!file.exists(f)) {
    stop("Missing demultiplex file: ", f)
  }
  x <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  barcode_col <- if ("Cell.barcode" %in% colnames(x)) "Cell.barcode" else "Cell barcode"
  x %>%
    transmute(
      Cell = paste0(.data[[barcode_col]], "-1"),
      prediction = ifelse(is.na(Sample_prediction) | Sample_prediction == "", "Undefined", as.character(Sample_prediction))
    )
}

zscore_positive_hla <- function(mat) {
  library_sizes <- rowSums(mat)
  keep <- library_sizes > 0
  mat <- mat[keep, , drop = FALSE]
  norm <- sweep(mat, 1, rowSums(mat), "/") * 100
  norm <- log1p(norm)
  norm[!is.finite(norm)] <- 0
  for (j in seq_len(ncol(norm))) {
    pos <- norm[, j] > 0
    if (sum(pos) > 1) {
      norm[pos, j] <- as.numeric(scale(norm[pos, j], center = FALSE, scale = TRUE))
    }
  }
  norm[!is.finite(norm)] <- 0
  norm
}

source_file <- file.path(result_dir, "HLA_TSNE_ONT_SBT_revision_source.csv")
force_tsne <- tolower(Sys.getenv("FORCE_HLA_TSNE_ONT_SBT", unset = "false")) %in% c("true", "1", "yes")
if (file.exists(source_file) && !force_tsne) {
  message("Using cached HLA t-SNE source: ", source_file)
  tsne_df <- read.csv(source_file, stringsAsFactors = FALSE, check.names = FALSE) %>% as_tibble()
} else {
  message("Reconstructing HLA t-SNE coordinates for ONT/SBT comparison...")
  hla_mat <- read_hla_matrix()
  sbt <- read_demux("SBT") %>% rename(HLA_res_SBT = prediction)
  ont <- read_demux("ONT") %>% rename(HLA_res_ONT = prediction)
  pred_df <- inner_join(sbt, ont, by = "Cell")
  common_cells <- intersect(pred_df$Cell, rownames(hla_mat))
  pred_df <- pred_df %>% filter(Cell %in% common_cells)
  hla_mat <- hla_mat[pred_df$Cell, , drop = FALSE]
  hla_norm <- zscore_positive_hla(hla_mat)
  pred_df <- pred_df %>% filter(Cell %in% rownames(hla_norm))
  hla_norm <- hla_norm[pred_df$Cell, , drop = FALSE]

  set.seed(321)
  perplexity_val <- min(max(20, floor((nrow(hla_norm) - 1) / 3)), 50)
  tsne_out <- Rtsne(
    as.matrix(hla_norm),
    dims = 2,
    pca = TRUE,
    max_iter = 1000,
    theta = 0.4,
    perplexity = perplexity_val,
    verbose = FALSE,
    check_duplicates = FALSE
  )
  tsne_df <- as_tibble(tsne_out$Y, .name_repair = "minimal") %>%
    setNames(c("tSNE1", "tSNE2")) %>%
    bind_cols(pred_df) %>%
    mutate(HLA_count = log10(rowSums(hla_mat[Cell, , drop = FALSE]) + 1))
# Export the analysis source data for figure reproduction.
  write.csv(tsne_df, source_file, row.names = FALSE, fileEncoding = "UTF-8")
}

tsne_df <- tsne_df %>%
  mutate(
    HLA_res_SBT = factor(ifelse(HLA_res_SBT %in% sample_levels, HLA_res_SBT, "Undefined"),
                         levels = c(sample_levels, "Undefined")),
    HLA_res_ONT = factor(ifelse(HLA_res_ONT %in% sample_levels, HLA_res_ONT, "Undefined"),
                         levels = c(sample_levels, "Undefined"))
  )

make_panel <- function(df, col, title) {
  label_df <- df %>%
    transmute(tSNE1, tSNE2, sample = .data[[col]]) %>%
    filter(!is.na(sample), sample != "Undefined") %>%
    group_by(sample) %>%
    summarise(
      n_cells = n(),
      tSNE1 = median(tSNE1, na.rm = TRUE),
      tSNE2 = median(tSNE2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_cells >= 20) %>%
    mutate(sample = factor(as.character(sample), levels = sample_levels))

  ggplot(df, aes(tSNE1, tSNE2, colour = .data[[col]])) +
    geom_point(size = 0.20, alpha = 0.62, stroke = 0) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(tSNE1, tSNE2, label = sample),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "bold",
      size = 2.1,
      colour = "black",
      min.segment.length = 0,
      segment.color = "grey50",
      segment.size = 0.14,
      box.padding = 0.16,
      point.padding = 0.08,
      max.overlaps = Inf,
      seed = 20260828
    ) +
    scale_colour_manual(
      values = sample_colors,
      breaks = sample_levels,
      limits = c(sample_levels, "Undefined"),
      drop = FALSE,
      name = "Sample"
    ) +
    labs(title = title, x = "t-SNE 1", y = "t-SNE 2") +
    guides(colour = guide_legend(
      ncol = 2,
      override.aes = list(size = 5.2, alpha = 1),
      byrow = TRUE
    )) +
    theme_classic(base_size = 8.8, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 10.6, face = "bold", hjust = 0.5, colour = "black"),
      axis.title = element_text(size = 9.2, face = "bold", colour = "black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.36),
      legend.title = element_text(size = 12.0, face = "bold", colour = "black"),
      legend.text = element_text(size = 10.2, colour = "black"),
      legend.key = element_blank(),
      legend.key.size = unit(0.48, "cm"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

p_sbt <- make_panel(tsne_df, "HLA_res_SBT", "SBT")
p_ont <- make_panel(tsne_df, "HLA_res_ONT", "ONT-Rapid")

p_combined <- (p_sbt | p_ont) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

save_plot_safe(p_combined, file.path(figure_dir, "HLA_TSNE_ONT_SBT"), width = 10.4, height = 4.2)

message("Saved HLA_TSNE_ONT_SBT with shared legend and Donor1-Donor30 order.")

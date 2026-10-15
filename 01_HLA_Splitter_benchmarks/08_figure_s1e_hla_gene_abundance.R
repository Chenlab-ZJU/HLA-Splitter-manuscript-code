#!/usr/bin/env Rscript
# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Final benchmark / HLA-signal figures.
# Original source: 14_revision_figure_s1e_hla_gene_abundance.R
# Figure/table scope: Fig1e; Fig1f; Fig1g; Fig1h; FigS2a; FigS2b; FigS2c
# Data cohort: synthetic. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# Rebuild Figure S1e: HLA gene/locus abundance statistics from Syn_data HLA counts.

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

project_dir <- Sys.getenv("SYN_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "synthetic")
}
project_dir <- gsub("\\\\", "/", project_dir)

hla_dir <- file.path(project_dir, "Data", "HLA_counts")
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures")
script_figure_dir <- file.path(getOption("hla.project_root"), "Script/Figures")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_figure_dir, recursive = TRUE, showWarnings = FALSE)

locus_levels <- c("A", "B", "C", "DRB1", "DQB1", "DPB1")
hla_colors <- c(
  "A" = "#F8766D",
  "B" = "#B79F00",
  "C" = "#00BA38",
  "DRB1" = "#00BFC4",
  "DQB1" = "#619CFF",
  "DPB1" = "#F564E3"
)

theme_s1e <- function(base_size = 10.5) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black", size = base_size + 1.5),
      axis.text = element_text(colour = "black", size = base_size),
      axis.line = element_line(colour = "black", linewidth = 0.38),
      axis.ticks = element_line(colour = "black", linewidth = 0.34),
      panel.border = element_rect(fill = NA, colour = "grey35", linewidth = 0.38),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.30),
      panel.grid.minor = element_line(colour = "grey94", linewidth = 0.22),
      legend.title = element_text(face = "bold", size = base_size + 0.5),
      legend.text = element_text(size = base_size),
      legend.key = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

read_one_pool <- function(pool_dir) {
  sample_id <- basename(pool_dir)
  genotype_file <- file.path(pool_dir, "HLA_genotype.tsv")
  barcode_file <- file.path(pool_dir, "barcodes.tsv")
  matrix_file <- file.path(pool_dir, "HLA_count_matrix.mtx")
  if (!all(file.exists(genotype_file, barcode_file, matrix_file))) {
    warning("Skipping incomplete HLA count folder: ", pool_dir)
    return(tibble())
  }

  alleles <- read.delim(genotype_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  barcodes <- read.delim(barcode_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  hla_mat <- Matrix::readMM(matrix_file)
  colnames(hla_mat) <- alleles
  rownames(hla_mat) <- barcodes

  bind_rows(lapply(locus_levels, function(locus) {
    allele_idx <- grepl(paste0("^", locus, "\\*"), alleles)
    if (!any(allele_idx)) return(tibble())
    locus_counts_per_cell <- Matrix::rowSums(hla_mat[, allele_idx, drop = FALSE])
    mean_count_per_cell <- mean(locus_counts_per_cell)
    median_count_per_cell <- median(locus_counts_per_cell)
    tibble(
      sample = sample_id,
      HLA_allele = locus,
      n_cells = nrow(hla_mat),
      n_alleles = sum(allele_idx),
      total_HLA_count = sum(locus_counts_per_cell),
      mean_count_per_cell = mean_count_per_cell,
      median_count_per_cell = median_count_per_cell,
      positive_cell_fraction = mean(locus_counts_per_cell > 0),
      log10_mean_count_per_cell = log10(mean_count_per_cell)
    )
  }))
}

abundance_source <- file.path(getOption("hla.data_root"), "synthetic", "Figure_S1e_HLA_gene_abundance_source.csv")
if (file.exists(abundance_source)) {
  # Reuse the final per-sample source table in the compact processed-data bundle.
  hla_abundance <- read.csv(abundance_source, stringsAsFactors = FALSE)
} else {
  pool_dirs <- list.dirs(hla_dir, full.names = TRUE, recursive = FALSE)
  if (length(pool_dirs) == 0) {
    stop("Provide the HLA abundance source table or HLA count folders in: ", hla_dir)
  }
  hla_abundance <- bind_rows(lapply(pool_dirs, read_one_pool))
}

hla_abundance <- hla_abundance %>%
  filter(HLA_allele %in% locus_levels, is.finite(log10_mean_count_per_cell)) %>%
  mutate(HLA_allele = factor(HLA_allele, levels = locus_levels)) %>%
  arrange(HLA_allele, sample)

summary_df <- hla_abundance %>%
  group_by(HLA_allele) %>%
  summarise(
    n_samples = n(),
    mean_log10_count = mean(log10_mean_count_per_cell, na.rm = TRUE),
    median_log10_count = median(log10_mean_count_per_cell, na.rm = TRUE),
    sd_log10_count = sd(log10_mean_count_per_cell, na.rm = TRUE),
    mean_positive_cell_fraction = mean(positive_cell_fraction, na.rm = TRUE),
    .groups = "drop"
  )

# Export the analysis source data for figure reproduction.
write.csv(hla_abundance, file.path(getOption("hla.data_root"), "synthetic", "Figure_S1e_HLA_gene_abundance_source.csv"), row.names = FALSE)
write.csv(summary_df, file.path(result_dir, "Figure_S1e_HLA_gene_abundance_summary.csv"), row.names = FALSE)

set.seed(20260910)
p_s1e <- ggplot(hla_abundance, aes(x = HLA_allele, y = log10_mean_count_per_cell, fill = HLA_allele, colour = HLA_allele)) +
  geom_violin(width = 0.74, alpha = 0.40, linewidth = 0.55, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", alpha = 0.82, colour = "grey20", linewidth = 0.40) +
  geom_point(
    position = position_jitter(width = 0.16, height = 0, seed = 20260910),
    size = 1.20,
    alpha = 0.92,
    stroke = 0
  ) +
  scale_fill_manual(values = hla_colors, name = "HLA allele", drop = FALSE) +
  scale_colour_manual(values = hla_colors, name = "HLA allele", drop = FALSE) +
  scale_y_continuous(
    limits = c(-2.2, 2.2),
    breaks = seq(-2, 2, by = 1),
    minor_breaks = seq(-2, 2, by = 0.5),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  labs(
    x = "HLA allele",
    y = expression(Log[10]~"(mean HLA count per cell)")
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.45, size = 3.2)),
    colour = "none"
  ) +
  theme_s1e(10.8) +
  theme(
    axis.text.x = element_text(size = 10.8),
    legend.position = "right",
    legend.key.height = unit(0.34, "in"),
    legend.key.width = unit(0.28, "in")
  )

save_plot <- function(plot, filename, width = 4.25, height = 3.45, dpi = 600) {
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(filename, ".png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf, bg = "white")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  }
  file.copy(tmp_pdf, pdf_path, overwrite = TRUE)
  file.copy(tmp_png, png_path, overwrite = TRUE)
}

save_plot(p_s1e, "Figure_S1e_HLA_gene_abundance")
message("Saved Figure S1e HLA gene abundance plot and source data.")

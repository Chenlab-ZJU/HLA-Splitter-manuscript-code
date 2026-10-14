#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Final urine T-cell, CXCR3 and expanded-clone figures.
# Original source: 33_revision_tcell_expanded_clone_fraction_by_group_celltype.R
# Figure/table scope: Fig5a; Fig5b; Fig5c; Fig5d; Fig5e; Fig5f; Fig5g; Fig5h; Fig5i; FigS9a; FigS9b; FigS9c; FigS9d; FigS9e; FigS9f; FigS9g
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Expanded T cells: clone_ratio > 1; distinguish cell counts from clonotype counts.
# Keep within-sample TCR denominators and clonotype identifiers as defined below.


# Expanded TCR clone fraction across urine T-cell subtypes and clinical groups.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "33_revision_tcell_expanded_clone_fraction_by_group_celltype.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

object_file <- file.path(project_dir, "Result", "Objects", "Urine_Tcell.rds")
figure_dir <- file.path(project_dir, "Figures", "TCR", "Revision_TCR_expanded_clone_fraction_celltype")
result_dir <- file.path(project_dir, "Result", "TCR", "Revision_TCR_expanded_clone_fraction_celltype")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "Mix" = "#FF8C00",
  "BK" = "#F898CB"
)

tcell_levels <- c(
  "Naive CD4", "CD4 Tem", "CD4 Tex", "Treg",
  "Naive CD8", "CXCR4+ CD8", "GZMB+ CD8", "ISG+ CD8", "Cycling T"
)

expanded_threshold <- 1
min_celltype_cells_for_fraction <- 5

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

save_plot_safe <- function(plot, filename, width, height) {
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(filename, ".png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
  file.copy(tmp_pdf, pdf_path, overwrite = TRUE)
  file.copy(tmp_png, png_path, overwrite = TRUE)
  invisible(c(pdf = pdf_path, png = png_path))
}

theme_pub <- function(base_size = 9.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.45),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      legend.title = element_text(face = "bold", colour = "black"),
      legend.text = element_text(colour = "black"),
      legend.key = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

if (!file.exists(object_file)) {
  stop("Cannot find urine T-cell object: ", object_file)
}

message("Loading urine T-cell object...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(object_file)

meta <- tcell@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
    Group = as.character(Group),
    T_celltype = as.character(revised_celltype),
    clone_ratio = suppressWarnings(as.numeric(clone_ratio)),
    has_tcr_clone = !is.na(clone_ratio) & clone_ratio > 0,
    expanded_clone = !is.na(clone_ratio) & clone_ratio > expanded_threshold
  ) %>%
  filter(!is.na(Sample), Sample != "", !is.na(Group), Group %in% group_levels, !is.na(T_celltype), T_celltype != "") %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    T_celltype = factor(T_celltype, levels = c(tcell_levels[tcell_levels %in% unique(T_celltype)], setdiff(sort(unique(T_celltype)), tcell_levels)))
  )

sample_celltype_fraction <- meta %>%
  count(Group, Sample, T_celltype, expanded_clone, name = "n_cells") %>%
  complete(
    nesting(Group, Sample),
    T_celltype = factor(levels(meta$T_celltype), levels = levels(meta$T_celltype)),
    expanded_clone = c(FALSE, TRUE),
    fill = list(n_cells = 0)
  ) %>%
  group_by(Group, Sample, T_celltype) %>%
  summarise(
    n_total_t_cells = sum(n_cells),
    n_expanded_cells = sum(n_cells[expanded_clone], na.rm = TRUE),
    expanded_clone_fraction = ifelse(n_total_t_cells > 0, n_expanded_cells / n_total_t_cells, 0),
    expanded_clone_percent = expanded_clone_fraction * 100,
    enough_cells_for_fraction = n_total_t_cells >= min_celltype_cells_for_fraction,
    .groups = "drop"
  )

tcr_detected_fraction <- meta %>%
  count(Group, Sample, T_celltype, has_tcr_clone, name = "n_cells") %>%
  complete(
    nesting(Group, Sample),
    T_celltype = factor(levels(meta$T_celltype), levels = levels(meta$T_celltype)),
    has_tcr_clone = c(FALSE, TRUE),
    fill = list(n_cells = 0)
  ) %>%
  group_by(Group, Sample, T_celltype) %>%
  summarise(
    n_total_t_cells = sum(n_cells),
    n_tcr_detected_cells = sum(n_cells[has_tcr_clone], na.rm = TRUE),
    tcr_detected_fraction = ifelse(n_total_t_cells > 0, n_tcr_detected_cells / n_total_t_cells, 0),
    .groups = "drop"
  )

plot_sample_celltype_fraction <- sample_celltype_fraction %>%
  filter(enough_cells_for_fraction)

group_celltype_summary <- plot_sample_celltype_fraction %>%
  group_by(Group, T_celltype) %>%
  summarise(
    n_samples = sum(!is.na(expanded_clone_percent)),
    mean_expanded_clone_percent = mean(expanded_clone_percent, na.rm = TRUE),
    median_expanded_clone_percent = median(expanded_clone_percent, na.rm = TRUE),
    sd_expanded_clone_percent = sd(expanded_clone_percent, na.rm = TRUE),
    sem_expanded_clone_percent = sd_expanded_clone_percent / sqrt(n_samples),
    total_t_cells = sum(n_total_t_cells, na.rm = TRUE),
    total_expanded_cells = sum(n_expanded_cells, na.rm = TRUE),
    pooled_expanded_clone_percent = ifelse(total_t_cells > 0, total_expanded_cells / total_t_cells * 100, NA_real_),
    .groups = "drop"
  ) %>%
  complete(
    Group = factor(group_levels, levels = group_levels),
    T_celltype = factor(levels(meta$T_celltype), levels = levels(meta$T_celltype)),
    fill = list(
      n_samples = 0,
      mean_expanded_clone_percent = 0,
      median_expanded_clone_percent = 0,
      sd_expanded_clone_percent = 0,
      sem_expanded_clone_percent = 0,
      total_t_cells = 0,
      total_expanded_cells = 0,
      pooled_expanded_clone_percent = 0
    )
  ) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    T_celltype = factor(as.character(T_celltype), levels = levels(meta$T_celltype)),
    celltype_index = as.numeric(T_celltype),
    group_index = as.numeric(Group),
    x_pos = celltype_index + (group_index - mean(seq_along(group_levels))) * 0.15
  )

# Export the analysis source data for figure reproduction.
write_csv_safe(meta, file.path(result_dir, "TCell_expanded_clone_fraction_cell_metadata_source.csv"))
write_csv_safe(sample_celltype_fraction, file.path(result_dir, "TCell_sample_celltype_expanded_clone_fraction.csv"))
write_csv_safe(tcr_detected_fraction, file.path(result_dir, "TCell_sample_celltype_TCR_detected_fraction.csv"))
write_csv_safe(group_celltype_summary, file.path(result_dir, "TCell_group_celltype_expanded_clone_fraction_summary.csv"))

set.seed(20260826)
plot_sample_celltype_fraction <- plot_sample_celltype_fraction %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    T_celltype = factor(as.character(T_celltype), levels = levels(meta$T_celltype)),
    celltype_index = as.numeric(T_celltype),
    group_index = as.numeric(Group),
    x_pos = celltype_index + (group_index - mean(seq_along(group_levels))) * 0.15 + runif(n(), -0.018, 0.018)
  )

p_bar <- ggplot(group_celltype_summary, aes(x = x_pos, y = mean_expanded_clone_percent, fill = Group)) +
  geom_col(
    width = 0.13,
    colour = "black",
    linewidth = 0.24,
    alpha = 0.90
  ) +
  geom_point(
    data = plot_sample_celltype_fraction,
    aes(x = x_pos, y = expanded_clone_percent, fill = Group),
    shape = 21,
    size = 1.75,
    stroke = 0.24,
    colour = "black"
  ,
    alpha = 0.92
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE, name = "Group") +
  scale_x_continuous(
    breaks = seq_along(levels(meta$T_celltype)),
    labels = levels(meta$T_celltype),
    expand = expansion(add = 0.42)
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    breaks = c(0, 2, 5, seq(10, 60, by = 10)),
    minor_breaks = NULL,
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    x = "T-cell subtype",
    y = paste0("Expanded clone fraction\n(clone ratio > ", expanded_threshold, ", %)")
  ) +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7.8),
    axis.text.y = element_text(size = 8.2),
    legend.position = "top",
    legend.title = element_text(size = 8.8),
    legend.text = element_text(size = 8.4),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_plot_safe(
  p_bar,
  "TCell_group_celltype_expanded_clone_fraction_barplot_revision",
  width = 8.4,
  height = 4.2
)

message("Saved urine T-cell expanded clone fraction barplot.")
message("Tables: ", result_dir)
message("Figures: ", figure_dir)

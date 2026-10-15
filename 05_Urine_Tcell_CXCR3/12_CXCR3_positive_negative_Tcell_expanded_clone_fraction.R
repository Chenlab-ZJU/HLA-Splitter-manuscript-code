#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Final urine T-cell, CXCR3 and expanded-clone figures.
# Original source: 34_revision_CXCR3_positive_negative_Tcell_expanded_clone_fraction.R
# Figure/table scope: Fig5a; Fig5b; Fig5c; Fig5d; Fig5e; Fig5f; Fig5g; Fig5h; Fig5i; FigS9a; FigS9b; FigS9c; FigS9d; FigS9e; FigS9f; FigS9g
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Expanded T cells: clone_ratio > 1; distinguish cell counts from clonotype counts.
# Keep within-sample TCR denominators and clonotype identifiers as defined below.


# Compare expanded clone fraction between CXCR3+ and CXCR3- urine T cells.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "34_revision_CXCR3_positive_negative_Tcell_expanded_clone_fraction.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

object_file <- file.path(getOption("hla.data_root"), "urine", "Urine_Tcell.rds")
figure_dir <- file.path(project_dir, "Figures", "TCR", "Revision_CXCR3_Tcell_expanded_clone_fraction")
result_dir <- file.path(project_dir, "Result", "TCR", "Revision_CXCR3_Tcell_expanded_clone_fraction")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
status_levels <- c("CXCR3- T", "CXCR3+ T")
status_colors <- c("CXCR3- T" = "#BDBDBD", "CXCR3+ T" = "#B2182B")

expanded_threshold <- 1
min_status_cells_for_fraction <- 5
status_offset <- 0.26

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
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

format_p <- function(p) {
  ifelse(
    is.na(p), "p=NA",
    ifelse(p < 0.001, "p<0.001", paste0("p=", formatC(p, format = "f", digits = 2)))
  )
}

safe_paired_wilcox <- function(df) {
  wide <- df %>%
    select(Sample, CXCR3_status, expanded_clone_percent) %>%
    pivot_wider(names_from = CXCR3_status, values_from = expanded_clone_percent) %>%
    filter(!is.na(.data[["CXCR3- T"]]), !is.na(.data[["CXCR3+ T"]]))
  if (nrow(wide) < 2) return(NA_real_)
  if (length(unique(c(wide[["CXCR3- T"]], wide[["CXCR3+ T"]]))) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(wide[["CXCR3+ T"]], wide[["CXCR3- T"]], paired = TRUE, exact = FALSE)$p.value)
}

if (!file.exists(object_file)) {
  stop("Cannot find urine T-cell object: ", object_file)
}

message("Loading urine T-cell object...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(object_file)
DefaultAssay(tcell) <- "RNA"
if (!"CXCR3" %in% rownames(tcell)) {
  stop("CXCR3 was not found in Urine_Tcell.rds.")
}

cxcr3_expr <- FetchData(tcell, vars = "CXCR3") %>%
  rownames_to_column("cell") %>%
  as_tibble()

cell_df <- tcell@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
    Group = as.character(Group),
    clone_ratio = suppressWarnings(as.numeric(clone_ratio))
  ) %>%
  left_join(cxcr3_expr, by = "cell") %>%
  filter(!is.na(Sample), Sample != "", !is.na(Group), Group %in% group_levels) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    CXCR3_status = factor(ifelse(CXCR3 > 0, "CXCR3+ T", "CXCR3- T"), levels = status_levels),
    expanded_clone = !is.na(clone_ratio) & clone_ratio > expanded_threshold
  )

sample_status_fraction <- cell_df %>%
  count(Group, Sample, CXCR3_status, expanded_clone, name = "n_cells") %>%
  complete(
    nesting(Group, Sample),
    CXCR3_status = factor(status_levels, levels = status_levels),
    expanded_clone = c(FALSE, TRUE),
    fill = list(n_cells = 0)
  ) %>%
  group_by(Group, Sample, CXCR3_status) %>%
  summarise(
    n_t_cells = sum(n_cells),
    n_expanded_cells = sum(n_cells[expanded_clone], na.rm = TRUE),
    expanded_clone_fraction = ifelse(n_t_cells > 0, n_expanded_cells / n_t_cells, 0),
    expanded_clone_percent = expanded_clone_fraction * 100,
    enough_cells_for_fraction = n_t_cells >= min_status_cells_for_fraction,
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    CXCR3_status = factor(as.character(CXCR3_status), levels = status_levels)
  )

plot_sample_status_fraction <- sample_status_fraction %>%
  filter(enough_cells_for_fraction)

group_status_summary <- plot_sample_status_fraction %>%
  group_by(Group, CXCR3_status) %>%
  summarise(
    n_samples = n_distinct(Sample),
    mean_expanded_clone_percent = mean(expanded_clone_percent, na.rm = TRUE),
    median_expanded_clone_percent = median(expanded_clone_percent, na.rm = TRUE),
    sd_expanded_clone_percent = sd(expanded_clone_percent, na.rm = TRUE),
    total_t_cells = sum(n_t_cells, na.rm = TRUE),
    total_expanded_cells = sum(n_expanded_cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    Group = factor(group_levels, levels = group_levels),
    CXCR3_status = factor(status_levels, levels = status_levels),
    fill = list(
      n_samples = 0,
      mean_expanded_clone_percent = 0,
      median_expanded_clone_percent = 0,
      sd_expanded_clone_percent = 0,
      total_t_cells = 0,
      total_expanded_cells = 0
    )
  ) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    CXCR3_status = factor(as.character(CXCR3_status), levels = status_levels),
    group_index = as.numeric(Group),
    status_index = as.numeric(CXCR3_status),
    x_pos = group_index + (status_index - mean(seq_along(status_levels))) * status_offset
  )

paired_stats <- plot_sample_status_fraction %>%
  group_by(Group) %>%
  summarise(
    group1 = "CXCR3- T",
    group2 = "CXCR3+ T",
    p.value = safe_paired_wilcox(pick(everything())),
    .groups = "drop"
  ) %>%
  mutate(
    p.adjust = p.adjust(p.value, method = "BH"),
    label = format_p(p.value),
    Group = factor(as.character(Group), levels = group_levels),
    group_index = as.numeric(Group),
    x_start = group_index - status_offset / 2,
    x_end = group_index + status_offset / 2
  ) %>%
  left_join(
    plot_sample_status_fraction %>%
      group_by(Group) %>%
      summarise(y_max = max(expanded_clone_percent, na.rm = TRUE), .groups = "drop"),
    by = "Group"
  ) %>%
  mutate(
    y_line = y_max + pmax(y_max * 0.08, 2),
    y_tick = y_line - pmax(y_max * 0.025, 0.8),
    y_text = y_line + pmax(y_max * 0.035, 1.0)
  )

set.seed(20260827)
plot_sample_status_fraction <- plot_sample_status_fraction %>%
  mutate(
    group_index = as.numeric(Group),
    status_index = as.numeric(CXCR3_status),
    x_pos = group_index + (status_index - mean(seq_along(status_levels))) * status_offset + runif(n(), -0.025, 0.025)
  )

# Export the analysis source data for figure reproduction.
write_csv_safe(cell_df, file.path(result_dir, "TCell_CXCR3_status_expanded_clone_cell_source.csv"))
write_csv_safe(sample_status_fraction, file.path(result_dir, "TCell_CXCR3_status_sample_expanded_clone_fraction.csv"))
write_csv_safe(group_status_summary, file.path(result_dir, "TCell_CXCR3_status_group_expanded_clone_fraction_summary.csv"))
write_csv_safe(paired_stats, file.path(result_dir, "TCell_CXCR3_positive_vs_negative_paired_wilcox_by_group.csv"))

p <- ggplot(group_status_summary, aes(x = x_pos, y = mean_expanded_clone_percent, fill = CXCR3_status)) +
  geom_col(width = 0.22, colour = "black", linewidth = 0.28, alpha = 0.88) +
  geom_point(
    data = plot_sample_status_fraction,
    aes(x = x_pos, y = expanded_clone_percent, fill = CXCR3_status),
    shape = 21,
    size = 1.9,
    stroke = 0.25,
    colour = "black",
    alpha = 0.92
  ) +
  geom_segment(
    data = paired_stats %>% filter(!is.na(p.value)),
    aes(x = x_start, xend = x_end, y = y_line, yend = y_line),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = paired_stats %>% filter(!is.na(p.value)),
    aes(x = x_start, xend = x_start, y = y_tick, yend = y_line),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = paired_stats %>% filter(!is.na(p.value)),
    aes(x = x_end, xend = x_end, y = y_tick, yend = y_line),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_text(
    data = paired_stats %>% filter(!is.na(p.value)),
    aes(x = group_index, y = y_text, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 2.55
  ) +
  scale_fill_manual(values = status_colors, drop = FALSE, name = NULL) +
  scale_x_continuous(
    breaks = seq_along(group_levels),
    labels = group_levels,
    expand = expansion(add = 0.55)
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    breaks = c(0, 2, 5, seq(10, 100, by = 10)),
    expand = expansion(mult = c(0, 0.16))
  ) +
  labs(
    x = "Group",
    y = paste0("Expanded clone fraction\n(clone ratio > ", expanded_threshold, ", %)")
  ) +
  theme_pub(9.4) +
  theme(
    axis.text.x = element_text(size = 8.8, face = "bold"),
    axis.text.y = element_text(size = 8.4),
    legend.position = "top",
    legend.text = element_text(size = 8.8),
    panel.grid = element_blank()
  )

save_plot_safe(
  p,
  "TCell_CXCR3_positive_negative_expanded_clone_fraction_barplot_revision",
  width = 5.8,
  height = 3.9
)

message("Saved CXCR3+ versus CXCR3- urine T-cell expanded clone fraction barplot.")
message("Tables: ", result_dir)
message("Figures: ", figure_dir)

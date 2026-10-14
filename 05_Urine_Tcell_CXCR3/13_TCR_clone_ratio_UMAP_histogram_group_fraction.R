#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Final urine T-cell, CXCR3 and expanded-clone figures.
# Original source: 35_revision_TCR_clone_ratio_UMAP_histogram_group_fraction.R
# Figure/table scope: Fig5a; Fig5b; Fig5c; Fig5d; Fig5e; Fig5f; Fig5g; Fig5h; Fig5i; FigS9a; FigS9b; FigS9c; FigS9d; FigS9e; FigS9f; FigS9g
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Expanded T cells: clone_ratio > 1; distinguish cell counts from clonotype counts.
# Keep within-sample TCR denominators and clonotype identifiers as defined below.


# TCR clone_ratio UMAP, per-sample histograms, and group-level expanded clone fraction.

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
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "35_revision_TCR_clone_ratio_UMAP_histogram_group_fraction.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

object_file <- file.path(project_dir, "Result", "Objects", "Urine_Tcell.rds")
figure_dir <- file.path(project_dir, "Figures", "TCR", "Revision_TCR_clone_ratio_distribution")
result_dir <- file.path(project_dir, "Result", "TCR", "Revision_TCR_clone_ratio_distribution")
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

expanded_threshold <- 1
histogram_xmax <- 10
top_expanded_clones_n <- 30

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

theme_pub <- function(base_size = 9.0) {
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
      strip.background = element_rect(fill = "grey93", colour = NA),
      strip.text = element_text(face = "bold", colour = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

format_p <- function(p, prefix = "p") {
  ifelse(
    is.na(p), paste0(prefix, "=NA"),
    ifelse(p < 0.001, paste0(prefix, "<0.001"), paste0(prefix, "=", formatC(p, format = "f", digits = 2)))
  )
}

safe_wilcox <- function(df, group_a, group_b, value_col = "expanded_clone_percent", group_col = "Group") {
  df <- df %>%
    filter(.data[[group_col]] %in% c(group_a, group_b), !is.na(.data[[value_col]])) %>%
    mutate(.group_tmp = droplevels(factor(.data[[group_col]], levels = c(group_a, group_b))))
  if (n_distinct(df$.group_tmp) < 2) return(NA_real_)
  if (length(unique(df[[value_col]])) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(df[[value_col]] ~ df$.group_tmp, exact = FALSE)$p.value)
}

if (!file.exists(object_file)) {
  stop("Cannot find urine T-cell object: ", object_file)
}

message("Loading urine T-cell object...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(object_file)
if (!"umap" %in% Reductions(tcell)) {
  stop("UMAP reduction was not found in Urine_Tcell.rds.")
}
if (!"clone_ratio" %in% colnames(tcell@meta.data)) {
  stop("clone_ratio was not found in Urine_Tcell.rds metadata.")
}

umap <- Embeddings(tcell, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell")
colnames(umap)[2:3] <- c("UMAP_1", "UMAP_2")

cell_df <- tcell@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
    Group = as.character(Group),
    T_celltype = as.character(revised_celltype),
    clonalFrequency = suppressWarnings(as.numeric(clonalFrequency)),
    clone_ratio = suppressWarnings(as.numeric(clone_ratio))
  ) %>%
  left_join(umap, by = "cell") %>%
  filter(!is.na(Sample), Sample != "", !is.na(Group), Group %in% group_levels, !is.na(UMAP_1), !is.na(UMAP_2)) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    clone_ratio = ifelse(is.na(clone_ratio), 0, clone_ratio),
    clone_ratio_log10p1 = log10(clone_ratio + 1),
    expanded_clone = clone_ratio > expanded_threshold,
    expanded_label = factor(ifelse(expanded_clone, "Expanded clone", "Non-expanded"), levels = c("Non-expanded", "Expanded clone"))
  )

sample_order <- cell_df %>%
  distinct(Group, Sample) %>%
  arrange(Group, Sample) %>%
  pull(Sample)

sample_fraction <- cell_df %>%
  group_by(Group, Sample) %>%
  summarise(
    n_t_cells = n(),
    n_tcr_detected = sum(clone_ratio > 0, na.rm = TRUE),
    n_expanded_cells = sum(expanded_clone, na.rm = TRUE),
    expanded_clone_fraction = n_expanded_cells / n_t_cells,
    expanded_clone_percent = expanded_clone_fraction * 100,
    mean_clone_ratio = mean(clone_ratio, na.rm = TRUE),
    median_clone_ratio = median(clone_ratio, na.rm = TRUE),
    max_clone_ratio = max(clone_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Group = factor(as.character(Group), levels = group_levels))

group_fraction_summary <- sample_fraction %>%
  group_by(Group) %>%
  summarise(
    n_samples = n_distinct(Sample),
    mean_expanded_clone_percent = mean(expanded_clone_percent, na.rm = TRUE),
    median_expanded_clone_percent = median(expanded_clone_percent, na.rm = TRUE),
    sd_expanded_clone_percent = sd(expanded_clone_percent, na.rm = TRUE),
    sem_expanded_clone_percent = sd_expanded_clone_percent / sqrt(n_samples),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    group_index = as.numeric(Group)
  )

pairwise_stats <- tidyr::crossing(
  group1 = group_levels,
  group2 = group_levels
) %>%
  filter(as.numeric(factor(group1, levels = group_levels)) < as.numeric(factor(group2, levels = group_levels))) %>%
  rowwise() %>%
  mutate(p.value = safe_wilcox(sample_fraction, group1, group2)) %>%
  ungroup() %>%
  mutate(
    p.adjust = p.adjust(p.value, method = "BH"),
    label = format_p(p.value)
  )

plot_stats <- pairwise_stats %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  mutate(
    x_start = as.numeric(factor(group1, levels = group_levels)),
    x_end = as.numeric(factor(group2, levels = group_levels))
  )

y_base <- max(sample_fraction$expanded_clone_percent, na.rm = TRUE)
y_step <- max(y_base * 0.10, 2.5)
plot_stats <- plot_stats %>%
  arrange(x_end) %>%
  mutate(
    y_line = y_base + row_number() * y_step,
    y_tick = y_line - y_step * 0.22,
    y_text = y_line + y_step * 0.34
  )

formula_text <- paste(
  "TCR clone_ratio calculation",
  "",
  "For each sample, a per-sample baseline clonal frequency is calculated as:",
  "  baseline_sample = mean(clonalFrequency of cells with clonalFrequency > 0 within that sample)",
  "",
  "For each T cell:",
  "  clone_ratio_cell = clonalFrequency_cell / baseline_sample",
  "",
  "If a sample has no positive clonalFrequency values, baseline_sample is set to 1.",
  paste0("Expanded clone definition used here: clone_ratio > ", expanded_threshold, "."),
  sep = "\n"
)
writeLines(formula_text, con = file.path(result_dir, "TCR_clone_ratio_formula_and_threshold.txt"), useBytes = TRUE)

# Export the analysis source data for figure reproduction.
write_csv_safe(cell_df, file.path(result_dir, "TCR_clone_ratio_cell_UMAP_histogram_source.csv"))
write_csv_safe(sample_fraction, file.path(result_dir, "TCR_clone_ratio_sample_expanded_fraction.csv"))
write_csv_safe(group_fraction_summary, file.path(result_dir, "TCR_clone_ratio_group_expanded_fraction_summary.csv"))
write_csv_safe(pairwise_stats, file.path(result_dir, "TCR_clone_ratio_group_expanded_fraction_pairwise_wilcox.csv"))

umap_bg <- cell_df %>% filter(!expanded_clone)
umap_fg <- cell_df %>% filter(expanded_clone) %>% arrange(clone_ratio)

p_umap <- ggplot() +
  geom_point(data = umap_bg, aes(UMAP_1, UMAP_2), colour = "grey68", size = 0.30, alpha = 0.68, stroke = 0) +
  geom_point(
    data = umap_fg,
    aes(UMAP_1, UMAP_2, colour = clone_ratio_log10p1, size = clone_ratio),
    alpha = 0.95,
    stroke = 0
  ) +
  scale_colour_gradientn(
    colours = c("#F7F7F7", "#FEE08B", "#F46D43", "#B2182B"),
    name = "log10(clone\nratio + 1)"
  ) +
  scale_size_continuous(range = c(0.55, 2.6), trans = "sqrt", name = "clone_ratio") +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.0) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )
save_plot_safe(p_umap, "TCR_clone_ratio_UMAP_revision", width = 4.8, height = 4.2)

p_expanded_umap <- ggplot() +
  geom_point(data = umap_bg, aes(UMAP_1, UMAP_2), colour = "grey68", size = 0.30, alpha = 0.68, stroke = 0) +
  geom_point(
    data = umap_fg,
    aes(UMAP_1, UMAP_2, size = clone_ratio),
    colour = "#B2182B",
    alpha = 0.94,
    stroke = 0
  ) +
  scale_size_continuous(range = c(0.60, 2.8), trans = "sqrt", name = "clone_ratio") +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.0) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )
save_plot_safe(p_expanded_umap, "TCR_expanded_T_cells_UMAP_revision", width = 4.8, height = 4.2)

expanded_clone_rank <- cell_df %>%
  filter(expanded_clone) %>%
  mutate(
    clonotype_id = dplyr::coalesce(as.character(tcell@meta.data[cell, "CTaa"]), as.character(tcell@meta.data[cell, "CTgene"]), cell),
    clonotype_id = ifelse(is.na(clonotype_id) | clonotype_id == "", cell, clonotype_id),
    clone_id = paste(Sample, clonotype_id, sep = "__")
  ) %>%
  group_by(Group, Sample, clone_id, clonotype_id) %>%
  summarise(
    expanded_cells = n(),
    mean_clone_ratio = mean(clone_ratio, na.rm = TRUE),
    max_clone_ratio = max(clone_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(expanded_cells), desc(max_clone_ratio)) %>%
  mutate(
    rank = row_number(),
    clone_label = paste0(Sample, " #", rank),
    Group = factor(as.character(Group), levels = group_levels)
  )

write_csv_safe(expanded_clone_rank, file.path(result_dir, "TCR_expanded_clone_cell_number_rank_source.csv"))

expanded_clone_rank_plot <- expanded_clone_rank %>%
  slice_head(n = top_expanded_clones_n) %>%
  mutate(clone_label = factor(clone_label, levels = rev(clone_label)))

p_clone_cells <- ggplot(expanded_clone_rank_plot, aes(expanded_cells, clone_label, fill = Group)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.24, alpha = 0.88) +
  scale_fill_manual(values = group_colors, drop = FALSE, name = "Group") +
  labs(x = "Expanded T cells per clone", y = NULL) +
  theme_pub(8.6) +
  theme(
    axis.text.y = element_text(size = 6.4),
    axis.text.x = element_text(size = 7.6),
    legend.position = "top"
  )
save_plot_safe(p_clone_cells, "TCR_expanded_clone_cell_number_rank_barplot_revision", width = 5.4, height = 6.4)

hist_df <- cell_df %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    clone_ratio_capped = pmin(clone_ratio, histogram_xmax),
    Group = factor(as.character(Group), levels = group_levels),
    clonotype_id = dplyr::coalesce(as.character(tcell@meta.data[cell, "CTaa"]), as.character(tcell@meta.data[cell, "CTgene"]), cell),
    clonotype_id = ifelse(is.na(clonotype_id) | clonotype_id == "", cell, clonotype_id),
    clone_id = paste(Sample, clonotype_id, sep = "__")
  )

hist_group_label <- hist_df %>%
  group_by(Group) %>%
  summarise(
    n_t_cells = n(),
    n_expanded_cells = sum(expanded_clone, na.rm = TRUE),
    n_expanded_clones = n_distinct(clone_id[expanded_clone]),
    expanded_cell_percent = 100 * n_expanded_cells / n_t_cells,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "TCR clones=", n_expanded_clones,
      "\nT cells=", n_expanded_cells, "/", n_t_cells,
      " (", sprintf("%.1f", expanded_cell_percent), "%)"
    )
  )

p_hist <- ggplot(hist_df, aes(x = clone_ratio_capped, fill = Group)) +
  geom_histogram(binwidth = 0.25, boundary = 0, colour = "white", linewidth = 0.08, alpha = 0.86) +
  geom_vline(xintercept = expanded_threshold, linetype = "dashed", linewidth = 0.38, colour = "black") +
  geom_text(
    data = hist_group_label,
    aes(x = histogram_xmax * 0.98, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1.16,
    family = "Arial",
    fontface = "bold",
    size = 2.35,
    lineheight = 0.90
  ) +
  facet_wrap(~Group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors, drop = FALSE, name = "Group") +
  scale_x_continuous(limits = c(0, histogram_xmax), breaks = c(0, 1, 2, 5, 10), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "TCR clone_ratio", y = "T cells") +
  theme_pub(8.6) +
  theme(
    strip.text = element_text(size = 8.6, face = "bold"),
    axis.text.x = element_text(size = 7.0),
    axis.text.y = element_text(size = 6.7),
    panel.spacing.y = unit(0.28, "lines"),
    legend.position = "none"
  )
save_plot_safe(p_hist, "TCR_clone_ratio_sample_histogram_revision", width = 3.8, height = 5.8)

p_hist_sample <- ggplot(hist_df, aes(x = clone_ratio_capped, fill = Group)) +
  geom_histogram(binwidth = 0.25, boundary = 0, colour = "white", linewidth = 0.10, alpha = 0.88) +
  geom_vline(xintercept = expanded_threshold, linetype = "dashed", linewidth = 0.32, colour = "black") +
  annotate(
    "text",
    x = expanded_threshold + 0.25,
    y = Inf,
    label = "threshold = 1\nexpanded clone",
    hjust = 0,
    vjust = 1.25,
    family = "Arial",
    fontface = "bold",
    size = 2.0
  ) +
  facet_wrap(~Sample, ncol = 6, scales = "free_y") +
  scale_fill_manual(values = group_colors, drop = FALSE, name = "Group") +
  scale_x_continuous(limits = c(0, histogram_xmax), breaks = c(0, 1, 2, 5, 10), expand = expansion(mult = c(0, 0.02))) +
  labs(x = paste0("TCR clone_ratio (values > ", histogram_xmax, " shown at ", histogram_xmax, ")"), y = "T cells") +
  theme_pub(8.0) +
  theme(
    strip.text = element_text(size = 6.5, face = "bold"),
    axis.text.x = element_text(size = 6.2),
    axis.text.y = element_text(size = 5.8),
    legend.position = "top",
    legend.text = element_text(size = 7.2)
  )
save_plot_safe(p_hist_sample, "TCR_clone_ratio_sample_histogram_by_group_color_revision", width = 10.8, height = 7.2)

group_expanded_count <- hist_df %>%
  group_by(Group) %>%
  summarise(
    total_t_cells = n(),
    expanded_t_cells = sum(expanded_clone, na.rm = TRUE),
    expanded_tcr_clones = n_distinct(clone_id[expanded_clone]),
    expanded_t_cell_percent = 100 * expanded_t_cells / total_t_cells,
    .groups = "drop"
  ) %>%
  mutate(Group = factor(as.character(Group), levels = group_levels))

write_csv_safe(group_expanded_count, file.path(result_dir, "TCR_expanded_clone_and_cell_count_by_group.csv"))

group_expanded_count_plot <- group_expanded_count %>%
  select(Group, total_t_cells, expanded_t_cells, expanded_tcr_clones, expanded_t_cell_percent) %>%
  tidyr::pivot_longer(
    cols = c(expanded_tcr_clones, expanded_t_cells),
    names_to = "metric",
    values_to = "count"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("expanded_tcr_clones", "expanded_t_cells"),
      labels = c("Expanded TCR clones", "Expanded T cells")
    ),
    label = ifelse(
      metric == "Expanded T cells",
      paste0(count, " (", sprintf("%.1f", expanded_t_cell_percent), "%)"),
      as.character(count)
    )
  )

p_group <- ggplot(group_expanded_count_plot, aes(Group, count, fill = Group)) +
  geom_col(width = 0.54, colour = "black", linewidth = 0.28, alpha = 0.86) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    family = "Arial",
    fontface = "bold",
    size = 2.55
  ) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors, drop = FALSE, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Count") +
  theme_pub(9.2) +
  theme(
    strip.text = element_text(size = 8.8, face = "bold"),
    axis.text.x = element_text(size = 8.5, face = "bold"),
    axis.text.y = element_text(size = 8.0),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.38, "lines")
  )
save_plot_safe(p_group, "TCR_expanded_clone_fraction_by_group_barplot_wilcox_revision", width = 4.2, height = 4.2)

message("Saved TCR clone_ratio UMAP, sample histograms, and group expanded clone fraction barplot.")
message("Formula: ", file.path(result_dir, "TCR_clone_ratio_formula_and_threshold.txt"))
message("Tables: ", result_dir)
message("Figures: ", figure_dir)

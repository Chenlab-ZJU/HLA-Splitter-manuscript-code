#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Final urine T-cell, CXCR3 and expanded-clone figures.
# Original source: 36_revision_expanded_TCR_clone_diversity_by_group.R
# Figure/table scope: Fig5a; Fig5b; Fig5c; Fig5d; Fig5e; Fig5f; Fig5g; Fig5h; Fig5i; FigS9a; FigS9b; FigS9c; FigS9d; FigS9e; FigS9f; FigS9g
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Expanded T cells: clone_ratio > 1; distinguish cell counts from clonotype counts.
# Keep within-sample TCR denominators and clonotype identifiers as defined below.


# Expanded TCR clone diversity comparison by clinical group.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

project_dir_env <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(project_dir_env)) {
  stopifnot(dir.exists(project_dir_env))
  project_dir <- gsub("\\\\", "/", project_dir_env)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "36_revision_expanded_TCR_clone_diversity_by_group.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

object_path <- file.path(project_dir, "Result", "Objects", "Urine_Tcell.rds")
figure_dir <- file.path(project_dir, "Figures", "TCR", "Revision_TCR_clone_diversity")
result_dir <- file.path(project_dir, "Result", "TCR", "Revision_TCR_clone_diversity")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#F3C969",
  "Mix" = "#5DB1DD",
  "BK" = "#E8A0C8"
)
expanded_threshold <- 1

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

save_plot_safe <- function(plot, name, width, height) {
  pdf_path <- file.path(figure_dir, paste0(name, ".pdf"))
  png_path <- file.path(figure_dir, paste0(name, ".png"))
  temp_png <- file.path(tempdir(), paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf)
  ggsave(temp_png, plot, width = width, height = height, dpi = 600, device = "png", bg = "white")
  file.copy(temp_png, png_path, overwrite = TRUE)
  unlink(temp_png)
}

theme_pub <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black", face = "bold"),
      strip.background = element_rect(fill = "grey94", colour = "grey82", linewidth = 0.35),
      strip.text = element_text(colour = "black", face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.key.size = unit(0.42, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

shannon_index <- function(counts) {
  counts <- counts[counts > 0]
  if (length(counts) == 0) return(0)
  p <- counts / sum(counts)
  -sum(p * log(p))
}

simpson_diversity <- function(counts) {
  counts <- counts[counts > 0]
  if (length(counts) == 0) return(0)
  p <- counts / sum(counts)
  1 - sum(p^2)
}

inverse_simpson <- function(counts) {
  counts <- counts[counts > 0]
  if (length(counts) == 0) return(0)
  p <- counts / sum(counts)
  1 / sum(p^2)
}

pielou_evenness <- function(counts) {
  counts <- counts[counts > 0]
  if (length(counts) <= 1) return(0)
  shannon_index(counts) / log(length(counts))
}

gini_index <- function(counts) {
  counts <- sort(counts[counts > 0])
  n <- length(counts)
  if (n <= 1 || sum(counts) == 0) return(0)
  (2 * sum(seq_len(n) * counts) / (n * sum(counts))) - ((n + 1) / n)
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "p<0.001", paste0("p=", sprintf("%.2g", p))))
}

safe_wilcox <- function(df, g1, g2) {
  dat <- df %>% filter(Group %in% c(g1, g2), is.finite(value))
  if (n_distinct(dat$Group) < 2 || any(table(dat$Group) < 2)) return(NA_real_)
  suppressWarnings(wilcox.test(value ~ Group, data = dat, exact = FALSE)$p.value)
}

safe_kruskal <- function(df) {
  dat <- df %>% filter(is.finite(value))
  if (n_distinct(dat$Group) < 2 || any(table(dat$Group) < 1)) return(NA_real_)
  suppressWarnings(kruskal.test(value ~ Group, data = dat)$p.value)
}

message("Loading urine T-cell object...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(object_path)
md <- tcell@meta.data %>% rownames_to_column("cell")

sample_col <- if ("HLA_Splitter" %in% colnames(md)) "HLA_Splitter" else "Sample"
clone_col <- intersect(c("CTstrict", "CTaa", "CTgene", "CTnt"), colnames(md))[1]
if (is.na(clone_col)) stop("No TCR clonotype column found: CTstrict, CTaa, CTgene, or CTnt.")
if (!"clone_ratio" %in% colnames(md)) stop("clone_ratio is missing from Urine_Tcell.rds.")

cell_df <- md %>%
  transmute(
    cell,
    Sample = as.character(.data[[sample_col]]),
    Group = factor(as.character(Group), levels = group_levels),
    clonotype = as.character(.data[[clone_col]]),
    clone_ratio = suppressWarnings(as.numeric(clone_ratio)),
    clonalFrequency = suppressWarnings(as.numeric(clonalFrequency))
  ) %>%
  filter(!is.na(Group), !is.na(Sample), Sample != "") %>%
  mutate(
    clonotype = ifelse(is.na(clonotype) | clonotype == "", cell, clonotype),
    clone_id = paste(Sample, clonotype, sep = "__"),
    clone_ratio = ifelse(is.na(clone_ratio), 0, clone_ratio),
    expanded_clone = clone_ratio > expanded_threshold
  )

sample_totals <- cell_df %>%
  group_by(Group, Sample) %>%
  summarise(
    total_t_cells = n(),
    tcr_detected_cells = sum(clonalFrequency > 0 | clone_ratio > 0, na.rm = TRUE),
    expanded_t_cells = sum(expanded_clone, na.rm = TRUE),
    .groups = "drop"
  )

expanded_clone_counts <- cell_df %>%
  filter(expanded_clone) %>%
  count(Group, Sample, clone_id, name = "clone_cells")

sample_diversity <- sample_totals %>%
  group_by(Group, Sample) %>%
  group_modify(function(.x, .y) {
    counts <- expanded_clone_counts %>%
      filter(Group == .y$Group, Sample == .y$Sample) %>%
      pull(clone_cells)
    n_clone <- length(counts)
    tibble(
      expanded_tcr_clones = n_clone,
      shannon = shannon_index(counts),
      inverse_simpson = inverse_simpson(counts),
      gini_simpson = simpson_diversity(counts),
      pielou_evenness = pielou_evenness(counts),
      clonality = 1 - pielou_evenness(counts),
      gini_index = gini_index(counts),
      dominant_clone_fraction = ifelse(sum(counts) > 0, max(counts) / sum(counts), 0)
    )
  }) %>%
  ungroup() %>%
  left_join(sample_totals, by = c("Group", "Sample")) %>%
  mutate(
    expanded_t_cell_fraction = expanded_t_cells / total_t_cells,
    expanded_t_cell_percent = expanded_t_cell_fraction * 100
  )

metric_labels <- c(
  expanded_tcr_clones = "Expanded TCR clones",
  expanded_t_cells = "Expanded T cells",
  expanded_t_cell_percent = "Expanded T cells (%)",
  shannon = "Shannon diversity",
  inverse_simpson = "Inverse Simpson",
  gini_simpson = "Gini-Simpson",
  pielou_evenness = "Pielou evenness",
  clonality = "Clonality",
  dominant_clone_fraction = "Dominant clone fraction",
  gini_index = "Gini index"
)

diversity_long <- sample_diversity %>%
  select(Group, Sample, all_of(names(metric_labels))) %>%
  pivot_longer(
    cols = all_of(names(metric_labels)),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels),
    Group = factor(as.character(Group), levels = group_levels)
  )

group_summary <- diversity_long %>%
  group_by(metric, Group) %>%
  summarise(
    n_samples = n_distinct(Sample),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    sem = sd / sqrt(n_samples),
    .groups = "drop"
  )

kw_stats <- diversity_long %>%
  group_by(metric) %>%
  summarise(
    p.value = safe_kruskal(pick(Group, value)),
    p.adjust = p.adjust(p.value, method = "BH"),
    label = format_p(p.value),
    .groups = "drop"
  )

pair_grid <- tidyr::crossing(group1 = group_levels, group2 = group_levels) %>%
  filter(as.numeric(factor(group1, levels = group_levels)) < as.numeric(factor(group2, levels = group_levels)))

pairwise_stats <- diversity_long %>%
  group_by(metric) %>%
  group_modify(function(.x, .y) {
    out <- pair_grid
    out$p.value <- mapply(function(g1, g2) safe_wilcox(.x, g1, g2), out$group1, out$group2)
    out$mean_group1 <- mapply(function(g1) mean(.x$value[.x$Group == g1], na.rm = TRUE), out$group1)
    out$mean_group2 <- mapply(function(g2) mean(.x$value[.x$Group == g2], na.rm = TRUE), out$group2)
    out
  }) %>%
  ungroup() %>%
  group_by(metric) %>%
  mutate(
    p.adjust = p.adjust(p.value, method = "BH"),
    label = format_p(p.value)
  ) %>%
  ungroup()

# Export the analysis source data for figure reproduction.
write_csv_safe(sample_diversity, file.path(result_dir, "Expanded_TCR_clone_diversity_by_sample.csv"))
write_csv_safe(diversity_long, file.path(result_dir, "Expanded_TCR_clone_diversity_long.csv"))
write_csv_safe(group_summary, file.path(result_dir, "Expanded_TCR_clone_diversity_group_summary.csv"))
write_csv_safe(kw_stats, file.path(result_dir, "Expanded_TCR_clone_diversity_kruskal_tests.csv"))
write_csv_safe(pairwise_stats, file.path(result_dir, "Expanded_TCR_clone_diversity_pairwise_wilcox_tests.csv"))

plot_metrics <- c(
  "Expanded TCR clones",
  "Expanded T cells (%)",
  "Shannon diversity",
  "Inverse Simpson",
  "Clonality",
  "Dominant clone fraction"
)

plot_df <- diversity_long %>%
  filter(metric %in% plot_metrics) %>%
  mutate(metric = factor(as.character(metric), levels = plot_metrics))

set.seed(20260827)
p_diversity <- ggplot(plot_df, aes(Group, value, fill = Group)) +
  geom_boxplot(width = 0.54, outlier.shape = NA, alpha = 0.28, linewidth = 0.42, colour = "black") +
  geom_point(
    position = position_jitter(width = 0.14, height = 0),
    shape = 21,
    size = 1.75,
    stroke = 0.24,
    colour = "black",
    alpha = 0.88
  ) +
  facet_wrap(~metric, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
  labs(x = NULL, y = "Metric value") +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(size = 8.0, face = "bold", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 7.8),
    strip.text = element_text(size = 8.4, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.55, "lines")
  )

save_plot_safe(p_diversity, "Expanded_TCR_clone_diversity_by_group_boxplot_revision", width = 7.2, height = 4.8)

calc_pooled_metrics <- function(dat) {
  totals <- dat %>%
    group_by(Group) %>%
    summarise(
      total_t_cells = n(),
      expanded_t_cells = sum(expanded_clone, na.rm = TRUE),
      .groups = "drop"
    )

  clone_counts <- dat %>%
    filter(expanded_clone) %>%
    count(Group, clone_id, name = "clone_cells")

  totals %>%
    group_by(Group) %>%
    group_modify(function(.x, .y) {
      counts <- clone_counts %>%
        filter(Group == .y$Group) %>%
        pull(clone_cells)
      tibble(
        expanded_tcr_clones = length(counts),
        expanded_t_cells = .x$expanded_t_cells[1],
        expanded_t_cell_percent = 100 * .x$expanded_t_cells[1] / .x$total_t_cells[1],
        shannon = shannon_index(counts),
        inverse_simpson = inverse_simpson(counts),
        gini_simpson = simpson_diversity(counts),
        pielou_evenness = pielou_evenness(counts),
        clonality = 1 - pielou_evenness(counts),
        dominant_clone_fraction = ifelse(sum(counts) > 0, max(counts) / sum(counts), 0),
        gini_index = gini_index(counts)
      )
    }) %>%
    ungroup()
}

stat_spread <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) < 2) return(NA_real_)
  sum((values - mean(values))^2)
}

pooled_long <- calc_pooled_metrics(cell_df) %>%
  pivot_longer(
    cols = all_of(names(metric_labels)),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels),
    Group = factor(as.character(Group), levels = group_levels)
  )

set.seed(20260827)
n_perm <- 199
observed_spread <- pooled_long %>%
  filter(metric %in% plot_metrics) %>%
  group_by(metric) %>%
  summarise(observed_spread = stat_spread(value), .groups = "drop")

perm_spread <- lapply(seq_len(n_perm), function(iter_id) {
  perm_df <- cell_df
  perm_df$Group <- sample(perm_df$Group)
  calc_pooled_metrics(perm_df) %>%
    pivot_longer(cols = all_of(names(metric_labels)), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = names(metric_labels), labels = metric_labels)) %>%
    filter(metric %in% plot_metrics) %>%
    group_by(metric) %>%
    summarise(perm_spread = stat_spread(value), .groups = "drop") %>%
    mutate(iteration = iter_id)
}) %>%
  bind_rows()

pooled_perm_stats <- lapply(seq_len(nrow(observed_spread)), function(i) {
  one_metric <- observed_spread$metric[i]
  obs_stat <- observed_spread$observed_spread[i]
  perm_values <- perm_spread %>%
    filter(metric == one_metric) %>%
    pull(perm_spread)
  p_val <- (sum(perm_values >= obs_stat, na.rm = TRUE) + 1) / (sum(is.finite(perm_values)) + 1)
  tibble(
    metric = one_metric,
    observed_spread = obs_stat,
    n_permutations = sum(is.finite(perm_values)),
    p.value = p_val,
    label = paste0("perm ", format_p(p_val))
  )
}) %>%
  bind_rows() %>%
  mutate(metric = factor(as.character(metric), levels = plot_metrics))

write_csv_safe(pooled_long, file.path(result_dir, "Expanded_TCR_clone_diversity_pooled_group_metrics.csv"))
write_csv_safe(pooled_perm_stats, file.path(result_dir, "Expanded_TCR_clone_diversity_pooled_group_permutation_tests.csv"))

rank_df <- pooled_long %>%
  filter(metric %in% plot_metrics) %>%
  mutate(metric = factor(as.character(metric), levels = plot_metrics))

p_summary <- ggplot(rank_df, aes(Group, value, fill = Group)) +
  geom_col(width = 0.54, colour = "black", linewidth = 0.25, alpha = 0.86) +
  geom_text(
    data = pooled_perm_stats,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.04,
    vjust = 1.20,
    family = "Arial",
    fontface = "bold",
    size = 2.25
  ) +
  facet_wrap(~metric, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  labs(x = NULL, y = "Pooled metric value") +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(size = 8.0, face = "bold", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 7.8),
    strip.text = element_text(size = 8.4, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.55, "lines")
  )

save_plot_safe(p_summary, "Expanded_TCR_clone_diversity_by_group_mean_sem_revision", width = 7.2, height = 4.8)

message("Saved expanded TCR clone diversity comparisons.")
message("Tables: ", result_dir)
message("Figures: ", figure_dir)

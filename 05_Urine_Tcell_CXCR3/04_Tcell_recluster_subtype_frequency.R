#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis.
# Original source: 19_revision_Tcell_recluster_subtype_frequency.R
# Figure/table scope: Fig5; FigS9; Supplementary Table 18; Methods CXCR3-axis analysis
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# Reclustered urine T-cell subtype frequency comparison between TCMR and ABMR.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
urine_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
result_dir <- file.path(urine_dir, "Result", "Reviewer_Tcell_ABMR_TCMR")
figure_dir <- file.path(urine_dir, "Figures", "Reviewer_Tcell_ABMR_TCMR", "Tcell_recluster_subtype_frequency")
table_dir <- file.path(result_dir, "Tcell_recluster_subtype_frequency")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

group_colors <- c("TCMR" = "#699ECA", "ABMR" = "#4DAF4A")
direction_colors <- c("ABMR enriched" = "#4DAF4A", "TCMR enriched" = "#699ECA")
min_tcells_for_test <- 10

theme_nc <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(colour = "black", size = base_size),
      axis.title = element_text(colour = "black", size = base_size + 1),
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.40, colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1.5),
      plot.subtitle = element_text(hjust = 0, size = base_size, colour = "grey30"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.35, "cm"),
      strip.background = element_rect(fill = "grey96", colour = "black", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

save_pub_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 600, bg = "white")
}

safe_wilcox <- function(x, g) {
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- droplevels(factor(g[keep], levels = c("TCMR", "ABMR")))
  if (length(unique(g)) != 2 || any(table(g) < 2)) return(NA_real_)
  wilcox.test(x ~ g, exact = FALSE)$p.value
}

message("Loading reclustered T-cell object...")
tcell_file <- file.path(result_dir, "Tcell_only_ABMR_TCMR_reclustered.rds")
if (!file.exists(tcell_file)) {
  stop("Cannot find reclustered T-cell object: ", tcell_file)
}

# Read the saved object or source table used by this analysis.
tcell <- readRDS(tcell_file)
metadata <- tcell@meta.data %>%
  tibble::rownames_to_column("cell_barcode")

required_cols <- c("Group", "Sample", "Reviewer_Celltype")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

metadata <- metadata %>%
  filter(Group %in% c("TCMR", "ABMR")) %>%
  mutate(
    Group = factor(as.character(Group), levels = c("TCMR", "ABMR")),
    Sample = as.character(Sample),
    Tcell_subtype = as.character(Reviewer_Celltype)
  ) %>%
  filter(!is.na(Sample), !is.na(Tcell_subtype), Sample != "", Tcell_subtype != "")

subtype_levels <- metadata %>%
  count(Tcell_subtype, sort = TRUE) %>%
  pull(Tcell_subtype)

sample_info <- metadata %>%
  distinct(Sample, Group) %>%
  arrange(Group, Sample)

sample_totals <- metadata %>%
  count(Group, Sample, name = "sample_tcell_count")

sample_subtype_frequency <- metadata %>%
  count(Group, Sample, Tcell_subtype, name = "cell_count") %>%
  complete(sample_info, Tcell_subtype = subtype_levels, fill = list(cell_count = 0)) %>%
  left_join(sample_totals, by = c("Group", "Sample")) %>%
  mutate(
    Group = factor(as.character(Group), levels = c("TCMR", "ABMR")),
    Tcell_subtype = factor(Tcell_subtype, levels = subtype_levels),
    frequency = cell_count / sample_tcell_count,
    frequency_percent = 100 * frequency,
    enough_cells_for_test = sample_tcell_count >= min_tcells_for_test
  ) %>%
  arrange(Group, Sample, Tcell_subtype)

group_subtype_summary <- sample_subtype_frequency %>%
  group_by(Group, Tcell_subtype) %>%
  summarise(
    n_samples = n_distinct(Sample),
    n_samples_for_test = n_distinct(Sample[enough_cells_for_test]),
    total_cells = sum(cell_count),
    mean_frequency = mean(frequency, na.rm = TRUE),
    median_frequency = median(frequency, na.rm = TRUE),
    sd_frequency = sd(frequency, na.rm = TRUE),
    sem_frequency = sd_frequency / sqrt(n_samples),
    .groups = "drop"
  )

frequency_stats <- sample_subtype_frequency %>%
  filter(enough_cells_for_test) %>%
  group_by(Tcell_subtype) %>%
  summarise(
    n_TCMR_samples = n_distinct(Sample[Group == "TCMR"]),
    n_ABMR_samples = n_distinct(Sample[Group == "ABMR"]),
    mean_TCMR = mean(frequency[Group == "TCMR"], na.rm = TRUE),
    mean_ABMR = mean(frequency[Group == "ABMR"], na.rm = TRUE),
    median_TCMR = median(frequency[Group == "TCMR"], na.rm = TRUE),
    median_ABMR = median(frequency[Group == "ABMR"], na.rm = TRUE),
    diff_ABMR_minus_TCMR = mean_ABMR - mean_TCMR,
    log2FC_ABMR_vs_TCMR = log2((mean_ABMR + 1e-4) / (mean_TCMR + 1e-4)),
    p_value = safe_wilcox(frequency, Group),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    direction = if_else(diff_ABMR_minus_TCMR >= 0, "ABMR enriched", "TCMR enriched"),
    abs_diff = abs(diff_ABMR_minus_TCMR),
    rank_score = abs_diff * (-log10(pmax(p_adj, .Machine$double.xmin)))
  ) %>%
  arrange(desc(abs_diff), p_adj)

top_subtype <- frequency_stats %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)

message("Most differential T-cell subtype by sample-level frequency:")
print(top_subtype)

# Export the analysis source data for figure reproduction.
write_csv(sample_subtype_frequency, file.path(table_dir, "Tcell_recluster_sample_subtype_frequency.csv"))
write_csv(group_subtype_summary, file.path(table_dir, "Tcell_recluster_group_subtype_frequency_summary.csv"))
write_csv(frequency_stats, file.path(table_dir, "Tcell_recluster_TCMR_ABMR_subtype_frequency_tests.csv"))
write_csv(top_subtype, file.path(table_dir, "Tcell_recluster_top_differential_subtype.csv"))

plot_data <- sample_subtype_frequency %>%
  mutate(
    Tcell_subtype = factor(Tcell_subtype, levels = frequency_stats$Tcell_subtype[order(frequency_stats$diff_ABMR_minus_TCMR)])
  )

p_frequency <- ggplot(plot_data, aes(x = Group, y = frequency, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.45, alpha = 0.82) +
  geom_point(
    aes(size = sample_tcell_count, alpha = enough_cells_for_test),
    position = position_jitter(width = 0.10, height = 0),
    shape = 21,
    stroke = 0.25,
    colour = "black"
  ) +
  facet_wrap(~ Tcell_subtype, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_alpha_manual(values = c("TRUE" = 0.88, "FALSE" = 0.35), guide = "none") +
  scale_size_continuous(range = c(1.8, 4.2), name = "T cells/sample") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.03, 0.18))) +
  labs(
    title = "Reclustered urine T-cell subtype frequencies",
    subtitle = "Each point represents one clinical sample; frequency is calculated within T cells of that sample",
    x = NULL,
    y = "Subtype frequency within T cells"
  ) +
  theme_nc(base_size = 9.5) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(face = "bold")
  )

diff_plot_data <- frequency_stats %>%
  mutate(
    Tcell_subtype = factor(Tcell_subtype, levels = Tcell_subtype[order(diff_ABMR_minus_TCMR)]),
    label_text = paste0(
      ifelse(diff_ABMR_minus_TCMR >= 0, "+", ""),
      percent(diff_ABMR_minus_TCMR, accuracy = 0.1)
    )
  )

p_difference <- ggplot(diff_plot_data, aes(x = diff_ABMR_minus_TCMR, y = Tcell_subtype, colour = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.45, linetype = "dashed", colour = "grey45") +
  geom_segment(aes(x = 0, xend = diff_ABMR_minus_TCMR, yend = Tcell_subtype), linewidth = 1.05, alpha = 0.72) +
  geom_point(size = 3.2) +
  geom_text(aes(label = label_text), hjust = ifelse(diff_plot_data$diff_ABMR_minus_TCMR >= 0, -0.12, 1.12),
            size = 3.0, colour = "black") +
  scale_colour_manual(values = direction_colors) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.18, 0.18))) +
  labs(
    title = "Subtype enrichment in ABMR versus TCMR",
    subtitle = paste0(
      "Top differential subtype: ",
      as.character(top_subtype$Tcell_subtype),
      " (",
      top_subtype$direction,
      ")"
    ),
    x = "Mean frequency difference (ABMR - TCMR)",
    y = NULL,
    colour = NULL
  ) +
  theme_nc(base_size = 9.5) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0)
  )

heatmap_data <- frequency_stats %>%
  select(Tcell_subtype, mean_TCMR, mean_ABMR) %>%
  pivot_longer(cols = c(mean_TCMR, mean_ABMR), names_to = "Group", values_to = "mean_frequency") %>%
  mutate(
    Group = recode(Group, mean_TCMR = "TCMR", mean_ABMR = "ABMR"),
    Group = factor(Group, levels = c("TCMR", "ABMR")),
    Tcell_subtype = factor(Tcell_subtype, levels = diff_plot_data$Tcell_subtype)
  )

p_heatmap <- ggplot(heatmap_data, aes(x = Group, y = Tcell_subtype, fill = mean_frequency)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = percent(mean_frequency, accuracy = 0.1)), size = 2.9, fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "Mean sample-level frequency",
    x = NULL,
    y = NULL,
    fill = "Mean\nfrequency"
  ) +
  theme_nc(base_size = 9.5) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 0),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )

p_ranked_summary <- p_difference + p_heatmap + plot_layout(widths = c(1.45, 0.82))
p_combined <- p_frequency / p_ranked_summary + plot_layout(heights = c(1.35, 1))

save_pub_plot(p_frequency, "Tcell_recluster_subtype_frequency_boxplot", width = 8.4, height = 7.2)
save_pub_plot(p_ranked_summary, "Tcell_recluster_subtype_enrichment_summary", width = 8.6, height = 4.4)
save_pub_plot(p_combined, "Tcell_recluster_subtype_frequency_TCMR_ABMR_combined", width = 8.8, height = 11.0)

message("Done.")
message("Tables: ", table_dir)
message("Figures: ", figure_dir)

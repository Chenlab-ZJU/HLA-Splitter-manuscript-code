#!/usr/bin/env Rscript
# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final urine immune-cell composition, module and demultiplexing figures.
# Original source: 41_revision_scRNA_celltype_propration_group_compare.R
# Figure/table scope: Fig4b; Fig4c; Fig4d; Fig4e; Fig4f; Fig4g; Fig4h; FigS7a; FigS7b; FigS7c; FigS7d; FigS7e
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# Redraw the legacy urine scRNA celltype proportion plot using current
# annotations and sample-level group comparisons.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(gtools)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "41_revision_scRNA_celltype_propration_group_compare.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures", "scRNA")
figure_comp_dir <- file.path(project_dir, "Figures", "scRNA", "Celltype_Composition")
result_dir <- file.path(project_dir, "Result", "scRNA", "Celltype_Composition")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_comp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "BK")
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "BK" = "#F898CB"
)

major_levels <- c("CD4+ T", "CD8+ T", "Cycling T", "NK", "B cell", "Monocyte", "Macrophage", "Neutrophil", "cDC", "pDC")
major_colors <- c(
  "CD4+ T" = "#D86A63",
  "CD8+ T" = "#7EA9CC",
  "Cycling T" = "#77B96E",
  "NK" = "#E9A24C",
  "B cell" = "#E7A2C9",
  "Monocyte" = "#C86A9B",
  "Macrophage" = "#87608C",
  "Neutrophil" = "#E7C46A",
  "cDC" = "#D79D4F",
  "pDC" = "#5A9BC3"
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    p < 0.01 ~ formatC(p, format = "f", digits = 3),
    TRUE ~ formatC(p, format = "f", digits = 2)
  )
}

theme_pub <- function(base_size = 10.2) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.text = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey94", colour = NA),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(4, 8, 8, 10)
    )
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

save_plot_to_dir <- function(plot, dir, filename, width, height) {
  pdf_path <- file.path(dir, paste0(filename, ".pdf"))
  png_path <- file.path(dir, paste0(filename, ".png"))
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

object_path <- file.path(project_dir, "Result", "Objects", "Urine_ALL.rds")
if (!file.exists(object_path)) {
  object_path <- file.path(project_dir, "Result", "Objects", "Urine_Allcells.rds")
}
if (!file.exists(object_path)) {
  stop("Urine scRNA object not found in Result/Objects.")
}

message("Loading urine all-cell object...")
# Read the saved object or source table used by this analysis.
obj <- readRDS(object_path)

meta <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
    Group = as.character(Group),
    Major_celltype = as.character(Major_celltype)
  ) %>%
  filter(
    !is.na(Sample), !Sample %in% c("Undefined", "Doublets", "Doublet", ""),
    !is.na(Group), Group %in% group_levels,
    !is.na(Major_celltype), Major_celltype %in% major_levels
  ) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Major_celltype = factor(Major_celltype, levels = major_levels)
  )

sample_order <- meta %>%
  distinct(Sample, Group) %>%
  arrange(Group, mixedorder(Sample)) %>%
  pull(Sample)

prop_df <- meta %>%
  count(Sample, Group, Major_celltype, name = "n_cells") %>%
  complete(
    nesting(Sample, Group),
    Major_celltype = factor(major_levels, levels = major_levels),
    fill = list(n_cells = 0)
  ) %>%
  group_by(Sample, Group) %>%
  mutate(
    total_cells = sum(n_cells),
    Percent = ifelse(total_cells > 0, n_cells / total_cells * 100, 0)
  ) %>%
  ungroup() %>%
  mutate(
    Sample = factor(as.character(Sample), levels = sample_order),
    Group = factor(as.character(Group), levels = group_levels),
    Major_celltype = factor(as.character(Major_celltype), levels = major_levels)
  )

summary_df <- prop_df %>%
  group_by(Group, Major_celltype) %>%
  summarise(
    mean_percent = mean(Percent, na.rm = TRUE),
    sem_percent = sd(Percent, na.rm = TRUE) / sqrt(dplyr::n()),
    n_samples = n_distinct(Sample),
    .groups = "drop"
  )

total_stack_df <- meta %>%
  count(Group, Major_celltype, name = "n_cells") %>%
  complete(
    Group = factor(group_levels, levels = group_levels),
    Major_celltype = factor(major_levels, levels = major_levels),
    fill = list(n_cells = 0)
  ) %>%
  group_by(Group) %>%
  mutate(
    group_total_cells = sum(n_cells),
    Percent = ifelse(group_total_cells > 0, n_cells / group_total_cells * 100, 0)
  ) %>%
  ungroup() %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    Major_celltype = factor(as.character(Major_celltype), levels = major_levels)
  )

pairwise_grid <- as_tibble(t(utils::combn(group_levels, 2)), .name_repair = "minimal") %>%
  setNames(c("group1", "group2"))

wilcox_df <- prop_df %>%
  group_by(Major_celltype) %>%
  group_modify(function(dat, key) {
    pairwise_grid %>%
      mutate(
        p.value = mapply(function(g1, g2) {
          v1 <- dat$Percent[dat$Group == g1]
          v2 <- dat$Percent[dat$Group == g2]
          if (length(v1) < 2 || length(v2) < 2 || dplyr::n_distinct(c(v1, v2)) < 2) return(NA_real_)
          suppressWarnings(stats::wilcox.test(v1, v2, exact = FALSE)$p.value)
        }, group1, group2)
      )
  }) %>%
  ungroup() %>%
  group_by(Major_celltype) %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH"), label = format_p(p.value)) %>%
  ungroup()

stat_plot <- wilcox_df %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  mutate(
    Major_celltype = factor(as.character(Major_celltype), levels = major_levels),
    group1 = factor(group1, levels = group_levels),
    group2 = factor(group2, levels = group_levels),
    x = as.numeric(group1),
    xend = as.numeric(group2)
  ) %>%
  left_join(
    prop_df %>%
      group_by(Major_celltype) %>%
      summarise(
        y_min = min(Percent, na.rm = TRUE),
        y_max = max(Percent, na.rm = TRUE),
        y_range = pmax(y_max - y_min, 1.0),
        .groups = "drop"
      ),
    by = "Major_celltype"
  ) %>%
  group_by(Major_celltype) %>%
  arrange(p.value, .by_group = TRUE) %>%
  mutate(
    stat_rank = row_number(),
    y = y_max + y_range * (0.18 + 0.25 * (stat_rank - 1)),
    tick = y_range * 0.055,
    y_text = y + y_range * 0.10
  ) %>%
  ungroup()

# Export the analysis source data for figure reproduction.
write_csv_safe(prop_df, file.path(result_dir, "scRNA_celltype_propration_major_celltype_sample_percent_revision.csv"))
write_csv_safe(summary_df, file.path(result_dir, "scRNA_celltype_propration_major_celltype_group_summary_revision.csv"))
write_csv_safe(total_stack_df, file.path(result_dir, "scRNA_celltype_propration_major_celltype_group_total_stack_revision.csv"))
write_csv_safe(wilcox_df, file.path(result_dir, "scRNA_celltype_propration_major_celltype_pairwise_wilcox_revision.csv"))
write_csv_safe(stat_plot, file.path(result_dir, "scRNA_celltype_propration_major_celltype_pairwise_wilcox_p_lt_0.1_revision.csv"))

p_stack <- ggplot(total_stack_df, aes(Group, Percent, fill = Major_celltype)) +
  geom_col(width = 0.66, colour = "white", linewidth = 0.28) +
  scale_fill_manual(values = major_colors, drop = FALSE, name = "Major celltype") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0.015))
  ) +
  labs(x = NULL, y = "Cells (%)") +
  theme_bw(base_size = 10.2, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(face = "bold", colour = "black"),
    axis.title.y = element_text(size = 12.0, margin = margin(r = 6)),
    axis.text.x = element_text(size = 10.5, face = "bold", angle = 25, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10.0),
    axis.line = element_line(colour = "black", linewidth = 0.38),
    axis.ticks = element_line(colour = "black", linewidth = 0.32),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.4),
    legend.key.size = unit(0.34, "cm"),
    plot.margin = margin(5, 8, 5, 8)
  ) +
  guides(fill = guide_legend(override.aes = list(linewidth = 0), ncol = 1))

save_plot_to_dir(p_stack, figure_comp_dir, "scRNA_celltype_propration", width = 4.9, height = 4.6)

p <- ggplot(prop_df, aes(Group, Percent, fill = Group, colour = Group)) +
  stat_summary(
    fun = mean,
    geom = "col",
    width = 0.58,
    colour = "black",
    linewidth = 0.42,
    alpha = 0.42
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.22,
    linewidth = 0.42,
    colour = "black"
  ) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 11),
    shape = 21,
    stroke = 0.52,
    size = 2.10,
    alpha = 0.92
  ) +
  geom_segment(
    data = stat_plot,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.38,
    colour = "black"
  ) +
  geom_segment(
    data = stat_plot,
    aes(x = x, xend = x, y = y, yend = y - tick),
    inherit.aes = FALSE,
    linewidth = 0.38,
    colour = "black"
  ) +
  geom_segment(
    data = stat_plot,
    aes(x = xend, xend = xend, y = y, yend = y - tick),
    inherit.aes = FALSE,
    linewidth = 0.38,
    colour = "black"
  ) +
  geom_text(
    data = stat_plot,
    aes(x = (x + xend) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 3.8,
    colour = "black"
  ) +
  facet_wrap(~Major_celltype, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = group_colors, drop = FALSE, guide = "none") +
  scale_colour_manual(values = group_colors, drop = FALSE, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.34))) +
  labs(x = NULL, y = "Frequency (%)") +
  theme_pub(10.2) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9.2),
    axis.text.y = element_text(size = 9.8),
    axis.title.y = element_text(size = 11.8, face = "bold"),
    strip.text = element_text(size = 9.2, face = "bold"),
    panel.spacing.x = unit(0.75, "lines"),
    panel.spacing.y = unit(0.72, "lines")
  )

save_plot_safe(p, "scRNA_celltype_propration", width = 11.0, height = 4.8)

message("Saved revised scRNA_celltype_propration to: ", figure_dir)

#!/usr/bin/env Rscript
# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 08b_revision_selected_celltype_compare_width.R
# Figure/table scope: Fig2b; Fig2c; Fig2d; Fig2e; Fig2f; FigS4c; FigS4d; FigS4e; FigS5a; FigS5b; FigS5c; FigS5d; FigS5e; FigS5f
# Data cohort: csf. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.


# Replot selected CSF celltype comparison using existing source tables.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

env_project_dir <- Sys.getenv("CSF_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "08b_revision_selected_celltype_compare_width.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
group_colors <- c(
  "Metastases" = "#E42320",
  "Hemorrhage" = "#699ECA",
  "Encephalitis" = "#4DAF4A"
)
selected_celltypes <- c("CD4+ T", "B", "Macrophages")

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", signif(p, 2)))
}

theme_pub <- function(base_size = 9.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black", face = "bold"),
      strip.background = element_rect(fill = "grey94", colour = "grey82", linewidth = 0.35),
      strip.text = element_text(colour = "black", face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.key.size = unit(0.38, "cm"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_plot_safe <- function(plot, filename, width, height) {
  target_pdf <- file.path(figure_dir, paste0(filename, ".pdf"))
  target_png <- file.path(figure_dir, paste0(filename, ".png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  file.copy(tmp_pdf, target_pdf, overwrite = TRUE)
  file.copy(tmp_png, target_png, overwrite = TRUE)
  unlink(c(tmp_pdf, tmp_png))
}

# Read the saved object or source table used by this analysis.
sample_qc <- readRDS(file.path(getOption("hla.data_root"), "csf", "CSF_scRNA.rds"))@meta.data %>%
  rownames_to_column("cell") %>%
  count(HLA_Splitter, Group, name = "n_cells") %>%
  rename(Sample = HLA_Splitter) %>%
  arrange(desc(n_cells))
# Export the analysis source data for figure reproduction.
write.csv(sample_qc, file.path(result_dir, "CSF_sample_HLA_Splitter_cell_count_qc.csv"), row.names = FALSE)
valid_samples <- sample_qc$Sample[sample_qc$n_cells > 100]

selected_comp <- read.csv(
  file.path(getOption("hla.data_root"), "csf", "scRNA_celltype_sample_percent_revision.csv"),
  check.names = FALSE
) %>%
  filter(Sample %in% valid_samples) %>%
  filter(Celltype %in% selected_celltypes) %>%
  mutate(
    Celltype = factor(Celltype, levels = selected_celltypes),
    Group = factor(Group, levels = group_levels),
    Percent = as.numeric(Percent)
  )

selected_wilcox_plot <- selected_comp %>%
  group_by(Celltype) %>%
  group_modify(function(.x, .y) {
    comparisons <- t(utils::combn(group_levels, 2))
    data.frame(
      group1 = comparisons[, 1],
      group2 = comparisons[, 2],
      x = match(comparisons[, 1], group_levels),
      xend = match(comparisons[, 2], group_levels),
      p.value = apply(comparisons, 1, function(pair) {
        values_1 <- .x$Percent[.x$Group == pair[1]]
        values_2 <- .x$Percent[.x$Group == pair[2]]
        if (length(values_1) < 2 || length(values_2) < 2 || length(unique(c(values_1, values_2))) < 2) return(NA_real_)
        suppressWarnings(stats::wilcox.test(values_1, values_2, exact = FALSE)$p.value)
      })
    )
  }) %>%
  ungroup() %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  left_join(
    selected_comp %>%
      group_by(Celltype) %>%
      summarise(y_max = max(Percent, na.rm = TRUE), y_min = min(Percent, na.rm = TRUE), .groups = "drop"),
    by = "Celltype"
  ) %>%
  group_by(Celltype) %>%
  arrange(p.value, .by_group = TRUE) %>%
  mutate(
    Celltype = factor(Celltype, levels = selected_celltypes),
    y_range = pmax(y_max - y_min, 1),
    y.position = y_max + row_number() * y_range * 0.18,
    tick.height = y_range * 0.035,
    label.y = y.position + y_range * 0.065,
    label = format_p(p.value)
  ) %>%
  ungroup()

write.csv(
  selected_wilcox_plot,
  file.path(result_dir, "scRNA_selected_celltype_wilcox_test_revision.csv"),
  row.names = FALSE
)

p_selected_celltype_box <- ggplot(selected_comp, aes(x = Group, y = Percent, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.38, alpha = 0.82) +
  geom_point(
    position = position_jitter(width = 0.11, height = 0, seed = 20260813),
    size = 1.35,
    alpha = 0.90,
    shape = 21,
    colour = "black",
    stroke = 0.20
  ) +
  geom_segment(
    data = selected_wilcox_plot,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = selected_wilcox_plot,
    aes(x = x, xend = x, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = selected_wilcox_plot,
    aes(x = xend, xend = xend, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_text(
    data = selected_wilcox_plot,
    aes(x = (x + xend) / 2, y = label.y, label = label),
    inherit.aes = FALSE,
    size = 2.45,
    family = "Arial",
    fontface = "bold",
    colour = "black"
  ) +
  facet_wrap(~Celltype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_x_discrete(limits = group_levels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "Cells per sample (%)") +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7.8),
    strip.text = element_text(size = 8.2),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

save_plot_safe(p_selected_celltype_box, "scRNA_selected_celltype_compare_wilcox", width = 5.2, height = 2.8)

message("Updated scRNA_selected_celltype_compare_wilcox width.")

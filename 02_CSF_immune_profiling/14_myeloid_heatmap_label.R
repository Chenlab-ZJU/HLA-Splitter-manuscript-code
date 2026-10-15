#!/usr/bin/env Rscript
# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 10c_revision_myeloid_heatmap_label.R
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


# Replot Myeloid group DEG key-gene heatmap from existing source table.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

env_project_dir <- Sys.getenv("CSF_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "10c_revision_myeloid_heatmap_label.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")

theme_pub <- function(base_size = 9.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black", face = "bold"),
      plot.title = element_text(colour = "black", face = "bold"),
      strip.text = element_text(colour = "black", face = "bold"),
      legend.title = element_text(face = "bold"),
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

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")

# Read the saved object or source table used by this analysis.
heat_df <- read.csv(
  file.path(getOption("hla.data_root"), "csf", "Myeloid_revision_group_DEG_heatmap_source.csv"),
  check.names = FALSE
) %>%
  mutate(
    module = factor(module, levels = unique(module)),
    module_gene = factor(module_gene, levels = rev(unique(module_gene))),
    Group = factor(Group, levels = group_levels),
    z_score = as.numeric(z_score)
  )

p_heat <- ggplot(heat_df, aes(Group, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#6F99C3", mid = "#F7F7F7", high = "#B7192B", midpoint = 0, name = "Z-score") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  labs(title = "Myeloid cells", x = NULL, y = NULL) +
  theme_pub(9.2) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 7.6),
    strip.text.y = element_text(size = 7.6, angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey93", colour = NA),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_plot_safe(p_heat, "Myeloid_group_DEG_key_gene_heatmap", width = 4.7, height = 7.05)

message("Updated Myeloid_group_DEG_key_gene_heatmap with Myeloid cells label.")

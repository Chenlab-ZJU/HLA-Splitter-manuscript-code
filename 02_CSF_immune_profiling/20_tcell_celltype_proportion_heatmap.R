# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 13b_revision_tcell_celltype_proportion_heatmap.R
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

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(gtools)
})

project_dir <- Sys.getenv("CSF_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "csf")
}
figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
tcell_levels <- c(
  "Naive CD4 T", "CD4 memory T", "GZMK+ CD4 T", "ISG+ CD4 T", "Treg",
  "Naive CD8 T", "GZMK+ CD8 T", "GZMB+ cytotoxic T", "NK-like T", "gd T"
)

theme_pub <- function(base_size = 9.6) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.border = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
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
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
  file.copy(tmp_pdf, target_pdf, overwrite = TRUE)
  file.copy(tmp_png, target_png, overwrite = TRUE)
}

message("Loading CSF_Tcell.rds and manual T cell reannotation...")
manual_annotation_path <- file.path(getOption("hla.data_root"), "csf", "Tcell_revision_manual_cluster_reannotation.csv")
if (!file.exists(manual_annotation_path)) {
  stop("Missing Tcell_revision_manual_cluster_reannotation.csv. Please run 11b_revision_redraw_tcell_umap_from_csf_tcell.R first.")
}
# Read the saved object or source table used by this analysis.
manual_annotation <- read.csv(manual_annotation_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
  transmute(Cluster = as.character(Cluster), T_celltype = as.character(T_celltype)) %>%
  filter(T_celltype %in% tcell_levels)

Tcell <- readRDS(file.path(getOption("hla.data_root"), "csf", "CSF_Tcell.rds"))
meta <- Tcell@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  transmute(
    cell,
    Sample = as.character(HLA_Splitter),
    Group = factor(as.character(Group), levels = group_levels),
    Cluster = as.character(seurat_clusters)
  ) %>%
  left_join(manual_annotation, by = "Cluster") %>%
  filter(!is.na(Sample), !is.na(Group), !is.na(T_celltype)) %>%
  mutate(T_celltype = factor(T_celltype, levels = tcell_levels))

prop_df <- meta %>%
  count(Sample, Group, T_celltype, name = "Count") %>%
  group_by(Sample, Group) %>%
  mutate(Total = sum(Count), Percent = Count / Total * 100) %>%
  ungroup() %>%
  complete(nesting(Sample, Group), T_celltype = tcell_levels, fill = list(Count = 0, Percent = 0)) %>%
  group_by(Sample, Group) %>%
  mutate(Total = sum(Count)) %>%
  ungroup()

sample_order <- prop_df %>%
  distinct(Sample, Group) %>%
  arrange(Group, gtools::mixedorder(Sample)) %>%
  pull(Sample)

heat_df <- prop_df %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    Group = factor(as.character(Group), levels = group_levels),
    T_celltype = factor(as.character(T_celltype), levels = rev(tcell_levels))
  ) %>%
  arrange(Group, Sample, T_celltype)

# Export the analysis source data for figure reproduction.
write.csv(prop_df, file.path(result_dir, "Tcell_revision_subtype_sample_percent.csv"), row.names = FALSE)
write.csv(heat_df, file.path(result_dir, "Tcell_revision_subtype_proportion_heatmap_source.csv"), row.names = FALSE)

p_heat <- ggplot(heat_df, aes(Sample, T_celltype, fill = Percent)) +
  geom_tile(colour = "white", linewidth = 0.28) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#D7E4EF", "#F7F7F7", "#E79A94", "#B7192B"),
    name = "Cells (%)",
    limits = c(0, max(heat_df$Percent, na.rm = TRUE))
  ) +
  labs(x = "Samples grouped by diagnosis", y = "T cell subset") +
  theme_pub(10.0) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.6),
    axis.text.y = element_text(size = 9.2),
    axis.title = element_text(size = 10.5),
    strip.text = element_text(size = 9.8, face = "bold"),
    legend.position = "right",
    legend.key.height = unit(0.42, "cm")
  )

save_plot_safe(p_heat, "Tcell_celltype_proportion_heatmap_revision", width = 6.5, height = 4.8)
message("Redrew Tcell_celltype_proportion_heatmap_revision with current manual T cell annotations.")

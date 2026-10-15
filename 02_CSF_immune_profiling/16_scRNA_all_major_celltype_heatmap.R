# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 16_revision_scRNA_all_major_celltype_heatmap.R
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
project_dir <- gsub("\\\\", "/", project_dir)
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
major_levels <- c("T cells", "B cells", "Myeloid", "DC", "Cycling")

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

theme_pub <- function(base_size = 10.5) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.text = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = "grey55", linewidth = 0.35),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey55", linewidth = 0.35),
      legend.title = element_text(face = "bold", colour = "black"),
      legend.text = element_text(colour = "black"),
      plot.margin = margin(4, 4, 4, 4)
    )
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

message("Loading all-sample CSF major celltype composition table...")
major_comp_raw <- read.csv(
  file.path(getOption("hla.data_root"), "csf", "scRNA_major_celltype_sample_percent_revision.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  filter(Sample %in% valid_samples) %>%
  filter(!is.na(Sample), !is.na(Group), Group %in% group_levels, Major_celltype %in% major_levels) %>%
  mutate(
    Sample = as.character(Sample),
    Group = factor(as.character(Group), levels = group_levels),
    Major_celltype = as.character(Major_celltype),
    Percent = as.numeric(Percent)
  )

sample_group <- major_comp_raw %>%
  distinct(Sample, Group)

sample_order <- sample_group %>%
  arrange(Group, mixedorder(Sample)) %>%
  pull(Sample)

major_comp <- major_comp_raw %>%
  complete(
    Sample = sample_order,
    Major_celltype = major_levels,
    fill = list(Count = 0, Percent = 0)
  ) %>%
  left_join(sample_group, by = "Sample") %>%
  mutate(
    Group = coalesce(Group.x, Group.y),
    Sample = factor(Sample, levels = sample_order),
    Group = factor(Group, levels = group_levels),
    Major_celltype = factor(Major_celltype, levels = rev(major_levels))
  ) %>%
  select(Sample, Group, Major_celltype, Count, Total_count, Percent)

write.csv(
  major_comp,
  file.path(result_dir, "scRNA_all_major_celltype_heatmap_source_revision.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

p_heatmap <- ggplot(major_comp, aes(x = Sample, y = Major_celltype, fill = Percent)) +
  geom_tile(colour = "white", linewidth = 0.34) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#6F99C3", "#D7E4EF", "#F7F7F7", "#E79A94", "#B7192B"),
    name = "Cells (%)",
    limits = c(0, max(major_comp$Percent, na.rm = TRUE))
  ) +
  labs(x = "Samples grouped by diagnosis", y = "Major celltype") +
  theme_pub(11.0) +
  theme(
    axis.title = element_text(size = 12.6, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10.2),
    axis.text.y = element_text(size = 11.0),
    strip.text = element_text(size = 11.5, face = "bold"),
    legend.title = element_text(size = 11.2, face = "bold"),
    legend.text = element_text(size = 10.2),
    legend.position = "right",
    legend.key.height = unit(0.55, "in"),
    legend.key.width = unit(0.16, "in")
  )

save_plot_safe(p_heatmap, "scRNA_all_celltype_heatmap", width = 6.7, height = 3.1)

message("Saved scRNA_all_celltype_heatmap using all samples from the major celltype composition table.")

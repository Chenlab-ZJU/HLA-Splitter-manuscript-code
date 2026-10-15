#!/usr/bin/env Rscript
# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final urine immune-cell composition, module and demultiplexing figures.
# Original source: 42_revision_scRNA_UMAP_Rejection_Major_celltype.R
# Figure/table scope: Fig4b; Fig4c; Fig4d; Fig4e; Fig4f; Fig4g; Fig4h; FigS7a; FigS7b; FigS7c; FigS7d; FigS7e
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# Redraw urine all-cell UMAP split by rejection group and colored by
# Major_celltype. Existing UMAP coordinates are reused; no reclustering.

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
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "42_revision_scRNA_UMAP_Rejection_Major_celltype.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures", "scRNA", "UMAP_Markers")
result_dir <- file.path(project_dir, "Result", "scRNA", "UMAP_Markers")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "BK")
major_levels <- c("CD4+ T", "CD8+ T", "Cycling T", "NK", "B cell", "Monocyte", "Macrophage", "Neutrophil", "cDC", "pDC")
major_colors <- c(
  "CD4+ T" = "#E42320",
  "CD8+ T" = "#699ECA",
  "Cycling T" = "#4DAF4A",
  "NK" = "#FF8C00",
  "B cell" = "#F898CB",
  "Monocyte" = "#D65190",
  "Macrophage" = "#731A73",
  "Neutrophil" = "#FFCB5B",
  "cDC" = "#E89B1E",
  "pDC" = "#0076B9"
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

theme_umap <- function(base_size = 10.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      strip.background = element_rect(fill = "grey94", colour = NA),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.margin = margin(5, 6, 5, 6)
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

object_path <- file.path(getOption("hla.data_root"), "urine", "Urine_ALL.rds")
if (!file.exists(object_path)) {
  object_path <- file.path(project_dir, "Result", "Objects", "Urine_Allcells.rds")
}
if (!file.exists(object_path)) {
  stop("Urine scRNA object not found in Result/Objects.")
}

message("Loading urine all-cell object...")
# Read the saved object or source table used by this analysis.
obj <- readRDS(object_path)
if (!"umap" %in% names(obj@reductions)) {
  stop("The urine scRNA object has no UMAP reduction.")
}

meta <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Group = as.character(Group),
    Major_celltype = as.character(Major_celltype)
  ) %>%
  filter(
    !is.na(Group), Group %in% group_levels,
    !is.na(Major_celltype), Major_celltype %in% major_levels
  ) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Major_celltype = factor(Major_celltype, levels = major_levels)
  )

umap_df <- Embeddings(obj, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  left_join(meta, by = "cell") %>%
  filter(!is.na(Group), !is.na(Major_celltype))
colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")

group_labels <- umap_df %>%
  count(Group, name = "n_cells") %>%
  mutate(Group_label = paste0(as.character(Group), "\n", scales::comma(n_cells), " cells"))

umap_df <- umap_df %>%
  left_join(group_labels %>% select(Group, Group_label), by = "Group") %>%
  mutate(Group_label = factor(Group_label, levels = group_labels$Group_label))

background_df <- umap_df %>%
  select(cell, UMAP_1, UMAP_2) %>%
  distinct() %>%
  tidyr::crossing(Group_label = factor(group_labels$Group_label, levels = group_labels$Group_label))

set.seed(20260831)
umap_df <- umap_df %>% arrange(Group_label) %>% group_by(Group_label) %>% sample_frac(1) %>% ungroup()

# Export the analysis source data for figure reproduction.
write_csv_safe(umap_df, file.path(result_dir, "scRNA_UMAP_Rejection_Major_celltype_source_no_Mix_revision.csv"))

p_umap <- ggplot() +
  geom_point(
    data = background_df,
    aes(UMAP_1, UMAP_2),
    colour = "grey82",
    size = 0.08,
    alpha = 0.45,
    stroke = 0
  ) +
  geom_point(
    data = umap_df,
    aes(UMAP_1, UMAP_2, colour = Major_celltype),
    size = 0.18,
    alpha = 0.88,
    stroke = 0
  ) +
  facet_wrap(~Group_label, ncol = 2) +
  scale_colour_manual(values = major_colors, drop = FALSE, name = "Major celltype") +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 3.6, alpha = 1), ncol = 1)) +
  theme_umap(10.2) +
  theme(
    axis.title = element_text(size = 11.5, face = "bold"),
    axis.text = element_text(size = 9.0),
    strip.text = element_text(size = 10.6, face = "bold", lineheight = 0.90),
    panel.spacing = unit(0.85, "lines"),
    legend.position = "right",
    legend.title = element_text(size = 10.8, face = "bold"),
    legend.text = element_text(size = 9.8),
    legend.key.height = unit(0.35, "cm"),
    plot.margin = margin(5, 7, 5, 7)
  )

save_plot_safe(p_umap, "scRNA_UMAP_Rejection_Major_celltype", width = 8.8, height = 7.2)
message("Saved revised scRNA_UMAP_Rejection_Major_celltype to: ", figure_dir)

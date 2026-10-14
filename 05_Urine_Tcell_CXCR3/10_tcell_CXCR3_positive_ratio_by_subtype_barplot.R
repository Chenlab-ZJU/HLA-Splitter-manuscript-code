#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Final urine T-cell, CXCR3 and expanded-clone figures.
# Original source: 39_revision_tcell_CXCR3_positive_ratio_by_subtype_barplot.R
# Figure/table scope: Fig5a; Fig5b; Fig5c; Fig5d; Fig5e; Fig5f; Fig5g; Fig5h; Fig5i; FigS9a; FigS9b; FigS9c; FigS9d; FigS9e; FigS9f; FigS9g
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# CXCR3-positive fraction across urine T-cell subtypes.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "39_revision_tcell_CXCR3_positive_ratio_by_subtype_barplot.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

object_file <- file.path(project_dir, "Result", "Objects", "Urine_Tcell.rds")
figure_dir <- file.path(project_dir, "Figures", "Tcell", "Marker_Frequency")
result_dir <- file.path(project_dir, "Result", "Tcell", "Marker_Frequency")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

tcell_levels <- c(
  "Naive CD4", "CD4 Tem", "CD4 Tex", "Treg",
  "Naive CD8", "CXCR4+ CD8", "GZMB+ CD8", "ISG+ CD8", "Cycling T"
)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")

tcell_colors <- c(
  "Naive CD4" = "#8FBBD9",
  "CD4 Tem" = "#5AA6B8",
  "CD4 Tex" = "#387AA3",
  "Treg" = "#78B77D",
  "Naive CD8" = "#C8A96A",
  "CXCR4+ CD8" = "#D6C15B",
  "GZMB+ CD8" = "#D98B6B",
  "ISG+ CD8" = "#B779A5",
  "Cycling T" = "#A9A9A9"
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

theme_pub <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black", size = base_size + 1),
      axis.text = element_text(colour = "black", size = base_size),
      axis.line = element_line(colour = "black", linewidth = 0.42),
      axis.ticks = element_line(colour = "black", linewidth = 0.32),
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey91", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 6, 4, 4)
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

if (!file.exists(object_file)) stop("Cannot find urine T-cell object: ", object_file)

message("Loading urine T-cell object and calculating CXCR3-positive fractions by subtype...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(object_file)
DefaultAssay(tcell) <- "RNA"
if (!"CXCR3" %in% rownames(tcell)) stop("CXCR3 was not found in Urine_Tcell.rds.")
if (!"revised_celltype" %in% colnames(tcell@meta.data)) stop("revised_celltype was not found in Urine_Tcell.rds metadata.")

cxcr3_df <- FetchData(tcell, vars = "CXCR3") %>%
  rownames_to_column("cell") %>%
  as_tibble() %>%
  left_join(
    tcell@meta.data %>%
      rownames_to_column("cell") %>%
      transmute(
        cell,
        Group = as.character(Group),
        Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
        T_celltype = as.character(revised_celltype)
      ),
    by = "cell"
  ) %>%
  filter(!is.na(T_celltype), T_celltype != "") %>%
  mutate(
    T_celltype = if_else(T_celltype %in% tcell_levels, T_celltype, T_celltype),
    CXCR3_positive = CXCR3 > 0
  )

subtype_summary <- cxcr3_df %>%
  group_by(T_celltype) %>%
  summarise(
    total_t_cells = n(),
    CXCR3_positive_cells = sum(CXCR3_positive, na.rm = TRUE),
    pct_CXCR3_positive = CXCR3_positive_cells / total_t_cells * 100,
    mean_CXCR3_expression = mean(CXCR3, na.rm = TRUE),
    median_CXCR3_expression = median(CXCR3, na.rm = TRUE),
    n_samples = n_distinct(Sample[!is.na(Sample) & Sample != ""]),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_CXCR3_positive), desc(total_t_cells)) %>%
  mutate(
    rank = row_number(),
    T_celltype = factor(T_celltype, levels = rev(T_celltype)),
    bar_label = paste0(sprintf("%.1f", pct_CXCR3_positive), "%")
  )

# Export the analysis source data for figure reproduction.
write_csv_safe(cxcr3_df, file.path(result_dir, "Tcell_CXCR3_positive_by_subtype_cell_source.csv"))
write_csv_safe(subtype_summary, file.path(result_dir, "Tcell_CXCR3_positive_ratio_by_subtype_summary.csv"))

group_subtype_summary <- cxcr3_df %>%
  filter(!is.na(Group), Group %in% group_levels) %>%
  group_by(Group, T_celltype) %>%
  summarise(
    total_t_cells = n(),
    CXCR3_positive_cells = sum(CXCR3_positive, na.rm = TRUE),
    pct_CXCR3_positive = CXCR3_positive_cells / total_t_cells * 100,
    mean_CXCR3_expression = mean(CXCR3, na.rm = TRUE),
    median_CXCR3_expression = median(CXCR3, na.rm = TRUE),
    n_samples = n_distinct(Sample[!is.na(Sample) & Sample != ""]),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    T_celltype = factor(T_celltype, levels = levels(subtype_summary$T_celltype)),
    pct_label = paste0(sprintf("%.1f", pct_CXCR3_positive), "%")
  )

write_csv_safe(group_subtype_summary, file.path(result_dir, "Tcell_CXCR3_positive_ratio_by_group_subtype_summary.csv"))

plot_colors <- tcell_colors[names(tcell_colors) %in% as.character(subtype_summary$T_celltype)]
missing_colors <- setdiff(as.character(subtype_summary$T_celltype), names(plot_colors))
if (length(missing_colors) > 0) {
  extra_cols <- setNames(hue_pal(l = 65, c = 55)(length(missing_colors)), missing_colors)
  plot_colors <- c(plot_colors, extra_cols)
}

p <- ggplot(subtype_summary, aes(x = 1, y = T_celltype)) +
  geom_point(
    aes(size = pct_CXCR3_positive, fill = mean_CXCR3_expression),
    shape = 21,
    colour = "black",
    stroke = 0.28,
    alpha = 0.96
  ) +
  geom_text(
    aes(x = 1.18, label = bar_label),
    hjust = 0,
    size = 3.05,
    family = "Arial",
    colour = "black"
  ) +
  scale_x_continuous(
    breaks = 1,
    labels = "CXCR3",
    limits = c(0.78, 1.48),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_size_area(
    name = "CXCR3+ cells",
    max_size = 9.5,
    limits = c(0, 80),
    breaks = c(20, 40, 60, 80),
    labels = function(x) paste0(x, "%")
  ) +
  scale_fill_gradientn(
    name = "Mean CXCR3",
    colours = c("#D9D9D9", "#FEE08B", "#D73027"),
    values = scales::rescale(c(0, 0.45, max(subtype_summary$mean_CXCR3_expression, na.rm = TRUE)))
  ) +
  guides(
    fill = guide_colourbar(order = 1, barheight = unit(1.9, "cm"), barwidth = unit(0.32, "cm")),
    size = guide_legend(order = 2, override.aes = list(fill = "#8DAFCB", colour = "black"))
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub(9.6) +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.ticks = element_blank(),
    axis.line.x = element_blank(),
    legend.title = element_text(face = "bold", size = 8.6),
    legend.text = element_text(size = 8.2),
    plot.margin = margin(4, 36, 4, 4)
  )

save_plot_safe(p, "Tcell_CXCR3_positive_ratio_by_subtype_barplot", width = 3.9, height = 3.55)
save_plot_safe(p, "Tcell_CXCR3_positive_ratio_by_subtype_dotplot", width = 3.9, height = 3.55)

p_group <- ggplot(group_subtype_summary, aes(x = Group, y = T_celltype)) +
  geom_point(
    aes(size = pct_CXCR3_positive, fill = mean_CXCR3_expression),
    shape = 21,
    colour = "black",
    stroke = 0.25,
    alpha = 0.96
  ) +
  scale_size_area(
    name = "CXCR3+ cells",
    max_size = 8.4,
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100),
    labels = function(x) paste0(x, "%")
  ) +
  scale_fill_gradientn(
    name = "Mean CXCR3",
    colours = c("#D9D9D9", "#FEE08B", "#D73027"),
    values = scales::rescale(c(0, 0.55, max(group_subtype_summary$mean_CXCR3_expression, na.rm = TRUE)))
  ) +
  guides(
    fill = guide_colourbar(order = 1, barheight = unit(1.8, "cm"), barwidth = unit(0.32, "cm")),
    size = guide_legend(order = 2, override.aes = list(fill = "#8DAFCB", colour = "black"))
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub(9.3) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
    legend.title = element_text(face = "bold", size = 8.4),
    legend.text = element_text(size = 8.0),
    plot.margin = margin(4, 8, 4, 4)
  )

save_plot_safe(p_group, "Tcell_CXCR3_positive_ratio_by_group_subtype_dotplot", width = 4.8, height = 3.75)

message("Done. Output saved to: ", figure_dir)

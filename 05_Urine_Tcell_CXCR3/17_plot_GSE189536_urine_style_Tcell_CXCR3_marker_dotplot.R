#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: External Harmony annotation and final CXCR3/module validation.
# Original source: plot_GSE189536_urine_style_Tcell_CXCR3_marker_dotplot.R
# Figure/table scope: FigS10a; FigS10b; FigS10c; FigS10d; FigS10e; FigS10f
# Data cohort: GSE189536. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# GSE189536 T-cell CXCR3 sample-level metrics and urine-style T-cell
# subtype marker dotplot.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

env_data_dir <- Sys.getenv("GSE189536_DATA_DIR", unset = "")
if (nzchar(env_data_dir)) {
  data_dir <- sub("/$", "", gsub("\\\\", "/", env_data_dir))
} else {
  data_dir <- file.path(getOption("hla.data_root"), "GSE189536")
}

annotation_dir <- file.path(data_dir, "Harmony_immune_cells_manual_annotation_results", "Urine_style_annotation")
input_rds <- file.path(data_dir, "Revision_rejection_axis_validation", "Objects", "final_objects", "GSE189536_Tcell_PCA_no_batch_res0.65_reclustered.rds")
fallback_immune_rds <- file.path(annotation_dir, "GSE189536_urine_style_annotated_immune_seurat.rds")
harmony_tcell_rds <- file.path(annotation_dir, "GSE189536_urine_style_Tcell_harmony_reclustered_seurat.rds")
figure_dir <- file.path(annotation_dir, "Figures", "Tcell_CXCR3_Marker")
result_dir <- file.path(annotation_dir, "Tables", "Tcell_CXCR3_Marker")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
group_colors <- c("NR" = "#E42320", "TCMR" = "#699ECA", "ABMR" = "#4DAF4A", "Mix" = "#FF8C00", "BK" = "#F898CB")
tcell_levels <- c("Naive CD4", "Activated CD4", "CCR6+ CD4", "Treg", "GZMK+ CD8", "CXCR4+ CD8", "GZMB+ CD8", "gd/NKT-like")
tcell_colors <- c(
  "Naive CD4" = "#E42320",
  "Activated CD4" = "#699ECA",
  "CCR6+ CD4" = "#D95F02",
  "Treg" = "#7570B3",
  "GZMK+ CD8" = "#4DAF4A",
  "CXCR4+ CD8" = "#1F78B4",
  "GZMB+ CD8" = "#FF8C00",
  "gd/NKT-like" = "#8A7355"
)
cluster_annotation_map <- c(
  "0" = "GZMK+ CD8",
  "1" = "Naive CD4",
  "2" = "Activated CD4",
  "3" = "GZMB+ CD8",
  "4" = "CXCR4+ CD8",
  "5" = "CCR6+ CD4",
  "6" = "Treg",
  "7" = "gd/NKT-like"
)

marker_sets <- list(
  "Naive CD4" = c("IL7R", "MAL"),
  "Activated CD4" = c("FOS", "CD69"),
  "CCR6+ CD4" = c("CCR6", "RORC"),
  "Treg" = c("FOXP3", "CTLA4"),
  "GZMK+ CD8" = c("CD8A", "GZMK"),
  "CXCR4+ CD8" = c("CD8A", "CXCR4"),
  "GZMB+ CD8" = c("GZMB", "NKG7"),
  "gd/NKT-like" = c("TRDC", "KLRD1")
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "p=NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
  )
}

theme_pub <- function(base_size = 8.8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      strip.background = element_rect(fill = "grey93", colour = NA),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 6, 4, 6)
    )
}

save_plot_safe <- function(plot, filename, width, height) {
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(filename, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(
    png_path, plot, width = width, height = height, units = "in", dpi = 600, bg = "white",
    device = function(filename, width, height, units, res, bg, ...) {
      grDevices::png(filename, width = width, height = height, units = units, res = res, bg = bg, type = "cairo")
    }
  )
  invisible(c(pdf = pdf_path, png = png_path))
}

safe_wilcox <- function(values, groups, g1, g2) {
  keep <- is.finite(values) & groups %in% c(g1, g2)
  values <- values[keep]
  groups <- factor(groups[keep], levels = c(g1, g2))
  if (nlevels(droplevels(groups)) < 2 || any(table(groups) < 2) || dplyr::n_distinct(values) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(values ~ groups, exact = FALSE)$p.value)
}

message("Loading GSE189536 T-cell object...")
tobj <- if (file.exists(harmony_tcell_rds)) {
  message("Using Harmony-reclustered T-cell object.")
  readRDS(harmony_tcell_rds)
} else if (file.exists(input_rds)) {
  readRDS(input_rds)
} else if (file.exists(fallback_immune_rds)) {
  message("T-cell-only object was not found; extracting annotated T cells from urine-style immune object.")
  immune_obj <- readRDS(fallback_immune_rds)
  t_meta <- immune_obj@meta.data
  t_keep <- !is.na(t_meta$Urine_Tcell_revised_celltype) & t_meta$Urine_Tcell_revised_celltype != ""
  subset(immune_obj, cells = rownames(t_meta)[t_keep])
} else {
  stop("Cannot find a usable GSE189536 T-cell object:\n", input_rds, "\n", fallback_immune_rds)
}
DefaultAssay(tobj) <- "RNA"

if ("seurat_clusters" %in% colnames(tobj@meta.data)) {
  cluster_label <- cluster_annotation_map[as.character(tobj$seurat_clusters)]
  cluster_label[is.na(cluster_label)] <- paste0("Cluster ", as.character(tobj$seurat_clusters)[is.na(cluster_label)])
  names(cluster_label) <- colnames(tobj)
  tobj$Tcell_cluster_annotation <- factor(cluster_label, levels = tcell_levels)
}

if (!"CXCR3" %in% rownames(tobj)) stop("CXCR3 was not found in the T-cell object.")
cxcr3_df <- FetchData(tobj, vars = "CXCR3")
cxcr3_expr <- cxcr3_df[, 1]
names(cxcr3_expr) <- rownames(cxcr3_df)
lineage_genes <- intersect(c("CD4", "CD8A", "CD8B"), rownames(tobj))
lineage_df <- FetchData(tobj, vars = lineage_genes) %>%
  rownames_to_column("cell")

subtype_col <- dplyr::case_when(
  "Tcell_cluster_annotation" %in% colnames(tobj@meta.data) ~ "Tcell_cluster_annotation",
  "tcell_type" %in% colnames(tobj@meta.data) ~ "tcell_type",
  "Urine_Tcell_revised_celltype" %in% colnames(tobj@meta.data) ~ "Urine_Tcell_revised_celltype",
  TRUE ~ NA_character_
)
if (is.na(subtype_col)) stop("No T-cell subtype annotation was found.")

meta <- tobj@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(lineage_df, by = "cell") %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    Biopsy = as.character(biopsy),
    Tcell_subtype_raw = as.character(.data[[subtype_col]]),
    CD4_expr = if ("CD4" %in% colnames(.)) as.numeric(CD4) else 0,
    CD8A_expr = if ("CD8A" %in% colnames(.)) as.numeric(CD8A) else 0,
    CD8B_expr = if ("CD8B" %in% colnames(.)) as.numeric(CD8B) else 0,
    Tcell_subtype = if_else(
      subtype_col != "Tcell_cluster_annotation" &
        Tcell_subtype_raw == "CD4 Tex" &
        (CD8A_expr + CD8B_expr) > CD4_expr &
        (CD8A_expr + CD8B_expr) > 0,
      "CD8 Tex-like",
      Tcell_subtype_raw
    ),
    Tcell_subtype = factor(Tcell_subtype, levels = tcell_levels),
    CXCR3_expr = as.numeric(cxcr3_expr[cell]),
    CXCR3_positive = CXCR3_expr > 0
  ) %>%
  filter(!is.na(Group), !is.na(Biopsy), Biopsy != "", !is.na(Tcell_subtype))

if ("Tcell_cluster_annotation" %in% colnames(tobj@meta.data)) {
# Export the analysis source data for figure reproduction.
  saveRDS(tobj, file.path(annotation_dir, "GSE189536_urine_style_Tcell_harmony_reclustered_cluster_annotated_seurat.rds"))
  write_csv_safe(
    meta %>%
      count(seurat_clusters, Tcell_subtype, Group, name = "cells") %>%
      arrange(as.numeric(as.character(seurat_clusters)), Tcell_subtype, Group),
    file.path(result_dir, "GSE189536_Tcell_harmony_cluster_reannotation_counts.csv")
  )
}

available_groups <- intersect(group_levels, unique(as.character(meta$Group)))

cxcr3_group_subtype_summary <- meta %>%
  group_by(Group, Tcell_subtype) %>%
  summarise(
    cells = n(),
    mean_CXCR3 = mean(CXCR3_expr, na.rm = TRUE),
    pct_CXCR3_positive = mean(CXCR3_positive, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(as.character(Group), levels = available_groups),
    Tcell_subtype = factor(as.character(Tcell_subtype), levels = rev(tcell_levels))
  )
write_csv_safe(
  cxcr3_group_subtype_summary,
  file.path(result_dir, "GSE189536_Tcell_CXCR3_by_group_subtype_dotplot_summary.csv")
)

cxcr3_subtype_order <- meta %>%
  group_by(Tcell_subtype) %>%
  summarise(
    cells = n(),
    mean_CXCR3 = mean(CXCR3_expr, na.rm = TRUE),
    median_CXCR3 = median(CXCR3_expr, na.rm = TRUE),
    pct_CXCR3_positive = mean(CXCR3_positive, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_CXCR3))
write_csv_safe(
  cxcr3_subtype_order,
  file.path(result_dir, "GSE189536_Tcell_CXCR3_subtype_expression_order.csv")
)

cxcr3_plot_groups <- intersect(c("NR", "TCMR", "ABMR"), available_groups)
cxcr3_sample_subtype_df <- meta %>%
  filter(as.character(Group) %in% cxcr3_plot_groups) %>%
  group_by(Group, Biopsy, diagnosis, Tcell_subtype) %>%
  summarise(
    n_T_cells = n(),
    mean_CXCR3 = mean(CXCR3_expr, na.rm = TRUE),
    median_CXCR3 = median(CXCR3_expr, na.rm = TRUE),
    pct_CXCR3_positive = mean(CXCR3_positive, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(as.character(Group), levels = cxcr3_plot_groups),
    Tcell_subtype = factor(as.character(Tcell_subtype), levels = rev(as.character(cxcr3_subtype_order$Tcell_subtype)))
  )
write_csv_safe(
  cxcr3_sample_subtype_df,
  file.path(result_dir, "GSE189536_Tcell_CXCR3_sample_level_by_group_subtype_boxpoint.csv")
)

p_cxcr3_group_dot <- ggplot(cxcr3_sample_subtype_df, aes(x = mean_CXCR3, y = Tcell_subtype, fill = Group)) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    linewidth = 0.30,
    colour = "grey25",
    alpha = 0.78
  ) +
  geom_point(
    aes(colour = Group),
    position = position_jitter(width = 0, height = 0.10, seed = 20260827),
    shape = 21,
    fill = "white",
    size = 1.65,
    stroke = 0.38,
    alpha = 0.92
  ) +
  facet_wrap(~Group, nrow = 1) +
  scale_fill_manual(values = group_colors[cxcr3_plot_groups], drop = TRUE, guide = "none") +
  scale_colour_manual(values = group_colors[cxcr3_plot_groups], drop = TRUE, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = "Mean CXCR3 expression per biopsy", y = NULL) +
  theme_pub(9.2) +
  theme(
    axis.title.x = element_text(face = "bold", size = 9.4),
    axis.text.x = element_text(size = 8.2, colour = "black"),
    axis.text.y = element_text(face = "bold", size = 8.4, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.34),
    axis.ticks = element_line(colour = "black", linewidth = 0.28),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text = element_text(face = "bold", size = 9.0, colour = "black"),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.22),
    panel.grid.major.y = element_blank(),
    panel.spacing.x = unit(0.14, "in"),
    plot.margin = margin(4, 5, 4, 4)
  )
save_plot_safe(p_cxcr3_group_dot, "GSE189536_Tcell_CXCR3_expression_by_group_subtype_dotplot", width = 5.9, height = 3.45)
save_plot_safe(p_cxcr3_group_dot, "GSE189536_Tcell_CXCR3_expression_by_group_subtype_boxpoint_NR_TCMR_ABMR", width = 5.9, height = 3.45)

sample_metrics <- meta %>%
  group_by(Group, diagnosis, Biopsy) %>%
  summarise(
    n_T_cells = n(),
    mean_CXCR3_T = mean(CXCR3_expr, na.rm = TRUE),
    median_CXCR3_T = median(CXCR3_expr, na.rm = TRUE),
    CXCR3_positive_T_cells = sum(CXCR3_positive, na.rm = TRUE),
    CXCR3_positive_T_percent = mean(CXCR3_positive, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(Group = factor(as.character(Group), levels = group_levels))

metric_long <- sample_metrics %>%
  select(Group, diagnosis, Biopsy, n_T_cells, mean_CXCR3_T, CXCR3_positive_T_percent) %>%
  pivot_longer(
    cols = c(mean_CXCR3_T, CXCR3_positive_T_percent),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("mean_CXCR3_T", "CXCR3_positive_T_percent"),
      labels = c("Mean CXCR3 in T cells", "CXCR3+ T cells (%)")
    )
  )

comparison_grid <- tibble(
  group1 = "TCMR",
  group2 = setdiff(available_groups, "TCMR")
) %>%
  filter(group1 %in% available_groups, group2 %in% available_groups)

wilcox_stats <- metric_long %>%
  group_by(Metric) %>%
  group_modify(function(dat, key) {
    comparison_grid %>%
      mutate(
        p_value = mapply(
          function(g1, g2) safe_wilcox(dat$Value, as.character(dat$Group), g1, g2),
          group1, group2
        )
      )
  }) %>%
  ungroup() %>%
  group_by(Metric) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"), label = format_p(p_value)) %>%
  ungroup()

plot_stats <- wilcox_stats %>%
  filter(!is.na(p_value)) %>%
  mutate(
    group1 = factor(group1, levels = group_levels),
    group2 = factor(group2, levels = group_levels),
    x = as.numeric(group1),
    xend = as.numeric(group2)
  ) %>%
  left_join(
    metric_long %>%
      group_by(Metric) %>%
      summarise(y_min = min(Value, na.rm = TRUE), y_max = max(Value, na.rm = TRUE), .groups = "drop") %>%
      mutate(y_span = pmax(y_max - y_min, 0.1), step_y = y_span * 0.19),
    by = "Metric"
  ) %>%
  group_by(Metric) %>%
  arrange(p_value, .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    y = y_max + step_y * (rank + 0.10),
    y_text = y + step_y * 0.18,
    tick = step_y * 0.12
  ) %>%
  ungroup()

write_csv_safe(sample_metrics, file.path(result_dir, "GSE189536_Tcell_CXCR3_sample_level_metrics.csv"))
write_csv_safe(wilcox_stats, file.path(result_dir, "GSE189536_Tcell_CXCR3_sample_level_wilcox_tests.csv"))

p_cxcr3 <- ggplot(metric_long, aes(Group, Value, fill = Group)) +
  geom_boxplot(
    width = 0.54,
    outlier.shape = NA,
    colour = "grey20",
    linewidth = 0.38,
    alpha = 0.52
  ) +
  geom_point(
    aes(colour = Group),
    position = position_jitter(width = 0.10, height = 0, seed = 20260827),
    shape = 16,
    stroke = 0,
    size = 2.35,
    alpha = 0.95
  ) +
  geom_segment(data = plot_stats, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.32, colour = "black") +
  geom_segment(data = plot_stats, aes(x = x, xend = x, y = y - tick, yend = y), inherit.aes = FALSE, linewidth = 0.32, colour = "black") +
  geom_segment(data = plot_stats, aes(x = xend, xend = xend, y = y - tick, yend = y), inherit.aes = FALSE, linewidth = 0.32, colour = "black") +
  geom_text(data = plot_stats, aes(x = (x + xend) / 2, y = y_text, label = label), inherit.aes = FALSE, family = "Arial", fontface = "bold", size = 3.15) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_colors[available_groups], drop = TRUE, name = "Group") +
  scale_colour_manual(values = group_colors[available_groups], drop = TRUE, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = NULL, y = NULL) +
  theme_pub(9.0) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.0),
    axis.text.y = element_text(size = 8.0),
    strip.text = element_text(size = 8.4),
    legend.position = "top",
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_blank()
  )
save_plot_safe(p_cxcr3, "GSE189536_Tcell_CXCR3_sample_level_bar_wilcox_TCMR_vs_others", width = 6.2, height = 3.6)
save_plot_safe(p_cxcr3, "GSE189536_Tcell_CXCR3_sample_level_scatter_bar_wilcox", width = 6.2, height = 3.6)
save_plot_safe(p_cxcr3, "GSE189536_Tcell_CXCR3_sample_level_scatter_box_wilcox_TCMR_vs_others", width = 6.2, height = 3.6)

marker_df <- tibble(
  Tcell_subtype = rep(names(marker_sets), lengths(marker_sets)),
  gene = unlist(marker_sets, use.names = FALSE)
) %>%
  filter(Tcell_subtype %in% unique(as.character(meta$Tcell_subtype)), gene %in% rownames(tobj)) %>%
  mutate(
    Tcell_subtype = factor(Tcell_subtype, levels = tcell_levels),
    gene = factor(gene, levels = unique(gene))
  )
write_csv_safe(marker_df, file.path(result_dir, "GSE189536_Tcell_subtype_marker_genes_used.csv"))

plot_subtypes <- intersect(tcell_levels, unique(as.character(meta$Tcell_subtype)))
tcell_subtype_plot <- setNames(as.character(meta$Tcell_subtype), meta$cell)
tobj$Tcell_subtype_plot <- factor(tcell_subtype_plot[colnames(tobj)], levels = plot_subtypes)
write_csv_safe(
  meta %>%
    count(Tcell_subtype_raw, Tcell_subtype, name = "n_cells") %>%
    arrange(Tcell_subtype),
  file.path(result_dir, "GSE189536_Tcell_subtype_plot_annotation_counts.csv")
)

p_cluster_anno <- NULL
if ("tcell_harmony_umap" %in% names(tobj@reductions)) {
  umap_df <- Embeddings(tobj, "tcell_harmony_umap") %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    left_join(
      meta %>% select(cell, seurat_clusters, Tcell_subtype, Group),
      by = "cell"
    )
  colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")
  p_cluster_anno <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Tcell_subtype)) +
    geom_point(size = 0.30, alpha = 0.86, stroke = 0) +
    coord_equal() +
    scale_colour_manual(values = tcell_colors, drop = TRUE, na.value = "grey82") +
    labs(x = "UMAP 1", y = "UMAP 2", colour = "T-cell subtype") +
    theme_classic(base_size = 9.0, base_family = "Arial") +
    theme(
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      legend.title = element_text(face = "bold", size = 8.4),
      legend.text = element_text(size = 7.8),
      plot.margin = margin(4, 6, 4, 6)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3.2, alpha = 1)))
  save_plot_safe(p_cluster_anno, "GSE189536_Tcell_harmony_cluster_reannotated_UMAP", width = 4.8, height = 4.1)

  p_group_facet <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Tcell_subtype)) +
    geom_point(size = 0.24, alpha = 0.86, stroke = 0) +
    facet_wrap(~Group, ncol = 2) +
    coord_equal() +
    scale_colour_manual(values = tcell_colors, drop = TRUE, na.value = "grey82") +
    labs(x = "UMAP 1", y = "UMAP 2", colour = "T-cell subtype") +
    theme_classic(base_size = 9.0, base_family = "Arial") +
    theme(
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.32),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold", size = 8.6, colour = "black"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 8.2),
      legend.text = element_text(size = 7.6),
      legend.key.height = unit(0.16, "in"),
      legend.key.width = unit(0.16, "in"),
      panel.spacing.x = unit(0.28, "in"),
      plot.margin = margin(4, 6, 4, 6)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3.0, alpha = 1)))
  save_plot_safe(p_group_facet, "GSE189536_Tcell_harmony_reannotated_UMAP_by_group_facet", width = 5.9, height = 4.85)
}

p_marker_base <- DotPlot(
  tobj,
  features = unique(as.character(marker_df$gene)),
  group.by = "Tcell_subtype_plot",
  cols = c("#E6E6E6", "#B7192B"),
  dot.scale = 4.2
) 
p_marker_data <- p_marker_base$data %>%
  mutate(
    id = factor(as.character(id), levels = plot_subtypes),
    features.plot = factor(as.character(features.plot), levels = rev(unique(as.character(marker_df$gene))))
  )
p_marker <- ggplot(p_marker_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, colour = avg.exp.scaled), alpha = 0.95) +
  scale_colour_gradientn(
    colours = c("#E6E6E6", "#E7AAA3", "#B7192B"),
    name = "Average\nexpression"
  ) +
  scale_size(
    range = c(0.2, 4.2),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75),
    name = "Percent\nexpressed"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10.0, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.8, colour = "black"),
    axis.text.y = element_text(size = 9.0, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.35),
    axis.ticks = element_line(colour = "black", linewidth = 0.30),
    legend.title = element_text(face = "bold", size = 8.4),
    legend.text = element_text(size = 7.6),
    legend.key.height = unit(0.18, "in"),
    legend.key.width = unit(0.18, "in"),
    plot.margin = margin(2, 3, 2, 4)
  )
save_plot_safe(p_marker, "GSE189536_Tcell_subtype_marker_gene_dotplot", width = 4.7, height = 5.4)

if (!is.null(p_cluster_anno)) {
  p_umap_compact <- p_cluster_anno +
    labs(title = NULL) +
    theme(
      axis.title = element_text(size = 8.4, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 7.2, face = "bold"),
      legend.text = element_text(size = 6.8),
      legend.key.height = unit(0.14, "in"),
      legend.key.width = unit(0.14, "in"),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    guides(colour = guide_legend(ncol = 4, override.aes = list(size = 2.4, alpha = 1)))

  p_marker_compact <- p_marker +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8.0, colour = "black"),
      axis.text.y = element_text(size = 8.4, colour = "black"),
      legend.title = element_text(face = "bold", size = 7.8),
      legend.text = element_text(size = 7.0),
      legend.key.height = unit(0.16, "in"),
      legend.key.width = unit(0.16, "in"),
      plot.margin = margin(2, 2, 2, 2)
    )

  p_umap_marker <- p_umap_compact + p_marker_compact + plot_layout(ncol = 2, widths = c(0.82, 1.18))
  save_plot_safe(p_umap_marker, "GSE189536_Tcell_reannotated_UMAP_marker_dotplot_combined", width = 9.1, height = 3.55)
}

message("Saved GSE189536 T-cell CXCR3 sample-level plots and subtype marker dotplot.")

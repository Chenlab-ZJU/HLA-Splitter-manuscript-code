# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 11b_revision_redraw_tcell_umap_from_csf_tcell.R
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
  library(ggplot2)
  library(ggrepel)
  library(gtools)
  library(cowplot)
  library(patchwork)
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
group_colors <- c("Metastases" = "#E42320", "Hemorrhage" = "#699ECA", "Encephalitis" = "#4DAF4A")
cluster_celltype_map <- c(
  "0" = "CD4 memory T",
  "1" = "Naive CD4 T",
  "2" = "GZMK+ CD8 T",
  "3" = "GZMK+ CD8 T",
  "4" = "GZMB+ cytotoxic T",
  "5" = "NK-like T",
  "6" = "Treg",
  "7" = "GZMK+ CD4 T",
  "8" = "Naive CD8 T",
  "9" = "NK-like T",
  "10" = "ISG+ CD4 T",
  "11" = "NK-like T",
  "12" = "NK-like T",
  "13" = "gd T"
)
tcell_colors <- c(
  "Naive CD4 T" = "#E42320",
  "CD4 memory T" = "#699ECA",
  "GZMK+ CD4 T" = "#4DAF4A",
  "ISG+ CD4 T" = "#FF8C00",
  "Treg" = "#731A73",
  "Naive CD8 T" = "#FFCB5B",
  "GZMK+ CD8 T" = "#0076B9",
  "GZMB+ cytotoxic T" = "#3D505A",
  "NK-like T" = "#0098B2",
  "gd T" = "#8A7355"
)

theme_pub <- function(base_size = 9.8) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
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

message("Loading CSF_Tcell.rds without reclustering...")
# Read the saved object or source table used by this analysis.
Tcell <- readRDS(file.path(getOption("hla.data_root"), "csf", "CSF_Tcell.rds"))

Tcell$T_celltype <- unname(cluster_celltype_map[as.character(Tcell$seurat_clusters)])
Tcell <- subset(Tcell, subset = !is.na(T_celltype) & T_celltype %in% names(tcell_colors))
Tcell$T_celltype <- factor(as.character(Tcell$T_celltype), levels = names(tcell_colors))
Tcell$Group <- factor(as.character(Tcell$Group), levels = group_levels)
Tcell$Cluster <- factor(as.character(Tcell$seurat_clusters), levels = gtools::mixedsort(unique(as.character(Tcell$seurat_clusters))))

marker_check_genes <- c(
  "CD4", "IL7R", "CCR7", "SELL", "TCF7", "LEF1", "MAL",
  "CD8A", "CD8B", "GZMK", "GZMB", "PRF1", "NKG7", "GNLY", "KLRD1",
  "FOXP3", "IL2RA", "CTLA4", "ISG15", "IFIT1", "IFIT3", "MX1", "TRDC",
  "HLA-DRA", "CCR6", "AQP3", "KLRB1"
)
marker_check_genes <- intersect(marker_check_genes, rownames(Tcell))
if (length(marker_check_genes) > 0) {
  avg_expr <- Seurat::AverageExpression(
    Tcell,
    features = marker_check_genes,
    group.by = "Cluster",
    assays = DefaultAssay(Tcell),
    slot = "data"
  )[[DefaultAssay(Tcell)]]
  avg_expr_df <- as.data.frame(as.matrix(avg_expr)) %>%
    tibble::rownames_to_column("gene")
# Export the analysis source data for figure reproduction.
  write.csv(avg_expr_df, file.path(result_dir, "Tcell_revision_manual_annotation_marker_average_expression.csv"), row.names = FALSE)
}

umap_mat <- Embeddings(Tcell, "umap")[, 1:2, drop = FALSE]
meta_df <- Tcell@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  transmute(
    cell,
    Cluster = Cluster,
    T_celltype = T_celltype,
    Group = Group
  )
umap_df <- as.data.frame(umap_mat) %>%
  setNames(c("UMAP_1", "UMAP_2")) %>%
  tibble::rownames_to_column("cell") %>%
  left_join(meta_df, by = "cell") %>%
  filter(!is.na(Cluster), !is.na(T_celltype), !is.na(Group))

write.csv(umap_df, file.path(result_dir, "Tcell_revision_umap_metadata_from_CSF_Tcell.csv"), row.names = FALSE)
write.csv(
  data.frame(
    Cluster = names(cluster_celltype_map),
    T_celltype = unname(cluster_celltype_map),
    color = unname(tcell_colors[unname(cluster_celltype_map)])
  ),
  file.path(getOption("hla.data_root"), "csf", "Tcell_revision_manual_cluster_reannotation.csv"),
  row.names = FALSE
)

umap_xlim <- range(umap_df$UMAP_1, finite = TRUE)
umap_ylim <- range(umap_df$UMAP_2, finite = TRUE)
coord_tcell_umap <- coord_cartesian(xlim = umap_xlim, ylim = umap_ylim, expand = FALSE)
tcell_colors_present <- tcell_colors[levels(droplevels(umap_df$T_celltype))]
tcell_colors_present <- tcell_colors_present[!is.na(tcell_colors_present)]
write.csv(
  data.frame(T_celltype = names(tcell_colors_present), color = unname(tcell_colors_present)),
  file.path(result_dir, "Tcell_revision_celltype_color_map.csv"),
  row.names = FALSE
)

cluster_centers <- umap_df %>%
  group_by(Cluster) %>%
  summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE), UMAP_2 = median(UMAP_2, na.rm = TRUE), .groups = "drop")

cluster_annotation <- umap_df %>%
  count(Cluster, T_celltype, name = "n") %>%
  group_by(Cluster) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()
cluster_palette <- setNames(
  unname(tcell_colors_present[as.character(cluster_annotation$T_celltype)]),
  as.character(cluster_annotation$Cluster)
)
write.csv(
  cluster_annotation %>%
    mutate(color = cluster_palette[as.character(Cluster)]) %>%
    arrange(gtools::mixedorder(as.character(Cluster))),
  file.path(result_dir, "Tcell_revision_cluster_celltype_color_map.csv"),
  row.names = FALSE
)

p_tcell_cluster <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Cluster)) +
  geom_point(size = 0.18, alpha = 0.82, stroke = 0) +
  geom_label_repel(
    data = cluster_centers,
    aes(x = UMAP_1, y = UMAP_2, label = Cluster),
    inherit.aes = FALSE,
    size = 3.0,
    family = "Arial",
    fontface = "bold",
    label.size = 0.15,
    label.padding = unit(0.10, "lines"),
    fill = scales::alpha("white", 0.82),
    colour = "black",
    seed = 20260813,
    max.overlaps = Inf
  ) +
  scale_colour_manual(values = cluster_palette, drop = TRUE) +
  coord_tcell_umap +
  labs(title = "T cell clusters", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(legend.position = "none", panel.grid = element_blank(), aspect.ratio = 1.08)

p_tcell_type <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = T_celltype)) +
  geom_point(size = 0.18, alpha = 0.82, stroke = 0) +
  scale_colour_manual(values = tcell_colors_present, limits = names(tcell_colors_present), drop = FALSE, name = "T cell subset") +
  guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1), ncol = 1)) +
  coord_tcell_umap +
  labs(title = "T cell subsets", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.44, "cm"),
    legend.key.width = unit(0.44, "cm"),
    panel.grid = element_blank(),
    aspect.ratio = 1.08
  )

p_tcell_group <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Group)) +
  geom_point(size = 0.18, alpha = 0.78, stroke = 0) +
  scale_colour_manual(values = group_colors, drop = FALSE, name = "Group") +
  guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1), ncol = 1)) +
  coord_tcell_umap +
  labs(title = "Clinical groups", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.44, "cm"),
    legend.key.width = unit(0.44, "cm"),
    panel.grid = element_blank(),
    aspect.ratio = 1.08
  )

tcell_subset_legend <- cowplot::get_legend(p_tcell_type)
tcell_group_legend <- cowplot::get_legend(p_tcell_group)
tcell_umap_panels <- (
  p_tcell_cluster +
    p_tcell_type + theme(legend.position = "none") +
    p_tcell_group + theme(legend.position = "none") +
    plot_layout(nrow = 1, widths = c(1, 1, 1))
)
tcell_umap_legends <- cowplot::plot_grid(
  tcell_subset_legend,
  tcell_group_legend,
  ncol = 1,
  rel_heights = c(1.7, 0.8),
  align = "v"
)
tcell_umap_combined <- cowplot::plot_grid(
  tcell_umap_panels,
  tcell_umap_legends,
  nrow = 1,
  rel_widths = c(3.0, 0.85),
  align = "h",
  axis = "tb"
)

save_plot_safe(tcell_umap_combined, "CSF_UMAP_Tcell_celltype_revision", width = 9.8, height = 4.2)
message("Redrew CSF_UMAP_Tcell_celltype_revision from CSF_Tcell.rds using manual cluster-level T cell annotations.")

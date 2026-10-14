# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 17_revision_tcell_celltype_heatmap_align_umap.R
# Figure/table scope: Fig2b; Fig2c; Fig2d; Fig2e; Fig2f; FigS4c; FigS4d; FigS4e; FigS5a; FigS5b; FigS5c; FigS5d; FigS5e; FigS5f
# Data cohort: 脑脊液样本. Inputs retain the original Data/Result hierarchy.
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
})

project_dir <- Sys.getenv("CSF_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "\u8111\u810a\u6db2\u6837\u672c")
}
project_dir <- gsub("\\\\", "/", project_dir)
figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

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

tcell_marker_modules <- list(
  "Naive/memory" = c("CCR7", "IL7R", "LEF1", "TCF7", "SELL", "MAL", "LTB"),
  "CD4 activation/Treg" = c("CD40LG", "CCR6", "AQP3", "FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "Cytotoxic" = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TOX", "ICOS", "MAF"),
  "IFN/chemokine" = c("IFNG", "CXCR3", "CCL3", "CCL4", "ISG15", "IFI6", "IFIT3", "MX1"),
  "Antigen presentation" = c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DPA1", "CD74")
)

theme_pub <- function(base_size = 9.5) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      strip.text.y = element_text(face = "bold", angle = 0, hjust = 0),
      strip.background = element_rect(fill = "grey93", colour = NA),
      panel.border = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
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

message("Loading CSF_Tcell.rds and applying the same manual annotations used for CSF_UMAP_Tcell_celltype_revision...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(file.path(result_dir, "CSF_Tcell.rds"))
DefaultAssay(tcell) <- "RNA"

tcell$T_celltype <- unname(cluster_celltype_map[as.character(tcell$seurat_clusters)])
tcell <- subset(tcell, subset = !is.na(T_celltype) & T_celltype %in% names(tcell_colors))
tcell$T_celltype <- factor(as.character(tcell$T_celltype), levels = names(tcell_colors))

present_celltypes <- levels(droplevels(tcell$T_celltype))
present_colors <- tcell_colors[present_celltypes]
marker_df <- tibble(
  module = rep(names(tcell_marker_modules), lengths(tcell_marker_modules)),
  gene = unlist(tcell_marker_modules, use.names = FALSE)
) %>%
  filter(gene %in% rownames(tcell)) %>%
  distinct(module, gene) %>%
  group_by(module) %>%
  mutate(module_gene = paste(module, gene, sep = "__")) %>%
  ungroup()

avg_expr <- AverageExpression(
  tcell,
  assays = "RNA",
  features = unique(marker_df$gene),
  group.by = "T_celltype",
  slot = "data"
)$RNA %>%
  as.matrix()

avg_expr <- avg_expr[, intersect(present_celltypes, colnames(avg_expr)), drop = FALSE]
z_expr <- t(scale(t(avg_expr)))
z_expr[is.na(z_expr)] <- 0
z_expr[z_expr > 2.2] <- 2.2
z_expr[z_expr < -2.2] <- -2.2

heatmap_df <- as.data.frame(z_expr) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "T_celltype", values_to = "z_score") %>%
  left_join(marker_df, by = "gene", relationship = "many-to-many") %>%
  mutate(
    module = factor(module, levels = names(tcell_marker_modules)),
    module_gene = factor(module_gene, levels = rev(marker_df$module_gene)),
    T_celltype = factor(T_celltype, levels = present_celltypes)
  )

annotation_summary <- tcell@meta.data %>%
  count(seurat_clusters, T_celltype, name = "n_cells") %>%
  arrange(as.numeric(as.character(seurat_clusters)))

# Export the analysis source data for figure reproduction.
write.csv(marker_df, file.path(result_dir, "Tcell_revision_heatmap_gene_modules_aligned_to_umap.csv"), row.names = FALSE)
write.csv(heatmap_df, file.path(result_dir, "Tcell_revision_celltype_heatmap_source.csv"), row.names = FALSE)
write.csv(annotation_summary, file.path(result_dir, "Tcell_revision_celltype_heatmap_annotation_check.csv"), row.names = FALSE)
write.csv(
  data.frame(T_celltype = present_celltypes, color = unname(present_colors)),
  file.path(result_dir, "Tcell_revision_celltype_heatmap_color_map.csv"),
  row.names = FALSE
)

p_heatmap <- ggplot(heatmap_df, aes(T_celltype, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.26) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_gradient2(
    low = "#6F99C3",
    mid = "#F7F7F7",
    high = "#B7192B",
    midpoint = 0,
    limits = c(-2.2, 2.2),
    name = "Z-score"
  ) +
  labs(x = "T cell subset", y = NULL) +
  theme_pub(9.4) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.6),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 5)),
    axis.text.y = element_text(size = 8.2),
    strip.text.y = element_text(size = 8.2),
    legend.position = "right",
    legend.key.height = unit(0.42, "cm"),
    legend.key.width = unit(0.35, "cm")
  )

save_plot_safe(p_heatmap, "Tcell_celltype_heatmap_revision", width = 6.6, height = 7.0)
message("Saved Tcell_celltype_heatmap_revision aligned to CSF_UMAP_Tcell_celltype_revision annotations.")

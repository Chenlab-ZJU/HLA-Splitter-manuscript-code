# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 10e_revision_myeloid_key_gene_feature_umap.R
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
  library(ggplot2)
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

key_plot_genes <- c("TREM2", "SPP1", "APOE", "C1QA", "CXCL9", "CXCL10", "IL1B", "HMOX1")

theme_feature_umap <- function(base_size = 9.8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 1.2),
      strip.text = element_text(face = "bold", colour = "black", size = base_size + 0.6),
      strip.background = element_rect(fill = "grey95", colour = NA),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 0.5),
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

message("Loading existing myeloid object...")
# Read the saved object or source table used by this analysis.
Myeloid <- readRDS(file.path(getOption("hla.data_root"), "csf", "CSF_Myeloid_revision.rds"))
DefaultAssay(Myeloid) <- "RNA"

genes_present <- intersect(key_plot_genes, rownames(Myeloid))
if (length(genes_present) == 0) {
  stop("None of the requested key genes were found in the myeloid object.")
}

umap_mat <- Embeddings(Myeloid, "umap")
umap_df <- as.data.frame(umap_mat) %>%
  setNames(c("UMAP_1", "UMAP_2")) %>%
  tibble::rownames_to_column("cell")

expr_df <- FetchData(Myeloid, vars = genes_present) %>%
  tibble::rownames_to_column("cell") %>%
  pivot_longer(cols = all_of(genes_present), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = key_plot_genes))

plot_df <- expr_df %>%
  left_join(umap_df, by = "cell") %>%
  group_by(gene) %>%
  mutate(
    q99 = stats::quantile(expression, 0.99, na.rm = TRUE, names = FALSE),
    expr_cap = pmin(expression, q99),
    expr_scaled = ifelse(q99 > 0, expr_cap / q99, 0),
    is_expressing = expression > 0
  ) %>%
  ungroup()

# Export the analysis source data for figure reproduction.
write.csv(
  plot_df %>% select(cell, gene, UMAP_1, UMAP_2, expression, expr_cap, expr_scaled, is_expressing),
  file.path(result_dir, "Myeloid_revision_key_gene_feature_UMAP_source.csv"),
  row.names = FALSE
)

background_df <- plot_df %>%
  distinct(gene, cell, UMAP_1, UMAP_2)

foreground_df <- plot_df %>%
  filter(is_expressing) %>%
  arrange(gene, expr_scaled)

p_feature <- ggplot() +
  geom_point(
    data = background_df,
    aes(UMAP_1, UMAP_2),
    colour = "grey82",
    size = 0.18,
    alpha = 0.62,
    stroke = 0
  ) +
  geom_point(
    data = foreground_df,
    aes(UMAP_1, UMAP_2, colour = expr_scaled),
    size = 0.28,
    alpha = 0.95,
    stroke = 0
  ) +
  facet_wrap(~gene, ncol = 2) +
  scale_colour_gradientn(
    colours = c("#D9D9D9", "#FEE08B", "#D73027", "#7F0000"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("Low", "Mid", "High"),
    name = "Relative\nexpression"
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_feature_umap(10.5) +
  theme(
    legend.position = "right",
    panel.spacing = unit(0.32, "lines")
  ) +
  guides(colour = guide_colourbar(barheight = unit(30, "mm"), barwidth = unit(4.5, "mm")))

save_plot_safe(p_feature, "Myeloid_key_gene_feature_UMAP", width = 5.8, height = 10.8)

message("Saved Myeloid_key_gene_feature_UMAP in 2 x 4 layout.")

# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Final BALF figures from saved objects and signatures.
# Original source: 23_revision_lavage_tcell_selected_module_score_umap.R
# Figure/table scope: Fig3b; Fig3c; Fig3d; Fig3e; Fig3f; Fig3g; Fig3h; Fig3i; FigS6a; FigS6b; FigS6c; FigS6d; FigS6e; FigS6f; FigS6g; FigS6h
# Data cohort: balf. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

project_dir <- Sys.getenv("LAVAGE_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "balf")
}
project_dir <- gsub("\\\\", "/", project_dir)
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures", "TCR")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

selected_modules <- list(
  "Antigen presentation" = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "HLA-A", "HLA-B", "HLA-C", "B2M"),
  "Cytotoxic" = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "CTSW", "CST7", "KLRD1", "KLRG1", "FGFBP2", "CCL5", "CD8A", "CD8B"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX", "ENTPD1", "LAYN", "CXCL13", "TNFRSF9", "MAF", "NR4A1", "NR4A2"),
  "Tissue resident" = c("ITGAE", "ITGA1", "CXCR6", "ZNF683", "CD69", "RGS1", "XCL1", "XCL2", "KLRB1", "KLRG1", "S1PR1")
)

theme_umap <- function(base_size = 9.6) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      strip.text = element_text(face = "bold", colour = "black", size = 10.5),
      strip.background = element_rect(fill = "grey94", colour = NA),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      panel.grid = element_blank(),
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
  invisible(c(pdf = target_pdf, png = target_png))
}

zscore_clip <- function(x, lower = -1.2, upper = 2.4) {
  z <- as.numeric(scale(x))
  z[is.na(z) | !is.finite(z)] <- 0
  pmax(pmin(z, upper), lower)
}

message("Loading lavage_Tcell.rds...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(file.path(getOption("hla.data_root"), "balf", "lavage_Tcell.rds"))
DefaultAssay(tcell) <- if ("RNA" %in% Assays(tcell)) "RNA" else DefaultAssay(tcell)
if (!"umap" %in% names(tcell@reductions)) {
  stop("lavage_Tcell.rds does not contain a UMAP reduction.")
}

umap <- as.data.frame(Embeddings(tcell, "umap")) %>%
  rownames_to_column("cell")
colnames(umap)[2:3] <- c("UMAP_1", "UMAP_2")

module_score <- lapply(names(selected_modules), function(module_name) {
  genes <- intersect(selected_modules[[module_name]], rownames(tcell))
  if (length(genes) < 2) {
    warning("Module has fewer than two detected genes: ", module_name)
    return(NULL)
  }
  expr <- FetchData(tcell, vars = genes)
  tibble(
    cell = rownames(expr),
    module = module_name,
    n_genes = length(genes),
    module_score = rowMeans(as.matrix(expr), na.rm = TRUE)
  )
}) %>%
  bind_rows() %>%
  group_by(module) %>%
  mutate(module_score_z = zscore_clip(module_score)) %>%
  ungroup() %>%
  left_join(umap, by = "cell") %>%
  mutate(
    module = factor(module, levels = names(selected_modules)),
    draw_order = module_score_z
  ) %>%
  arrange(module, draw_order)

# Export the analysis source data for figure reproduction.
write.csv(module_score, file.path(result_dir, "lavage_Tcell_selected_module_score_umap_source_revision.csv"), row.names = FALSE)

message("Drawing selected module score UMAP...")
p_umap <- ggplot(module_score, aes(UMAP_1, UMAP_2, colour = module_score_z)) +
  geom_point(size = 0.34, alpha = 0.90, stroke = 0) +
  facet_wrap(~module, ncol = 2) +
  coord_fixed() +
  scale_colour_gradientn(
    colours = c("#D0D0D0", "#FEE08B", "#F46D43", "#B2182B"),
    limits = c(-1.2, 2.4),
    oob = scales::squish,
    name = "Module score\n(z-score)",
    guide = guide_colourbar(
      title.position = "top",
      barheight = unit(38, "mm"),
      barwidth = unit(6.2, "mm"),
      ticks.linewidth = 0.45,
      frame.linewidth = 0.35
    )
  ) +
  labs(title = "BALF T cells: selected functional module scores", x = "UMAP 1", y = "UMAP 2") +
  theme_umap(9.6) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11.6, face = "bold"),
    legend.text = element_text(size = 11.0),
    legend.key.height = unit(8, "mm"),
    legend.key.width = unit(6, "mm"),
    axis.text = element_text(size = 8.4),
    axis.title = element_text(size = 10.5, face = "bold")
  )

save_plot_safe(p_umap, "Tcell_selected_functional_module_score_UMAP_revision", width = 6.8, height = 6.1)
message("Saved Tcell_selected_functional_module_score_UMAP_revision.pdf and .png")

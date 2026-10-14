# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Final BALF figures from saved objects and signatures.
# Original source: 19_revision_lavage_tcell_subtype_deg_module_heatmap.R
# Figure/table scope: Fig3b; Fig3c; Fig3d; Fig3e; Fig3f; Fig3g; Fig3h; Fig3i; FigS6a; FigS6b; FigS6c; FigS6d; FigS6e; FigS6f; FigS6g; FigS6h
# Data cohort: 肺泡灌洗液. Inputs retain the original Data/Result hierarchy.
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
  library(pheatmap)
  library(grid)
})

project_dir <- Sys.getenv("LAVAGE_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "\u80ba\u6ce1\u704c\u6d17\u6db2")
}
project_dir <- gsub("\\\\", "/", project_dir)
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

tcell_type_levels <- c(
  "Naive CD4", "Th1-like", "Tfh", "Th17", "Treg", "GZMK+ CD8",
  "ITGAE+ CD8", "NKT", "ISG+ T", "Cycling T", "gdT", "MAIT"
)
global_colors <- c(
  "#E42320", "#699ECA", "#4DAF4A", "#FF8C00", "#F898CB", "#D65190",
  "#731A73", "#FFCB5B", "#E89B1E", "#0076B9", "#3D505A", "#0098B2",
  "#8A7355", "#CC5B45", "#008B8B", "#4E5689"
)
tcell_type_colors <- setNames(global_colors[seq_along(tcell_type_levels)], tcell_type_levels)
module_colors <- c(
  "Cytotoxic" = "#C74546",
  "Exhaustion" = "#8A7355",
  "IFN response" = "#0076B9",
  "Chemokine/cytokine" = "#E89B1E",
  "Activation" = "#D65190",
  "Naive/memory" = "#699ECA",
  "Tfh" = "#4DAF4A",
  "Treg" = "#731A73",
  "Tissue resident" = "#008B8B",
  "Cycling" = "#FFCB5B",
  "gdT" = "#3D505A",
  "Antigen presentation" = "#4E5689"
)
gene_modules <- list(
  "Cytotoxic" = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "KLRD1", "KLRG1", "CTSW", "CST7", "FGFBP2", "CCL5", "CD8A", "CD8B"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX", "ENTPD1", "LAYN", "CXCL13", "BTLA", "TNFRSF9", "MAF", "CD200"),
  "IFN response" = c("ISG15", "IFI6", "IFI27", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "RSAD2", "STAT1", "IRF7", "IRF8", "GBP1", "CXCL10"),
  "Chemokine/cytokine" = c("CXCR3", "CXCR4", "CCR7", "CCR5", "CCR6", "CXCR6", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2", "IFNG", "TNF", "LTB", "IL7R", "IL32", "GZMK"),
  "Activation" = c("IL2RA", "CD69", "ICOS", "TNFRSF9", "CD40LG", "FOS", "JUN", "JUNB", "DUSP1", "DUSP2", "NFKBIA", "CD69", "HLA-DRA", "HLA-DQA1"),
  "Naive/memory" = c("CCR7", "IL7R", "SELL", "TCF7", "LEF1", "MAL", "LTB", "S100A4", "AQP3", "PLAC8", "NOSIP", "LDHB", "TXK"),
  "Tfh" = c("CXCR5", "BCL6", "SH2D1A", "TOX2", "IL21", "ICOS", "PDCD1", "MAF", "CXCL13", "CD40LG"),
  "Treg" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "TNFRSF18", "RTKN2", "LAYN", "TNFRSF4"),
  "Tissue resident" = c("ITGAE", "ITGA1", "CXCR6", "ZNF683", "CD69", "RGS1", "XCL1", "XCL2", "KLRB1", "KLRG1"),
  "Cycling" = c("MKI67", "TOP2A", "STMN1", "TUBA1B", "TYMS", "HMGB2", "UBE2C", "MCM5", "PCNA"),
  "gdT" = c("TRDC", "TRGC1", "TRDV2", "TRGV9", "KLRD1", "NKG7", "GZMB"),
  "Antigen presentation" = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "HLA-A", "HLA-B", "HLA-C")
)
exclude_genes <- c("TRGC2", "TRDV3", "TRDV1")

assign_module <- function(genes) {
  vapply(genes, function(g) {
    hits <- names(gene_modules)[vapply(gene_modules, function(x) g %in% x, logical(1))]
    if (length(hits) > 0) hits[1] else "Other"
  }, character(1))
}

zscore_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  pmax(pmin(z, 2.5), -2.5)
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

theme_heatmap <- function(base_size = 9.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      strip.text.y = element_text(face = "bold", angle = 0, hjust = 0, colour = "black"),
      strip.background = element_rect(fill = "grey93", colour = NA),
      panel.grid = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(4, 5, 4, 5)
    )
}

message("Loading lavage_Tcell.rds...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(file.path(result_dir, "lavage_Tcell.rds"))
DefaultAssay(tcell) <- if ("RNA" %in% Assays(tcell)) "RNA" else DefaultAssay(tcell)
if (!"revised_celltype" %in% colnames(tcell@meta.data)) {
  stop("lavage_Tcell.rds does not contain revised_celltype.")
}

tcell$T_celltype <- factor(as.character(tcell$revised_celltype), levels = tcell_type_levels)
tcell <- subset(tcell, subset = !is.na(T_celltype))
present_types <- levels(droplevels(tcell$T_celltype))
Idents(tcell) <- "T_celltype"

deg_path <- file.path(result_dir, "lavage_Tcell_subtype_DEG_wilcox_revision.csv")
if (file.exists(deg_path)) {
  message("Reusing existing T-cell subtype DEG table...")
  deg_all <- read.csv(deg_path, stringsAsFactors = FALSE, check.names = FALSE) %>%
    as_tibble() %>%
    mutate(
      T_celltype_up = factor(as.character(T_celltype_up), levels = present_types),
      module = assign_module(gene)
    ) %>%
    arrange(T_celltype_up, p_val_adj, desc(avg_log2FC))
} else {
  message("Calculating T-cell subtype-specific DEGs with Wilcoxon test...")
  deg_all <- FindAllMarkers(
    tcell,
    assay = DefaultAssay(tcell),
    test.use = "wilcox",
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.25,
    return.thresh = 1,
    verbose = FALSE
  ) %>%
    as_tibble() %>%
    rename(T_celltype_up = cluster) %>%
    mutate(
      T_celltype_up = factor(as.character(T_celltype_up), levels = present_types),
      module = assign_module(gene)
    ) %>%
    arrange(T_celltype_up, p_val_adj, desc(avg_log2FC))
# Export the analysis source data for figure reproduction.
  write.csv(deg_all, deg_path, row.names = FALSE)
}

deg_sig <- deg_all %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.05, avg_log2FC > 0.25, module != "Other")

typical_gene_df <- tibble(
  module = rep(names(gene_modules), lengths(gene_modules)),
  gene = unlist(gene_modules, use.names = FALSE)
) %>%
  filter(gene %in% rownames(tcell)) %>%
  filter(!gene %in% exclude_genes) %>%
  distinct(module, gene) %>%
  left_join(
    deg_all %>%
      filter(module != "Other") %>%
      group_by(gene, module) %>%
      arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(gene, module, T_celltype_up, p_val_adj, avg_log2FC, pct.1, pct.2),
    by = c("gene", "module")
  ) %>%
  mutate(
    T_celltype_up = ifelse(is.na(T_celltype_up), "Module marker", as.character(T_celltype_up)),
    p_val_adj = ifelse(is.na(p_val_adj), NA_real_, p_val_adj),
    avg_log2FC = ifelse(is.na(avg_log2FC), NA_real_, avg_log2FC),
    source = ifelse(is.na(p_val_adj), "canonical_marker", "subtype_DEG")
  )

deg_selected <- deg_sig %>%
  group_by(T_celltype_up, module) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup() %>%
  group_by(module) %>%
  slice_head(n = 8) %>%
  ungroup()

selected_genes <- bind_rows(typical_gene_df, deg_selected) %>%
  filter(module != "Other", gene %in% rownames(tcell)) %>%
  filter(!gene %in% exclude_genes) %>%
  arrange(match(module, names(module_colors)), p_val_adj, desc(avg_log2FC), gene) %>%
  distinct(gene, .keep_all = TRUE) %>%
  group_by(module) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(module = factor(module, levels = names(module_colors))) %>%
  arrange(module, gene)

if (nrow(selected_genes) < 30) {
  selected_genes <- bind_rows(selected_genes, deg_sig) %>%
    filter(module != "Other", gene %in% rownames(tcell)) %>%
    filter(!gene %in% exclude_genes) %>%
    group_by(module) %>%
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 8) %>%
    ungroup() %>%
    arrange(match(module, names(module_colors)), p_val_adj, desc(avg_log2FC), gene) %>%
    distinct(gene, .keep_all = TRUE) %>%
    mutate(module = factor(module, levels = names(module_colors))) %>%
    arrange(module, gene)
}
write.csv(selected_genes, file.path(result_dir, "lavage_Tcell_subtype_DEG_module_heatmap_selected_genes_revision.csv"), row.names = FALSE)

message("Building subtype-average expression matrix...")
genes_to_plot <- intersect(selected_genes$gene, rownames(tcell))
avg <- AverageExpression(
  tcell,
  assays = DefaultAssay(tcell),
  group.by = "T_celltype",
  slot = "data",
  return.seurat = FALSE,
  verbose = FALSE
)[[DefaultAssay(tcell)]]
avg <- as.matrix(avg)
avg <- avg[genes_to_plot, present_types, drop = FALSE]
mat_z <- zscore_rows(avg)

row_anno <- selected_genes %>%
  filter(gene %in% rownames(mat_z)) %>%
  distinct(gene, module, T_celltype_up, .keep_all = TRUE) %>%
  transmute(gene, Module = factor(module, levels = names(module_colors)))

peak_subtype <- as.data.frame(mat_z) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = all_of(present_types), names_to = "peak_T_celltype", values_to = "z_score") %>%
  group_by(gene) %>%
  arrange(desc(z_score), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    peak_T_celltype = factor(peak_T_celltype, levels = present_types),
    peak_rank = match(as.character(peak_T_celltype), present_types)
  ) %>%
  select(gene, peak_T_celltype, peak_rank, peak_z_score = z_score)

row_anno <- row_anno %>%
  left_join(peak_subtype, by = "gene") %>%
  arrange(Module, peak_rank, desc(peak_z_score), gene)

module_order <- row_anno %>%
  group_by(Module) %>%
  summarise(
    median_peak_rank = median(peak_rank, na.rm = TRUE),
    min_peak_rank = min(peak_rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_peak_rank, min_peak_rank, match(as.character(Module), names(module_colors))) %>%
  pull(Module) %>%
  as.character()

row_anno <- row_anno %>%
  mutate(Module = factor(as.character(Module), levels = module_order)) %>%
  arrange(Module, peak_rank, desc(peak_z_score), gene) %>%
  column_to_rownames("gene")
write.csv(
  row_anno %>% rownames_to_column("gene"),
  file.path(result_dir, "lavage_Tcell_subtype_DEG_module_heatmap_gene_peak_order_revision.csv"),
  row.names = FALSE
)
mat_z <- mat_z[rownames(row_anno), , drop = FALSE]
col_anno <- data.frame(`T cell subtype` = factor(colnames(mat_z), levels = present_types), row.names = colnames(mat_z), check.names = FALSE)

row_gaps <- which(row_anno$Module[-1] != row_anno$Module[-nrow(row_anno)])
annotation_colors <- list(
  `T cell subtype` = tcell_type_colors[present_types],
  Module = module_colors[module_order]
)
heat_colors <- colorRampPalette(c("#376A9E", "#F7F7F7", "#B7192B"))(101)

message("Drawing module-faceted T-cell subtype DEG heatmap...")
module_gene_levels <- row_anno %>%
  rownames_to_column("gene") %>%
  mutate(
    Module = factor(Module, levels = module_order),
    module_gene = paste(Module, gene, sep = "__")
  ) %>%
  arrange(Module, peak_rank, desc(peak_z_score), gene) %>%
  pull(module_gene)

heat_df <- as.data.frame(mat_z) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = all_of(present_types), names_to = "T_celltype", values_to = "z_score") %>%
  left_join(row_anno %>% rownames_to_column("gene"), by = "gene") %>%
  mutate(
    Module = factor(Module, levels = module_order),
    module_gene = factor(paste(Module, gene, sep = "__"), levels = rev(module_gene_levels)),
    T_celltype = factor(T_celltype, levels = present_types)
  )
write.csv(heat_df, file.path(result_dir, "lavage_Tcell_subtype_DEG_module_heatmap_source_revision.csv"), row.names = FALSE)

p_heat <- ggplot(heat_df, aes(T_celltype, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(Module ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_gradient2(
    low = "#6F99C3",
    mid = "#F7F7F7",
    high = "#B7192B",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Z-score"
  ) +
  labs(title = "BALF T cells: subtype-specific DEGs") +
  theme_heatmap(9.2) +
  theme(
    axis.text.x = element_text(size = 8.0),
    axis.text.y = element_text(size = 6.8),
    strip.text.y = element_text(size = 8.1),
    legend.position = "right",
    legend.title = element_text(size = 10.2, face = "bold"),
    legend.text = element_text(size = 9.4)
  )

save_plot_safe(p_heat, "Tcell_subtype_DEG_module_heatmap_revision", width = 6.8, height = 9.2)

message("Saved Tcell_subtype_DEG_module_heatmap_revision.pdf and .png")

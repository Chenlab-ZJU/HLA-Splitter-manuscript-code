# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Final BALF figures from saved objects and signatures.
# Original source: 22_revision_lavage_tcell_CAP_LC_DEG_volcano.R
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
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

project_dir <- Sys.getenv("LAVAGE_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "\u80ba\u6ce1\u704c\u6d17\u6db2")
}
project_dir <- gsub("\\\\", "/", project_dir)
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures")
volcano_dir <- file.path(figure_dir, "TCell_CAP_LC_DEG")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

group_colors <- c("CAP" = "#699ECA", "LC" = "#E42320", "Not significant" = "#B8B8B8")
lfc_cutoff <- 0.25
padj_cutoff <- 0.05
y_axis_cap <- 120

save_plot_safe <- function(plot, path_no_ext, width, height) {
  tmp_pdf <- file.path(tempdir(), paste0(basename(path_no_ext), "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(basename(path_no_ext), "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
  file.copy(tmp_pdf, paste0(path_no_ext, ".pdf"), overwrite = TRUE)
  file.copy(tmp_png, paste0(path_no_ext, ".png"), overwrite = TRUE)
}

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.45),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.position = "right",
      plot.margin = margin(5, 5, 5, 5)
    )
}

message("Loading lavage_Tcell.rds...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(file.path(result_dir, "lavage_Tcell.rds"))
DefaultAssay(tcell) <- if ("RNA" %in% Assays(tcell)) "RNA" else DefaultAssay(tcell)
required_cols <- c("Group", "HLA_Splitter")
missing_cols <- setdiff(required_cols, colnames(tcell@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

tcell$Group <- factor(as.character(tcell$Group), levels = c("CAP", "LC", "ILD"))
tcell <- subset(tcell, subset = Group %in% c("CAP", "LC"))
Idents(tcell) <- "Group"

group_count <- tcell@meta.data %>%
  as_tibble() %>%
  count(Group, HLA_Splitter, name = "n_cells") %>%
  arrange(Group, HLA_Splitter)
# Export the analysis source data for figure reproduction.
write.csv(group_count, file.path(result_dir, "lavage_Tcell_CAP_LC_DEG_volcano_sample_cell_counts_revision.csv"), row.names = FALSE)

message("Running CAP vs LC differential expression...")
# Discover marker genes with the original assay and comparison settings.
deg <- FindMarkers(
  tcell,
  ident.1 = "CAP",
  ident.2 = "LC",
  test.use = "wilcox",
  only.pos = FALSE,
  min.pct = 0.10,
  logfc.threshold = 0,
  verbose = FALSE
) %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

write.csv(deg, file.path(result_dir, "lavage_Tcell_CAP_vs_LC_DEG_wilcox_revision.csv"), row.names = FALSE)

min_nonzero_padj <- suppressWarnings(min(deg$p_val_adj[deg$p_val_adj > 0], na.rm = TRUE))
if (!is.finite(min_nonzero_padj)) {
  min_nonzero_padj <- .Machine$double.xmin
}

plot_df <- deg %>%
  mutate(
    p_val_adj_plot = ifelse(is.na(p_val_adj), 1, ifelse(p_val_adj <= 0, min_nonzero_padj * 0.1, p_val_adj)),
    neg_log10_padj_raw = -log10(p_val_adj_plot),
    neg_log10_padj = pmin(neg_log10_padj_raw, y_axis_cap),
    direction = case_when(
      p_val_adj < padj_cutoff & avg_log2FC >= lfc_cutoff ~ "CAP",
      p_val_adj < padj_cutoff & avg_log2FC <= -lfc_cutoff ~ "LC",
      TRUE ~ "Not significant"
    ),
    direction = factor(direction, levels = c("LC", "Not significant", "CAP"))
  )

label_exclude_pattern <- "^(RPS|RPL|MT-|AC[0-9]|AL[0-9]|LINC|MTRNR|TMSB|TOMM|EEF)"
label_genes <- plot_df %>%
  filter(
    direction != "Not significant",
    abs(avg_log2FC) >= 0.50,
    !grepl(label_exclude_pattern, gene)
  ) %>%
  group_by(direction) %>%
  arrange(desc(abs(avg_log2FC)), p_val_adj_plot, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

key_genes <- c("IFNG", "GZMB", "PRF1", "NKG7", "CXCR3", "GZMK", "HLA-DRA", "HLA-DRB1", "HLA-DPB1", "LGALS1")
label_genes <- bind_rows(
  label_genes,
  plot_df %>% filter(gene %in% key_genes, direction != "Not significant")
) %>%
  distinct(gene, .keep_all = TRUE)

summary_tbl <- plot_df %>%
  count(direction, name = "n_genes") %>%
  mutate(direction = as.character(direction))
write.csv(summary_tbl, file.path(result_dir, "lavage_Tcell_CAP_vs_LC_DEG_volcano_summary_revision.csv"), row.names = FALSE)

p_volcano <- ggplot(plot_df, aes(avg_log2FC, neg_log10_padj)) +
  geom_point(aes(colour = direction), size = 1.25, alpha = 0.70, stroke = 0) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = gene),
    family = "Arial",
    colour = "black",
    size = 3.2,
    box.padding = unit(0.24, "lines"),
    point.padding = unit(0.18, "lines"),
    segment.color = "grey65",
    segment.size = 0.22,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, y_axis_cap),
    expand = expansion(mult = c(0.01, 0.05)),
    breaks = seq(0, y_axis_cap, by = 30)
  ) +
  labs(
    x = expression(Log[2]~fold~change~"(CAP / LC)"),
    y = expression(-Log[10]~adjusted~italic(P))
  ) +
  theme_pub(11.5) +
  theme(
    axis.title = element_text(size = 13.2, face = "bold"),
    axis.text = element_text(size = 11.2),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 9.8),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, -3, 0),
    plot.margin = margin(2, 4, 4, 4)
  ) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.2, alpha = 0.90)))

save_plot_safe(p_volcano, file.path(volcano_dir, "Tcell_CAP_vs_LC_DEG_volcano_revision"), width = 4.6, height = 4.4)
save_plot_safe(p_volcano, file.path(figure_dir, "Tcell_CAP_vs_LC_DEG_volcano_revision"), width = 4.6, height = 4.4)

message("Saved CAP vs LC DEG volcano plot and source tables.")

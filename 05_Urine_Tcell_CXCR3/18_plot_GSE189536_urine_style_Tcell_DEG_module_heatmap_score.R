#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: External Harmony annotation and final CXCR3/module validation.
# Original source: plot_GSE189536_urine_style_Tcell_DEG_module_heatmap_score.R
# Figure/table scope: FigS10a; FigS10b; FigS10c; FigS10d; FigS10e; FigS10f
# Data cohort: gse189536. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# GSE189536 biopsy validation: urine-style T-cell module heatmap and
# sample-level module score boxplots.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

env_data_dir <- Sys.getenv("GSE189536_DATA_DIR", unset = "")
if (nzchar(env_data_dir)) {
  data_dir <- sub("/$", "", gsub("\\\\", "/", env_data_dir))
} else {
  data_dir <- file.path(getOption("hla.data_root"), "gse189536")
}

annotation_dir <- file.path(data_dir, "Harmony_immune_cells_manual_annotation_results", "Urine_style_annotation")
input_rds <- file.path(getOption("hla.data_root"), "gse189536", "Immune_annotated.rds")
figure_dir <- file.path(annotation_dir, "Figures", "Tcell_Group_DEG_Function")
result_dir <- file.path(annotation_dir, "Tables", "Tcell_Group_DEG_Function")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "Mix" = "#FF8C00",
  "BK" = "#F898CB"
)

tcell_marker_modules <- list(
  "Naive/memory" = c("CCR7", "IL7R", "LEF1", "TCF7", "SELL", "MAL", "LTB"),
  "Activation/Treg" = c("CD40LG", "CCR6", "AQP3", "FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "Cytotoxic" = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TOX", "ICOS", "MAF"),
  "IFN/chemokine" = c("IFNG", "CXCR3", "CCL3", "CCL4", "ISG15", "IFI6", "IFIT3", "MX1"),
  "Antigen presentation" = c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DPA1", "CD74")
)
tcell_score_modules <- tcell_marker_modules

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

get_data_matrix <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "data")
  )
}

format_p <- function(p, prefix = "Wilcoxon p") {
  dplyr::case_when(
    is.na(p) ~ paste0(prefix, "=NA"),
    p < 0.001 ~ paste0(prefix, "<0.001"),
    p < 0.01 ~ paste0(prefix, "=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0(prefix, "=", formatC(p, format = "f", digits = 2))
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

safe_wilcox_two_groups <- function(df, value_col = "module_score", group_col = "Group", g1 = "TCMR", g2 = "ABMR") {
  vals <- df[[value_col]]
  grp <- as.character(df[[group_col]])
  keep <- is.finite(vals) & grp %in% c(g1, g2)
  vals <- vals[keep]
  grp <- factor(grp[keep], levels = c(g1, g2))
  if (nlevels(droplevels(grp)) < 2 || any(table(grp) < 2) || dplyr::n_distinct(vals) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(vals ~ grp, exact = FALSE)$p.value)
}

message("Loading urine-style annotated GSE189536 immune object...")
obj <- tryCatch(
  readRDS(input_rds),
  error = function(e) stop("Cannot read input object: ", input_rds, "\n", conditionMessage(e))
)
DefaultAssay(obj) <- "RNA"

meta_all <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Sample = as.character(biopsy),
    Group = as.character(Group),
    Urine_major_celltype = as.character(Urine_major_celltype),
    Urine_revised_celltype = as.character(Urine_revised_celltype),
    Urine_Tcell_revised_celltype = as.character(Urine_Tcell_revised_celltype)
  )

tcell_cells <- meta_all %>%
  filter(Urine_major_celltype %in% c("CD4+ T", "CD8+ T", "Cycling T"), Group %in% group_levels, !is.na(Sample), Sample != "") %>%
  pull(cell)
if (length(tcell_cells) < 100) stop("Too few urine-style T cells found in GSE189536 object.")

Tcell <- subset(obj, cells = tcell_cells)
expr <- get_data_matrix(Tcell, assay = "RNA")

meta <- Tcell@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Sample = as.character(biopsy),
    Group = as.character(Group),
    T_celltype = as.character(Urine_Tcell_revised_celltype),
    Urine_revised_celltype = as.character(Urine_revised_celltype)
  ) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    T_celltype = ifelse(is.na(T_celltype) | T_celltype == "", Urine_revised_celltype, T_celltype)
  ) %>%
  filter(!is.na(Group), !is.na(Sample), Sample != "")

available_groups <- intersect(group_levels, unique(as.character(meta$Group)))
if (length(available_groups) < 2) stop("Need at least two groups for GSE189536 T-cell validation.")

tcell_marker_df <- tibble(
  module = rep(names(tcell_marker_modules), lengths(tcell_marker_modules)),
  gene = unlist(tcell_marker_modules, use.names = FALSE)
) %>%
  filter(gene %in% rownames(Tcell)) %>%
  distinct(module, gene) %>%
  mutate(module = factor(module, levels = names(tcell_marker_modules))) %>%
  arrange(module, gene) %>%
  group_by(module) %>%
  mutate(module_gene = paste(module, gene, sep = "__")) %>%
  ungroup()
if (nrow(tcell_marker_df) == 0) stop("None of the selected T-cell module genes were found.")

message("Testing module genes across clinical groups...")
gene_group_tests <- bind_rows(lapply(tcell_marker_df$gene, function(gene) {
  values <- as.numeric(expr[gene, meta$cell, drop = TRUE])
  bind_rows(lapply(available_groups, function(grp) {
    in_group <- as.character(meta$Group) == grp
    if (sum(in_group) < 10 || sum(!in_group) < 10 || dplyr::n_distinct(values) < 2) {
      return(tibble(gene = gene, enriched_group = grp, avg_log2FC = NA_real_, p_value = NA_real_))
    }
    avg_in <- mean(values[in_group], na.rm = TRUE)
    avg_out <- mean(values[!in_group], na.rm = TRUE)
    tibble(
      gene = gene,
      enriched_group = grp,
      avg_log2FC = avg_in - avg_out,
      p_value = suppressWarnings(wilcox.test(values[in_group], values[!in_group], exact = FALSE)$p.value)
    )
  }))
})) %>%
  group_by(enriched_group) %>%
  mutate(p_val_adj_group = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  group_by(gene) %>%
  slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(p_val_adj = p.adjust(p_value, method = "BH"))

tcell_marker_df <- tcell_marker_df %>%
  left_join(gene_group_tests, by = "gene")
# Export the analysis source data for figure reproduction.
write_csv_safe(tcell_marker_df, file.path(result_dir, "GSE189536_Tcell_urine_style_heatmap_gene_modules.csv"))
write_csv_safe(gene_group_tests, file.path(result_dir, "GSE189536_Tcell_urine_style_module_gene_group_tests.csv"))

heatmap_genes <- unique(tcell_marker_df$gene)
group_avg <- AverageExpression(
  Tcell,
  assays = "RNA",
  features = heatmap_genes,
  group.by = "Group",
  slot = "data",
  verbose = FALSE
)$RNA %>%
  as.matrix()
group_avg <- group_avg[, intersect(group_levels, colnames(group_avg)), drop = FALSE]
group_z <- t(scale(t(group_avg)))
group_z[is.na(group_z)] <- 0
group_z <- pmax(pmin(group_z, 2.2), -2.2)

group_heat_df <- as.data.frame(group_z) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Group", values_to = "z_score") %>%
  left_join(tcell_marker_df, by = "gene", relationship = "many-to-many") %>%
  mutate(
    module = factor(as.character(module), levels = names(tcell_marker_modules)),
    module_gene = factor(module_gene, levels = rev(tcell_marker_df$module_gene)),
    Group = factor(Group, levels = group_levels)
  )

deg_gene_rank <- tcell_marker_df %>%
  mutate(
    rank_p = ifelse(is.na(p_val_adj), 1, p_val_adj),
    rank_fc = ifelse(is.na(avg_log2FC), 0, abs(avg_log2FC))
  ) %>%
  group_by(module) %>%
  arrange(rank_p, desc(rank_fc), gene, .by_group = TRUE) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  mutate(module_gene = paste(module, gene, sep = "__"))

deg_group_heat_df <- group_heat_df %>%
  select(-module_gene) %>%
  inner_join(deg_gene_rank %>% select(module, gene, module_gene), by = c("module", "gene")) %>%
  mutate(module_gene = factor(module_gene, levels = rev(deg_gene_rank$module_gene)))

write_csv_safe(group_heat_df, file.path(result_dir, "GSE189536_Tcell_group_key_gene_module_heatmap_source.csv"))
write_csv_safe(deg_group_heat_df, file.path(result_dir, "GSE189536_Tcell_group_DEG_gene_module_heatmap_source.csv"))

p_group_deg_module_heat <- ggplot(deg_group_heat_df, aes(Group, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_gradient2(
    low = "#6F99C3",
    mid = "#F7F7F7",
    high = "#B7192B",
    midpoint = 0,
    limits = c(-2.2, 2.2),
    oob = scales::squish,
    name = "Z-score"
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub(9.4) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.5, face = "bold"),
    axis.text.y = element_text(size = 8.0),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text.y = element_text(size = 8.0, angle = 0, hjust = 0, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )
save_plot_safe(p_group_deg_module_heat, "GSE189536_Tcell_group_DEG_gene_module_heatmap_revision", width = 4.9, height = 6.8)

message("Calculating sample-level T-cell module scores...")
module_score_cell <- lapply(names(tcell_score_modules), function(module_name) {
  genes <- intersect(tcell_score_modules[[module_name]], rownames(Tcell))
  if (length(genes) < 2) return(NULL)
  tibble(
    cell = colnames(expr),
    module = module_name,
    n_genes = length(genes),
    module_score = as.numeric(Matrix::colMeans(expr[genes, , drop = FALSE]))
  )
}) %>%
  bind_rows() %>%
  left_join(meta %>% select(cell, Sample, Group), by = "cell") %>%
  filter(!is.na(Sample), !is.na(Group))

module_score_sample <- module_score_cell %>%
  group_by(Sample, Group, module) %>%
  summarise(
    n_cells = n(),
    module_score = mean(module_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    module = factor(module, levels = names(tcell_score_modules))
  )

module_tcmr_abmr_wilcox <- module_score_sample %>%
  group_by(module) %>%
  summarise(
    group1 = "TCMR",
    group2 = "ABMR",
    p.value = safe_wilcox_two_groups(pick(everything()), "module_score", "Group", "TCMR", "ABMR"),
    .groups = "drop"
  ) %>%
  mutate(
    p.adjust = p.adjust(p.value, method = "BH"),
    label = format_p(p.value, prefix = "Wilcoxon p")
  )

write_csv_safe(module_score_cell, file.path(result_dir, "GSE189536_Tcell_functional_module_cell_scores.csv"))
write_csv_safe(module_score_sample, file.path(result_dir, "GSE189536_Tcell_functional_module_sample_scores.csv"))
write_csv_safe(module_tcmr_abmr_wilcox, file.path(result_dir, "GSE189536_Tcell_functional_module_TCMR_vs_ABMR_wilcox_tests.csv"))

wilcox_plot <- module_tcmr_abmr_wilcox %>%
  left_join(
    module_score_sample %>%
      group_by(module) %>%
      summarise(
        y_min = min(module_score, na.rm = TRUE),
        y_max = max(module_score, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "module"
  ) %>%
  mutate(
    y_span = pmax(y_max - y_min, 0.1),
    x_start = 2,
    x_end = 3,
    y_line = y_max + pmax(y_span * 0.08, 0.03),
    y_tick = y_line - pmax(y_span * 0.025, 0.01),
    y_text = y_line + pmax(y_span * 0.045, 0.015),
    module = factor(as.character(module), levels = names(tcell_score_modules))
  )

p_module <- ggplot(module_score_sample, aes(Group, module_score, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.36, alpha = 0.82) +
  geom_point(
    position = position_jitter(width = 0.11, height = 0, seed = 20260827),
    shape = 21,
    colour = "black",
    stroke = 0.18,
    size = 1.25,
    alpha = 0.88
  ) +
  geom_segment(data = wilcox_plot, aes(x = x_start, xend = x_end, y = y_line, yend = y_line), inherit.aes = FALSE, linewidth = 0.34, colour = "black", lineend = "round") +
  geom_segment(data = wilcox_plot, aes(x = x_start, xend = x_start, y = y_tick, yend = y_line), inherit.aes = FALSE, linewidth = 0.34, colour = "black", lineend = "round") +
  geom_segment(data = wilcox_plot, aes(x = x_end, xend = x_end, y = y_tick, yend = y_line), inherit.aes = FALSE, linewidth = 0.34, colour = "black", lineend = "round") +
  geom_text(
    data = wilcox_plot,
    aes(x = (x_start + x_end) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 2.25
  ) +
  facet_wrap(~module, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_colors[available_groups], drop = TRUE, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "Module score per biopsy") +
  theme_pub(8.8) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7.0),
    strip.text = element_text(size = 7.8),
    legend.position = "top",
    legend.text = element_text(size = 7.8),
    legend.title = element_text(size = 8.2, face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_blank()
  )
save_plot_safe(p_module, "GSE189536_Tcell_functional_module_score_compare_revision", width = 4.8, height = 6.0)

message("Saved GSE189536 urine-style T-cell DEG module heatmap and module score comparison.")

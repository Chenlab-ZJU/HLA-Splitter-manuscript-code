#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: External Harmony annotation and final CXCR3/module validation.
# Original source: run_GSE189536_revision_rejection_axis_validation.R
# Figure/table scope: FigS10a; FigS10b; FigS10c; FigS10d; FigS10e; FigS10f
# Data cohort: gse189536. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# Independent validation of urine rejection-associated findings in GSE189536.
# The analysis uses manually annotated immune cells and performs statistics at
# biopsy/sample level to avoid cell-number driven bias.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Matrix)
  library(stringr)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

get_data_matrix <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "data")
  )
}

score_gene_set <- function(expr, gene_set) {
  genes <- intersect(gene_set, rownames(expr))
  if (!length(genes)) return(rep(NA_real_, ncol(expr)))
  Matrix::colMeans(expr[genes, , drop = FALSE])
}

safe_kruskal <- function(values, groups) {
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- droplevels(factor(groups[keep]))
  if (length(unique(groups)) < 2 || any(table(groups) < 2)) return(NA_real_)
  kruskal.test(values ~ groups)$p.value
}

safe_wilcox_pair <- function(values, groups, g1, g2) {
  keep <- is.finite(values) & groups %in% c(g1, g2)
  values <- values[keep]
  groups <- droplevels(factor(groups[keep], levels = c(g1, g2)))
  if (length(unique(groups)) != 2 || any(table(groups) < 2)) return(NA_real_)
  wilcox.test(values ~ groups, exact = FALSE)$p.value
}

theme_nc <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(colour = "black", size = base_size),
      axis.title = element_text(colour = "black", size = base_size + 1),
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.40, colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1.5),
      plot.subtitle = element_text(hjust = 0, size = base_size, colour = "grey30"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.35, "cm"),
      strip.background = element_rect(fill = "grey96", colour = "black", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

theme_umap <- function(base_size = 9) {
  theme_void(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1.5),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.35, "cm"),
      strip.background = element_rect(fill = "grey96", colour = "black", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(4, 5, 4, 4)
    )
}

save_pub_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 600, bg = "white")
}

script_dir <- get_script_dir()
data_dir <- getOption("hla.cohort_dir")
manual_dir <- file.path(data_dir, "Harmony_immune_cells_manual_annotation_results")
input_rds <- file.path(getOption("hla.data_root"), "gse189536", "Immune_filtered.rds")
result_dir <- file.path(data_dir, "Revision_rejection_axis_validation")
figure_dir <- file.path(result_dir, "Figures")
table_dir <- file.path(result_dir, "Tables")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix")
group_colors <- c("NR" = "#E42320", "TCMR" = "#699ECA", "ABMR" = "#4DAF4A", "Mix" = "#FF8C00")
heat_colors <- c("#F7F7F7", "#FDBE85", "#FB6A4A", "#DE2D26", "#A50F15")
diverging_colors <- c("#2C7BB6", "#F7F7F7", "#B2182B")

tcell_types <- c("CD4 T", "CD8 T", "gamma-delta T", "Treg")
source_types <- c("Myeloid", "ISG+ immune", "B cell", "Plasma", "NK", "Neutrophil")
cytotoxic_genes <- c("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH", "GZMK", "CST7", "CTSW", "IFNG", "FASLG", "CCL5")
cxcr3_axis_genes <- c("CXCR3", "CXCL9", "CXCL10", "CXCL11", "IFNG")
needed_genes <- unique(c(cytotoxic_genes, cxcr3_axis_genes))

message("Loading manually annotated GSE189536 immune object...")
if (!file.exists(input_rds)) stop("Cannot find input Seurat object: ", input_rds)
# Read the saved object or source table used by this analysis.
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

if (!"manual_cell_type" %in% colnames(obj@meta.data)) stop("manual_cell_type metadata column was not found.")
if (!"sample_type" %in% colnames(obj@meta.data)) stop("sample_type metadata column was not found.")
if (!"biopsy" %in% colnames(obj@meta.data)) stop("biopsy metadata column was not found.")

valid_genes <- intersect(needed_genes, rownames(obj))
missing_genes <- setdiff(needed_genes, valid_genes)
message("Detected genes: ", paste(valid_genes, collapse = ", "))
if (length(missing_genes)) message("Missing genes skipped: ", paste(missing_genes, collapse = ", "))

expr <- get_data_matrix(obj, assay = "RNA")
meta <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    Group = case_when(
      sample_type %in% group_levels ~ as.character(sample_type),
      diagnosis == "No rejection" ~ "NR",
      diagnosis == "TCMR" ~ "TCMR",
      diagnosis == "ABMR" ~ "ABMR",
      diagnosis == "Mixed" ~ "Mix",
      TRUE ~ NA_character_
    ),
    Group = factor(Group, levels = group_levels),
    Biopsy = as.character(biopsy),
    Manual_celltype = as.character(manual_cell_type),
    is_T_cell = Manual_celltype %in% tcell_types,
    is_source_cell = Manual_celltype %in% source_types
  )

expr <- expr[, meta$cell, drop = FALSE]
for (gene in needed_genes) {
  meta[[paste0(gene, "_expr")]] <- if (gene %in% rownames(expr)) as.numeric(expr[gene, ]) else 0
}
meta$CXCR3_positive <- meta$is_T_cell & meta$CXCR3_expr > 0
meta$CXCL9_CXCL10_ligand_score <- rowMeans(meta[, c("CXCL9_expr", "CXCL10_expr"), drop = FALSE], na.rm = TRUE)
meta$CXCL9_CXCL10_CXCL11_ligand_score <- rowMeans(meta[, c("CXCL9_expr", "CXCL10_expr", "CXCL11_expr"), drop = FALSE], na.rm = TRUE)
meta$cytotoxic_score <- score_gene_set(expr, cytotoxic_genes)

# Export the analysis source data for figure reproduction.
write_csv(
  meta %>%
    count(Group, diagnosis, Biopsy, library_id, fraction, Manual_celltype, name = "n_cells") %>%
    arrange(Group, Biopsy, desc(n_cells)),
  file.path(table_dir, "GSE189536_cell_counts_by_sample_celltype.csv")
)

message("Summarising T-cell CXCR3 and cytotoxic signatures at sample level...")
tcell_sample <- meta %>%
  filter(is_T_cell, !is.na(Group)) %>%
  group_by(Group, diagnosis, Biopsy) %>%
  summarise(
    n_T_cells = n(),
    CXCR3_positive_T_cells = sum(CXCR3_positive, na.rm = TRUE),
    CXCR3_positive_T_fraction = CXCR3_positive_T_cells / n_T_cells,
    mean_CXCR3_T = mean(CXCR3_expr, na.rm = TRUE),
    median_CXCR3_T = median(CXCR3_expr, na.rm = TRUE),
    mean_IFNG_T = mean(IFNG_expr, na.rm = TRUE),
    mean_cytotoxic_score_T = mean(cytotoxic_score, na.rm = TRUE),
    median_cytotoxic_score_T = median(cytotoxic_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(CXCR3_positive_T_percent = 100 * CXCR3_positive_T_fraction)

tcell_celltype_sample <- meta %>%
  filter(is_T_cell, !is.na(Group)) %>%
  group_by(Group, diagnosis, Biopsy, Manual_celltype) %>%
  summarise(
    n_T_cells = n(),
    CXCR3_positive_T_fraction = mean(CXCR3_expr > 0, na.rm = TRUE),
    mean_CXCR3_T = mean(CXCR3_expr, na.rm = TRUE),
    mean_cytotoxic_score_T = mean(cytotoxic_score, na.rm = TRUE),
    mean_IFNG_T = mean(IFNG_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(CXCR3_positive_T_percent = 100 * CXCR3_positive_T_fraction)

cytotoxic_gene_sample <- meta %>%
  filter(is_T_cell, !is.na(Group)) %>%
  select(Group, diagnosis, Biopsy, Manual_celltype, ends_with("_expr")) %>%
  pivot_longer(cols = all_of(paste0(intersect(cytotoxic_genes, rownames(expr)), "_expr")),
               names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = sub("_expr$", "", Gene)) %>%
  group_by(Group, diagnosis, Biopsy, Gene) %>%
  summarise(
    n_T_cells = n(),
    mean_expression = mean(Expression, na.rm = TRUE),
    pct_positive = mean(Expression > 0, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(tcell_sample, file.path(table_dir, "GSE189536_Tcell_CXCR3_cytotoxic_sample_summary.csv"))
write_csv(tcell_celltype_sample, file.path(table_dir, "GSE189536_Tcell_CXCR3_cytotoxic_by_Tcelltype.csv"))
write_csv(cytotoxic_gene_sample, file.path(table_dir, "GSE189536_Tcell_cytotoxic_gene_sample_expression.csv"))

message("Estimating CXCR3-CXCL9/CXCL10 ligand-receptor axis strength...")
source_sample <- meta %>%
  filter(is_source_cell, !is.na(Group)) %>%
  group_by(Group, diagnosis, Biopsy, Source_celltype = Manual_celltype) %>%
  summarise(
    n_source_cells = n(),
    mean_CXCL9 = mean(CXCL9_expr, na.rm = TRUE),
    mean_CXCL10 = mean(CXCL10_expr, na.rm = TRUE),
    mean_CXCL11 = mean(CXCL11_expr, na.rm = TRUE),
    pct_CXCL9_positive = mean(CXCL9_expr > 0, na.rm = TRUE),
    pct_CXCL10_positive = mean(CXCL10_expr > 0, na.rm = TRUE),
    pct_CXCL11_positive = mean(CXCL11_expr > 0, na.rm = TRUE),
    mean_CXCL9_CXCL10 = mean(CXCL9_CXCL10_ligand_score, na.rm = TRUE),
    mean_CXCL9_CXCL10_CXCL11 = mean(CXCL9_CXCL10_CXCL11_ligand_score, na.rm = TRUE),
    .groups = "drop"
  )

target_sample <- meta %>%
  filter(is_T_cell, !is.na(Group)) %>%
  group_by(Group, diagnosis, Biopsy, Target_Tcelltype = Manual_celltype) %>%
  summarise(
    n_target_T_cells = n(),
    mean_CXCR3 = mean(CXCR3_expr, na.rm = TRUE),
    pct_CXCR3_positive = mean(CXCR3_expr > 0, na.rm = TRUE),
    .groups = "drop"
  )

axis_pair_sample <- source_sample %>%
  inner_join(target_sample, by = c("Group", "diagnosis", "Biopsy"), relationship = "many-to-many") %>%
  mutate(
    CXCL9_CXCR3_score = mean_CXCL9 * mean_CXCR3 * sqrt(pct_CXCL9_positive * pct_CXCR3_positive),
    CXCL10_CXCR3_score = mean_CXCL10 * mean_CXCR3 * sqrt(pct_CXCL10_positive * pct_CXCR3_positive),
    CXCL9_CXCL10_CXCR3_score = mean_CXCL9_CXCL10 * mean_CXCR3 *
      sqrt(pmax((pct_CXCL9_positive + pct_CXCL10_positive) / 2, 0) * pct_CXCR3_positive),
    CXCL9_CXCL10_CXCL11_CXCR3_score = mean_CXCL9_CXCL10_CXCL11 * mean_CXCR3 *
      sqrt(pmax((pct_CXCL9_positive + pct_CXCL10_positive + pct_CXCL11_positive) / 3, 0) * pct_CXCR3_positive)
  )

axis_sample_overall <- axis_pair_sample %>%
  group_by(Group, diagnosis, Biopsy) %>%
  summarise(
    max_CXCL9_CXCR3_score = max(CXCL9_CXCR3_score, na.rm = TRUE),
    max_CXCL10_CXCR3_score = max(CXCL10_CXCR3_score, na.rm = TRUE),
    max_CXCL9_CXCL10_CXCR3_score = max(CXCL9_CXCL10_CXCR3_score, na.rm = TRUE),
    mean_CXCL9_CXCL10_CXCR3_score = mean(CXCL9_CXCL10_CXCR3_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x) | is.nan(.x), NA_real_, .x)))

write_csv(source_sample, file.path(table_dir, "GSE189536_source_cell_CXCL9_CXCL10_sample_summary.csv"))
write_csv(target_sample, file.path(table_dir, "GSE189536_target_Tcell_CXCR3_sample_summary.csv"))
write_csv(axis_pair_sample, file.path(table_dir, "GSE189536_CXCR3_CXCL9_CXCL10_source_target_LR_scores.csv"))
write_csv(axis_sample_overall, file.path(table_dir, "GSE189536_CXCR3_CXCL9_CXCL10_axis_sample_summary.csv"))

message("Running sample-level statistics...")
core_metrics <- bind_rows(
  tcell_sample %>%
    select(Group, diagnosis, Biopsy, CXCR3_positive_T_percent, mean_CXCR3_T, mean_cytotoxic_score_T) %>%
    pivot_longer(cols = c(CXCR3_positive_T_percent, mean_CXCR3_T, mean_cytotoxic_score_T),
                 names_to = "Metric", values_to = "Value"),
  axis_sample_overall %>%
    select(Group, diagnosis, Biopsy, max_CXCL9_CXCR3_score, max_CXCL10_CXCR3_score,
           max_CXCL9_CXCL10_CXCR3_score, mean_CXCL9_CXCL10_CXCR3_score) %>%
    pivot_longer(cols = starts_with(c("max_", "mean_")), names_to = "Metric", values_to = "Value")
)

metric_tests <- core_metrics %>%
  group_by(Metric) %>%
  summarise(
    n_samples = n_distinct(Biopsy),
    p_kruskal = safe_kruskal(Value, Group),
    p_ABMR_vs_NR = safe_wilcox_pair(Value, Group, "NR", "ABMR"),
    p_TCMR_vs_NR = safe_wilcox_pair(Value, Group, "NR", "TCMR"),
    p_ABMR_vs_TCMR = safe_wilcox_pair(Value, Group, "TCMR", "ABMR"),
    median_NR = median(Value[Group == "NR"], na.rm = TRUE),
    median_TCMR = median(Value[Group == "TCMR"], na.rm = TRUE),
    median_ABMR = median(Value[Group == "ABMR"], na.rm = TRUE),
    median_Mix = median(Value[Group == "Mix"], na.rm = TRUE),
    ABMR_minus_NR = median_ABMR - median_NR,
    ABMR_minus_TCMR = median_ABMR - median_TCMR,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_kruskal = p.adjust(p_kruskal, method = "BH"),
    p_adj_ABMR_vs_NR = p.adjust(p_ABMR_vs_NR, method = "BH"),
    p_adj_TCMR_vs_NR = p.adjust(p_TCMR_vs_NR, method = "BH"),
    p_adj_ABMR_vs_TCMR = p.adjust(p_ABMR_vs_TCMR, method = "BH")
  ) %>%
  arrange(p_adj_kruskal, desc(abs(ABMR_minus_NR)))
write_csv(core_metrics, file.path(table_dir, "GSE189536_revision_core_metrics_long.csv"))
write_csv(metric_tests, file.path(table_dir, "GSE189536_revision_sample_level_statistical_tests.csv"))

cytotoxic_gene_tests <- cytotoxic_gene_sample %>%
  group_by(Gene) %>%
  summarise(
    median_NR = median(mean_expression[Group == "NR"], na.rm = TRUE),
    median_TCMR = median(mean_expression[Group == "TCMR"], na.rm = TRUE),
    median_ABMR = median(mean_expression[Group == "ABMR"], na.rm = TRUE),
    ABMR_minus_NR = median_ABMR - median_NR,
    ABMR_minus_TCMR = median_ABMR - median_TCMR,
    p_ABMR_vs_NR = safe_wilcox_pair(mean_expression, Group, "NR", "ABMR"),
    p_ABMR_vs_TCMR = safe_wilcox_pair(mean_expression, Group, "TCMR", "ABMR"),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_ABMR_vs_NR = p.adjust(p_ABMR_vs_NR, method = "BH"),
    p_adj_ABMR_vs_TCMR = p.adjust(p_ABMR_vs_TCMR, method = "BH")
  ) %>%
  arrange(desc(ABMR_minus_NR), p_ABMR_vs_NR)
write_csv(cytotoxic_gene_tests, file.path(table_dir, "GSE189536_ABMR_Tcell_cytotoxic_gene_tests.csv"))

message("Drawing publication-style validation plots...")
metric_labels <- c(
  CXCR3_positive_T_percent = "CXCR3+ cells in T cells (%)",
  mean_CXCR3_T = "Mean CXCR3 in T cells only",
  mean_cytotoxic_score_T = "Cytotoxic score in T cells only",
  max_CXCL9_CXCR3_score = "Max CXCL9-CXCR3 score",
  max_CXCL10_CXCR3_score = "Max CXCL10-CXCR3 score",
  max_CXCL9_CXCL10_CXCR3_score = "Max CXCL9/10-CXCR3 score\nsource cells to T cells",
  mean_CXCL9_CXCL10_CXCR3_score = "Mean CXCL9/10-CXCR3 score\nsource cells to T cells"
)

p_core <- core_metrics %>%
  filter(Metric %in% c("CXCR3_positive_T_percent", "mean_CXCR3_T", "mean_cytotoxic_score_T",
                       "max_CXCL9_CXCL10_CXCR3_score")) %>%
  mutate(Metric = factor(Metric, levels = c("CXCR3_positive_T_percent", "mean_CXCR3_T",
                                            "mean_cytotoxic_score_T", "max_CXCL9_CXCL10_CXCR3_score"),
                         labels = metric_labels[c("CXCR3_positive_T_percent", "mean_CXCR3_T",
                                                  "mean_cytotoxic_score_T", "max_CXCL9_CXCL10_CXCR3_score")])) %>%
  ggplot(aes(Group, Value, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.42, alpha = 0.78) +
  geom_point(size = 1.8, position = position_jitter(width = 0.11, height = 0),
             shape = 21, colour = "black", stroke = 0.25, alpha = 0.9) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    title = "Core rejection-axis validation metrics",
    subtitle = "Sample-level analysis; each panel states the cell compartment used.",
    x = NULL,
    y = "Sample-level value"
  ) +
  theme_nc(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
        legend.position = "none")

cyto_gene_plot <- cytotoxic_gene_tests %>%
  mutate(Gene = factor(Gene, levels = rev(Gene))) %>%
  ggplot(aes(ABMR_minus_NR, Gene)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = ABMR_minus_NR, yend = Gene), linewidth = 0.65, colour = "#4DAF4A") +
  geom_point(aes(size = -log10(pmax(p_ABMR_vs_NR, .Machine$double.xmin))), colour = "#4DAF4A") +
  scale_size_continuous(range = c(1.4, 3.4), name = expression(-log[10]("P"))) +
  labs(x = "Median expression difference (ABMR - NR)", y = NULL) +
  theme_nc(9)

axis_heat_data <- axis_pair_sample %>%
  group_by(Group, Source_celltype, Target_Tcelltype) %>%
  summarise(
    mean_LR_score = mean(CXCL9_CXCL10_CXCR3_score, na.rm = TRUE),
    n_samples = n_distinct(Biopsy),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Pair = paste(Source_celltype, Target_Tcelltype, sep = " -> ")
  ) %>%
  group_by(Pair) %>%
  mutate(max_pair_score = max(mean_LR_score, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(max_pair_score > 0) %>%
  arrange(desc(max_pair_score)) %>%
  mutate(Pair = factor(Pair, levels = rev(unique(Pair))))
write_csv(axis_heat_data, file.path(table_dir, "GSE189536_CXCR3_CXCL9_CXCL10_axis_group_source_target_heatmap.csv"))

p_axis_heat <- ggplot(axis_heat_data, aes(Group, Pair, fill = mean_LR_score)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = signif(mean_LR_score, 2)), size = 2.2, colour = "black") +
  scale_fill_gradientn(colours = heat_colors, name = "Mean\nLR score") +
  labs(x = NULL, y = NULL) +
  theme_nc(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 6.8),
        panel.grid = element_blank())

cyto_heat_data <- cytotoxic_gene_sample %>%
  group_by(Group, Gene) %>%
  summarise(mean_expression = mean(mean_expression, na.rm = TRUE), .groups = "drop") %>%
  group_by(Gene) %>%
  mutate(z = as.numeric(scale(mean_expression))) %>%
  ungroup() %>%
  mutate(
    z = ifelse(is.na(z), 0, z),
    Group = factor(Group, levels = group_levels),
    Gene = factor(Gene, levels = rev(cytotoxic_genes[cytotoxic_genes %in% unique(Gene)]))
  )

p_cyto_heat <- ggplot(cyto_heat_data, aes(Group, Gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = diverging_colors[1], mid = diverging_colors[2], high = diverging_colors[3],
                       midpoint = 0, name = "Scaled\nexpression") +
  labs(x = NULL, y = NULL) +
  theme_nc(8.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
        panel.grid = element_blank())

umap_df <- as.data.frame(Embeddings(obj, reduction = "umap")) %>%
  tibble::rownames_to_column("cell") %>%
  left_join(
    meta %>%
      select(cell, Group, Manual_celltype, is_T_cell, CXCR3_expr, CXCR3_positive, cytotoxic_score,
             CXCL9_CXCL10_ligand_score),
    by = "cell"
  ) %>%
  mutate(
    CXCR3_T_status = case_when(
      CXCR3_positive ~ "CXCR3+ T cells",
      is_T_cell ~ "CXCR3- T cells",
      Manual_celltype %in% source_types & CXCL9_CXCL10_ligand_score > 0 ~ "CXCL9/10+ source cells",
      TRUE ~ "Other immune cells"
    ),
    CXCR3_T_status = factor(CXCR3_T_status,
                            levels = c("Other immune cells", "CXCR3- T cells", "CXCL9/10+ source cells", "CXCR3+ T cells"))
  )
write_csv(umap_df, file.path(table_dir, "GSE189536_revision_validation_umap_metadata.csv"))

p_umap <- ggplot() +
  geom_point(data = filter(umap_df, CXCR3_T_status == "Other immune cells"),
             aes(UMAP_1, UMAP_2), colour = "#D9D9D9", size = 0.16, alpha = 0.35, stroke = 0) +
  geom_point(data = filter(umap_df, CXCR3_T_status == "CXCR3- T cells"),
             aes(UMAP_1, UMAP_2), colour = "#808080", size = 0.20, alpha = 0.55, stroke = 0) +
  geom_point(data = filter(umap_df, CXCR3_T_status == "CXCL9/10+ source cells"),
             aes(UMAP_1, UMAP_2), colour = "#FF8C00", size = 0.34, alpha = 0.80, stroke = 0) +
  geom_point(data = filter(umap_df, CXCR3_T_status == "CXCR3+ T cells"),
             aes(UMAP_1, UMAP_2), colour = "#B2182B", size = 0.40, alpha = 0.92, stroke = 0) +
  coord_equal() +
  labs(title = "CXCR3+ T cells and CXCL9/10+ source cells") +
  theme_umap(9)

combined <- (p_umap | p_core) / (p_cyto_heat | p_axis_heat) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_pub_plot(p_core, "GSE189536_revision_core_metric_boxplots", 7.2, 5.3)
save_pub_plot(cyto_gene_plot, "GSE189536_ABMR_Tcell_cytotoxic_gene_lollipop", 4.4, 4.0)
save_pub_plot(p_cyto_heat, "GSE189536_Tcell_cytotoxic_gene_group_heatmap", 3.8, 4.4)
save_pub_plot(p_axis_heat, "GSE189536_CXCR3_CXCL9_CXCL10_axis_source_target_heatmap", 7.8, 5.2)
save_pub_plot(p_umap, "GSE189536_CXCR3_Tcell_CXCL9_CXCL10_source_UMAP", 4.6, 4.2)
save_pub_plot(combined, "GSE189536_revision_rejection_axis_validation_combined", 12.4, 10.0)

summary_table <- tibble(
  validation_question = c(
    "T-cell CXCR3 expression",
    "CXCR3+ T-cell proportion",
    "ABMR T-cell cytotoxic gene program",
    "CXCR3-CXCL9/CXCL10 ligand-receptor axis"
  ),
  primary_metric = c(
    "mean_CXCR3_T",
    "CXCR3_positive_T_percent",
    "mean_cytotoxic_score_T and individual cytotoxic genes",
    "max_CXCL9_CXCL10_CXCR3_score"
  ),
  statistical_unit = "biopsy/sample",
  key_table = c(
    "GSE189536_Tcell_CXCR3_cytotoxic_sample_summary.csv",
    "GSE189536_Tcell_CXCR3_cytotoxic_sample_summary.csv",
    "GSE189536_ABMR_Tcell_cytotoxic_gene_tests.csv",
    "GSE189536_CXCR3_CXCL9_CXCL10_axis_sample_summary.csv"
  )
)
write_csv(summary_table, file.path(table_dir, "GSE189536_revision_validation_output_summary.csv"))

message("Done. Results: ", result_dir)
message("Figures: ", figure_dir)
message("Tables: ", table_dir)

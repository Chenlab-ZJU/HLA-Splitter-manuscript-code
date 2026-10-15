# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final robust signature construction, scoring and cross-cohort evaluation.
# Original source: 34_quick_cross_cohort_robust_signature_distribution_plots.R
# Figure/table scope: Fig4i; Fig4j; FigS8a; FigS8b; FigS8c; FigS8d; FigS8e
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  project_dir <- file.path(getOption("hla.data_root"), "urine")
}

figure_dir <- file.path(project_dir, "Figures", "Signature_score_analysis", "Cross_cohort_robust_selection")
result_dir <- file.path(project_dir, "Result", "Signature_score_analysis", "Cross_cohort_robust_selection")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
signature_levels <- c("TCMR", "ABMR", "BK")
signature_colors <- c(
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "BK" = "#F898CB"
)

curated_rejection_gene_sets <- list(
  TCMR = c(
    "CXCR3", "CXCL9", "CXCL10", "CXCL11", "CCL5", "IFNG", "TNF",
    "GZMK", "GZMB", "GZMA", "PRF1", "GNLY", "NKG7", "CST7",
    "CD3D", "CD3E", "CD8A", "CD8B", "TRAC", "STAT1", "GBP1",
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "TAP1", "TAP2"
  ),
  ABMR = c(
    "CXCL9", "CXCL10", "CXCL11", "CCL2", "CCL5", "IFNG", "TNF",
    "STAT1", "GBP1", "GBP5", "IRF1", "ISG15", "IFI44L", "MX1", "OAS1",
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "TAP1", "TAP2",
    "MS4A1", "CD79A", "MZB1", "JCHAIN", "IGHG1", "FCGR3A", "NKG7", "GZMB",
    "CXCR3", "SLC7A11", "CHI3L1", "IL6"
  ),
  BK = c(
    "CXCL9", "CXCL10", "CXCL11", "CCL5", "IFNG", "STAT1", "STAT2", "IRF7",
    "ISG15", "IFI6", "IFI27", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "MX1",
    "OAS1", "OAS2", "OAS3", "RSAD2", "HLA-DRA", "CD74", "TAP1", "TAP2",
    "NKG7", "GZMB", "PRF1"
  )
)

curated_gene_modules <- list(
  chemokine_CXCR3_axis = c("CXCR3", "CXCL9", "CXCL10", "CXCL11", "CCL2", "CCL5"),
  T_NK_cytotoxicity = c("GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7", "CST7", "FCGR3A"),
  IFN_response = c("IFNG", "TNF", "STAT1", "STAT2", "GBP1", "GBP5", "IRF1", "IRF7", "ISG15", "IFI6", "IFI27", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2", "OAS3", "RSAD2"),
  antigen_presentation_MHC = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "TAP1", "TAP2"),
  T_cell_identity = c("CD3D", "CD3E", "CD8A", "CD8B", "TRAC"),
  B_plasma_response = c("MS4A1", "CD79A", "MZB1", "JCHAIN", "IGHG1"),
  inflammation_remodeling = c("SLC7A11", "CHI3L1", "IL6")
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

theme_pub <- function(base_size = 9.2) {
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

get_assay_matrix <- function(obj, assay = "RNA") {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    return(GetAssayData(obj, assay = assay, layer = "data"))
  }
  GetAssayData(obj, assay = assay, slot = "data")
}

score_auc <- function(expr, gene_sets) {
  rankings <- AUCell_buildRankings(expr, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE)
  auc <- AUCell_calcAUC(gene_sets, rankings, aucMaxRank = ceiling(nrow(rankings) * 0.10), verbose = FALSE)
  as.data.frame(t(getAUC(auc))) %>%
    rownames_to_column("cell")
}

celltype_contribution <- function(cell_score_long, type_col, target_group_only = TRUE) {
  group_mean <- cell_score_long %>%
    group_by(signature, celltype = .data[[type_col]], Group) %>%
    summarise(group_mean_auc = mean(auc_score, na.rm = TRUE), .groups = "drop")

  target_mean <- group_mean %>%
    filter(as.character(Group) == as.character(signature)) %>%
    select(signature, celltype, target_mean_auc = group_mean_auc)

  other_mean <- group_mean %>%
    filter(as.character(Group) != as.character(signature)) %>%
    group_by(signature, celltype) %>%
    summarise(other_mean_auc = mean(group_mean_auc, na.rm = TRUE), .groups = "drop")

  base_df <- if (target_group_only) {
    cell_score_long %>% filter(as.character(Group) == as.character(signature))
  } else {
    cell_score_long
  }

  base_df %>%
    group_by(signature, celltype = .data[[type_col]]) %>%
    summarise(
      n_cells = n(),
      mean_auc = mean(auc_score, na.rm = TRUE),
      total_auc = sum(auc_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(target_mean, by = c("signature", "celltype")) %>%
    left_join(other_mean, by = c("signature", "celltype")) %>%
    group_by(signature) %>%
    mutate(
      contribution_pct = 100 * total_auc / sum(total_auc, na.rm = TRUE),
      target_enrichment = target_mean_auc - other_mean_auc,
      target_ratio = (target_mean_auc + 1e-4) / (other_mean_auc + 1e-4),
      contribution_range = max(contribution_pct, na.rm = TRUE) - min(contribution_pct, na.rm = TRUE),
      enrichment_range = max(target_enrichment, na.rm = TRUE) - min(target_enrichment, na.rm = TRUE),
      scaled_contribution = (contribution_pct - min(contribution_pct, na.rm = TRUE)) / pmax(contribution_range, .Machine$double.eps),
      scaled_enrichment = (target_enrichment - min(target_enrichment, na.rm = TRUE)) / pmax(enrichment_range, .Machine$double.eps),
      relative_target_enrichment = scaled_enrichment,
      celltype_specificity_score = scaled_contribution + scaled_enrichment,
      rank_by_contribution = rank(-celltype_specificity_score, ties.method = "first")
    ) %>%
    ungroup() %>%
    select(-contribution_range, -enrichment_range, -scaled_contribution, -scaled_enrichment) %>%
    filter(!is.na(celltype), celltype != "")
}

plot_contribution <- function(df, filename, width, height, n_top = NULL) {
  plot_df <- df
  if (!is.null(n_top)) {
    plot_df <- df %>%
      group_by(signature) %>%
      slice_min(rank_by_contribution, n = n_top, with_ties = FALSE) %>%
      ungroup()
  }
  type_levels <- plot_df %>%
    group_by(celltype) %>%
    summarise(max_contribution = max(contribution_pct, na.rm = TRUE), .groups = "drop") %>%
    arrange(max_contribution) %>%
    pull(celltype)
  plot_df <- plot_df %>%
    mutate(
      celltype = factor(celltype, levels = type_levels),
      signature = factor(signature, levels = signature_levels)
    )

  p <- ggplot(plot_df, aes(x = signature, y = celltype)) +
    geom_point(aes(size = contribution_pct, fill = relative_target_enrichment), shape = 21, colour = "grey25", stroke = 0.22, alpha = 0.94) +
    scale_size_continuous(range = c(1.5, 8.5), name = "Contribution (%)") +
    scale_fill_gradient2(
      low = "#3F86C5",
      mid = "#F7F7F7",
      high = "#9E0142",
      midpoint = 0.5,
      limits = c(0, 1),
      name = "Relative target\nenrichment"
    ) +
    labs(x = "Signature", y = NULL) +
    theme_pub(9.2) +
    theme(
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right"
    )
  save_plot_safe(p, filename, width = width, height = height)
  invisible(p)
}

read_signature_deg_table <- function(project_dir) {
  deg_dir <- file.path(project_dir, "Result", "Signature_score_analysis", "Differential_Expression")
  files <- c(
    TCMR = file.path(getOption("hla.data_root"), "urine", "scRNA_TCMR_DEG.csv"),
    ABMR = file.path(getOption("hla.data_root"), "urine", "scRNA_ABMR_DEG.csv"),
    BK = file.path(getOption("hla.data_root"), "urine", "scRNA_BK_DEG.csv")
  )
  bind_rows(lapply(names(files), function(sig) {
    if (!file.exists(files[[sig]])) return(tibble())
    x <- read.csv(files[[sig]], stringsAsFactors = FALSE, check.names = FALSE)
    gene_col <- if ("gene" %in% names(x)) "gene" else tail(names(x), 1)
    tibble(
      signature = sig,
      gene = as.character(x[[gene_col]]),
      avg_log2FC = suppressWarnings(as.numeric(x[["avg_log2FC"]])),
      p_val_adj = suppressWarnings(as.numeric(x[["p_val_adj"]])),
      pct.1 = if ("pct.1" %in% names(x)) suppressWarnings(as.numeric(x[["pct.1"]])) else NA_real_,
      pct.2 = if ("pct.2" %in% names(x)) suppressWarnings(as.numeric(x[["pct.2"]])) else NA_real_
    )
  })) %>%
    filter(signature %in% signature_levels, !is.na(gene), gene != "") %>%
    group_by(signature, gene) %>%
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
}

gene_module_for <- function(gene) {
  matched <- names(curated_gene_modules)[vapply(curated_gene_modules, function(x) gene %in% x, logical(1))]
  if (length(matched) == 0) "robust_data_driven" else matched[1]
}

message("Reading robust selected genes and urine object...")
selected_genes_file <- file.path(getOption("hla.data_root"), "urine", "Cross_cohort_robust_selected_genes.csv")
if (!file.exists(selected_genes_file)) {
  stop("Missing selected gene table: ", selected_genes_file)
}
# Read the saved object or source table used by this analysis.
selected_genes <- read.csv(selected_genes_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(signature %in% signature_levels, !is.na(gene), gene != "") %>%
  arrange(factor(signature, levels = signature_levels), gene_rank)

object_file <- file.path(getOption("hla.data_root"), "urine", "Urine_ALL.rds")
if (!file.exists(object_file)) {
  object_file <- file.path(project_dir, "Result", "Objects", "Urine_ALL_Intergrate.rds")
}
scRNA <- readRDS(object_file)
DefaultAssay(scRNA) <- "RNA"

meta_cols <- colnames(scRNA@meta.data)
celltype_col <- if ("revised_celltype" %in% meta_cols) "revised_celltype" else if ("celltype" %in% meta_cols) "celltype" else NA_character_
major_col <- if ("Major_celltype" %in% meta_cols) "Major_celltype" else celltype_col
sample_col <- if ("Sample" %in% meta_cols) "Sample" else if ("HLA_Splitter" %in% meta_cols) "HLA_Splitter" else "orig.ident"
if (is.na(celltype_col) || !"Group" %in% meta_cols) {
  stop("The urine object needs Group and celltype metadata.")
}

expr <- get_assay_matrix(scRNA, assay = "RNA")
features <- rownames(expr)
deg_lookup <- read_signature_deg_table(project_dir)
curated_genes <- bind_rows(lapply(names(curated_rejection_gene_sets), function(sig) {
  tibble(
    signature = sig,
    gene = unique(curated_rejection_gene_sets[[sig]]),
    curated_source = "curated_rejection_literature_and_project_markers",
    functional_module = vapply(unique(curated_rejection_gene_sets[[sig]]), gene_module_for, character(1))
  )
}))

selected_gene_core <- selected_genes %>%
  transmute(
    signature,
    gene,
    avg_log2FC,
    p_val_adj,
    pct.1,
    pct.2,
    sample_gene_auc,
    selected_candidate_set,
    selected_ranking_method,
    selected_n_genes,
    urine_auc,
    gse_auc,
    robust_score,
    original_gene_rank = gene_rank,
    gene_set_source = "cross_cohort_robust_selected",
    functional_module = vapply(gene, gene_module_for, character(1))
  )

curated_gene_aug <- curated_genes %>%
  anti_join(selected_gene_core %>% select(signature, gene), by = c("signature", "gene")) %>%
  left_join(deg_lookup, by = c("signature", "gene")) %>%
  mutate(
    sample_gene_auc = NA_real_,
    selected_candidate_set = NA_character_,
    selected_ranking_method = NA_character_,
    selected_n_genes = NA_real_,
    urine_auc = NA_real_,
    gse_auc = NA_real_,
    robust_score = NA_real_,
    original_gene_rank = NA_integer_,
    gene_set_source = curated_source
  ) %>%
  select(-curated_source)

augmented_gene_table <- bind_rows(selected_gene_core, curated_gene_aug) %>%
  filter(gene %in% features) %>%
  group_by(signature, gene) %>%
  arrange(match(gene_set_source, c("cross_cohort_robust_selected", "curated_rejection_literature_and_project_markers")), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(signature) %>%
  mutate(
    gene_rank = row_number(),
    avg_log2FC = ifelse(is.na(avg_log2FC), 1, avg_log2FC),
    p_val_adj = ifelse(is.na(p_val_adj), NA_real_, p_val_adj),
    functional_module = ifelse(is.na(functional_module) | functional_module == "", "robust_data_driven", functional_module)
  ) %>%
  ungroup() %>%
  arrange(factor(signature, levels = signature_levels), gene_rank)

selected_gene_sets <- split(augmented_gene_table$gene, augmented_gene_table$signature)
selected_gene_sets <- lapply(selected_gene_sets[signature_levels], function(x) unique(intersect(x, features)))
selected_gene_sets <- selected_gene_sets[lengths(selected_gene_sets) > 0]
if (length(selected_gene_sets) == 0) {
  stop("No robust signature genes are present in the urine expression matrix.")
}

# Export the analysis source data for figure reproduction.
write_csv_safe(augmented_gene_table, file.path(getOption("hla.data_root"), "urine", "Cross_cohort_robust_augmented_rejection_gene_sets_used_revision.csv"))

curated_requested_table <- curated_genes %>%
  mutate(present_in_urine_object = gene %in% features)
write_csv_safe(curated_requested_table, file.path(result_dir, "Cross_cohort_robust_curated_rejection_genes_requested_revision.csv"))

gene_coverage <- augmented_gene_table %>%
  mutate(present_in_urine_object = gene %in% features) %>%
  group_by(signature) %>%
  summarise(
    final_gene_set_size = n_distinct(gene),
    present_genes = sum(present_in_urine_object),
    missing_genes = final_gene_set_size - present_genes,
    coverage_pct = 100 * present_genes / final_gene_set_size,
    robust_genes = sum(gene_set_source == "cross_cohort_robust_selected" & present_in_urine_object),
    curated_added_genes = sum(gene_set_source == "curated_rejection_literature_and_project_markers" & present_in_urine_object),
    .groups = "drop"
  )
write_csv_safe(gene_coverage, file.path(result_dir, "Cross_cohort_robust_signature_gene_coverage.csv"))

message("Calculating AUCell scores for robust signatures...")
auc_df <- score_auc(expr, selected_gene_sets)
signature_cols <- names(selected_gene_sets)
colnames(auc_df)[match(signature_cols, colnames(auc_df))] <- paste0(signature_cols, "_Signature_AUC")
signature_auc_cols <- paste0(names(selected_gene_sets), "_Signature_AUC")

score_df <- scRNA@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(auc_df, by = "cell") %>%
  transmute(
    cell,
    Sample = as.character(.data[[sample_col]]),
    Group = as.character(Group),
    Major_celltype = as.character(.data[[major_col]]),
    Celltype = as.character(.data[[celltype_col]]),
    across(all_of(signature_auc_cols), as.numeric)
  ) %>%
  filter(!is.na(Group), Group %in% group_levels, !is.na(Sample), Sample != "")
write_csv_safe(score_df, file.path(result_dir, "Cross_cohort_robust_cell_signature_AUCell_scores.csv"))

cell_score_long <- score_df %>%
  pivot_longer(cols = all_of(signature_auc_cols), names_to = "signature_auc", values_to = "auc_score") %>%
  mutate(
    signature = sub("_Signature_AUC$", "", signature_auc),
    signature = factor(signature, levels = signature_levels),
    Group = factor(Group, levels = group_levels)
  ) %>%
  filter(!is.na(signature))

major_contrib <- celltype_contribution(cell_score_long, "Major_celltype")
fine_contrib <- celltype_contribution(cell_score_long, "Celltype")
write_csv_safe(major_contrib, file.path(result_dir, "Cross_cohort_robust_signature_score_major_celltype_contribution_revision.csv"))
write_csv_safe(fine_contrib, file.path(result_dir, "Cross_cohort_robust_signature_score_revised_celltype_contribution_revision.csv"))

message("Estimating gene-level contribution in high-score cells...")
top_cutoffs <- cell_score_long %>%
  group_by(signature) %>%
  summarise(top10_auc_cutoff = quantile(auc_score, 0.90, na.rm = TRUE), .groups = "drop")

gene_contribution <- bind_rows(lapply(names(selected_gene_sets), function(sig) {
  high_cells <- cell_score_long %>%
    filter(signature == sig) %>%
    left_join(top_cutoffs, by = "signature") %>%
    filter(auc_score >= top10_auc_cutoff) %>%
    pull(cell)
  genes <- selected_gene_sets[[sig]]
  if (length(high_cells) == 0 || length(genes) == 0) return(tibble())
  gene_sum <- Matrix::rowSums(expr[genes, high_cells, drop = FALSE])
  total_sum <- sum(gene_sum, na.rm = TRUE)
  tibble(
    signature = sig,
    gene = names(gene_sum),
    expression_sum_in_high_score_cells = as.numeric(gene_sum),
    gene_score_contribution_pct = ifelse(total_sum > 0, 100 * as.numeric(gene_sum) / total_sum, 0)
  )
})) %>%
  left_join(augmented_gene_table, by = c("signature", "gene")) %>%
  group_by(signature) %>%
  mutate(
    contribution_denominator = sum(expression_sum_in_high_score_cells, na.rm = TRUE),
    gene_score_contribution_pct = if_else(
      contribution_denominator > 0,
      100 * expression_sum_in_high_score_cells / contribution_denominator,
      0
    )
  ) %>%
  ungroup() %>%
  select(-contribution_denominator) %>%
  arrange(factor(signature, levels = signature_levels), desc(gene_score_contribution_pct))
write_csv_safe(gene_contribution, file.path(result_dir, "Cross_cohort_robust_gene_score_contribution.csv"))

signature_n_label <- augmented_gene_table %>%
  count(signature, name = "n_genes") %>%
  mutate(signature_label = paste0(signature, "\n", "n=", n_genes, " genes"))
signature_label_map <- setNames(signature_n_label$signature_label, signature_n_label$signature)

gene_plot_df <- gene_contribution %>%
  group_by(signature) %>%
  slice_max(gene_score_contribution_pct, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    signature = factor(signature, levels = signature_levels),
    display_gene = paste(signature, gene, sep = "__")
  ) %>%
  arrange(signature, gene_score_contribution_pct) %>%
  mutate(display_gene = factor(display_gene, levels = display_gene))

p_gene_sets <- ggplot(gene_plot_df, aes(x = gene_score_contribution_pct, y = display_gene, colour = signature)) +
  geom_segment(aes(x = 0, xend = gene_score_contribution_pct, yend = display_gene), linewidth = 0.35, alpha = 0.55) +
  geom_point(aes(size = avg_log2FC), alpha = 0.92) +
  facet_wrap(~ signature, scales = "free_y", nrow = 1, labeller = as_labeller(signature_label_map)) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_colour_manual(values = signature_colors, guide = "none") +
  scale_size_continuous(range = c(1.8, 4.3), name = "Average log2FC") +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4),
    labels = function(x) paste0(formatC(x, format = "f", digits = 1), "%"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Gene contribution to high-score cells", y = "Genes used for AUCell scoring") +
  theme_pub(8.6) +
  theme(
    axis.text.y = element_text(size = 6.8),
    axis.text.x = element_text(size = 7.8),
    strip.text = element_text(size = 9.2),
    legend.position = "right",
    panel.spacing.x = unit(0.25, "lines"),
    plot.margin = margin(3, 4, 3, 4)
  )
save_plot_safe(p_gene_sets, "Signature_gene_sets_used_top10_revision", width = 6.0, height = 3.2)

plot_contribution(major_contrib, "Signature_major_celltype_score_contribution_revision", width = 4.8, height = 4.2)
plot_contribution(fine_contrib, "Signature_revised_celltype_score_contribution_top12_revision", width = 5.2, height = 5.4, n_top = 12)

top_celltype_stack <- cell_score_long %>%
  left_join(top_cutoffs, by = "signature") %>%
  filter(auc_score >= top10_auc_cutoff) %>%
  count(signature, Major_celltype, name = "high_score_cells") %>%
  group_by(signature) %>%
  mutate(high_score_pct = 100 * high_score_cells / sum(high_score_cells)) %>%
  ungroup()
write_csv_safe(top_celltype_stack, file.path(result_dir, "Cross_cohort_robust_top10_percent_high_score_celltype_composition_revision.csv"))

if ("umap" %in% names(scRNA@reductions)) {
  umap_df <- Embeddings(scRNA, "umap") %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")
  umap_df <- umap_df %>%
    left_join(score_df %>% select(cell, Group, Major_celltype, all_of(signature_auc_cols)), by = "cell") %>%
    filter(!is.na(Group))

  umap_long <- umap_df %>%
    pivot_longer(cols = all_of(signature_auc_cols), names_to = "signature_auc", values_to = "auc_score") %>%
    mutate(
      signature = factor(sub("_Signature_AUC$", "", signature_auc), levels = signature_levels),
      auc_score = as.numeric(auc_score)
    ) %>%
    group_by(signature) %>%
    mutate(
      q995 = quantile(auc_score, 0.995, na.rm = TRUE),
      signature_score = ifelse(q995 > 0, pmin(auc_score / q995, 1), 0)
    ) %>%
    ungroup() %>%
    arrange(signature, signature_score)

  p_umap <- ggplot(umap_long, aes(UMAP_1, UMAP_2)) +
    geom_point(colour = "grey88", size = 0.08, alpha = 0.42) +
    geom_point(aes(colour = signature_score), size = 0.11, alpha = 0.92) +
    facet_wrap(~ signature, nrow = 1) +
    scale_colour_gradientn(
      colours = c("#D9D9D9", "#FFF7BC", "#FEC44F", "#FC4E2A", "#99000D"),
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      oob = scales::squish,
      name = "Signature\nscore"
    ) +
    coord_equal() +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_pub(8.5) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      strip.text = element_text(size = 9),
      legend.position = "right",
      panel.spacing.x = unit(0.50, "lines")
    )
  save_plot_safe(p_umap, "Signature_AUCell_score_UMAP_revision", width = 9.2, height = 2.4)
}

message("Done. Outputs written to: ", figure_dir)

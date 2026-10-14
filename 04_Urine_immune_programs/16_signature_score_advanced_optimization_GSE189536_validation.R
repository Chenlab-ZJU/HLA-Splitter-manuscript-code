# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final robust signature construction, scoring and cross-cohort evaluation.
# Original source: 31_revision_signature_score_advanced_optimization_GSE189536_validation.R
# Figure/table scope: Fig4i; Fig4j; FigS8a; FigS8b; FigS8c; FigS8d; FigS8e
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
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
  library(pROC)
})

urine_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (!nzchar(urine_dir)) urine_dir <- file.path(getOption("hla.data_root"), "\u5c3f\u6db2\u6837\u672c")
urine_dir <- gsub("\\\\", "/", urine_dir)
if (!dir.exists(urine_dir)) stop("Cannot find URINE_PROJECT_DIR: ", urine_dir)

gse_dir <- Sys.getenv("GSE189536_PROJECT_DIR", unset = "")
if (!nzchar(gse_dir)) gse_dir <- file.path(getOption("hla.data_root"), "GSE189536")
gse_dir <- gsub("\\\\", "/", gse_dir)
if (!dir.exists(gse_dir)) stop("Cannot find GSE189536_PROJECT_DIR: ", gse_dir)

figure_dir <- file.path(urine_dir, "Figures", "Signature_score_analysis", "Advanced_optimization_GSE189536_validation")
result_dir <- file.path(urine_dir, "Result", "Signature_score_analysis", "Advanced_optimization_GSE189536_validation")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

revision_result_dir <- file.path(urine_dir, "Result", "Signature_score_analysis", "Revision_gene_contribution")
deg_dir <- file.path(urine_dir, "Result", "Signature_score_analysis", "Differential_Expression")

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
signature_levels <- c("TCMR", "ABMR", "BK")
signature_colors <- c("TCMR" = "#699ECA", "ABMR" = "#4DAF4A", "BK" = "#F898CB")
group_colors <- c("NR" = "#E42320", "TCMR" = "#699ECA", "ABMR" = "#4DAF4A", "Mix" = "#FF8C00", "BK" = "#F898CB")

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

read_deg <- function(file, signature_name, cluster_keep) {
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[, nzchar(names(x)), drop = FALSE]
  x %>%
    transmute(
      source_file = basename(file),
      signature = signature_name,
      cluster = as.character(cluster),
      gene = as.character(gene),
      avg_log2FC = suppressWarnings(as.numeric(avg_log2FC)),
      p_val_adj = suppressWarnings(as.numeric(p_val_adj)),
      pct.1 = if ("pct.1" %in% names(.)) suppressWarnings(as.numeric(`pct.1`)) else NA_real_,
      pct.2 = if ("pct.2" %in% names(.)) suppressWarnings(as.numeric(`pct.2`)) else NA_real_
    ) %>%
    filter(cluster %in% cluster_keep, !is.na(gene), gene != "", !is.na(avg_log2FC), avg_log2FC > 0)
}

average_sparse_by_group <- function(mat, group) {
  group <- factor(group)
  design <- sparseMatrix(
    i = seq_along(group),
    j = as.integer(group),
    x = 1,
    dims = c(length(group), nlevels(group)),
    dimnames = list(NULL, levels(group))
  )
  avg <- mat %*% design
  avg <- t(t(avg) / as.numeric(tabulate(as.integer(group), nbins = nlevels(group))))
  avg
}

safe_auc <- function(truth, score) {
  keep <- is.finite(score) & !is.na(truth)
  truth <- truth[keep]
  score <- score[keep]
  if (length(unique(truth)) < 2 || length(unique(score)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(truth, score, quiet = TRUE, direction = "<")))
}

build_candidate_sets <- function(deg_table, sample_gene_auc) {
  deg_table <- deg_table %>%
    left_join(sample_gene_auc, by = c("signature", "gene")) %>%
    mutate(
      pct_diff = pmax(pct.1 - pct.2, 0),
      neg_log10_padj = -log10(pmax(p_val_adj, 1e-300)),
      sample_gene_auc = ifelse(is.na(sample_gene_auc), 0.5, sample_gene_auc),
      rank_pval_logfc = rank(p_val_adj, ties.method = "first") + rank(-avg_log2FC, ties.method = "first") / 100000,
      rank_logfc = rank(-avg_log2FC, ties.method = "first"),
      rank_pctdiff = rank(-pct_diff, ties.method = "first"),
      rank_combined = rank(-(neg_log10_padj * avg_log2FC * pmax(pct_diff, 0.01)), ties.method = "first"),
      rank_sample_auc = rank(-sample_gene_auc, ties.method = "first"),
      rank_sample_auc_logfc = rank(-(sample_gene_auc * avg_log2FC), ties.method = "first"),
      rank_sample_auc_pctdiff = rank(-(sample_gene_auc * pmax(pct_diff, 0.01)), ties.method = "first")
    )

  ranking_methods <- c(
    "pval_logfc" = "rank_pval_logfc",
    "logFC" = "rank_logfc",
    "pct_diff" = "rank_pctdiff",
    "combined" = "rank_combined",
    "sample_auc" = "rank_sample_auc",
    "sample_auc_logFC" = "rank_sample_auc_logfc",
    "sample_auc_pctdiff" = "rank_sample_auc_pctdiff"
  )
  top_n_grid <- c(3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 60, 80, 100, 120, 150, 200)

  candidate_sets <- list()
  candidate_info <- list()
  for (sig in signature_levels) {
    sig_table <- deg_table %>% filter(signature == sig)
    for (method_name in names(ranking_methods)) {
      rank_col <- ranking_methods[[method_name]]
      ordered_genes <- sig_table %>%
        arrange(.data[[rank_col]], p_val_adj, desc(avg_log2FC)) %>%
        pull(gene) %>%
        unique()
      for (n_genes in top_n_grid) {
        genes <- head(ordered_genes, min(n_genes, length(ordered_genes)))
        if (length(genes) < 3) next
        candidate_set <- paste(sig, method_name, length(genes), sep = "__")
        candidate_sets[[candidate_set]] <- genes
        candidate_info[[candidate_set]] <- tibble(
          candidate_set = candidate_set,
          signature = sig,
          ranking_method = method_name,
          n_genes = length(genes)
        )
      }
    }
  }
  list(
    deg_table = deg_table,
    candidate_sets = candidate_sets,
    candidate_info = bind_rows(candidate_info)
  )
}

score_auc_candidates <- function(candidate_sets, candidate_info, rankings, meta_df, prefix) {
  message("Scoring ", prefix, " candidate gene sets: ", length(candidate_sets))
  auc_obj <- AUCell_calcAUC(candidate_sets, rankings, aucMaxRank = ceiling(nrow(rankings) * 0.1), verbose = FALSE)
  auc_mat <- t(as.matrix(getAUC(auc_obj)))
  candidate_sample <- as.data.frame(auc_mat) %>%
    rownames_to_column("cell") %>%
    left_join(meta_df %>% select(cell, Sample, Group), by = "cell") %>%
    filter(!is.na(Sample), !is.na(Group), Group %in% group_levels) %>%
    pivot_longer(cols = all_of(names(candidate_sets)), names_to = "candidate_set", values_to = "auc_score") %>%
    group_by(candidate_set, Sample, Group) %>%
    summarise(mean_auc = mean(auc_score, na.rm = TRUE), n_cells = n(), .groups = "drop") %>%
    left_join(candidate_info, by = "candidate_set")

  perf <- candidate_sample %>%
    group_by(candidate_set, signature, ranking_method, n_genes) %>%
    summarise(
      sample_auc = safe_auc(as.integer(as.character(Group) == signature[1]), mean_auc),
      target_median = median(mean_auc[as.character(Group) == signature[1]], na.rm = TRUE),
      other_median = median(mean_auc[as.character(Group) != signature[1]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(median_delta = target_median - other_median) %>%
    arrange(signature, desc(sample_auc), desc(median_delta), n_genes)
  list(sample_scores = candidate_sample, performance = perf)
}

plot_roc <- function(roc_df, roc_summary, color_values, color_col, filename, width = 3.8, height = 3.4) {
  p <- ggplot(roc_df, aes(x = specificity, y = sensitivity, colour = .data[[color_col]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
    geom_path(linewidth = 1.0) +
    scale_colour_manual(values = color_values, name = NULL, labels = setNames(roc_summary$ROC_label, roc_summary[[color_col]])) +
    scale_x_reverse(limits = c(100, 0), breaks = c(100, 75, 50, 25, 0), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Specificity", y = "Sensitivity") +
    theme_classic(base_size = 10.5, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      legend.position = c(0.58, 0.22),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 9.0),
      plot.margin = margin(5, 6, 5, 6)
    )
  save_plot_safe(p, filename, width, height)
  invisible(p)
}

message("Reading urine object.")
urine_object <- c(
  file.path(urine_dir, "Result", "Objects", "Urine_scRNA.rds"),
  file.path(urine_dir, "Result", "Objects", "Urine_ALL.rds")
)
urine_object <- urine_object[file.exists(urine_object)][1]
# Read the saved object or source table used by this analysis.
urine <- readRDS(urine_object)
DefaultAssay(urine) <- "RNA"
urine_expr <- GetAssayData(urine, assay = "RNA", layer = "data")

sample_col <- if ("Sample" %in% colnames(urine@meta.data)) "Sample" else "HLA_Splitter"
urine_meta <- urine@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(cell, Sample = as.character(.data[[sample_col]]), Group = as.character(Group)) %>%
  filter(!is.na(Sample), Sample != "", !is.na(Group), Group %in% group_levels)

deg_all <- bind_rows(
  read_deg(file.path(deg_dir, "scRNA_TCMR_DEG.csv"), "TCMR", "TCMR"),
  read_deg(file.path(deg_dir, "scRNA_ABMR_DEG.csv"), "ABMR", "ABMR"),
  read_deg(file.path(deg_dir, "scRNA_BK_DEG.csv"), "BK", "BK")
) %>%
  filter(gene %in% rownames(urine_expr), is.na(p_val_adj) | p_val_adj < 0.05, avg_log2FC >= 0.5) %>%
  group_by(signature, gene) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

candidate_genes <- unique(deg_all$gene)
urine_sample_avg <- average_sparse_by_group(urine_expr[candidate_genes, urine_meta$cell, drop = FALSE], urine_meta$Sample)
urine_sample_meta <- urine_meta %>%
  distinct(Sample, Group) %>%
  filter(Sample %in% colnames(urine_sample_avg)) %>%
  arrange(match(Sample, colnames(urine_sample_avg)))

sample_gene_auc <- bind_rows(lapply(signature_levels, function(sig) {
  genes <- deg_all$gene[deg_all$signature == sig]
  genes <- intersect(genes, rownames(urine_sample_avg))
  tibble(
    signature = sig,
    gene = genes,
    sample_gene_auc = vapply(genes, function(g) {
      safe_auc(as.integer(urine_sample_meta$Group == sig), as.numeric(urine_sample_avg[g, urine_sample_meta$Sample]))
    }, numeric(1))
  )
}))
# Export the analysis source data for figure reproduction.
write_csv_safe(sample_gene_auc, file.path(result_dir, "Advanced_signature_single_gene_sample_AUC_urine.csv"))

candidate_pack <- build_candidate_sets(deg_all, sample_gene_auc)
urine_rankings <- AUCell_buildRankings(urine_expr, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE)
urine_candidate <- score_auc_candidates(candidate_pack$candidate_sets, candidate_pack$candidate_info, urine_rankings, urine_meta, "urine")
write_csv_safe(urine_candidate$performance, file.path(result_dir, "Advanced_signature_gene_set_candidates_urine_AUC.csv"))

selected <- urine_candidate$performance %>%
  group_by(signature) %>%
  arrange(desc(sample_auc), desc(median_delta), n_genes, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()
write_csv_safe(selected, file.path(result_dir, "Advanced_signature_gene_set_selected_urine.csv"))

selected_gene_table <- bind_rows(lapply(seq_len(nrow(selected)), function(i) {
  sig <- selected$signature[i]
  candidate_name <- selected$candidate_set[i]
  genes <- candidate_pack$candidate_sets[[candidate_name]]
  candidate_pack$deg_table %>%
    filter(signature == sig, gene %in% genes) %>%
    mutate(
      selected_candidate_set = candidate_name,
      selected_ranking_method = selected$ranking_method[i],
      selected_n_genes = selected$n_genes[i],
      selected_urine_auc = selected$sample_auc[i],
      gene_rank = match(gene, genes)
    ) %>%
    arrange(gene_rank)
})) %>%
  mutate(signature = factor(signature, levels = signature_levels)) %>%
  arrange(signature, gene_rank)
write_csv_safe(selected_gene_table, file.path(result_dir, "Advanced_signature_gene_sets_selected_genes.csv"))

selected_gene_summary <- selected_gene_table %>%
  group_by(signature, selected_ranking_method, selected_n_genes, selected_urine_auc) %>%
  summarise(genes = paste(gene, collapse = "; "), top10_genes = paste(head(gene, 10), collapse = "; "), .groups = "drop")
write_csv_safe(selected_gene_summary, file.path(result_dir, "Advanced_signature_gene_sets_selected_summary.csv"))

selected_gene_sets <- split(selected_gene_table$gene, as.character(selected_gene_table$signature))
selected_gene_sets <- selected_gene_sets[signature_levels]
names(selected_gene_sets) <- paste0(names(selected_gene_sets), "_Advanced_Signature_AUC")

urine_final_auc <- AUCell_calcAUC(selected_gene_sets, urine_rankings, aucMaxRank = ceiling(nrow(urine_rankings) * 0.1), verbose = FALSE)
urine_final_mat <- t(as.matrix(getAUC(urine_final_auc)))
urine_score <- as.data.frame(urine_final_mat) %>%
  rownames_to_column("cell") %>%
  left_join(urine_meta, by = "cell") %>%
  pivot_longer(cols = all_of(colnames(urine_final_mat)), names_to = "signature_auc", values_to = "score") %>%
  mutate(signature = sub("_Advanced_Signature_AUC$", "", signature_auc)) %>%
  group_by(signature, Sample, Group) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), n_cells = n(), .groups = "drop") %>%
  mutate(signature = factor(signature, levels = signature_levels), Group = factor(Group, levels = group_levels))
write_csv_safe(urine_score, file.path(result_dir, "Advanced_signature_sample_scores_urine.csv"))

urine_roc_df <- bind_rows(lapply(signature_levels, function(sig) {
  tmp <- urine_score %>% filter(signature == sig)
  roc_obj <- pROC::roc(as.integer(as.character(tmp$Group) == sig), tmp$mean_score, quiet = TRUE, direction = "<")
  tibble(signature = sig, auc = as.numeric(pROC::auc(roc_obj)), specificity = roc_obj$specificities * 100, sensitivity = roc_obj$sensitivities * 100)
})) %>%
  mutate(signature = factor(signature, levels = signature_levels)) %>%
  arrange(signature, desc(specificity), sensitivity)
urine_roc_summary <- urine_roc_df %>%
  group_by(signature) %>%
  summarise(AUC = first(auc), .groups = "drop") %>%
  mutate(ROC_label = paste0(signature, " AUC = ", sprintf("%.3f", AUC)))
write_csv_safe(urine_roc_df, file.path(result_dir, "Advanced_signature_ROC_curve_urine.csv"))
write_csv_safe(urine_roc_summary, file.path(result_dir, "Advanced_signature_ROC_AUC_urine.csv"))
plot_roc(urine_roc_df, urine_roc_summary, signature_colors, "signature", "Advanced_signature_ROC_AUC_urine", 3.8, 3.4)

p_urine_box <- ggplot(urine_score, aes(Group, mean_score, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.38, linewidth = 0.32) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.5, alpha = 0.86) +
  facet_wrap(~ signature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = alpha(group_colors, 0.45), drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(x = NULL, y = "Mean AUCell score per sample") +
  theme_pub(8.6) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
save_plot_safe(p_urine_box, "Advanced_signature_sample_scores_urine_boxplot", 7.6, 3.2)

message("Reading GSE189536 biopsy immune object.")
gse_object <- file.path(gse_dir, "Harmony_immune_cells_manual_annotation_results", "GSE189536_manual_annotation_immune_filtered_seurat.rds")
gse <- readRDS(gse_object)
DefaultAssay(gse) <- "RNA"
gse_expr <- GetAssayData(gse, assay = "RNA", layer = "data")
gse_meta <- gse@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Sample = as.character(biopsy),
    Group = as.character(sample_type),
    manual_cell_type = as.character(manual_cell_type)
  ) %>%
  filter(!is.na(Sample), Sample != "", !is.na(Group), Group %in% group_levels)

gse_gene_sets <- lapply(selected_gene_sets, function(genes) intersect(genes, rownames(gse_expr)))
gse_gene_coverage <- bind_rows(lapply(names(selected_gene_sets), function(set_name) {
  tibble(
    signature = sub("_Advanced_Signature_AUC$", "", set_name),
    total_genes = length(selected_gene_sets[[set_name]]),
    genes_found = length(gse_gene_sets[[set_name]]),
    genes_missing = length(setdiff(selected_gene_sets[[set_name]], rownames(gse_expr))),
    found_gene_list = paste(gse_gene_sets[[set_name]], collapse = "; "),
    missing_gene_list = paste(setdiff(selected_gene_sets[[set_name]], rownames(gse_expr)), collapse = "; ")
  )
}))
write_csv_safe(gse_gene_coverage, file.path(result_dir, "Advanced_signature_gene_coverage_GSE189536.csv"))
gse_gene_sets <- gse_gene_sets[lengths(gse_gene_sets) >= 3]

gse_rankings <- AUCell_buildRankings(gse_expr, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE)
gse_auc <- AUCell_calcAUC(gse_gene_sets, gse_rankings, aucMaxRank = ceiling(nrow(gse_rankings) * 0.1), verbose = FALSE)
gse_mat <- t(as.matrix(getAUC(gse_auc)))
gse_score <- as.data.frame(gse_mat) %>%
  rownames_to_column("cell") %>%
  left_join(gse_meta, by = "cell") %>%
  pivot_longer(cols = all_of(colnames(gse_mat)), names_to = "signature_auc", values_to = "score") %>%
  mutate(signature = sub("_Advanced_Signature_AUC$", "", signature_auc)) %>%
  group_by(signature, Sample, Group) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), n_cells = n(), .groups = "drop") %>%
  mutate(signature = factor(signature, levels = signature_levels), Group = factor(Group, levels = group_levels))
write_csv_safe(gse_score, file.path(result_dir, "Advanced_signature_sample_scores_GSE189536.csv"))

gse_roc_df <- bind_rows(lapply(intersect(c("TCMR", "ABMR"), unique(as.character(gse_score$signature))), function(sig) {
  tmp <- gse_score %>% filter(signature == sig)
  roc_obj <- pROC::roc(as.integer(as.character(tmp$Group) == sig), tmp$mean_score, quiet = TRUE, direction = "<")
  tibble(signature = sig, auc = as.numeric(pROC::auc(roc_obj)), specificity = roc_obj$specificities * 100, sensitivity = roc_obj$sensitivities * 100)
})) %>%
  mutate(signature = factor(signature, levels = signature_levels)) %>%
  arrange(signature, desc(specificity), sensitivity)
gse_roc_summary <- gse_roc_df %>%
  group_by(signature) %>%
  summarise(AUC = first(auc), .groups = "drop") %>%
  mutate(ROC_label = paste0(signature, " AUC = ", sprintf("%.3f", AUC)))
write_csv_safe(gse_roc_df, file.path(result_dir, "Advanced_signature_ROC_curve_GSE189536.csv"))
write_csv_safe(gse_roc_summary, file.path(result_dir, "Advanced_signature_ROC_AUC_GSE189536.csv"))
plot_roc(gse_roc_df, gse_roc_summary, signature_colors, "signature", "Advanced_signature_ROC_AUC_GSE189536", 3.8, 3.4)

p_gse_box <- ggplot(gse_score, aes(Group, mean_score, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.38, linewidth = 0.32) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.8, alpha = 0.90) +
  facet_wrap(~ signature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = alpha(group_colors, 0.45), drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(x = NULL, y = "Mean AUCell score per biopsy") +
  theme_pub(8.6) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
save_plot_safe(p_gse_box, "Advanced_signature_sample_scores_GSE189536_boxplot", 7.2, 3.2)

validation_note <- tibble(
  note = c(
    "Advanced signatures were optimized only in urine samples by scanning DEG ranking strategies and top-N gene counts.",
    "GSE189536 kidney-biopsy immune cells were used only for fixed-gene-set validation.",
    "GSE189536 contains NR, TCMR, ABMR and Mix groups but no BK-positive biopsy samples; BK ROC validation is therefore not applicable in this dataset."
  )
)
write_csv_safe(validation_note, file.path(result_dir, "Advanced_signature_validation_notes.csv"))

message("Done. Advanced optimization and GSE189536 validation outputs: ", result_dir)

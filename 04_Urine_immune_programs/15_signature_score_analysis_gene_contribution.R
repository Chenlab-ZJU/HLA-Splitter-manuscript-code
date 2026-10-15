# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final robust signature construction, scoring and cross-cohort evaluation.
# Original source: 30_revision_signature_score_analysis_gene_contribution.R
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
  library(patchwork)
  library(scales)
  library(pROC)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "30_revision_signature_score_analysis_gene_contribution.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures", "Signature_score_analysis")
result_dir <- file.path(project_dir, "Result", "Signature_score_analysis")
deg_dir <- file.path(result_dir, "Differential_Expression")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "Revision_gene_contribution"), recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
signature_levels <- c("TCMR", "ABMR", "BK")
signature_colors <- c(
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "BK" = "#F898CB"
)
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "Mix" = "#FF8C00",
  "BK" = "#F898CB"
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

format_p_value <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    p < 0.01 ~ formatC(p, format = "f", digits = 3),
    TRUE ~ formatC(p, format = "f", digits = 2)
  )
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

read_deg <- function(file, signature_name = NULL, cluster_keep = NULL) {
  if (!file.exists(file)) return(tibble())
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[, nzchar(names(x)), drop = FALSE]
  required <- c("gene", "cluster", "avg_log2FC", "p_val_adj")
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    stop("DEG file lacks required columns: ", basename(file), " -> ", paste(missing_cols, collapse = ", "))
  }
  x <- x %>%
    transmute(
      source_file = basename(file),
      signature = if (is.null(signature_name)) as.character(cluster) else signature_name,
      cluster = as.character(cluster),
      gene = as.character(gene),
      avg_log2FC = suppressWarnings(as.numeric(avg_log2FC)),
      p_val_adj = suppressWarnings(as.numeric(p_val_adj)),
      pct.1 = if ("pct.1" %in% names(.)) suppressWarnings(as.numeric(`pct.1`)) else NA_real_,
      pct.2 = if ("pct.2" %in% names(.)) suppressWarnings(as.numeric(`pct.2`)) else NA_real_
    ) %>%
    filter(!is.na(gene), gene != "", !is.na(avg_log2FC), avg_log2FC > 0)
  if (!is.null(cluster_keep)) x <- x %>% filter(cluster %in% cluster_keep)
  x
}

make_gene_sets <- function(features, max_genes = 150, min_genes = 10) {
  deg_all <- bind_rows(
    read_deg(file.path(getOption("hla.data_root"), "urine", "scRNA_TCMR_DEG.csv"), signature_name = "TCMR", cluster_keep = "TCMR"),
    read_deg(file.path(getOption("hla.data_root"), "urine", "scRNA_ABMR_DEG.csv"), signature_name = "ABMR", cluster_keep = "ABMR"),
    read_deg(file.path(getOption("hla.data_root"), "urine", "scRNA_BK_DEG.csv"), signature_name = "BK", cluster_keep = "BK")
  ) %>%
    filter(signature %in% signature_levels, gene %in% features) %>%
    group_by(signature, gene) %>%
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    filter(is.na(p_val_adj) | p_val_adj < 0.05, avg_log2FC >= 1) %>%
    arrange(factor(signature, levels = signature_levels), p_val_adj, desc(avg_log2FC))

  gene_table <- deg_all %>%
    mutate(
      pct_diff = pmax(pct.1 - pct.2, 0),
      neg_log10_padj = -log10(pmax(p_val_adj, 1e-300)),
      ranking_pval_logfc = rank(p_val_adj, ties.method = "first") + rank(-avg_log2FC, ties.method = "first") / 100000,
      ranking_logfc = rank(-avg_log2FC, ties.method = "first"),
      ranking_pctdiff = rank(-pct_diff, ties.method = "first"),
      ranking_combined = rank(-(neg_log10_padj * avg_log2FC * pmax(pct_diff, 0.01)), ties.method = "first")
    ) %>%
    group_by(signature) %>%
    slice_head(n = max_genes) %>%
    mutate(
      signature = factor(signature, levels = signature_levels),
      gene_rank = row_number()
    ) %>%
    ungroup() %>%
    arrange(signature, gene_rank)

  counts <- table(gene_table$signature)
  too_few <- names(counts[counts < min_genes])
  if (length(too_few) > 0) {
    stop("Too few genes for signature(s): ", paste(too_few, collapse = ", "))
  }

  gene_sets <- split(gene_table$gene, as.character(gene_table$signature))
  gene_sets <- gene_sets[signature_levels[signature_levels %in% names(gene_sets)]]
  names(gene_sets) <- paste0(names(gene_sets), "_Signature_AUC")
  attr(gene_sets, "gene_table") <- gene_table
  gene_sets
}

optimize_signature_gene_sets <- function(gene_table, cells_rankings, meta_df, optimize_signatures = c("TCMR", "BK")) {
  top_n_grid <- c(5, 8, 10, 12, 15, 20, 30, 40, 60, 80, 100, 120, 150)
  ranking_methods <- c(
    "pval_logfc" = "ranking_pval_logfc",
    "logFC" = "ranking_logfc",
    "pct_diff" = "ranking_pctdiff",
    "combined" = "ranking_combined"
  )

  candidate_sets <- list()
  candidate_info <- list()
  for (sig in signature_levels) {
    sig_table <- gene_table %>% filter(as.character(signature) == sig)
    if (nrow(sig_table) == 0) next
    methods_to_use <- if (sig %in% optimize_signatures) names(ranking_methods) else "pval_logfc"
    top_n_to_use <- if (sig %in% optimize_signatures) top_n_grid else 80
    for (method_name in methods_to_use) {
      rank_col <- ranking_methods[[method_name]]
      ordered_genes <- sig_table %>%
        arrange(.data[[rank_col]], p_val_adj, desc(avg_log2FC)) %>%
        pull(gene) %>%
        unique()
      for (n_genes in top_n_to_use) {
        genes <- head(ordered_genes, min(n_genes, length(ordered_genes)))
        if (length(genes) < 5) next
        set_name <- paste(sig, method_name, length(genes), sep = "__")
        candidate_sets[[set_name]] <- genes
        candidate_info[[set_name]] <- tibble(
          candidate_set = set_name,
          signature = sig,
          ranking_method = method_name,
          n_genes = length(genes)
        )
      }
    }
  }

  message("Optimizing candidate gene sets: ", length(candidate_sets), " candidates.")
  candidate_auc <- AUCell_calcAUC(candidate_sets, cells_rankings, aucMaxRank = ceiling(nrow(cells_rankings) * 0.1), verbose = FALSE)
  candidate_mat <- t(as.matrix(getAUC(candidate_auc)))
  candidate_score <- as.data.frame(candidate_mat) %>%
    rownames_to_column("cell") %>%
    left_join(meta_df %>% select(cell, Sample, Group), by = "cell") %>%
    filter(!is.na(Sample), !is.na(Group), Group %in% group_levels)

  candidate_long <- candidate_score %>%
    pivot_longer(cols = all_of(names(candidate_sets)), names_to = "candidate_set", values_to = "auc_score") %>%
    group_by(candidate_set, Sample, Group) %>%
    summarise(mean_auc = mean(auc_score, na.rm = TRUE), n_cells = n(), .groups = "drop") %>%
    left_join(bind_rows(candidate_info), by = "candidate_set")

  candidate_perf <- candidate_long %>%
    group_by(candidate_set, signature, ranking_method, n_genes) %>%
    summarise(
      sample_auc = {
        truth <- as.integer(as.character(Group) == signature[1])
        if (length(unique(truth)) < 2 || length(unique(mean_auc)) < 2) {
          NA_real_
        } else {
          as.numeric(pROC::auc(pROC::roc(truth, mean_auc, quiet = TRUE, direction = "<")))
        }
      },
      target_median = median(mean_auc[as.character(Group) == signature[1]], na.rm = TRUE),
      other_median = median(mean_auc[as.character(Group) != signature[1]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      median_delta = target_median - other_median,
      optimized = signature %in% optimize_signatures
    ) %>%
    arrange(signature, desc(sample_auc), desc(median_delta), n_genes)

  selected <- candidate_perf %>%
    group_by(signature) %>%
    arrange(desc(sample_auc), desc(median_delta), n_genes, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()

  final_gene_table <- bind_rows(lapply(seq_len(nrow(selected)), function(i) {
    sig <- selected$signature[i]
    method_name <- selected$ranking_method[i]
    n_genes <- selected$n_genes[i]
    rank_col <- ranking_methods[[method_name]]
    gene_table %>%
      filter(as.character(signature) == sig) %>%
      arrange(.data[[rank_col]], p_val_adj, desc(avg_log2FC)) %>%
      slice_head(n = n_genes) %>%
      mutate(
        optimized_ranking_method = method_name,
        optimized_n_genes = n_genes,
        optimized_sample_auc = selected$sample_auc[i],
        gene_rank = row_number()
      )
  })) %>%
    mutate(signature = factor(as.character(signature), levels = signature_levels)) %>%
    arrange(signature, gene_rank)

  final_gene_sets <- split(final_gene_table$gene, as.character(final_gene_table$signature))
  final_gene_sets <- final_gene_sets[signature_levels[signature_levels %in% names(final_gene_sets)]]
  names(final_gene_sets) <- paste0(names(final_gene_sets), "_Signature_AUC")
  attr(final_gene_sets, "gene_table") <- final_gene_table
  attr(final_gene_sets, "candidate_performance") <- candidate_perf
  attr(final_gene_sets, "selected_performance") <- selected
  final_gene_sets
}

object_candidates <- c(
  file.path(project_dir, "Result", "Objects", "Urine_scRNA.rds"),
  file.path(getOption("hla.data_root"), "urine", "Urine_ALL.rds")
)
object_file <- object_candidates[file.exists(object_candidates)][1]
if (is.na(object_file)) stop("Cannot find Urine scRNA object in Result/Objects.")

message("Reading object: ", object_file)
# Read the saved object or source table used by this analysis.
scRNA <- readRDS(object_file)
DefaultAssay(scRNA) <- "RNA"

meta_cols <- colnames(scRNA@meta.data)
celltype_col <- if ("revised_celltype" %in% meta_cols) "revised_celltype" else if ("celltype" %in% meta_cols) "celltype" else NA_character_
major_col <- if ("Major_celltype" %in% meta_cols) "Major_celltype" else celltype_col
sample_col <- if ("Sample" %in% meta_cols) "Sample" else if ("HLA_Splitter" %in% meta_cols) "HLA_Splitter" else "orig.ident"
if (is.na(celltype_col) || !"Group" %in% meta_cols) stop("The object needs Group and celltype metadata.")

expr <- GetAssayData(scRNA, assay = "RNA", layer = "data")
message("Building AUCell rankings.")
cells_rankings <- AUCell_buildRankings(expr, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE)

optimization_meta <- scRNA@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Sample = as.character(.data[[sample_col]]),
    Group = as.character(Group)
  )
gene_set_source <- file.path(getOption("hla.data_root"), "urine", "Signature_gene_sets_used_revision.csv")
recompute_sets <- tolower(Sys.getenv("HLA_RECOMPUTE_SIGNATURES", "false")) %in% c("true", "1", "yes")
if (file.exists(gene_set_source) && !recompute_sets) {
  # Reproduce the existing selected signatures without repeating the candidate search.
  cached_genes <- read.csv(gene_set_source, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("signature", "gene", "gene_rank") %in% names(cached_genes)))
  cached_genes <- cached_genes %>% arrange(signature, gene_rank)
  gene_sets <- split(cached_genes$gene, cached_genes$signature)
  gene_sets <- gene_sets[signature_levels[signature_levels %in% names(gene_sets)]]
  names(gene_sets) <- paste0(names(gene_sets), "_Signature_AUC")
  attr(gene_sets, "gene_table") <- cached_genes
  message("Using saved selected signature genes: ", gene_set_source)
} else {
  initial_gene_sets <- make_gene_sets(rownames(expr), max_genes = 150)
  initial_gene_table <- attr(initial_gene_sets, "gene_table")
  gene_sets <- optimize_signature_gene_sets(initial_gene_table, cells_rankings, optimization_meta, optimize_signatures = c("TCMR", "BK"))
  write_csv_safe(attr(gene_sets, "candidate_performance"), file.path(result_dir, "Revision_gene_contribution", "Signature_gene_set_optimization_candidates_revision.csv"))
  write_csv_safe(attr(gene_sets, "selected_performance"), file.path(result_dir, "Revision_gene_contribution", "Signature_gene_set_optimization_selected_revision.csv"))
}

gene_table <- attr(gene_sets, "gene_table") %>%
  mutate(signature = as.character(signature), signature_auc = paste0(signature, "_Signature_AUC"))
if (!file.exists(gene_set_source) || recompute_sets) {
  write_csv_safe(gene_table, gene_set_source)
}

gene_set_summary <- gene_table %>%
  group_by(signature, signature_auc, optimized_ranking_method, optimized_n_genes, optimized_sample_auc) %>%
  summarise(
    n_genes = n(),
    genes = paste(gene, collapse = "; "),
    top20_genes = paste(head(gene, 20), collapse = "; "),
    .groups = "drop"
  )
write_csv_safe(gene_set_summary, file.path(result_dir, "Revision_gene_contribution", "Signature_gene_sets_used_summary_revision.csv"))

message("Calculating AUCell scores for optimized signatures: ", length(gene_sets), " signatures.")
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, aucMaxRank = ceiling(nrow(cells_rankings) * 0.1), verbose = FALSE)
auc_mat <- t(as.matrix(getAUC(cells_AUC)))
for (sig in colnames(auc_mat)) {
  scRNA[[sig]] <- as.numeric(auc_mat[, sig])
}

signature_cols <- colnames(auc_mat)
signature_labels <- sub("_Signature_AUC$", "", signature_cols)

score_df <- scRNA@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Sample = as.character(.data[[sample_col]]),
    Group = as.character(Group),
    Major_celltype = as.character(.data[[major_col]]),
    Celltype = as.character(.data[[celltype_col]]),
    across(all_of(signature_cols), as.numeric)
  ) %>%
  filter(!is.na(Group), Group %in% group_levels, !is.na(Sample), Sample != "") %>%
  mutate(Group = factor(Group, levels = group_levels))
write_csv_safe(score_df, file.path(result_dir, "Revision_gene_contribution", "Cell_signature_AUCell_scores_revision.csv"))

sample_score <- score_df %>%
  group_by(Sample, Group) %>%
  summarise(across(all_of(signature_cols), \(x) mean(x, na.rm = TRUE)), n_cells = n(), .groups = "drop") %>%
  pivot_longer(cols = all_of(signature_cols), names_to = "signature_auc", values_to = "mean_auc") %>%
  mutate(signature = factor(sub("_Signature_AUC$", "", signature_auc), levels = signature_levels))
write_csv_safe(sample_score, file.path(result_dir, "Revision_gene_contribution", "Sample_signature_AUCell_scores_revision.csv"))

cell_score_long <- score_df %>%
  pivot_longer(cols = all_of(signature_cols), names_to = "signature_auc", values_to = "auc_score") %>%
  mutate(signature = factor(sub("_Signature_AUC$", "", signature_auc), levels = signature_levels))

celltype_contribution <- function(df, type_col) {
  df %>%
    group_by(signature, celltype = .data[[type_col]]) %>%
    summarise(
      n_cells = n(),
      mean_auc = mean(auc_score, na.rm = TRUE),
      median_auc = median(auc_score, na.rm = TRUE),
      total_auc = sum(auc_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(signature) %>%
    mutate(
      contribution_pct = 100 * total_auc / sum(total_auc, na.rm = TRUE),
      rank_by_contribution = rank(-contribution_pct, ties.method = "first"),
      signature = factor(signature, levels = signature_levels)
    ) %>%
    ungroup() %>%
    arrange(signature, rank_by_contribution)
}

major_contrib <- celltype_contribution(cell_score_long, "Major_celltype")
fine_contrib <- celltype_contribution(cell_score_long, "Celltype")
write_csv_safe(major_contrib, file.path(result_dir, "Revision_gene_contribution", "Signature_score_major_celltype_contribution_revision.csv"))
write_csv_safe(fine_contrib, file.path(result_dir, "Revision_gene_contribution", "Signature_score_revised_celltype_contribution_revision.csv"))

top_cell_thresholds <- cell_score_long %>%
  group_by(signature) %>%
  summarise(top10_auc_cutoff = quantile(auc_score, 0.90, na.rm = TRUE), .groups = "drop")

top_score_cells <- cell_score_long %>%
  left_join(top_cell_thresholds, by = "signature") %>%
  filter(auc_score >= top10_auc_cutoff) %>%
  group_by(signature, cell) %>%
  slice_max(order_by = auc_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(signature, desc(auc_score))
write_csv_safe(top_score_cells, file.path(result_dir, "Revision_gene_contribution", "Signature_top10_percent_high_score_cells_revision.csv"))

top_celltype_stack <- top_score_cells %>%
  count(signature, Major_celltype, name = "high_score_cells") %>%
  group_by(signature) %>%
  mutate(high_score_pct = 100 * high_score_cells / sum(high_score_cells)) %>%
  ungroup()
write_csv_safe(top_celltype_stack, file.path(result_dir, "Revision_gene_contribution", "Signature_top10_percent_high_score_celltype_composition_revision.csv"))

top_cells_ranked <- top_score_cells %>%
  group_by(signature) %>%
  slice_max(order_by = auc_score, n = 200, with_ties = FALSE) %>%
  ungroup()
write_csv_safe(top_cells_ranked, file.path(result_dir, "Revision_gene_contribution", "Signature_top200_cells_by_AUCell_score_revision.csv"))

gene_contribution <- lapply(signature_levels, function(sig) {
  sig_auc <- paste0(sig, "_Signature_AUC")
  sig_genes <- gene_table$gene[as.character(gene_table$signature) == sig]
  sig_genes <- intersect(sig_genes, rownames(expr))
  high_cells <- top_score_cells$cell[as.character(top_score_cells$signature) == sig]
  high_cells <- intersect(high_cells, colnames(expr))
  if (length(sig_genes) == 0 || length(high_cells) == 0) return(tibble())
  gene_total <- Matrix::rowSums(expr[sig_genes, high_cells, drop = FALSE])
  tibble(
    signature = sig,
    signature_auc = sig_auc,
    gene = names(gene_total),
    high_score_gene_expression = as.numeric(gene_total),
    gene_score_contribution_pct = 100 * as.numeric(gene_total) / sum(as.numeric(gene_total), na.rm = TRUE)
  )
}) %>%
  bind_rows() %>%
  left_join(gene_table %>% select(signature, gene, avg_log2FC, p_val_adj, pct.1, pct.2), by = c("signature", "gene")) %>%
  group_by(signature) %>%
  arrange(desc(gene_score_contribution_pct), p_val_adj, .by_group = TRUE) %>%
  mutate(contribution_rank = row_number()) %>%
  ungroup()
write_csv_safe(gene_contribution, file.path(result_dir, "Revision_gene_contribution", "Signature_gene_score_contribution_revision.csv"))

gene_plot_df <- gene_contribution %>%
  group_by(signature) %>%
  slice_head(n = 10) %>%
  mutate(
    display_gene = paste(as.character(signature), gene, sep = "__"),
    signature = factor(signature, levels = signature_levels)
  ) %>%
  ungroup() %>%
  mutate(display_gene = factor(display_gene, levels = rev(display_gene)))

p_gene_sets <- ggplot(gene_plot_df, aes(x = gene_score_contribution_pct, y = display_gene, colour = signature)) +
  geom_segment(aes(x = 0, xend = gene_score_contribution_pct, yend = display_gene), linewidth = 0.35, alpha = 0.55) +
  geom_point(aes(size = avg_log2FC), alpha = 0.92) +
  facet_wrap(~ signature, scales = "free_y", nrow = 1) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_colour_manual(values = signature_colors, guide = "none") +
  scale_size_continuous(range = c(1.7, 4.2), name = "Average log2FC") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.05))) +
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

plot_contribution <- function(df, filename, width, height, n_top = NULL) {
  plot_df <- df
  if (!is.null(n_top)) {
    plot_df <- df %>% group_by(signature) %>% slice_min(rank_by_contribution, n = n_top, with_ties = FALSE) %>% ungroup()
  }
  type_levels <- plot_df %>%
    group_by(celltype) %>%
    summarise(max_contribution = max(contribution_pct, na.rm = TRUE), .groups = "drop") %>%
    arrange(max_contribution) %>%
    pull(celltype)
  plot_df <- plot_df %>%
    mutate(celltype = factor(celltype, levels = type_levels), signature = factor(signature, levels = signature_levels))

  p <- ggplot(plot_df, aes(x = signature, y = celltype)) +
    geom_point(aes(size = contribution_pct, fill = mean_auc), shape = 21, colour = "grey25", stroke = 0.22, alpha = 0.94) +
    scale_size_continuous(range = c(1.5, 8.5), name = "Contribution (%)") +
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#D0E1F2", "#8EBAE5", "#3F86C5", "#F6D55C", "#E76F51", "#9E0142"),
      name = "Mean AUCell"
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

plot_contribution(major_contrib, "Signature_major_celltype_score_contribution_revision", width = 4.8, height = 4.2)
plot_contribution(fine_contrib, "Signature_revised_celltype_score_contribution_top12_revision", width = 5.2, height = 5.4, n_top = 12)

stack_df <- top_celltype_stack %>%
  mutate(signature = factor(signature, levels = signature_levels)) %>%
  group_by(Major_celltype) %>%
  mutate(total = sum(high_score_pct)) %>%
  ungroup()
major_types <- stack_df %>% distinct(Major_celltype, total) %>% arrange(desc(total)) %>% pull(Major_celltype)
major_palette <- setNames(hue_pal(l = 62, c = 65)(length(major_types)), major_types)

p_top_stack <- ggplot(stack_df, aes(x = signature, y = high_score_pct, fill = factor(Major_celltype, levels = major_types))) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.18) +
  scale_fill_manual(values = major_palette, name = "Major celltype") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Signature", y = "High-score cells (%)") +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )
save_plot_safe(p_top_stack, "Signature_top10_percent_high_score_celltype_stack_revision", width = 5.4, height = 4.2)

if ("umap" %in% names(scRNA@reductions)) {
  umap_df <- Embeddings(scRNA, "umap") %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")
  umap_df <- umap_df %>%
    left_join(score_df %>% select(cell, Group, Major_celltype, all_of(signature_cols)), by = "cell") %>%
    filter(!is.na(Group))

  umap_long <- umap_df %>%
    pivot_longer(cols = all_of(signature_cols), names_to = "signature_auc", values_to = "auc_score") %>%
    mutate(signature = factor(sub("_Signature_AUC$", "", signature_auc), levels = signature_levels)) %>%
    group_by(signature) %>%
    mutate(
      normalized_auc_score = ifelse(
        max(auc_score, na.rm = TRUE) > min(auc_score, na.rm = TRUE),
        (auc_score - min(auc_score, na.rm = TRUE)) / (max(auc_score, na.rm = TRUE) - min(auc_score, na.rm = TRUE)),
        0
      )
    ) %>%
    ungroup()

  p_umap <- ggplot(umap_long, aes(UMAP_1, UMAP_2)) +
    geom_point(colour = "grey88", size = 0.08, alpha = 0.45) +
    geom_point(aes(colour = auc_score), size = 0.10, alpha = 0.90) +
    facet_wrap(~ signature, nrow = 1) +
    scale_colour_gradientn(
      colours = c("#D9D9D9", "#FFF7BC", "#FEC44F", "#FC4E2A", "#99000D"),
      limits = c(0, quantile(umap_long$auc_score, 0.995, na.rm = TRUE)),
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

sample_wilcox <- lapply(signature_levels, function(sig) {
  target_groups <- switch(
    sig,
    "TCMR" = list(c("NR", "TCMR"), c("TCMR", "ABMR")),
    "ABMR" = list(c("NR", "ABMR"), c("TCMR", "ABMR")),
    "BK" = list(c("NR", "BK"), c("TCMR", "BK")),
    list(c("NR", sig))
  )
  tmp <- sample_score %>% filter(signature == sig)
  max_y <- max(tmp$mean_auc, na.rm = TRUE)
  min_y <- min(tmp$mean_auc, na.rm = TRUE)
  step_y <- max((max_y - min_y) * 0.12, max_y * 0.04, 0.004)
  bind_rows(lapply(seq_along(target_groups), function(i) {
    g1 <- target_groups[[i]][1]
    g2 <- target_groups[[i]][2]
    v1 <- tmp$mean_auc[as.character(tmp$Group) == g1]
    v2 <- tmp$mean_auc[as.character(tmp$Group) == g2]
    p <- if (length(v1) > 0 && length(v2) > 0 && length(unique(c(v1, v2))) > 1) {
      suppressWarnings(wilcox.test(v1, v2, exact = FALSE)$p.value)
    } else {
      NA_real_
    }
    tibble(
      signature = factor(sig, levels = signature_levels),
      group1 = g1,
      group2 = g2,
      x = g1,
      xend = g2,
      y = max_y + step_y * (i + 0.25),
      y_text = max_y + step_y * (i + 0.55),
      p_value = p,
      p_label = format_p_value(p)
    )
  }))
}) %>%
  bind_rows()
write_csv_safe(sample_wilcox, file.path(result_dir, "Revision_gene_contribution", "Sample_signature_AUCell_wilcox_revision.csv"))

sample_wilcox_plot <- sample_wilcox %>%
  filter(!is.na(p_value), p_value < 0.1)

p_sample <- ggplot(sample_score, aes(x = Group, y = mean_auc, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.38, linewidth = 0.32) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.5, alpha = 0.86) +
  geom_segment(
    data = sample_wilcox_plot,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = sample_wilcox_plot,
    aes(x = x, xend = x, y = y, yend = y - 0.003),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = sample_wilcox_plot,
    aes(x = xend, xend = xend, y = y, yend = y - 0.003),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_text(
    data = sample_wilcox_plot,
    aes(x = (as.numeric(factor(x, levels = group_levels)) + as.numeric(factor(xend, levels = group_levels))) / 2, y = y_text, label = p_label),
    inherit.aes = FALSE,
    size = 3.4,
    family = "Arial",
    colour = "black"
  ) +
  facet_wrap(~ signature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = alpha(group_colors, 0.45), drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "Mean AUCell score per sample") +
  theme_pub(8.6) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none",
    panel.spacing.x = unit(0.65, "lines")
  )
save_plot_safe(p_sample, "Signature_sample_level_AUCell_score_revision", width = 7.6, height = 3.2)

roc_df <- lapply(signature_levels, function(sig) {
  sig_auc <- paste0(sig, "_Signature_AUC")
  tmp <- sample_score %>% filter(signature == sig, !is.na(mean_auc))
  truth <- as.integer(as.character(tmp$Group) == sig)
  if (length(unique(truth)) < 2 || length(unique(tmp$mean_auc)) < 2) return(tibble())
  roc_obj <- pROC::roc(truth, tmp$mean_auc, quiet = TRUE, direction = "<")
  tibble(
    signature = sig,
    auc = as.numeric(pROC::auc(roc_obj)),
    specificity = roc_obj$specificities * 100,
    sensitivity = roc_obj$sensitivities * 100
  )
}) %>%
  bind_rows() %>%
  mutate(signature = factor(signature, levels = signature_levels)) %>%
  arrange(signature, desc(specificity), sensitivity)
roc_summary <- roc_df %>%
  group_by(signature) %>%
  summarise(AUC = first(auc), .groups = "drop") %>%
  mutate(ROC_label = paste0(signature, " AUC = ", sprintf("%.3f", AUC)))
write_csv_safe(roc_df, file.path(result_dir, "Revision_gene_contribution", "Signature_sample_level_ROC_curve_revision.csv"))
write_csv_safe(roc_summary, file.path(result_dir, "Revision_gene_contribution", "Signature_sample_level_ROC_AUC_revision.csv"))

p_roc <- ggplot(roc_df, aes(x = specificity, y = sensitivity, colour = signature)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
  geom_path(linewidth = 1.0) +
  scale_colour_manual(values = signature_colors, name = NULL, labels = setNames(roc_summary$ROC_label, roc_summary$signature)) +
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
    legend.text = element_text(size = 9.5),
    plot.margin = margin(5, 6, 5, 6)
  )
save_plot_safe(p_roc, "Signature_sample_level_ROC_AUC_revision", width = 3.7, height = 3.4)

combo_definitions <- list(
  "TCMR+ABMR" = c("TCMR", "ABMR"),
  "TCMR+BK" = c("TCMR", "BK"),
  "ABMR+BK" = c("ABMR", "BK")
)
combo_colors <- c(
  "TCMR+ABMR" = signature_colors[["TCMR"]],
  "TCMR+BK" = signature_colors[["BK"]],
  "ABMR+BK" = signature_colors[["ABMR"]]
)

sample_score_wide <- sample_score %>%
  select(Sample, Group, signature, mean_auc) %>%
  mutate(signature = as.character(signature)) %>%
  pivot_wider(names_from = signature, values_from = mean_auc) %>%
  mutate(across(all_of(signature_levels), \(x) as.numeric(scale(x))))

combo_roc_df <- lapply(names(combo_definitions), function(combo_name) {
  combo_groups <- combo_definitions[[combo_name]]
  tmp <- sample_score_wide %>%
    mutate(
      combo_score = rowMeans(across(all_of(combo_groups)), na.rm = TRUE),
      positive = as.integer(as.character(Group) %in% combo_groups)
    )
  if (length(unique(tmp$positive)) < 2 || length(unique(tmp$combo_score)) < 2) return(tibble())
  roc_obj <- pROC::roc(tmp$positive, tmp$combo_score, quiet = TRUE, direction = "<")
  tibble(
    combination = combo_name,
    auc = as.numeric(pROC::auc(roc_obj)),
    specificity = roc_obj$specificities * 100,
    sensitivity = roc_obj$sensitivities * 100
  )
}) %>%
  bind_rows() %>%
  mutate(combination = factor(combination, levels = names(combo_definitions))) %>%
  arrange(combination, desc(specificity), sensitivity)

combo_roc_summary <- combo_roc_df %>%
  group_by(combination) %>%
  summarise(AUC = first(auc), .groups = "drop") %>%
  mutate(ROC_label = paste0(combination, " AUC = ", sprintf("%.3f", AUC)))
write_csv_safe(combo_roc_df, file.path(result_dir, "Revision_gene_contribution", "Signature_score_combination_ROC_curve_revision.csv"))
write_csv_safe(combo_roc_summary, file.path(result_dir, "Revision_gene_contribution", "Signature_score_combination_ROC_AUC_revision.csv"))

p_combo_roc <- ggplot(combo_roc_df, aes(x = specificity, y = sensitivity, colour = combination)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
  geom_path(linewidth = 1.0) +
  scale_colour_manual(values = combo_colors, name = NULL, labels = setNames(combo_roc_summary$ROC_label, combo_roc_summary$combination)) +
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
    legend.text = element_text(size = 8.8),
    plot.margin = margin(5, 6, 5, 6)
  )
save_plot_safe(p_combo_roc, "Signature_score_combination_ROC_AUC_revision", width = 4.2, height = 3.4)

message("Done. Outputs written to: ", figure_dir)

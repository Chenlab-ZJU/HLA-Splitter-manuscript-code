# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final robust signature construction, scoring and cross-cohort evaluation.
# Original source: 32_revision_signature_cross_cohort_robust_selection.R
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
  library(pROC)
})

urine_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (!nzchar(urine_dir)) urine_dir <- file.path(getOption("hla.data_root"), "urine")
urine_dir <- gsub("\\\\", "/", urine_dir)
gse_dir <- Sys.getenv("GSE189536_PROJECT_DIR", unset = "")
if (!nzchar(gse_dir)) gse_dir <- file.path(getOption("hla.data_root"), "gse189536")
gse_dir <- gsub("\\\\", "/", gse_dir)

out_dir <- file.path(urine_dir, "Result", "Signature_score_analysis", "Cross_cohort_robust_selection")
fig_dir <- file.path(urine_dir, "Figures", "Signature_score_analysis", "Cross_cohort_robust_selection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

advanced_dir <- file.path(urine_dir, "Result", "Signature_score_analysis", "Advanced_optimization_GSE189536_validation")
deg_dir <- file.path(urine_dir, "Result", "Signature_score_analysis", "Differential_Expression")

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
signature_levels <- c("TCMR", "ABMR", "BK")
signature_colors <- c("TCMR" = "#699ECA", "ABMR" = "#4DAF4A", "BK" = "#F898CB")
group_colors <- c("NR" = "#E42320", "TCMR" = "#699ECA", "ABMR" = "#4DAF4A", "Mix" = "#FF8C00", "BK" = "#F898CB")

write_csv_safe <- function(x, path) utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")

read_deg <- function(file, signature_name, cluster_keep) {
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[, nzchar(names(x)), drop = FALSE]
  x %>%
    transmute(
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

safe_auc <- function(truth, score) {
  keep <- is.finite(score) & !is.na(truth)
  truth <- truth[keep]
  score <- score[keep]
  if (length(unique(truth)) < 2 || length(unique(score)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(truth, score, quiet = TRUE, direction = "<")))
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

save_plot <- function(plot, filename, width, height) {
  pdf_path <- file.path(fig_dir, paste0(filename, ".pdf"))
  png_path <- file.path(fig_dir, paste0(filename, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
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
      ordered_genes <- sig_table %>% arrange(.data[[rank_col]], p_val_adj, desc(avg_log2FC)) %>% pull(gene) %>% unique()
      for (n_genes in top_n_grid) {
        genes <- head(ordered_genes, min(n_genes, length(ordered_genes)))
        if (length(genes) < 3) next
        candidate_set <- paste(sig, method_name, length(genes), sep = "__")
        candidate_sets[[candidate_set]] <- genes
        candidate_info[[candidate_set]] <- tibble(candidate_set, signature = sig, ranking_method = method_name, n_genes = length(genes))
      }
    }
  }
  list(deg_table = deg_table, candidate_sets = candidate_sets, candidate_info = bind_rows(candidate_info))
}

message("Reconstructing candidate gene sets.")
deg_all <- bind_rows(
  read_deg(file.path(getOption("hla.data_root"), "urine", "scRNA_TCMR_DEG.csv"), "TCMR", "TCMR"),
  read_deg(file.path(getOption("hla.data_root"), "urine", "scRNA_ABMR_DEG.csv"), "ABMR", "ABMR"),
  read_deg(file.path(getOption("hla.data_root"), "urine", "scRNA_BK_DEG.csv"), "BK", "BK")
) %>%
  filter(is.na(p_val_adj) | p_val_adj < 0.05, avg_log2FC >= 0.5) %>%
  group_by(signature, gene) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Read the saved object or source table used by this analysis.
sample_gene_auc <- read.csv(file.path(advanced_dir, "Advanced_signature_single_gene_sample_AUC_urine.csv"), stringsAsFactors = FALSE)
candidate_pack <- build_candidate_sets(deg_all, sample_gene_auc)
urine_perf <- read.csv(file.path(advanced_dir, "Advanced_signature_gene_set_candidates_urine_AUC.csv"), stringsAsFactors = FALSE)

message("Scoring candidate gene sets in GSE189536.")
gse_object <- file.path(getOption("hla.data_root"), "gse189536", "Immune_filtered.rds")
gse <- readRDS(gse_object)
DefaultAssay(gse) <- "RNA"
gse_expr <- GetAssayData(gse, assay = "RNA", layer = "data")
gse_meta <- gse@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(cell, Sample = as.character(biopsy), Group = as.character(sample_type)) %>%
  filter(!is.na(Sample), Sample != "", !is.na(Group), Group %in% group_levels)

gse_sets <- lapply(candidate_pack$candidate_sets, function(genes) intersect(genes, rownames(gse_expr)))
gse_sets <- gse_sets[lengths(gse_sets) >= 3]
gse_info <- candidate_pack$candidate_info %>% filter(candidate_set %in% names(gse_sets))
gse_rankings <- AUCell_buildRankings(gse_expr, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE)
gse_auc <- AUCell_calcAUC(gse_sets, gse_rankings, aucMaxRank = ceiling(nrow(gse_rankings) * 0.1), verbose = FALSE)
gse_mat <- t(as.matrix(getAUC(gse_auc)))
gse_candidate_sample <- as.data.frame(gse_mat) %>%
  rownames_to_column("cell") %>%
  left_join(gse_meta, by = "cell") %>%
  pivot_longer(cols = all_of(colnames(gse_mat)), names_to = "candidate_set", values_to = "auc_score") %>%
  group_by(candidate_set, Sample, Group) %>%
  summarise(mean_auc = mean(auc_score, na.rm = TRUE), n_cells = n(), .groups = "drop") %>%
  left_join(gse_info, by = "candidate_set")

gse_perf <- gse_candidate_sample %>%
  filter(signature %in% c("TCMR", "ABMR")) %>%
  group_by(candidate_set, signature, ranking_method, n_genes) %>%
  summarise(
    gse_auc = safe_auc(as.integer(as.character(Group) == signature[1]), mean_auc),
    gse_target_median = median(mean_auc[as.character(Group) == signature[1]], na.rm = TRUE),
    gse_other_median = median(mean_auc[as.character(Group) != signature[1]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(gse_median_delta = gse_target_median - gse_other_median)
# Export the analysis source data for figure reproduction.
write_csv_safe(gse_perf, file.path(out_dir, "Cross_cohort_candidate_GSE189536_AUC.csv"))

joint_perf <- urine_perf %>%
  left_join(gse_perf %>% select(candidate_set, gse_auc, gse_median_delta), by = "candidate_set") %>%
  mutate(
    gse_auc = ifelse(signature == "BK", NA_real_, gse_auc),
    gse_median_delta = ifelse(signature == "BK", NA_real_, gse_median_delta),
    robust_score = dplyr::case_when(
      signature == "BK" ~ sample_auc,
      !is.na(gse_auc) & sample_auc > 0.5 & gse_auc > 0.5 ~ 2 / (1 / sample_auc + 1 / gse_auc),
      !is.na(gse_auc) ~ pmin(sample_auc, gse_auc),
      TRUE ~ NA_real_
    )
  ) %>%
  arrange(signature, desc(robust_score), desc(sample_auc), desc(gse_auc), n_genes)
write_csv_safe(joint_perf, file.path(out_dir, "Cross_cohort_candidate_joint_AUC.csv"))

robust_selected <- bind_rows(lapply(signature_levels, function(sig) {
  x <- joint_perf %>% filter(signature == sig)
  if (sig %in% c("TCMR", "ABMR")) {
    filtered <- x %>% filter(sample_auc >= 0.65, gse_auc >= 0.55, median_delta > 0, gse_median_delta > 0)
    if (nrow(filtered) == 0) filtered <- x %>% filter(!is.na(gse_auc), median_delta > 0, gse_median_delta > 0)
    if (nrow(filtered) == 0) filtered <- x %>% filter(!is.na(gse_auc))
    filtered %>% arrange(desc(robust_score), desc(sample_auc), desc(gse_auc), n_genes) %>% slice(1)
  } else {
    x %>% filter(median_delta > 0) %>% arrange(desc(sample_auc), n_genes) %>% slice(1)
  }
}))
write_csv_safe(robust_selected, file.path(out_dir, "Cross_cohort_robust_selected_gene_sets.csv"))

selected_gene_table <- bind_rows(lapply(seq_len(nrow(robust_selected)), function(i) {
  sig <- robust_selected$signature[i]
  candidate_name <- robust_selected$candidate_set[i]
  genes <- candidate_pack$candidate_sets[[candidate_name]]
  candidate_pack$deg_table %>%
    filter(signature == sig, gene %in% genes) %>%
    mutate(
      selected_candidate_set = candidate_name,
      selected_ranking_method = robust_selected$ranking_method[i],
      selected_n_genes = robust_selected$n_genes[i],
      urine_auc = robust_selected$sample_auc[i],
      gse_auc = robust_selected$gse_auc[i],
      robust_score = robust_selected$robust_score[i],
      gene_rank = match(gene, genes)
    ) %>%
    arrange(gene_rank)
})) %>%
  mutate(signature = factor(signature, levels = signature_levels)) %>%
  arrange(signature, gene_rank)
write_csv_safe(selected_gene_table, file.path(getOption("hla.data_root"), "urine", "Cross_cohort_robust_selected_genes.csv"))

selected_sets <- split(selected_gene_table$gene, as.character(selected_gene_table$signature))
selected_sets <- selected_sets[signature_levels]
names(selected_sets) <- paste0(names(selected_sets), "_Robust_Signature_AUC")

score_selected <- function(expr, meta, rankings, dataset_label) {
  sets <- lapply(selected_sets, function(genes) intersect(genes, rownames(expr)))
  sets <- sets[lengths(sets) >= 3]
  auc_obj <- AUCell_calcAUC(sets, rankings, aucMaxRank = ceiling(nrow(rankings) * 0.1), verbose = FALSE)
  mat <- t(as.matrix(getAUC(auc_obj)))
  as.data.frame(mat) %>%
    rownames_to_column("cell") %>%
    left_join(meta, by = "cell") %>%
    pivot_longer(cols = all_of(colnames(mat)), names_to = "signature_auc", values_to = "score") %>%
    mutate(signature = sub("_Robust_Signature_AUC$", "", signature_auc), dataset = dataset_label) %>%
    group_by(dataset, signature, Sample, Group) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), n_cells = n(), .groups = "drop")
}

message("Scoring robust selected sets in urine and GSE189536.")
urine_object <- c(file.path(urine_dir, "Result", "Objects", "Urine_scRNA.rds"), file.path(getOption("hla.data_root"), "urine", "Urine_ALL.rds"))
urine_object <- urine_object[file.exists(urine_object)][1]
urine <- readRDS(urine_object)
DefaultAssay(urine) <- "RNA"
urine_expr <- GetAssayData(urine, assay = "RNA", layer = "data")
sample_col <- if ("Sample" %in% colnames(urine@meta.data)) "Sample" else "HLA_Splitter"
urine_meta <- urine@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(cell, Sample = as.character(.data[[sample_col]]), Group = as.character(Group)) %>%
  filter(!is.na(Sample), Sample != "", !is.na(Group), Group %in% group_levels)
urine_rankings <- AUCell_buildRankings(urine_expr, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE)
urine_score <- score_selected(urine_expr, urine_meta, urine_rankings, "Urine")
gse_score <- score_selected(gse_expr, gse_meta, gse_rankings, "GSE189536")
all_scores <- bind_rows(urine_score, gse_score) %>%
  mutate(signature = factor(signature, levels = signature_levels), Group = factor(Group, levels = group_levels))
write_csv_safe(all_scores, file.path(getOption("hla.data_root"), "urine", "Cross_cohort_robust_sample_scores.csv"))

roc_from_scores <- function(df, dataset_label, signatures) {
  bind_rows(lapply(signatures, function(sig) {
    tmp <- df %>% filter(dataset == dataset_label, signature == sig)
    roc_obj <- pROC::roc(as.integer(as.character(tmp$Group) == sig), tmp$mean_score, quiet = TRUE, direction = "<")
    tibble(dataset = dataset_label, signature = sig, auc = as.numeric(pROC::auc(roc_obj)), specificity = roc_obj$specificities * 100, sensitivity = roc_obj$sensitivities * 100)
  })) %>%
    mutate(signature = factor(signature, levels = signature_levels)) %>%
    arrange(signature, desc(specificity), sensitivity)
}
urine_roc <- roc_from_scores(all_scores, "Urine", signature_levels)
gse_roc <- roc_from_scores(all_scores, "GSE189536", c("TCMR", "ABMR"))
roc_df <- bind_rows(urine_roc, gse_roc)
roc_summary <- roc_df %>%
  group_by(dataset, signature) %>%
  summarise(AUC = first(auc), .groups = "drop") %>%
  mutate(ROC_label = paste0(signature, " AUC = ", sprintf("%.3f", AUC)))
write_csv_safe(roc_df, file.path(out_dir, "Cross_cohort_robust_ROC_curve.csv"))
write_csv_safe(roc_summary, file.path(out_dir, "Cross_cohort_robust_ROC_AUC.csv"))

plot_roc_one <- function(dataset_label, signatures, filename) {
  rdf <- roc_df %>% filter(dataset == dataset_label, signature %in% signatures)
  rsum <- roc_summary %>% filter(dataset == dataset_label, signature %in% signatures)
  p <- ggplot(rdf, aes(specificity, sensitivity, colour = signature)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
    geom_path(linewidth = 1.0) +
    scale_colour_manual(values = signature_colors, name = NULL, labels = setNames(rsum$ROC_label, rsum$signature)) +
    scale_x_reverse(limits = c(100, 0), breaks = c(100, 75, 50, 25, 0), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Specificity", y = "Sensitivity") +
    theme_pub(10.5) +
    theme(legend.position = c(0.58, 0.22), legend.background = element_blank(), legend.text = element_text(size = 9))
  save_plot(p, filename, 3.8, 3.4)
}
plot_roc_one("Urine", signature_levels, "Cross_cohort_robust_ROC_AUC_urine")
plot_roc_one("GSE189536", c("TCMR", "ABMR"), "Cross_cohort_robust_ROC_AUC_GSE189536")

format_p_value <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    p < 0.01 ~ formatC(p, format = "f", digits = 3),
    TRUE ~ formatC(p, format = "f", digits = 2)
  )
}

make_pairwise_target_stats <- function(df, signatures, groups) {
  bind_rows(lapply(signatures, function(sig) {
    tmp <- df %>%
      filter(as.character(signature) == sig, Group %in% groups) %>%
      mutate(Group = factor(as.character(Group), levels = groups))
    target_group <- sig
    other_groups <- setdiff(groups[groups %in% as.character(tmp$Group)], target_group)
    other_groups <- other_groups[vapply(other_groups, function(g) any(as.character(tmp$Group) == g), logical(1))]
    target_values <- tmp$mean_score[as.character(tmp$Group) == target_group]
    max_y <- max(tmp$mean_score, na.rm = TRUE)
    min_y <- min(tmp$mean_score, na.rm = TRUE)
    y_range <- max(max_y - min_y, 0.001)
    step_y <- max(0.135 * y_range, 0.00135)
    tick_height <- max(0.030 * y_range, 0.00035)
    bind_rows(lapply(seq_along(other_groups), function(i) {
      other_group <- other_groups[i]
      other_values <- tmp$mean_score[as.character(tmp$Group) == other_group]
      p <- if (length(target_values) > 0 && length(other_values) > 0 && length(unique(c(target_values, other_values))) > 1) {
        suppressWarnings(wilcox.test(target_values, other_values, exact = FALSE)$p.value)
      } else {
        NA_real_
      }
      x_target <- as.numeric(factor(target_group, levels = groups))
      x_other <- as.numeric(factor(other_group, levels = groups))
      tibble(
        dataset = unique(as.character(tmp$dataset)),
        signature = factor(sig, levels = signatures),
        target_group = target_group,
        compared_group = other_group,
        group1 = target_group,
        group2 = other_group,
        x = min(x_target, x_other),
        xend = max(x_target, x_other),
        y = max_y + step_y * i,
        y_text = max_y + step_y * i + tick_height * 0.75,
        tick_height = tick_height,
        p_value = p,
        p_label = format_p_value(p)
      )
    }))
  }))
}

plot_score_boxplot <- function(df, signatures, groups, filename, width, height, y_label) {
  df <- df %>%
    filter(as.character(signature) %in% signatures, Group %in% groups, is.finite(mean_score)) %>%
    mutate(
      signature = factor(as.character(signature), levels = signatures),
      Group = factor(as.character(Group), levels = groups)
    )
  stat_df <- make_pairwise_target_stats(df, signatures, groups)
  write_csv_safe(stat_df, file.path(out_dir, paste0(filename, "_wilcox.csv")))
  stat_df_plot <- stat_df %>%
    filter(!is.na(p_value), p_value < 0.1)
  p <- ggplot(df, aes(Group, mean_score, fill = Group, colour = Group)) +
    geom_boxplot(width = 0.56, outlier.shape = NA, alpha = 0.38, linewidth = 0.34) +
    geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.55, alpha = 0.88) +
    geom_segment(data = stat_df_plot, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.30, colour = "black") +
    geom_segment(data = stat_df_plot, aes(x = x, xend = x, y = y, yend = y - tick_height), inherit.aes = FALSE, linewidth = 0.30, colour = "black") +
    geom_segment(data = stat_df_plot, aes(x = xend, xend = xend, y = y, yend = y - tick_height), inherit.aes = FALSE, linewidth = 0.30, colour = "black") +
    geom_text(
      data = stat_df_plot,
      aes(x = (x + xend) / 2, y = y_text, label = p_label),
      inherit.aes = FALSE,
      size = 3.45,
      family = "Arial",
      colour = "black"
    ) +
    facet_wrap(~ signature, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = group_colors, guide = "none", drop = FALSE) +
    scale_colour_manual(values = group_colors, guide = "none", drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.035, 0.10))) +
    labs(x = NULL, y = y_label) +
    theme_pub(9.1) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 8.8),
      axis.title.y = element_text(size = 10.8),
      strip.text = element_text(size = 10.2),
      panel.spacing.x = unit(0.72, "lines")
    )
  save_plot(p, filename, width, height)
}

plot_score_boxplot(
  all_scores %>% filter(dataset == "Urine") %>% mutate(signature = droplevels(signature)),
  c("TCMR", "ABMR", "BK"),
  c("NR", "TCMR", "ABMR", "Mix", "BK"),
  "Cross_cohort_robust_sample_score_boxplots_urine",
  7.8,
  3.7,
  "Mean AUCell score per urine sample"
)

plot_score_boxplot(
  all_scores %>% filter(dataset == "GSE189536", signature %in% c("TCMR", "ABMR")) %>% mutate(signature = droplevels(signature)),
  c("TCMR", "ABMR"),
  c("NR", "TCMR", "ABMR", "Mix"),
  "Cross_cohort_robust_sample_score_boxplots_GSE189536",
  5.7,
  3.7,
  "Mean AUCell score per biopsy"
)

note <- tibble(note = c(
  "This analysis uses GSE189536 for cross-cohort robust selection among urine-derived DEG candidate signatures.",
  "Because GSE189536 participates in candidate selection, the GSE AUC here is not a fully independent validation estimate.",
  "GSE189536 contains no BK-positive biopsy samples; BK can only be optimized/evaluated in the urine cohort."
))
write_csv_safe(note, file.path(out_dir, "Cross_cohort_robust_selection_notes.csv"))

message("Done. Cross-cohort robust selection outputs: ", out_dir)

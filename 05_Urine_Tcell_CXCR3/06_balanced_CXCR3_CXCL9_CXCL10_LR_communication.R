#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis.
# Original source: 17_revision_balanced_CXCR3_CXCL9_CXCL10_LR_communication.R
# Figure/table scope: Fig5; FigS9; Supplementary Table 18; Methods CXCR3-axis analysis
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

pick_col <- function(meta, candidates) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (!length(hit)) stop("None of these metadata columns were found: ", paste(candidates, collapse = ", "))
  hit[1]
}

read_first_existing_rds <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop("No object file found. Checked: ", paste(paths, collapse = "; "))
  readRDS(hit[1])
}

safe_wilcox_p <- function(values, groups, min_groups = 2) {
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- droplevels(factor(groups[keep]))
  if (length(unique(groups)) != 2) return(NA_real_)
  if (any(table(groups) < min_groups)) return(NA_real_)
  wilcox.test(values ~ groups, exact = FALSE)$p.value
}

safe_cor <- function(x, y, min_n = 5) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < min_n) {
    return(tibble(n_samples = sum(keep), spearman_rho = NA_real_, p_value = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  tibble(n_samples = sum(keep), spearman_rho = unname(ct$estimate), p_value = ct$p.value)
}

theme_nc <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      axis.text = element_text(colour = "black", size = base_size),
      axis.title = element_text(colour = "black", size = base_size + 1),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.30, colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1.5),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.28, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(3, 3, 3, 3)
    )
}

save_panel <- function(plot, stem, width, height) {
  pdf_file <- file.path(panel_dir, paste0(stem, ".pdf"))
  png_file <- file.path(panel_dir, paste0(stem, ".png"))
  ggsave(pdf_file, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(png_file, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  invisible(c(pdf_file, png_file))
}

script_dir <- get_script_dir()
urine_dir <- getOption("hla.cohort_dir")
project_dir <- getOption("hla.project_root")
revision_dir <- getOption("hla.revision_root")
out_dir <- file.path(revision_dir, "revised figure")
table_dir <- file.path(out_dir, "tables")
panel_dir <- file.path(out_dir, "panels", "Balanced_CXCR3_CXCL9_CXCL10_LR_communication")
result_dir <- file.path(urine_dir, "Result", "Reviewer_Tcell_ABMR_TCMR")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

GROUP_LEVELS <- c("NR", "TCMR", "ABMR", "Mix", "BK")
GROUP_COLORS <- c(
  "NR" = "#E42320", "TCMR" = "#699ECA", "ABMR" = "#4DAF4A",
  "Mix" = "#FF8C00", "BK" = "#F898CB"
)
MIN_SOURCE_CELLS <- 10
MIN_TARGET_CELLS <- 5

message("Loading urine all-cell and T-cell objects...")
all_cells <- read_first_existing_rds(c(
  file.path(urine_dir, "Result", "Objects", "Urine_scRNA.rds"),
  file.path(urine_dir, "Result", "Objects", "Urine_Allcells.rds"),
  file.path(getOption("hla.data_root"), "urine", "Urine_ALL.rds"),
  file.path(urine_dir, "Result", "Objects", "Urine_ALL_Intergrate.rds")
))
tcell <- read_first_existing_rds(c(
  file.path(getOption("hla.data_root"), "urine", "Urine_Tcell.rds"),
  file.path(getOption("hla.data_root"), "urine", "Urine_Tcell.rds")
))
DefaultAssay(all_cells) <- "RNA"
DefaultAssay(tcell) <- "RNA"

all_major_col <- pick_col(all_cells@meta.data, c("Major_celltype", "Celltype", "revised_celltype", "celltype"))
all_sample_col <- pick_col(all_cells@meta.data, c("HLA_Splitter", "Sample.ID", "Sample", "orig.ident"))
tcell_type_col <- pick_col(tcell@meta.data, c("revised_celltype", "Celltype", "celltype", "Major_celltype"))
tcell_sample_col <- pick_col(tcell@meta.data, c("HLA_Splitter", "Sample.ID", "Sample", "orig.ident"))

if (!all(c("CXCL9", "CXCL10") %in% rownames(all_cells))) {
  stop("CXCL9/CXCL10 were not both found in the all-cell object.")
}
if (!"CXCR3" %in% rownames(tcell)) stop("CXCR3 was not found in the T-cell object.")

all_cells$Source_celltype <- as.character(all_cells@meta.data[[all_major_col]])
all_cells$Sample <- as.character(all_cells@meta.data[[all_sample_col]])
all_cells$Group <- factor(as.character(all_cells$Group), levels = GROUP_LEVELS)
all_cells <- subset(all_cells, subset = Group %in% GROUP_LEVELS)

tcell$Target_Tcelltype <- as.character(tcell@meta.data[[tcell_type_col]])
tcell$Sample <- as.character(tcell@meta.data[[tcell_sample_col]])
tcell$Group <- factor(as.character(tcell$Group), levels = GROUP_LEVELS)
tcell <- subset(tcell, subset = Group %in% GROUP_LEVELS)

message("Calculating sample-level ligand expression by source cell type...")
source_expr <- FetchData(all_cells, vars = c("Group", "Sample", "Source_celltype", "CXCL9", "CXCL10")) %>%
  tibble::rownames_to_column("barcode") %>%
  as_tibble() %>%
  pivot_longer(cols = c("CXCL9", "CXCL10"), names_to = "Ligand", values_to = "Ligand_expression") %>%
  group_by(Group, Sample, Source_celltype, Ligand) %>%
  summarise(
    n_source_cells = n(),
    mean_ligand_expression = mean(Ligand_expression, na.rm = TRUE),
    pct_ligand_positive = mean(Ligand_expression > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  filter(n_source_cells >= MIN_SOURCE_CELLS)

message("Calculating sample-level CXCR3 expression by T-cell subtype...")
target_expr <- FetchData(tcell, vars = c("Group", "Sample", "Target_Tcelltype", "CXCR3")) %>%
  tibble::rownames_to_column("barcode") %>%
  as_tibble() %>%
  group_by(Group, Sample, Target_Tcelltype) %>%
  summarise(
    n_target_Tcells = n(),
    mean_CXCR3_expression = mean(CXCR3, na.rm = TRUE),
    pct_CXCR3_positive = mean(CXCR3 > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  filter(n_target_Tcells >= MIN_TARGET_CELLS)

message("Building composition-balanced ligand-receptor communication scores...")
lr_scores <- source_expr %>%
  inner_join(target_expr, by = c("Group", "Sample"), relationship = "many-to-many") %>%
  mutate(
    Rejection_status = if_else(Group %in% c("TCMR", "ABMR", "Mix"), "Rejection", "Non-rejection"),
    LR_pair = paste(Ligand, "-> CXCR3"),
    expression_product_score = mean_ligand_expression * mean_CXCR3_expression,
    positive_fraction_product_score = (pct_ligand_positive / 100) * (pct_CXCR3_positive / 100),
    balanced_LR_score = expression_product_score * positive_fraction_product_score
  ) %>%
  select(
    Group, Rejection_status, Sample, Source_celltype, Target_Tcelltype, LR_pair,
    n_source_cells, n_target_Tcells,
    mean_ligand_expression, pct_ligand_positive,
    mean_CXCR3_expression, pct_CXCR3_positive,
    expression_product_score, positive_fraction_product_score, balanced_LR_score
  )

# Export the analysis source data for figure reproduction.
write_csv(lr_scores, file.path(result_dir, "Balanced_CXCR3_CXCL9_CXCL10_LR_scores_by_sample_source_target.csv"))
write_csv(lr_scores, file.path(table_dir, "Balanced_CXCR3_CXCL9_CXCL10_LR_scores_by_sample_source_target.csv"))

axis_sample_summary <- lr_scores %>%
  group_by(Group, Rejection_status, Sample, LR_pair) %>%
  summarise(
    n_source_target_combinations = n(),
    mean_balanced_LR_score = mean(balanced_LR_score, na.rm = TRUE),
    median_balanced_LR_score = median(balanced_LR_score, na.rm = TRUE),
    mean_expression_product_score = mean(expression_product_score, na.rm = TRUE),
    mean_positive_fraction_product_score = mean(positive_fraction_product_score, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(axis_sample_summary, file.path(table_dir, "Balanced_CXCR3_CXCL9_CXCL10_axis_sample_summary.csv"))

axis_group_tests <- axis_sample_summary %>%
  group_by(LR_pair) %>%
  summarise(
    n_non_rejection = n_distinct(Sample[Rejection_status == "Non-rejection"]),
    n_rejection = n_distinct(Sample[Rejection_status == "Rejection"]),
    median_non_rejection = median(mean_balanced_LR_score[Rejection_status == "Non-rejection"], na.rm = TRUE),
    median_rejection = median(mean_balanced_LR_score[Rejection_status == "Rejection"], na.rm = TRUE),
    rejection_minus_non_rejection = median_rejection - median_non_rejection,
    p_rejection_vs_nonrejection = safe_wilcox_p(mean_balanced_LR_score, Rejection_status),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_rejection_vs_nonrejection, method = "BH")) %>%
  arrange(p_adj, desc(rejection_minus_non_rejection))
write_csv(axis_group_tests, file.path(table_dir, "Balanced_CXCR3_CXCL9_CXCL10_rejection_vs_nonrejection_tests.csv"))

group_pair_summary <- axis_sample_summary %>%
  group_by(Group, LR_pair) %>%
  summarise(
    n_samples = n_distinct(Sample),
    median_balanced_LR_score = median(mean_balanced_LR_score, na.rm = TRUE),
    mean_balanced_LR_score = mean(mean_balanced_LR_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(Group = factor(GROUP_LEVELS, levels = GROUP_LEVELS), LR_pair, fill = list(
    n_samples = 0, median_balanced_LR_score = NA_real_, mean_balanced_LR_score = NA_real_
  )) %>%
  arrange(LR_pair, Group)
write_csv(group_pair_summary, file.path(table_dir, "Balanced_CXCR3_CXCL9_CXCL10_group_pair_summary.csv"))

group_pair_tests <- axis_sample_summary %>%
  group_by(LR_pair) %>%
  group_modify(~{
    dat <- .x
    groups <- levels(droplevels(dat$Group))
    pairs <- t(combn(groups, 2))
    bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
      g1 <- pairs[i, 1]
      g2 <- pairs[i, 2]
      tmp <- dat %>% filter(Group %in% c(g1, g2))
      tibble(
        Comparison = paste(g2, "vs", g1),
        Group_1 = g1,
        Group_2 = g2,
        n_Group_1 = n_distinct(tmp$Sample[tmp$Group == g1]),
        n_Group_2 = n_distinct(tmp$Sample[tmp$Group == g2]),
        median_Group_1 = median(tmp$mean_balanced_LR_score[tmp$Group == g1], na.rm = TRUE),
        median_Group_2 = median(tmp$mean_balanced_LR_score[tmp$Group == g2], na.rm = TRUE),
        Group_2_minus_Group_1 = median_Group_2 - median_Group_1,
        p_value = safe_wilcox_p(tmp$mean_balanced_LR_score, tmp$Group)
      )
    }))
  }) %>%
  ungroup() %>%
  group_by(LR_pair) %>%
  mutate(p_adj_within_pair = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(LR_pair, p_adj_within_pair, desc(abs(Group_2_minus_Group_1)))
write_csv(group_pair_tests, file.path(table_dir, "Balanced_CXCR3_CXCL9_CXCL10_pairwise_group_tests.csv"))

source_target_tests <- lr_scores %>%
  group_by(Source_celltype, Target_Tcelltype, LR_pair) %>%
  summarise(
    n_non_rejection = n_distinct(Sample[Rejection_status == "Non-rejection"]),
    n_rejection = n_distinct(Sample[Rejection_status == "Rejection"]),
    median_non_rejection = median(balanced_LR_score[Rejection_status == "Non-rejection"], na.rm = TRUE),
    median_rejection = median(balanced_LR_score[Rejection_status == "Rejection"], na.rm = TRUE),
    rejection_minus_non_rejection = median_rejection - median_non_rejection,
    p_rejection_vs_nonrejection = safe_wilcox_p(balanced_LR_score, Rejection_status),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_rejection_vs_nonrejection, method = "BH")) %>%
  arrange(p_adj, desc(rejection_minus_non_rejection))
write_csv(source_target_tests, file.path(table_dir, "Balanced_CXCR3_CXCL9_CXCL10_source_target_rejection_tests.csv"))

key_rejection_edges <- source_target_tests %>%
  filter(rejection_minus_non_rejection > 0) %>%
  arrange(p_adj, desc(rejection_minus_non_rejection)) %>%
  slice_head(n = 40)
write_csv(key_rejection_edges, file.path(table_dir, "Balanced_CXCR3_CXCL9_CXCL10_key_rejection_enhanced_source_target_edges.csv"))

source_target_cor <- lr_scores %>%
  group_by(Source_celltype, Target_Tcelltype, LR_pair) %>%
  group_modify(~{
    sample_df <- .x %>%
      group_by(Sample, Group, Rejection_status) %>%
      summarise(
        ligand_signal = mean(mean_ligand_expression, na.rm = TRUE),
        target_CXCR3_positive = mean(pct_CXCR3_positive, na.rm = TRUE),
        balanced_LR_score = mean(balanced_LR_score, na.rm = TRUE),
        .groups = "drop"
      )
    safe_cor(sample_df$ligand_signal, sample_df$target_CXCR3_positive)
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj, desc(spearman_rho))
write_csv(source_target_cor, file.path(table_dir, "Balanced_CXCR3_CXCL9_CXCL10_source_target_sample_correlations.csv"))

message("Drawing balanced communication panels...")
p_group <- axis_sample_summary %>%
  mutate(Group = factor(Group, levels = GROUP_LEVELS)) %>%
  ggplot(aes(Group, mean_balanced_LR_score, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 1.0, alpha = 0.75, colour = "black") +
  facet_wrap(~LR_pair, scales = "free_y") +
  scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
  labs(x = NULL, y = "Composition-balanced LR score") +
  theme_nc(8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p_group, "Balanced_CXCR3_CXCL9_CXCL10_group_boxplots", 5.2, 3.2)

top_sources <- lr_scores %>%
  group_by(Source_celltype) %>%
  summarise(total = sum(balanced_LR_score, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 8) %>%
  pull(Source_celltype)
top_targets <- lr_scores %>%
  group_by(Target_Tcelltype) %>%
  summarise(total = sum(balanced_LR_score, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 8) %>%
  pull(Target_Tcelltype)

heat_tab <- lr_scores %>%
  filter(Source_celltype %in% top_sources, Target_Tcelltype %in% top_targets) %>%
  group_by(Rejection_status, Source_celltype, Target_Tcelltype) %>%
  summarise(mean_balanced_LR_score = mean(balanced_LR_score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Rejection_status, values_from = mean_balanced_LR_score, values_fill = 0) %>%
  mutate(rejection_minus_nonrejection = Rejection - `Non-rejection`)

p_heat <- heat_tab %>%
  mutate(
    Source_celltype = factor(Source_celltype, levels = rev(top_sources)),
    Target_Tcelltype = factor(Target_Tcelltype, levels = top_targets)
  ) %>%
  ggplot(aes(Target_Tcelltype, Source_celltype, fill = rejection_minus_nonrejection)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#699ECA", mid = "white", high = "#E6550D", midpoint = 0, name = "Rejection -\nnon-rejection") +
  labs(x = "CXCR3 target T-cell subtype", y = "CXCL9/CXCL10 source cell type") +
  theme_nc(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(p_heat, "Balanced_CXCR3_CXCL9_CXCL10_rejection_minus_nonrejection_heatmap", 6.3, 4.6)

plot_edges <- key_rejection_edges %>%
  slice_head(n = 20) %>%
  mutate(edge = paste(Source_celltype, "->", Target_Tcelltype), edge = factor(edge, levels = rev(unique(edge))))
p_edges <- ggplot(plot_edges, aes(rejection_minus_non_rejection, edge, fill = LR_pair)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("CXCL9 -> CXCR3" = "#E6550D", "CXCL10 -> CXCR3" = "#756BB1")) +
  labs(x = "Increase of balanced LR score in rejection", y = NULL, fill = NULL) +
  theme_nc(7) +
  theme(legend.position = "top")
save_panel(p_edges, "Balanced_CXCR3_CXCL9_CXCL10_key_rejection_edges", 5.5, 4.2)

message("Done. Tables were written to: ", table_dir)
message("Panels were written to: ", panel_dir)

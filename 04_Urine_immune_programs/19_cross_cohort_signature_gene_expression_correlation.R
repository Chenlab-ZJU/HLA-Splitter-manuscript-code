# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final robust signature construction, scoring and cross-cohort evaluation.
# Original source: 35_cross_cohort_signature_gene_expression_correlation.R
# Figure/table scope: Fig4i; Fig4j; FigS8a; FigS8b; FigS8c; FigS8d; FigS8e
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
project_dir <- if (nzchar(env_project_dir)) gsub("\\\\", "/", env_project_dir) else file.path(getOption("hla.data_root"), "urine")
env_gse_dir <- Sys.getenv("GSE189536_PROJECT_DIR", unset = "")
gse_dir <- if (nzchar(env_gse_dir)) gsub("\\\\", "/", env_gse_dir) else file.path(getOption("hla.data_root"), "gse189536")

figure_dir <- file.path(project_dir, "Figures", "Signature_score_analysis")
result_dir <- file.path(project_dir, "Result", "Signature_score_analysis", "Cross_cohort_robust_selection")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

signature_levels <- c("TCMR", "ABMR", "BK")
validation_signatures <- c("TCMR", "ABMR")
signature_colors <- c(
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "BK" = "#F898CB"
)

priority_genes <- c(
  "CXCR3", "CXCL9", "CXCL10", "CXCL11", "CCL5", "IFNG", "TNF",
  "GZMB", "GZMK", "PRF1", "GNLY", "NKG7",
  "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "TAP1", "TAP2",
  "STAT1", "GBP1", "GBP5", "ISG15", "IFI44L", "MX1",
  "MS4A1", "CD79A", "MZB1", "JCHAIN", "IGHG1"
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
  DefaultAssay(obj) <- assay
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    return(GetAssayData(obj, assay = assay, layer = "data"))
  }
  GetAssayData(obj, assay = assay, slot = "data")
}

mean_expression_by_group <- function(obj, group_col, groups, genes, dataset_name) {
  expr <- get_assay_matrix(obj, assay = "RNA")
  present_genes <- intersect(genes, rownames(expr))
  group_vec <- as.character(obj@meta.data[[group_col]])
  bind_rows(lapply(groups, function(grp) {
    cells <- colnames(obj)[group_vec == grp]
    if (length(cells) == 0 || length(present_genes) == 0) return(tibble())
    tibble(
      dataset = dataset_name,
      Group = grp,
      gene = present_genes,
      mean_log_expression = as.numeric(Matrix::rowMeans(expr[present_genes, cells, drop = FALSE])),
      pct_expressing = as.numeric(Matrix::rowMeans(expr[present_genes, cells, drop = FALSE] > 0) * 100),
      n_cells = length(cells)
    )
  }))
}

cor_summary <- function(df, x_col, y_col) {
  df %>%
    group_by(signature) %>%
    summarise(
      n_genes = n(),
      pearson_r = if (n() >= 3) suppressWarnings(cor(.data[[x_col]], .data[[y_col]], method = "pearson", use = "complete.obs")) else NA_real_,
      pearson_p = if (n() >= 3) suppressWarnings(cor.test(.data[[x_col]], .data[[y_col]], method = "pearson")$p.value) else NA_real_,
      spearman_rho = if (n() >= 3) suppressWarnings(cor(.data[[x_col]], .data[[y_col]], method = "spearman", use = "complete.obs")) else NA_real_,
      spearman_p = if (n() >= 3) suppressWarnings(cor.test(.data[[x_col]], .data[[y_col]], method = "spearman", exact = FALSE)$p.value) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0(
        "n=", n_genes,
        "\nPearson r=", sprintf("%.2f", pearson_r),
        "\nSpearman rho=", sprintf("%.2f", spearman_rho)
      )
    )
}

add_gene_labels <- function(plot_df) {
  plot_df %>%
    group_by(signature) %>%
    mutate(
      rank_abs_delta = rank(-abs(urine_delta_vs_NR - gse_delta_vs_NR), ties.method = "first"),
      label_gene = gene %in% priority_genes | rank_abs_delta <= 3
    ) %>%
    ungroup()
}

gene_label_layer <- function(label_df, size = 2.35) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    return(ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = size,
      family = "Arial",
      min.segment.length = 0,
      segment.size = 0.18,
      segment.colour = "grey55",
      box.padding = 0.22,
      point.padding = 0.12,
      max.overlaps = Inf,
      seed = 20260826
    ))
  }
  geom_text(
    data = label_df,
    aes(label = gene),
    size = size,
    family = "Arial",
    check_overlap = TRUE,
    vjust = -0.45
  )
}

message("Reading augmented robust signature genes...")
gene_file <- file.path(getOption("hla.data_root"), "urine", "Cross_cohort_robust_augmented_rejection_gene_sets_used_revision.csv")
if (!file.exists(gene_file)) {
  gene_file <- file.path(getOption("hla.data_root"), "urine", "Cross_cohort_robust_selected_genes.csv")
}
# Read the saved object or source table used by this analysis.
gene_table <- read.csv(gene_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(signature %in% validation_signatures, !is.na(gene), gene != "") %>%
  distinct(signature, gene, .keep_all = TRUE)
signature_score_genes <- unique(gene_table$gene)

message("Reading urine and GSE189536 immune objects...")
urine_file <- file.path(getOption("hla.data_root"), "urine", "Urine_ALL.rds")
urine <- readRDS(urine_file)
gse_file <- file.path(getOption("hla.data_root"), "gse189536", "Immune_filtered.rds")
gse <- readRDS(gse_file)

urine_group_col <- "Group"
gse_group_col <- if ("sample_type" %in% colnames(gse@meta.data)) {
  "sample_type"
} else if ("Group" %in% colnames(gse@meta.data)) {
  "Group"
} else {
  stop("Cannot find group column in GSE object.")
}

message("Calculating group-level mean expression for all shared genes...")
groups_to_use <- c("NR", validation_signatures)
genes <- intersect(rownames(get_assay_matrix(urine, assay = "RNA")), rownames(get_assay_matrix(gse, assay = "RNA")))
urine_mean <- mean_expression_by_group(urine, urine_group_col, groups_to_use, genes, "Urine")
gse_mean <- mean_expression_by_group(gse, gse_group_col, groups_to_use, genes, "GSE189536")
# Export the analysis source data for figure reproduction.
write_csv_safe(bind_rows(urine_mean, gse_mean), file.path(result_dir, "Cross_cohort_signature_gene_group_mean_expression.csv"))

mean_joined <- urine_mean %>%
  select(Group, gene, urine_mean_log_expression = mean_log_expression, urine_pct_expressing = pct_expressing) %>%
  inner_join(
    gse_mean %>% select(Group, gene, gse_mean_log_expression = mean_log_expression, gse_pct_expressing = pct_expressing),
    by = c("Group", "gene")
  )

expr_wide <- tidyr::expand_grid(signature = validation_signatures, gene = genes) %>%
  left_join(mean_joined, by = c("signature" = "Group", "gene")) %>%
  left_join(
    gene_table %>% select(signature, gene, gene_set_source, functional_module),
    by = c("signature", "gene")
  ) %>%
  mutate(
    is_signature_score_gene = !is.na(gene_set_source),
    gene_set_source = dplyr::case_when(
      gene_set_source == "cross_cohort_robust_selected" ~ "Robust selected",
      gene_set_source == "curated_rejection_literature_and_project_markers" ~ "Curated rejection",
      TRUE ~ "All other genes"
    ),
    functional_module = ifelse(is_signature_score_gene, functional_module, "All other genes")
  ) %>%
  mutate(signature = factor(signature, levels = validation_signatures))

delta_wide <- expr_wide %>%
  select(signature, gene, gene_set_source, functional_module, is_signature_score_gene, urine_target = urine_mean_log_expression, gse_target = gse_mean_log_expression) %>%
  left_join(
    urine_mean %>% filter(Group == "NR") %>% select(gene, urine_NR = mean_log_expression),
    by = "gene"
  ) %>%
  left_join(
    gse_mean %>% filter(Group == "NR") %>% select(gene, gse_NR = mean_log_expression),
    by = "gene"
  ) %>%
  mutate(
    urine_delta_vs_NR = urine_target - urine_NR,
    gse_delta_vs_NR = gse_target - gse_NR,
    signature = factor(signature, levels = validation_signatures)
  ) %>%
  filter(is.finite(urine_delta_vs_NR), is.finite(gse_delta_vs_NR))

write_csv_safe(expr_wide, file.path(result_dir, "Cross_cohort_signature_gene_target_group_expression_correlation_table.csv"))
write_csv_safe(delta_wide, file.path(result_dir, "Cross_cohort_signature_gene_delta_vs_NR_correlation_table.csv"))

target_cor <- cor_summary(expr_wide, "urine_mean_log_expression", "gse_mean_log_expression")
delta_cor <- cor_summary(delta_wide, "urine_delta_vs_NR", "gse_delta_vs_NR")
write_csv_safe(target_cor, file.path(result_dir, "Cross_cohort_signature_gene_target_group_expression_correlation_summary.csv"))
write_csv_safe(delta_cor, file.path(result_dir, "Cross_cohort_signature_gene_delta_vs_NR_correlation_summary.csv"))

cor_diagnostics <- bind_rows(
  delta_wide %>%
    mutate(gene_subset = "All shared genes"),
  delta_wide %>%
    filter(is_signature_score_gene) %>%
    mutate(gene_subset = "Highlighted rejection genes"),
  delta_wide %>%
    filter(is_signature_score_gene, abs(urine_delta_vs_NR) > 0.02 | abs(gse_delta_vs_NR) > 0.02) %>%
    mutate(gene_subset = "Highlighted genes with detectable delta")
) %>%
  group_by(signature, gene_subset) %>%
  summarise(
    n_genes = n(),
    median_abs_urine_delta = median(abs(urine_delta_vs_NR), na.rm = TRUE),
    median_abs_gse_delta = median(abs(gse_delta_vs_NR), na.rm = TRUE),
    pearson_r = if (n() >= 3) suppressWarnings(cor(urine_delta_vs_NR, gse_delta_vs_NR, method = "pearson", use = "complete.obs")) else NA_real_,
    spearman_rho = if (n() >= 3) suppressWarnings(cor(urine_delta_vs_NR, gse_delta_vs_NR, method = "spearman", use = "complete.obs")) else NA_real_,
    concordant_direction_pct = mean(sign(urine_delta_vs_NR) == sign(gse_delta_vs_NR), na.rm = TRUE) * 100,
    .groups = "drop"
  )
write_csv_safe(cor_diagnostics, file.path(result_dir, "Cross_cohort_signature_gene_delta_vs_NR_correlation_diagnostics.csv"))

target_label_pos <- expr_wide %>%
  group_by(signature) %>%
  summarise(
    x = min(urine_mean_log_expression, na.rm = TRUE) + 0.04 * diff(range(urine_mean_log_expression, na.rm = TRUE)),
    y = max(gse_mean_log_expression, na.rm = TRUE) - 0.04 * diff(range(gse_mean_log_expression, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(target_cor, by = "signature")

delta_label_pos <- delta_wide %>%
  group_by(signature) %>%
  summarise(
    x = min(urine_delta_vs_NR, na.rm = TRUE) + 0.04 * diff(range(urine_delta_vs_NR, na.rm = TRUE)),
    y = max(gse_delta_vs_NR, na.rm = TRUE) - 0.04 * diff(range(gse_delta_vs_NR, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(delta_cor, by = "signature")

plot_target <- add_gene_labels(delta_wide %>%
  select(signature, gene, urine_delta_vs_NR, gse_delta_vs_NR) %>%
  right_join(expr_wide, by = c("signature", "gene")))

target_label_genes <- plot_target %>%
  mutate(mean_expr_rank_value = pmax(urine_mean_log_expression, gse_mean_log_expression, na.rm = TRUE)) %>%
  filter(is_signature_score_gene) %>%
  group_by(signature) %>%
  arrange(desc(gene %in% priority_genes), desc(mean_expr_rank_value), .by_group = TRUE) %>%
  ungroup()

p_target <- ggplot(plot_target, aes(urine_mean_log_expression, gse_mean_log_expression)) +
  geom_point(
    data = plot_target %>% filter(!is_signature_score_gene),
    colour = "grey82",
    size = 0.58,
    alpha = 0.42,
    stroke = 0
  ) +
  geom_point(
    data = plot_target %>% filter(is_signature_score_gene),
    fill = "#D73027",
    shape = 21,
    size = 2.65,
    alpha = 0.94,
    colour = "grey15",
    stroke = 0.22
  ) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.55, colour = "grey25") +
  gene_label_layer(target_label_genes, size = 2.05) +
  geom_text(
    data = target_label_pos,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.75,
    fontface = "bold",
    family = "Arial"
  ) +
  facet_wrap(~ signature, nrow = 1, scales = "free") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 2, sigma = 0.01)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 2, sigma = 0.01)) +
  labs(x = "Urine mean expression (log2 scale)", y = "GSE189536 mean expression (log2 scale)") +
  theme_pub(8.8) +
  theme(legend.position = "right", panel.spacing.x = unit(0.8, "lines"))
save_plot_safe(p_target, "Cross_cohort_signature_gene_target_group_expression_scatter", width = 7.6, height = 3.8)

plot_delta <- add_gene_labels(delta_wide)

delta_label_genes <- plot_delta %>%
  mutate(delta_rank_value = abs(urine_delta_vs_NR) + abs(gse_delta_vs_NR)) %>%
  filter(is_signature_score_gene) %>%
  group_by(signature) %>%
  arrange(desc(gene %in% priority_genes), desc(delta_rank_value), .by_group = TRUE) %>%
  ungroup()

p_delta <- ggplot(plot_delta, aes(urine_delta_vs_NR, gse_delta_vs_NR)) +
  geom_hline(yintercept = 0, colour = "grey82", linewidth = 0.28) +
  geom_vline(xintercept = 0, colour = "grey82", linewidth = 0.28) +
  geom_point(
    data = plot_delta %>% filter(!is_signature_score_gene),
    colour = "grey82",
    size = 0.58,
    alpha = 0.42,
    stroke = 0
  ) +
  geom_point(
    data = plot_delta %>% filter(is_signature_score_gene),
    fill = "#D73027",
    shape = 21,
    size = 2.65,
    alpha = 0.94,
    colour = "grey15",
    stroke = 0.22
  ) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.55, colour = "grey25") +
  gene_label_layer(delta_label_genes, size = 2.05) +
  geom_text(
    data = delta_label_pos,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.75,
    fontface = "bold",
    family = "Arial"
  ) +
  facet_wrap(~ signature, nrow = 1, scales = "free") +
  labs(x = "Urine disease vs NR expression", y = "GSE189536 disease vs NR expression") +
  theme_pub(8.8) +
  theme(legend.position = "right", panel.spacing.x = unit(0.8, "lines"))
save_plot_safe(p_delta, "Cross_cohort_signature_gene_delta_vs_NR_scatter", width = 7.6, height = 3.8)

message("Done. Correlation plots written to: ", figure_dir)

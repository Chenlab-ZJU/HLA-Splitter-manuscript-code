# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 10d_revision_myeloid_key_gene_expression_wilcox.R
# Figure/table scope: Fig2b; Fig2c; Fig2d; Fig2e; Fig2f; FigS4c; FigS4d; FigS4e; FigS5a; FigS5b; FigS5c; FigS5d; FigS5e; FigS5f
# Data cohort: 脑脊液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

project_dir <- Sys.getenv("CSF_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "\u8111\u810a\u6db2\u6837\u672c")
}
figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
group_colors <- c(
  "Metastases" = "#E42320",
  "Hemorrhage" = "#699ECA",
  "Encephalitis" = "#4DAF4A"
)

format_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "<0.001", signif(p, 2))
  )
}

theme_pub <- function(base_size = 9.5) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
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
}

# Read the saved object or source table used by this analysis.
sample_qc <- readRDS(file.path(result_dir, "CSF_scRNA.rds"))@meta.data %>%
  rownames_to_column("cell") %>%
  count(HLA_Splitter, Group, name = "n_cells") %>%
  rename(Sample = HLA_Splitter) %>%
  arrange(desc(n_cells))
# Export the analysis source data for figure reproduction.
write.csv(sample_qc, file.path(result_dir, "CSF_sample_HLA_Splitter_cell_count_qc.csv"), row.names = FALSE)
valid_samples <- sample_qc$Sample[sample_qc$n_cells > 100]

sample_expr <- read.csv(
  file.path(result_dir, "Myeloid_revision_key_gene_sample_expression.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  filter(Sample %in% valid_samples) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    gene = as.character(gene),
    mean_expression = as.numeric(mean_expression)
  ) %>%
  filter(!is.na(Group), !is.na(mean_expression))

key_plot_genes <- c("TREM2", "SPP1", "APOE", "C1QA", "CXCL9", "CXCL10", "IL1B", "HMOX1")
key_plot_genes <- intersect(key_plot_genes, unique(sample_expr$gene))

pairwise_wilcox <- sample_expr %>%
  filter(gene %in% key_plot_genes) %>%
  group_by(gene) %>%
  group_modify(function(.x, .y) {
    comparisons <- t(utils::combn(group_levels, 2))
    tibble(
      group1 = comparisons[, 1],
      group2 = comparisons[, 2],
      x = match(comparisons[, 1], group_levels),
      xend = match(comparisons[, 2], group_levels),
      p.value = apply(comparisons, 1, function(pair) {
        values_1 <- .x$mean_expression[.x$Group == pair[1]]
        values_2 <- .x$mean_expression[.x$Group == pair[2]]
        if (length(values_1) < 2 || length(values_2) < 2) return(NA_real_)
        stats::wilcox.test(values_1, values_2, exact = FALSE)$p.value
      })
    )
  }) %>%
  ungroup() %>%
  mutate(
    p.adjust = p.adjust(p.value, method = "BH"),
    label = format_p(p.value)
  )

write.csv(
  pairwise_wilcox,
  file.path(result_dir, "Myeloid_revision_key_gene_pairwise_wilcox_tests.csv"),
  row.names = FALSE
)

pairwise_wilcox_plot <- pairwise_wilcox %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  left_join(
    sample_expr %>%
      filter(gene %in% key_plot_genes) %>%
      group_by(gene) %>%
      summarise(
        y_max = max(mean_expression, na.rm = TRUE),
        y_min = min(mean_expression, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "gene"
  ) %>%
  group_by(gene) %>%
  arrange(p.value, .by_group = TRUE) %>%
  mutate(
    y_range = pmax(y_max - y_min, 0.08),
    y.position = y_max + row_number() * y_range * 0.18,
    tick.height = y_range * 0.035,
    label.y = y.position + y_range * 0.065
  ) %>%
  ungroup()

plot_data <- sample_expr %>%
  filter(gene %in% key_plot_genes) %>%
  mutate(
    gene = factor(gene, levels = key_plot_genes),
    Group = factor(as.character(Group), levels = group_levels)
  )

p_key_expr <- ggplot(plot_data, aes(Group, mean_expression, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.82, linewidth = 0.38) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 20260813),
    shape = 21,
    colour = "black",
    stroke = 0.20,
    size = 1.45,
    alpha = 0.90
  ) +
  geom_segment(
    data = pairwise_wilcox_plot,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.30
  ) +
  geom_segment(
    data = pairwise_wilcox_plot,
    aes(x = x, xend = x, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.30
  ) +
  geom_segment(
    data = pairwise_wilcox_plot,
    aes(x = xend, xend = xend, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.30
  ) +
  geom_text(
    data = pairwise_wilcox_plot,
    aes(x = (x + xend) / 2, y = label.y, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 3.25
  ) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  labs(x = NULL, y = "Mean expression per sample") +
  theme_pub(9.6) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.0),
    strip.text = element_text(size = 9.4),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

save_plot_safe(p_key_expr, "Myeloid_key_gene_expression_by_group", width = 8.2, height = 5.8)

message("Updated Myeloid_key_gene_expression_by_group with pairwise Wilcoxon labels (p.value < 0.1).")

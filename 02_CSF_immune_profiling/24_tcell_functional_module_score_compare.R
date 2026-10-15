# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 11c_revision_tcell_functional_module_score_compare.R
# Figure/table scope: Fig2b; Fig2c; Fig2d; Fig2e; Fig2f; FigS4c; FigS4d; FigS4e; FigS5a; FigS5b; FigS5c; FigS5d; FigS5e; FigS5f
# Data cohort: csf. Inputs use the compact English data layout.
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
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

project_dir <- Sys.getenv("CSF_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "csf")
}
figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
group_colors <- c("Metastases" = "#E42320", "Hemorrhage" = "#699ECA", "Encephalitis" = "#4DAF4A")

tcell_modules <- list(
  "Cytotoxic" = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TOX", "ICOS", "MAF"),
  "IFN/chemokine" = c("IFNG", "CXCR3", "CCL3", "CCL4", "ISG15", "IFI6", "IFIT3", "MX1"),
  "CD4 activation/Treg" = c("CD40LG", "CCR6", "AQP3", "FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "Naive/memory" = c("CCR7", "IL7R", "LEF1", "TCF7", "SELL", "MAL", "LTB"),
  "Antigen presentation" = c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DPA1", "CD74")
)

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "p<0.001", paste0("p=", signif(p, 2))))
}

theme_pub <- function(base_size = 9.2) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey94", colour = "grey82", linewidth = 0.35),
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

message("Loading CSF_Tcell.rds and calculating sample-level module scores...")
# Read the saved object or source table used by this analysis.
Tcell <- readRDS(file.path(getOption("hla.data_root"), "csf", "CSF_Tcell.rds"))
DefaultAssay(Tcell) <- "RNA"

meta <- Tcell@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  transmute(
    cell,
    Sample = as.character(HLA_Splitter),
    Group = factor(as.character(Group), levels = group_levels)
  ) %>%
  filter(!is.na(Sample), !is.na(Group))

module_score_cell <- lapply(names(tcell_modules), function(module_name) {
  genes <- intersect(tcell_modules[[module_name]], rownames(Tcell))
  if (length(genes) < 2) return(NULL)
  mat <- FetchData(Tcell, vars = genes)
  tibble(
    cell = rownames(mat),
    module = module_name,
    n_genes = length(genes),
    module_score = rowMeans(as.matrix(mat), na.rm = TRUE)
  )
}) %>%
  bind_rows() %>%
  left_join(meta, by = "cell") %>%
  filter(!is.na(Sample), !is.na(Group)) %>%
  mutate(module = factor(module, levels = names(tcell_modules)))

module_score_sample <- module_score_cell %>%
  group_by(Sample, Group, module) %>%
  summarise(
    n_cells = n(),
    n_genes = max(n_genes, na.rm = TRUE),
    module_score = mean(module_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Group = factor(as.character(Group), levels = group_levels), module = factor(module, levels = names(tcell_modules)))

# Export the analysis source data for figure reproduction.
write.csv(module_score_cell, file.path(result_dir, "Tcell_revision_functional_module_cell_scores.csv"), row.names = FALSE)
write.csv(module_score_sample, file.path(result_dir, "Tcell_revision_functional_module_sample_scores.csv"), row.names = FALSE)
write.csv(
  tibble(module = names(tcell_modules), genes = vapply(tcell_modules, paste, collapse = ", ", FUN.VALUE = character(1))),
  file.path(result_dir, "Tcell_revision_functional_module_gene_sets.csv"),
  row.names = FALSE
)

comparisons <- t(utils::combn(group_levels, 2))
module_wilcox <- module_score_sample %>%
  group_by(module) %>%
  group_modify(function(.x, .y) {
    y_max <- max(.x$module_score, na.rm = TRUE)
    y_min <- min(.x$module_score, na.rm = TRUE)
    y_range <- max(y_max - y_min, 0.08)
    tibble(
      group1 = comparisons[, 1],
      group2 = comparisons[, 2],
      x = match(comparisons[, 1], group_levels),
      xend = match(comparisons[, 2], group_levels),
      p.value = apply(comparisons, 1, function(pair) {
        v1 <- .x$module_score[.x$Group == pair[1]]
        v2 <- .x$module_score[.x$Group == pair[2]]
        if (length(v1) < 2 || length(v2) < 2 || dplyr::n_distinct(c(v1, v2)) < 2) return(NA_real_)
        stats::wilcox.test(v1, v2, exact = FALSE)$p.value
      }),
      y_max = y_max,
      y_range = y_range
    )
  }) %>%
  ungroup() %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH"), label = format_p(p.value))

module_wilcox_plot <- module_wilcox %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  group_by(module) %>%
  arrange(p.value, .by_group = TRUE) %>%
  mutate(
    y.position = y_max + row_number() * y_range * 0.16 + y_range * 0.04,
    tick.height = y_range * 0.035,
    label.y = y.position + y_range * 0.055
  ) %>%
  ungroup()

write.csv(module_wilcox, file.path(result_dir, "Tcell_revision_functional_module_wilcox_tests.csv"), row.names = FALSE)
write.csv(module_wilcox_plot, file.path(result_dir, "Tcell_revision_functional_module_wilcox_tests_p_lt_0.1.csv"), row.names = FALSE)

p_module_score <- ggplot(module_score_sample, aes(Group, module_score, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, colour = "black", linewidth = 0.42, alpha = 0.34) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 11),
    shape = 21,
    stroke = 0.52,
    size = 2.10,
    alpha = 0.92
  ) +
  geom_segment(
    data = module_wilcox_plot,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = module_wilcox_plot,
    aes(x = x, xend = x, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = module_wilcox_plot,
    aes(x = xend, xend = xend, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_text(
    data = module_wilcox_plot,
    aes(x = (x + xend) / 2, y = label.y, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 3.0,
    colour = "black"
  ) +
  facet_wrap(~module, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.24))) +
  labs(x = NULL, y = "Module score per sample") +
  theme_pub(10.2) +
  theme(
    strip.background = element_rect(fill = "grey94", colour = NA),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9.2),
    axis.text.y = element_text(size = 9.8),
    axis.title.y = element_text(size = 11.8, face = "bold"),
    strip.text = element_text(size = 9.2, face = "bold"),
    panel.spacing.x = unit(0.75, "lines"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

save_plot_safe(p_module_score, "Tcell_functional_module_score_compare_revision", width = 11.0, height = 3.4)
message("Updated Tcell_functional_module_score_compare_revision with heatmap-matched gene modules.")

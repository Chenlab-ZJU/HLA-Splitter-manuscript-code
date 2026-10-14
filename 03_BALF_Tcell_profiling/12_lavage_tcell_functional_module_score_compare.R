# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Final BALF figures from saved objects and signatures.
# Original source: 20_revision_lavage_tcell_functional_module_score_compare.R
# Figure/table scope: Fig3b; Fig3c; Fig3d; Fig3e; Fig3f; Fig3g; Fig3h; Fig3i; FigS6a; FigS6b; FigS6c; FigS6d; FigS6e; FigS6f; FigS6g; FigS6h
# Data cohort: 肺泡灌洗液. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

project_dir <- Sys.getenv("LAVAGE_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "\u80ba\u6ce1\u704c\u6d17\u6db2")
}
project_dir <- gsub("\\\\", "/", project_dir)
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("CAP", "LC", "ILD")
group_colors <- c("CAP" = "#699ECA", "LC" = "#E42320", "ILD" = "#4DAF4A")

tcell_modules <- list(
  "Antigen presentation" = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "HLA-A", "HLA-B", "HLA-C", "B2M"),
  "Cytotoxic" = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "CTSW", "CST7", "KLRD1", "KLRG1", "FGFBP2", "CCL5", "CD8A", "CD8B"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX", "ENTPD1", "LAYN", "CXCL13", "TNFRSF9", "MAF", "NR4A1", "NR4A2"),
  "Tissue resident" = c("ITGAE", "ITGA1", "CXCR6", "ZNF683", "CD69", "RGS1", "XCL1", "XCL2", "KLRB1", "KLRG1", "S1PR1")
)

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
  )
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
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
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
  invisible(c(pdf = target_pdf, png = target_png))
}

message("Loading lavage_Tcell.rds and calculating module scores...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(file.path(result_dir, "lavage_Tcell.rds"))
DefaultAssay(tcell) <- if ("RNA" %in% Assays(tcell)) "RNA" else DefaultAssay(tcell)
required_cols <- c("Group", "HLA_Splitter")
missing_cols <- setdiff(required_cols, colnames(tcell@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

meta <- tcell@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Sample = as.character(HLA_Splitter),
    Group = factor(as.character(Group), levels = group_levels)
  ) %>%
  filter(!is.na(Sample), !is.na(Group))

module_gene_table <- tibble(
  module = rep(names(tcell_modules), lengths(tcell_modules)),
  gene = unlist(tcell_modules, use.names = FALSE)
) %>%
  filter(gene %in% rownames(tcell)) %>%
  distinct(module, gene)
# Export the analysis source data for figure reproduction.
write.csv(module_gene_table, file.path(result_dir, "lavage_Tcell_functional_module_gene_sets_revision.csv"), row.names = FALSE)

module_score_cell <- lapply(names(tcell_modules), function(module_name) {
  genes <- intersect(tcell_modules[[module_name]], rownames(tcell))
  if (length(genes) < 2) return(NULL)
  expr <- FetchData(tcell, vars = genes)
  tibble(
    cell = rownames(expr),
    module = module_name,
    n_genes = length(genes),
    module_score = rowMeans(as.matrix(expr), na.rm = TRUE)
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
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    module = factor(as.character(module), levels = names(tcell_modules))
  )

write.csv(module_score_cell, file.path(result_dir, "lavage_Tcell_functional_module_cell_scores_revision.csv"), row.names = FALSE)
write.csv(module_score_sample, file.path(result_dir, "lavage_Tcell_functional_module_sample_scores_revision.csv"), row.names = FALSE)

comparisons <- as_tibble(t(utils::combn(group_levels, 2)), .name_repair = "minimal") %>%
  setNames(c("group1", "group2"))

wilcox_stats <- module_score_sample %>%
  group_by(module) %>%
  group_modify(function(.x, .y) {
    y_max <- max(.x$module_score, na.rm = TRUE)
    y_min <- min(.x$module_score, na.rm = TRUE)
    y_range <- max(y_max - y_min, 0.05)
    comparisons %>%
      mutate(
        x = match(group1, group_levels),
        xend = match(group2, group_levels),
        p.value = purrr::map2_dbl(group1, group2, function(g1, g2) {
          v1 <- .x$module_score[.x$Group == g1]
          v2 <- .x$module_score[.x$Group == g2]
          if (length(v1) < 2 || length(v2) < 2 || dplyr::n_distinct(c(v1, v2)) < 2) return(NA_real_)
          stats::wilcox.test(v1, v2, exact = FALSE)$p.value
        }),
        y_max = y_max,
        y_range = y_range
      )
  }) %>%
  ungroup() %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH"), label = format_p(p.value))

wilcox_plot <- wilcox_stats %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  group_by(module) %>%
  arrange(p.value, .by_group = TRUE) %>%
  mutate(
    y.position = y_max + row_number() * y_range * 0.16 + y_range * 0.05,
    tick.height = y_range * 0.035,
    label.y = y.position + y_range * 0.05
  ) %>%
  ungroup()

write.csv(wilcox_stats, file.path(result_dir, "lavage_Tcell_functional_module_wilcox_tests_revision.csv"), row.names = FALSE)
write.csv(wilcox_plot, file.path(result_dir, "lavage_Tcell_functional_module_wilcox_tests_p_lt_0.1_revision.csv"), row.names = FALSE)

p_module <- ggplot(module_score_sample, aes(Group, module_score, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.34, linewidth = 0.42, colour = "black") +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 11),
    shape = 21,
    stroke = 0.52,
    size = 2.10,
    alpha = 0.92
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = x, xend = x, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = xend, xend = xend, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_text(
    data = wilcox_plot,
    aes(x = (x + xend) / 2, y = label.y, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    size = 3.0,
    vjust = -0.30,
    colour = "black"
  ) +
  facet_wrap(~module, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "Module score per sample") +
  theme_pub(10.2) +
  theme(
    strip.background = element_rect(fill = "grey94", colour = NA),
    strip.text = element_text(size = 10.8, face = "bold"),
    axis.title.y = element_text(size = 11.8, face = "bold"),
    axis.text.x = element_text(size = 10.2, colour = "black"),
    axis.text.y = element_text(size = 9.8, colour = "black"),
    panel.spacing.x = unit(0.75, "lines"),
    panel.spacing.y = unit(0.75, "lines"),
    legend.position = "none",
    plot.margin = margin(3, 4, 3, 4)
  )

save_plot_safe(p_module, "Tcell_functional_module_score_compare_revision", width = 5.6, height = 5.5)
message("Saved Tcell_functional_module_score_compare_revision.pdf and .png")

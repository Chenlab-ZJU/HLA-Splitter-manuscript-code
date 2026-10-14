# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 15_revision_tcell_key_gene_violin_compare.R
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
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

env_project_dir <- Sys.getenv("CSF_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "15_revision_tcell_key_gene_violin_compare.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
group_colors <- c("Metastases" = "#E42320", "Hemorrhage" = "#699ECA", "Encephalitis" = "#4DAF4A")
tcell_subtype_levels <- c("Naive CD4", "Naive CD8", "ISG+ CD4", "Treg", "GZMK+ CD8", "GZMB+ CD8", "NKT-like")
tcell_subtype_colors <- c(
  "Naive CD4" = "#699ECA",
  "Naive CD8" = "#9FDAF7",
  "ISG+ CD4" = "#FF8C00",
  "Treg" = "#731A73",
  "GZMK+ CD8" = "#0076B9",
  "GZMB+ CD8" = "#3D505A",
  "NKT-like" = "#0098B2"
)
plot_genes <- c(
  "CXCR3", "GZMK", "NKG7", "PRF1", "IFNG", "PDCD1",
  "ISG15", "GZMB", "GNLY", "LAG3", "CTLA4", "HLA-DRA"
)

read_csv_safe <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    p < 0.1 ~ paste0("p=", formatC(p, format = "f", digits = 2)),
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
  )
}

theme_pub <- function(base_size = 8.8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.margin = margin(4, 6, 4, 8)
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

tcell_rds <- file.path(result_dir, "CSF_Tcell.rds")
# Read the saved object or source table used by this analysis.
Tcell <- readRDS(tcell_rds)
DefaultAssay(Tcell) <- if ("RNA" %in% Assays(Tcell)) "RNA" else DefaultAssay(Tcell)
if (!"revised_celltype" %in% colnames(Tcell@meta.data)) {
  stop("CSF_Tcell.rds does not contain revised_celltype. Please provide the manually revised T cell subtype annotation field.")
}
Tcell$T_celltype <- factor(as.character(Tcell$revised_celltype), levels = tcell_subtype_levels)
Tcell <- subset(Tcell, subset = !is.na(T_celltype))
tcell_subtype_colors_present <- tcell_subtype_colors[levels(droplevels(Tcell$T_celltype))]
# Export the analysis source data for figure reproduction.
write_csv_safe(
  data.frame(T_celltype = names(tcell_subtype_colors_present), color = unname(tcell_subtype_colors_present)),
  file.path(result_dir, "Tcell_revision_key_gene_subtype_color_map.csv")
)
plot_genes <- intersect(plot_genes, rownames(Tcell))
if (length(plot_genes) == 0) {
  stop("None of the requested genes were found in the T cell object.")
}

cell_meta <- Tcell@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Sample = as.character(HLA_Splitter),
    Group = factor(as.character(Group), levels = group_levels),
    T_celltype = factor(as.character(T_celltype), levels = tcell_subtype_levels)
  ) %>%
  filter(!is.na(Sample), !is.na(Group), !is.na(T_celltype))

expr_df <- FetchData(Tcell, vars = plot_genes) %>%
  as_tibble(rownames = "cell") %>%
  inner_join(cell_meta, by = "cell") %>%
  pivot_longer(cols = all_of(plot_genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = plot_genes))

sample_expr <- expr_df %>%
  group_by(Sample, Group, gene) %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE),
    pct_expressing = mean(expression > 0, na.rm = TRUE) * 100,
    n_cells = dplyr::n(),
    .groups = "drop"
  )
write_csv_safe(sample_expr, file.path(result_dir, "Tcell_revision_key_gene_violin_sample_expression.csv"))

subtype_sample_expr <- expr_df %>%
  group_by(Sample, Group, T_celltype, gene) %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE),
    pct_expressing = mean(expression > 0, na.rm = TRUE) * 100,
    n_cells = dplyr::n(),
    .groups = "drop"
  )
write_csv_safe(subtype_sample_expr, file.path(result_dir, "Tcell_revision_key_gene_violin_subtype_sample_expression.csv"))

kw_stats <- sample_expr %>%
  group_by(gene) %>%
  summarise(
    p.value = ifelse(n_distinct(Group) >= 2 && n_distinct(mean_expression) > 1, kruskal.test(mean_expression ~ Group)$p.value, NA_real_),
    y_max = max(mean_expression, na.rm = TRUE),
    y_min = min(mean_expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    p.adjust = p.adjust(p.value, method = "BH"),
    label = paste0("KW ", format_p(p.value)),
    y_range = pmax(y_max - y_min, 0.08),
    label.y = y_max + y_range * 0.16
  )
write_csv_safe(kw_stats, file.path(result_dir, "Tcell_revision_key_gene_violin_KW_tests.csv"))

selected_plot_genes <- kw_stats %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  arrange(p.value, match(as.character(gene), plot_genes)) %>%
  pull(gene) %>%
  as.character()
if (length(selected_plot_genes) == 0) {
  stop("No T cell key gene has KW p.value < 0.1; no differential gene panel to plot.")
}
if (length(selected_plot_genes) < 6) {
  selected_plot_genes <- c(
    selected_plot_genes,
    setdiff(plot_genes, selected_plot_genes)
  )
}
selected_plot_genes <- selected_plot_genes[seq_len(min(6, length(selected_plot_genes)))]
write_csv_safe(
  tibble(gene = selected_plot_genes),
  file.path(result_dir, "Tcell_revision_key_gene_boxplot_selected_genes.csv")
)

sample_expr <- sample_expr %>%
  filter(as.character(gene) %in% selected_plot_genes) %>%
  mutate(gene = factor(as.character(gene), levels = selected_plot_genes))
expr_df <- expr_df %>%
  filter(as.character(gene) %in% selected_plot_genes) %>%
  mutate(gene = factor(as.character(gene), levels = selected_plot_genes))
kw_stats <- kw_stats %>%
  filter(as.character(gene) %in% selected_plot_genes) %>%
  mutate(gene = factor(as.character(gene), levels = selected_plot_genes))

comparisons <- as_tibble(t(utils::combn(group_levels, 2)), .name_repair = "minimal") %>%
  setNames(c("group1", "group2"))

wilcox_stats <- sample_expr %>%
  group_by(gene) %>%
  group_modify(function(.x, .y) {
    y_max <- max(.x$mean_expression, na.rm = TRUE)
    y_min <- min(.x$mean_expression, na.rm = TRUE)
    y_range <- max(y_max - y_min, 0.08)
    comparisons %>%
      mutate(
        x = match(group1, group_levels),
        xend = match(group2, group_levels),
        p.value = purrr::map2_dbl(group1, group2, function(g1, g2) {
          v1 <- .x$mean_expression[.x$Group == g1]
          v2 <- .x$mean_expression[.x$Group == g2]
          if (length(v1) < 2 || length(v2) < 2 || dplyr::n_distinct(c(v1, v2)) < 2) return(NA_real_)
          stats::wilcox.test(v1, v2, exact = FALSE)$p.value
        }),
        y_max = y_max,
        y_range = y_range
      )
  }) %>%
  ungroup() %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH"), label = format_p(p.value))
write_csv_safe(wilcox_stats, file.path(result_dir, "Tcell_revision_key_gene_violin_wilcox_tests.csv"))

wilcox_plot <- wilcox_stats %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  group_by(gene) %>%
  arrange(p.value, .by_group = TRUE) %>%
  mutate(
    y.position = y_max + row_number() * y_range * 0.18 + y_range * 0.09,
    tick.height = y_range * 0.04,
    label.y = y.position + y_range * 0.05
  ) %>%
  ungroup()

p_violin <- ggplot(sample_expr, aes(Group, mean_expression, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.36, alpha = 0.82) +
  geom_point(
    data = sample_expr,
    aes(Group, mean_expression, fill = Group),
    position = position_jitter(width = 0.08, height = 0, seed = 20260816),
    shape = 21,
    colour = "black",
    stroke = 0.18,
    size = 1.35,
    alpha = 0.88
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.28
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = x, xend = x, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.28
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = xend, xend = xend, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.28
  ) +
  geom_text(
    data = wilcox_plot,
    aes(x = (x + xend) / 2, y = label.y, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 2.35
  ) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.28))) +
  labs(x = NULL, y = "Mean expression per sample") +
  theme_pub(8.8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 7.0),
    axis.text.y = element_text(size = 7.8),
    axis.title.y = element_text(size = 9.0),
    strip.text = element_text(size = 8.4),
    legend.position = "top",
    legend.text = element_text(size = 8.0),
    legend.title = element_text(size = 8.4, face = "bold"),
    plot.margin = margin(4, 6, 6, 14)
  )

save_plot_safe(p_violin, "Tcell_key_gene_expression_violin_revision", width = 6.4, height = 4.9)
message("Saved T cell key gene expression boxplot.")

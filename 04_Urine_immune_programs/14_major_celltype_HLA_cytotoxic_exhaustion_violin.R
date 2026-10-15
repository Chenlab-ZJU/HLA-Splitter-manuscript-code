# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final urine immune-cell composition, module and demultiplexing figures.
# Original source: 29_revision_major_celltype_HLA_cytotoxic_exhaustion_violin.R
# Figure/table scope: Fig4b; Fig4c; Fig4d; Fig4e; Fig4f; Fig4g; Fig4h; FigS7a; FigS7b; FigS7c; FigS7d; FigS7e
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
  library(tibble)
  library(ggplot2)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "29_revision_major_celltype_HLA_cytotoxic_exhaustion_violin.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures", "scRNA", "Major_celltype_Module_Score")
result_dir <- file.path(project_dir, "Result", "scRNA", "Major_celltype_Module_Score")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
plot_groups <- c("TCMR", "ABMR")
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "Mix" = "#FF8C00",
  "BK" = "#F898CB"
)

module_colors <- c(
  "MHC I HLA" = "#8FB8D8",
  "MHC II HLA" = "#B7D8A5",
  "Cytotoxic" = "#D9A0A6",
  "Exhaustion" = "#B7A6D8"
)

major_levels <- c("CD4+ T", "CD8+ T", "Cycling T", "NK", "B cell", "Monocyte", "Macrophage", "Neutrophil", "cDC", "pDC")

score_modules <- list(
  "MHC I HLA" = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M"),
  "MHC II HLA" = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "CD74"),
  "Cytotoxic" = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TOX", "ICOS", "MAF")
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

theme_pub <- function(base_size = 8.8) {
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
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 6, 4, 6)
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

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "p=NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
  )
}

safe_wilcox <- function(v1, v2) {
  if (length(v1) < 2 || length(v2) < 2 || dplyr::n_distinct(c(v1, v2)) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(v1, v2, exact = FALSE)$p.value)
}

object_path <- file.path(getOption("hla.data_root"), "urine", "Urine_ALL.rds")
if (!file.exists(object_path)) {
  object_path <- file.path(project_dir, "Result", "Objects", "Urine_Allcells.rds")
}
if (!file.exists(object_path)) {
  stop("Urine all-cell object not found in Result/Objects.")
}

message("Loading urine all-cell object...")
# Read the saved object or source table used by this analysis.
scRNA <- readRDS(object_path)
DefaultAssay(scRNA) <- "RNA"

meta <- scRNA@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
    Group = as.character(Group),
    Major_celltype = trimws(as.character(Major_celltype))
  ) %>%
  filter(
    !is.na(Sample), !Sample %in% c("Undefined", "Doublets", "Doublet", ""),
    !is.na(Group), Group %in% plot_groups,
    !is.na(Major_celltype), Major_celltype != ""
  ) %>%
  mutate(
    Group = factor(Group, levels = plot_groups),
    Major_celltype = factor(Major_celltype, levels = c(major_levels[major_levels %in% unique(Major_celltype)], setdiff(sort(unique(Major_celltype)), major_levels)))
  )

available_modules <- tibble(
  module = rep(names(score_modules), lengths(score_modules)),
  gene = unlist(score_modules, use.names = FALSE)
) %>%
  filter(gene %in% rownames(scRNA)) %>%
  distinct(module, gene) %>%
  mutate(module = factor(module, levels = names(score_modules))) %>%
  arrange(module, gene)

if (nrow(available_modules) == 0) {
  stop("None of the selected HLA/function genes were found.")
}
# Export the analysis source data for figure reproduction.
write_csv_safe(available_modules, file.path(result_dir, "Urine_major_celltype_module_score_genes.csv"))

message("Calculating per-cell module scores...")
score_cell <- lapply(names(score_modules), function(module_name) {
  genes <- intersect(score_modules[[module_name]], rownames(scRNA))
  if (length(genes) < 2) return(NULL)
  mat <- FetchData(scRNA, vars = genes, cells = meta$cell)
  tibble(
    cell = rownames(mat),
    module = module_name,
    n_genes = length(genes),
    module_score = rowMeans(as.matrix(mat), na.rm = TRUE),
    detected_gene_count = rowSums(as.matrix(mat) > 0, na.rm = TRUE),
    positive = detected_gene_count > 0
  )
}) %>%
  bind_rows() %>%
  left_join(meta, by = "cell") %>%
  filter(!is.na(Group), !is.na(Major_celltype)) %>%
  mutate(module = factor(module, levels = names(score_modules)))

score_sample <- score_cell %>%
  group_by(Group, Sample, Major_celltype, module) %>%
  summarise(
    n_cells = n(),
    mean_module_score = mean(module_score, na.rm = TRUE),
    median_module_score = median(module_score, na.rm = TRUE),
    mean_detected_gene_count = mean(detected_gene_count, na.rm = TRUE),
    pct_positive = mean(positive, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(as.character(Group), levels = plot_groups),
    module = factor(as.character(module), levels = names(score_modules))
  )

celltype_summary <- score_cell %>%
  group_by(Group, Major_celltype, module) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(Sample),
    mean_module_score = mean(module_score, na.rm = TRUE),
    median_module_score = median(module_score, na.rm = TRUE),
    mean_detected_gene_count = mean(detected_gene_count, na.rm = TRUE),
    pct_positive = mean(positive, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(module, Group, desc(mean_module_score))

sample_stats <- score_sample %>%
  group_by(Major_celltype, module) %>%
  summarise(
    n_TCMR_samples = sum(Group == "TCMR"),
    n_ABMR_samples = sum(Group == "ABMR"),
    mean_TCMR = mean(mean_module_score[Group == "TCMR"], na.rm = TRUE),
    mean_ABMR = mean(mean_module_score[Group == "ABMR"], na.rm = TRUE),
    delta_ABMR_vs_TCMR = mean_ABMR - mean_TCMR,
    p_value = safe_wilcox(mean_module_score[Group == "TCMR"], mean_module_score[Group == "ABMR"]),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(module, p_adj, p_value)

write_csv_safe(score_cell, file.path(result_dir, "Urine_TCMR_ABMR_major_celltype_module_cell_scores.csv"))
write_csv_safe(score_sample, file.path(result_dir, "Urine_TCMR_ABMR_major_celltype_module_sample_scores.csv"))
write_csv_safe(celltype_summary, file.path(result_dir, "Urine_TCMR_ABMR_major_celltype_module_summary.csv"))
write_csv_safe(sample_stats, file.path(result_dir, "Urine_TCMR_ABMR_major_celltype_module_sample_wilcox.csv"))

rank_df <- score_cell %>%
  group_by(module, Major_celltype) %>%
  summarise(rank_score = mean(module_score, na.rm = TRUE), .groups = "drop") %>%
  arrange(module, desc(rank_score)) %>%
  group_by(module) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(y_key = paste0(sprintf("%02d", rank), "__", Major_celltype, "__", module))

y_levels <- rank_df %>%
  arrange(module, desc(rank)) %>%
  pull(y_key)

plot_df <- score_cell %>%
  left_join(rank_df, by = c("module", "Major_celltype")) %>%
  mutate(
    y_key = factor(y_key, levels = y_levels)
  )

p_violin <- ggplot(plot_df, aes(module_score, y_key, fill = module)) +
  geom_violin(scale = "width", trim = TRUE, colour = "grey28", linewidth = 0.25, alpha = 0.78) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", colour = "grey25", linewidth = 0.25, alpha = 0.82) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 1.05, stroke = 0.20, fill = "black", colour = "black") +
  facet_grid(module ~ Group, scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("__[^_]+$", "", sub("^\\d+__", "", x))) +
  scale_fill_manual(values = module_colors, drop = FALSE, name = "Module") +
  labs(x = "Module score per cell", y = NULL) +
  theme_pub(8.8) +
  theme(
    axis.text.y = element_text(size = 7.3),
    axis.text.x = element_text(size = 7.4),
    axis.title.x = element_text(size = 9.2),
    strip.text = element_text(size = 8.6),
    legend.position = "top",
    legend.text = element_text(size = 7.8),
    legend.title = element_text(size = 8.0, face = "bold"),
    panel.spacing.x = unit(0.60, "lines"),
    panel.spacing.y = unit(0.45, "lines")
  )
save_plot_safe(p_violin, "Urine_TCMR_ABMR_major_celltype_HLA_cytotoxic_exhaustion_score_violin", width = 7.6, height = 7.2)

sample_plot_df <- score_sample %>%
  left_join(rank_df, by = c("module", "Major_celltype")) %>%
  mutate(y_key = factor(y_key, levels = y_levels))

write_csv_safe(sample_plot_df, file.path(result_dir, "Urine_TCMR_ABMR_major_celltype_module_sample_score_boxplot_source.csv"))

p_sample <- ggplot(sample_plot_df, aes(y_key, mean_module_score, fill = Group)) +
  geom_boxplot(
    width = 0.56,
    outlier.shape = NA,
    colour = "grey24",
    linewidth = 0.32,
    alpha = 0.78
  ) +
  geom_point(
    aes(size = n_cells),
    position = position_jitter(width = 0.08, height = 0, seed = 20260826),
    shape = 21,
    colour = "grey18",
    stroke = 0.18,
    alpha = 0.86
  ) +
  facet_grid(Group ~ module, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) sub("__[^_]+$", "", sub("^\\d+__", "", x))) +
  scale_fill_manual(values = group_colors[plot_groups], drop = FALSE, name = "Group") +
  scale_size_continuous(range = c(0.8, 2.4), name = "Cells/sample") +
  labs(x = NULL, y = "Mean module score per sample") +
  theme_pub(9.0) +
  theme(
    axis.text.x = element_text(size = 7.2, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7.8),
    axis.title.y = element_text(size = 9.3, face = "bold", margin = margin(r = 5)),
    strip.text = element_text(size = 8.4, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 8.0),
    legend.title = element_text(size = 8.4, face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.spacing.x = unit(0.62, "lines"),
    panel.spacing.y = unit(0.48, "lines")
  )
save_plot_safe(p_sample, "Urine_TCMR_ABMR_major_celltype_HLA_cytotoxic_exhaustion_sample_score_boxplot", width = 9.2, height = 5.2)

message("Saved TCMR/ABMR major-celltype HLA, cytotoxic and exhaustion module score figures.")

# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Final urine T-cell, CXCR3 and expanded-clone figures.
# Original source: 28_revision_tcell_group_DEG_module_heatmap_score.R
# Figure/table scope: Fig5a; Fig5b; Fig5c; Fig5d; Fig5e; Fig5f; Fig5g; Fig5h; Fig5i; FigS9a; FigS9b; FigS9c; FigS9d; FigS9e; FigS9f; FigS9g
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
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "28_revision_tcell_group_DEG_module_heatmap_score.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures", "Tcell", "Group_DEG_Function")
result_dir <- file.path(project_dir, "Result", "Tcell", "Group_DEG_Function")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "Mix" = "#FF8C00",
  "BK" = "#F898CB"
)

tcell_marker_modules <- list(
  "Naive/memory" = c("CCR7", "IL7R", "LEF1", "TCF7", "SELL", "MAL", "LTB"),
  "Activation/Treg" = c("CD40LG", "CCR6", "AQP3", "FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "Cytotoxic" = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TOX", "ICOS", "MAF"),
  "IFN/chemokine" = c("IFNG", "CXCR3", "CCL3", "CCL4", "ISG15", "IFI6", "IFIT3", "MX1"),
  "Antigen presentation" = c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DPA1", "CD74")
)

tcell_score_modules <- tcell_marker_modules

tcell_celltype_levels <- c(
  "Naive CD4", "CD4 Tem", "CD4 Tex", "Treg",
  "Naive CD8", "CXCR4+ CD8", "GZMB+ CD8", "ISG+ CD8", "Cycling T"
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

format_p <- function(p, prefix = "KW p") {
  dplyr::case_when(
    is.na(p) ~ paste0(prefix, "=NA"),
    p < 0.001 ~ paste0(prefix, "<0.001"),
    p < 0.01 ~ paste0(prefix, "=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0(prefix, "=", formatC(p, format = "f", digits = 2))
  )
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

safe_kw <- function(df, value_col = "module_score", group_col = "Group") {
  vals <- df[[value_col]]
  grp <- droplevels(factor(df[[group_col]]))
  keep <- is.finite(vals) & !is.na(grp)
  vals <- vals[keep]
  grp <- droplevels(grp[keep])
  if (nlevels(grp) < 2 || dplyr::n_distinct(vals) < 2) return(NA_real_)
  suppressWarnings(kruskal.test(vals ~ grp)$p.value)
}

safe_wilcox_two_groups <- function(df, value_col = "module_score", group_col = "Group", g1 = "TCMR", g2 = "ABMR") {
  vals <- df[[value_col]]
  grp <- as.character(df[[group_col]])
  keep <- is.finite(vals) & grp %in% c(g1, g2)
  vals <- vals[keep]
  grp <- factor(grp[keep], levels = c(g1, g2))
  if (nlevels(droplevels(grp)) < 2 || any(table(grp) < 2) || dplyr::n_distinct(vals) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(vals ~ grp, exact = FALSE)$p.value)
}

pairwise_wilcox <- function(df, value_col = "module_score", facet_col = "module") {
  comparisons <- as_tibble(t(utils::combn(group_levels, 2)), .name_repair = "minimal") %>%
    setNames(c("group1", "group2"))
  df %>%
    group_by(.data[[facet_col]]) %>%
    group_modify(function(.x, .y) {
      comparisons %>%
        mutate(
          p.value = mapply(function(g1, g2) {
            v1 <- .x[[value_col]][.x$Group == g1]
            v2 <- .x[[value_col]][.x$Group == g2]
            if (length(v1) < 2 || length(v2) < 2 || dplyr::n_distinct(c(v1, v2)) < 2) return(NA_real_)
            suppressWarnings(wilcox.test(v1, v2, exact = FALSE)$p.value)
          }, group1, group2)
        )
    }) %>%
    ungroup() %>%
    mutate(p.adjust = p.adjust(p.value, method = "BH"))
}

object_path <- file.path(getOption("hla.data_root"), "urine", "Urine_Tcell.rds")
if (!file.exists(object_path)) {
  stop("Cannot find urine T-cell object: ", object_path)
}

# Read the saved object or source table used by this analysis.
Tcell <- readRDS(object_path)
DefaultAssay(Tcell) <- "RNA"

meta <- Tcell@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
    Group = as.character(Group),
    T_celltype = as.character(revised_celltype)
  ) %>%
  filter(!is.na(Sample), Sample != "", !is.na(Group), Group %in% group_levels) %>%
  mutate(Group = factor(Group, levels = group_levels))

available_groups <- intersect(group_levels, unique(as.character(meta$Group)))
if (length(available_groups) < 2) {
  stop("Need at least two clinical groups in the urine T-cell object.")
}

tcell_marker_df <- tibble(
  module = rep(names(tcell_marker_modules), lengths(tcell_marker_modules)),
  gene = unlist(tcell_marker_modules, use.names = FALSE)
) %>%
  filter(gene %in% rownames(Tcell)) %>%
  distinct(module, gene) %>%
  mutate(module = factor(module, levels = names(tcell_marker_modules))) %>%
  arrange(module, gene) %>%
  group_by(module) %>%
  mutate(module_gene = paste(module, gene, sep = "__")) %>%
  ungroup()

if (nrow(tcell_marker_df) == 0) {
  stop("None of the selected T-cell key genes were found in the urine T-cell object.")
}

deg_file <- file.path(getOption("hla.data_root"), "urine", "Tcell_Group_DEG.csv")
if (file.exists(deg_file)) {
  old_deg_raw <- read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
  old_deg_raw <- old_deg_raw[, nzchar(names(old_deg_raw)), drop = FALSE]
  old_deg <- old_deg_raw %>%
    transmute(
      gene = as.character(gene),
      enriched_group = as.character(cluster),
      avg_log2FC = suppressWarnings(as.numeric(avg_log2FC)),
      p_val_adj = suppressWarnings(as.numeric(p_val_adj))
    )
  tcell_marker_df <- tcell_marker_df %>%
    left_join(
      old_deg %>%
        filter(gene %in% tcell_marker_df$gene) %>%
        group_by(gene) %>%
        slice_min(order_by = p_val_adj, n = 1, with_ties = FALSE) %>%
        ungroup(),
      by = "gene"
    )
} else {
  tcell_marker_df <- tcell_marker_df %>%
    mutate(enriched_group = NA_character_, avg_log2FC = NA_real_, p_val_adj = NA_real_)
}

# Export the analysis source data for figure reproduction.
write_csv_safe(tcell_marker_df, file.path(result_dir, "Tcell_revision_heatmap_gene_modules.csv"))

heatmap_genes <- unique(tcell_marker_df$gene)
celltype_levels <- meta %>%
  filter(!is.na(T_celltype), T_celltype != "") %>%
  count(T_celltype, sort = TRUE) %>%
  pull(T_celltype)
celltype_levels <- c(tcell_celltype_levels[tcell_celltype_levels %in% celltype_levels], setdiff(celltype_levels, tcell_celltype_levels))

Tcell$T_celltype_revision <- factor(as.character(Tcell$revised_celltype), levels = celltype_levels)

celltype_avg <- AverageExpression(
  Tcell,
  assays = "RNA",
  features = heatmap_genes,
  group.by = "T_celltype_revision",
  slot = "data",
  verbose = FALSE
)$RNA %>%
  as.matrix()

celltype_avg <- celltype_avg[, intersect(celltype_levels, colnames(celltype_avg)), drop = FALSE]
celltype_z <- t(scale(t(celltype_avg)))
celltype_z[is.na(celltype_z)] <- 0
celltype_z <- pmax(pmin(celltype_z, 2.2), -2.2)

celltype_heat_df <- as.data.frame(celltype_z) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "T_celltype", values_to = "z_score") %>%
  left_join(tcell_marker_df, by = "gene", relationship = "many-to-many") %>%
  mutate(
    module = factor(as.character(module), levels = names(tcell_marker_modules)),
    module_gene = factor(module_gene, levels = rev(tcell_marker_df$module_gene)),
    T_celltype = factor(T_celltype, levels = celltype_levels)
  )
write_csv_safe(celltype_heat_df, file.path(result_dir, "Tcell_revision_celltype_heatmap_source.csv"))

p_celltype_heat <- ggplot(celltype_heat_df, aes(T_celltype, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_gradient2(
    low = "#6F99C3",
    mid = "#F7F7F7",
    high = "#B7192B",
    midpoint = 0,
    limits = c(-2.2, 2.2),
    oob = scales::squish,
    name = "Z-score"
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.3),
    axis.text.y = element_text(size = 7.8),
    strip.text.y = element_text(size = 7.8, angle = 0, hjust = 0),
    panel.grid = element_blank(),
    legend.position = "right"
  )
save_plot_safe(p_celltype_heat, "Tcell_celltype_heatmap_revision", width = 6.5, height = 7.0)

group_avg <- AverageExpression(
  Tcell,
  assays = "RNA",
  features = heatmap_genes,
  group.by = "Group",
  slot = "data",
  verbose = FALSE
)$RNA %>%
  as.matrix()

group_avg <- group_avg[, intersect(group_levels, colnames(group_avg)), drop = FALSE]
group_z <- t(scale(t(group_avg)))
group_z[is.na(group_z)] <- 0
group_z <- pmax(pmin(group_z, 2.2), -2.2)

group_heat_df <- as.data.frame(group_z) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Group", values_to = "z_score") %>%
  left_join(tcell_marker_df, by = "gene", relationship = "many-to-many") %>%
  mutate(
    module = factor(as.character(module), levels = names(tcell_marker_modules)),
    module_gene = factor(module_gene, levels = rev(tcell_marker_df$module_gene)),
    Group = factor(Group, levels = group_levels)
  )
write_csv_safe(group_heat_df, file.path(result_dir, "Tcell_revision_group_DEG_key_gene_heatmap_source.csv"))

p_group_heat <- ggplot(group_heat_df, aes(Group, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_gradient2(
    low = "#6F99C3",
    mid = "#F7F7F7",
    high = "#B7192B",
    midpoint = 0,
    limits = c(-2.2, 2.2),
    oob = scales::squish,
    name = "Z-score"
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.2),
    axis.text.y = element_text(size = 7.8),
    strip.text.y = element_text(size = 7.8, angle = 0, hjust = 0),
    panel.grid = element_blank(),
    legend.position = "right"
  )
save_plot_safe(p_group_heat, "Tcell_group_DEG_key_gene_heatmap_revision", width = 4.8, height = 7.0)

deg_group_heat_df <- group_heat_df %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
  mutate(
    enriched_group = factor(enriched_group, levels = group_levels),
    module = factor(as.character(module), levels = names(tcell_marker_modules))
  ) %>%
  arrange(module, enriched_group, desc(abs(avg_log2FC)), gene) %>%
  mutate(module_gene = paste(module, gene, sep = "__")) %>%
  mutate(module_gene = factor(module_gene, levels = rev(unique(module_gene))))

write_csv_safe(deg_group_heat_df, file.path(result_dir, "Tcell_revision_group_DEG_gene_module_heatmap_source.csv"))

p_group_deg_module_heat <- ggplot(deg_group_heat_df, aes(Group, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_gradient2(
    low = "#6F99C3",
    mid = "#F7F7F7",
    high = "#B7192B",
    midpoint = 0,
    limits = c(-2.2, 2.2),
    oob = scales::squish,
    name = "Z-score"
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub(9.4) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.5, face = "bold"),
    axis.text.y = element_text(size = 8.0),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text.y = element_text(size = 8.0, angle = 0, hjust = 0, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_plot_safe(p_group_deg_module_heat, "Tcell_group_DEG_gene_module_heatmap_revision", width = 4.9, height = 6.8)

module_score_cell <- lapply(names(tcell_score_modules), function(module_name) {
  genes <- intersect(tcell_score_modules[[module_name]], rownames(Tcell))
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
  left_join(meta %>% dplyr::select(cell, Sample, Group), by = "cell") %>%
  filter(!is.na(Sample), !is.na(Group))

module_score_sample <- module_score_cell %>%
  group_by(Sample, Group, module) %>%
  summarise(
    n_cells = n(),
    module_score = mean(module_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    module = factor(module, levels = names(tcell_score_modules))
  )

module_kw <- module_score_sample %>%
  group_by(module) %>%
  summarise(p.value = safe_kw(pick(everything())), .groups = "drop") %>%
  mutate(
    p.adjust = p.adjust(p.value, method = "BH"),
    label = format_p(p.value, prefix = "KW p")
  )

module_wilcox <- pairwise_wilcox(module_score_sample, "module_score", "module")
module_tcmr_abmr_wilcox <- module_score_sample %>%
  group_by(module) %>%
  summarise(
    group1 = "TCMR",
    group2 = "ABMR",
    p.value = safe_wilcox_two_groups(pick(everything()), "module_score", "Group", "TCMR", "ABMR"),
    .groups = "drop"
  ) %>%
  mutate(
    p.adjust = p.adjust(p.value, method = "BH"),
    label = format_p(p.value, prefix = "Wilcoxon p")
  )

write_csv_safe(module_score_cell, file.path(result_dir, "Tcell_revision_functional_module_cell_scores.csv"))
write_csv_safe(module_score_sample, file.path(result_dir, "Tcell_revision_functional_module_sample_scores.csv"))
write_csv_safe(module_kw, file.path(result_dir, "Tcell_revision_functional_module_KW_tests.csv"))
write_csv_safe(module_wilcox, file.path(result_dir, "Tcell_revision_functional_module_wilcox_tests.csv"))
write_csv_safe(module_tcmr_abmr_wilcox, file.path(result_dir, "Tcell_revision_functional_module_TCMR_vs_ABMR_wilcox_tests.csv"))

wilcox_plot <- module_tcmr_abmr_wilcox %>%
  left_join(
    module_score_sample %>%
      group_by(module) %>%
      summarise(
        y_min = min(module_score, na.rm = TRUE),
        y_max = max(module_score, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "module"
  ) %>%
  mutate(
    y_span = pmax(y_max - y_min, 0.1),
    x_start = 2,
    x_end = 3,
    y_line = y_max + pmax(y_span * 0.08, 0.03),
    y_tick = y_line - pmax(y_span * 0.025, 0.01),
    y_text = y_line + pmax(y_span * 0.045, 0.015),
    module = factor(as.character(module), levels = names(tcell_score_modules))
  )

p_module <- ggplot(module_score_sample, aes(Group, module_score, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.36, alpha = 0.82) +
  geom_point(
    position = position_jitter(width = 0.11, height = 0, seed = 20260826),
    shape = 21,
    colour = "black",
    stroke = 0.18,
    size = 1.25,
    alpha = 0.88
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = x_start, xend = x_end, y = y_line, yend = y_line),
    inherit.aes = FALSE,
    linewidth = 0.34,
    colour = "black",
    lineend = "round"
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = x_start, xend = x_start, y = y_tick, yend = y_line),
    inherit.aes = FALSE,
    linewidth = 0.34,
    colour = "black",
    lineend = "round"
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = x_end, xend = x_end, y = y_tick, yend = y_line),
    inherit.aes = FALSE,
    linewidth = 0.34,
    colour = "black",
    lineend = "round"
  ) +
  geom_text(
    data = wilcox_plot,
    aes(x = (x_start + x_end) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 2.25
  ) +
  facet_wrap(~module, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_colors, drop = FALSE, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "Module score per sample") +
  theme_pub(8.8) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7.0),
    strip.text = element_text(size = 7.8),
    legend.position = "top",
    legend.text = element_text(size = 7.8),
    legend.title = element_text(size = 8.2, face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_blank()
  )
save_plot_safe(p_module, "Tcell_functional_module_score_compare_revision", width = 4.8, height = 6.0)

message("Saved urine T-cell group key-gene heatmap and module score comparison.")

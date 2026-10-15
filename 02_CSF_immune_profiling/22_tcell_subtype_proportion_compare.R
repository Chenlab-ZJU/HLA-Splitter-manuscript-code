# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 13_revision_tcell_subtype_proportion_compare.R
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
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

env_project_dir <- Sys.getenv("CSF_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "13_revision_tcell_subtype_proportion_compare.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
group_colors <- c("Metastases" = "#E42320", "Hemorrhage" = "#699ECA", "Encephalitis" = "#4DAF4A")
tcell_levels <- c(
  "Naive CD4 T", "CD4 memory T", "GZMK+ CD4 T", "ISG+ CD4 T", "Treg",
  "Naive CD8 T", "GZMK+ CD8 T", "GZMB+ cytotoxic T", "NK-like T", "gd T"
)

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    p < 0.1 ~ paste0("p=", formatC(p, format = "f", digits = 2)),
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
  )
}

theme_pub <- function(base_size = 9.2) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
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

read_csv_safe <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

suppressPackageStartupMessages(library(Seurat))
manual_annotation_path <- file.path(getOption("hla.data_root"), "csf", "Tcell_revision_manual_cluster_reannotation.csv")
if (!file.exists(manual_annotation_path)) {
  stop("Missing Tcell_revision_manual_cluster_reannotation.csv. Please run 11b_revision_redraw_tcell_umap_from_csf_tcell.R first.")
}
manual_annotation <- read_csv_safe(manual_annotation_path) %>%
  transmute(Cluster = as.character(Cluster), T_celltype = as.character(T_celltype)) %>%
  filter(!is.na(Cluster), !is.na(T_celltype), T_celltype %in% tcell_levels)

# Read the saved object or source table used by this analysis.
Tcell <- readRDS(file.path(getOption("hla.data_root"), "csf", "CSF_Tcell.rds"))
meta <- Tcell@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    Sample = as.character(HLA_Splitter),
    Group = as.character(Group),
    Cluster = as.character(seurat_clusters)
  ) %>%
  left_join(manual_annotation, by = "Cluster") %>%
  filter(!is.na(Sample), !is.na(Group), !is.na(T_celltype), T_celltype %in% tcell_levels)

prop_df <- meta %>%
  count(Sample, Group, T_celltype, name = "Count") %>%
  group_by(Sample, Group) %>%
  mutate(Total = sum(Count), Percent = Count / Total * 100) %>%
  ungroup() %>%
  complete(nesting(Sample, Group), T_celltype = tcell_levels, fill = list(Count = 0, Percent = 0)) %>%
  group_by(Sample, Group) %>%
  mutate(Total = sum(Count)) %>%
  ungroup()

prop_df <- prop_df %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    T_celltype = factor(as.character(T_celltype), levels = tcell_levels)
  ) %>%
  filter(!is.na(Group), !is.na(T_celltype)) %>%
  arrange(T_celltype, Group, Sample)

# Export the analysis source data for figure reproduction.
write_csv_safe(prop_df, file.path(result_dir, "Tcell_revision_subtype_sample_percent.csv"))

kw_stats <- prop_df %>%
  group_by(T_celltype) %>%
  summarise(
    n_samples = n_distinct(Sample),
    n_groups = n_distinct(Group),
    p.value = ifelse(n_groups >= 2 && n_distinct(Percent) > 1, kruskal.test(Percent ~ Group)$p.value, NA_real_),
    y_max = max(Percent, na.rm = TRUE),
    y_min = min(Percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    p.adjust = p.adjust(p.value, method = "BH"),
    label = paste0("KW ", format_p(p.value)),
    y_range = pmax(y_max - y_min, 0.8),
    label.y = y_max + y_range * 0.11
  )
write_csv_safe(kw_stats, file.path(result_dir, "Tcell_revision_subtype_proportion_KW_tests.csv"))

selected_tcell_levels <- c("GZMK+ CD4 T", "ISG+ CD4 T", "GZMK+ CD8 T")
selected_tcell_levels <- selected_tcell_levels[selected_tcell_levels %in% as.character(unique(prop_df$T_celltype))]
if (length(selected_tcell_levels) < 3) {
  stop("Not all requested T cell subtypes were found in the current annotation.")
}
write_csv_safe(
  tibble(T_celltype = selected_tcell_levels),
  file.path(result_dir, "Tcell_revision_subtype_proportion_selected_diff_subtypes.csv")
)

prop_df <- prop_df %>%
  filter(as.character(T_celltype) %in% selected_tcell_levels) %>%
  mutate(T_celltype = factor(as.character(T_celltype), levels = selected_tcell_levels))
kw_stats <- kw_stats %>%
  filter(as.character(T_celltype) %in% selected_tcell_levels) %>%
  mutate(T_celltype = factor(as.character(T_celltype), levels = selected_tcell_levels))

comparisons <- as_tibble(t(utils::combn(group_levels, 2)), .name_repair = "minimal") %>%
  setNames(c("group1", "group2"))

wilcox_stats <- prop_df %>%
  group_by(T_celltype) %>%
  group_modify(function(.x, .y) {
    y_max <- max(.x$Percent, na.rm = TRUE)
    y_min <- min(.x$Percent, na.rm = TRUE)
    y_range <- max(y_max - y_min, 0.8)
    comparisons %>%
      mutate(
        x = match(group1, group_levels),
        xend = match(group2, group_levels),
        p.value = purrr::map2_dbl(group1, group2, function(g1, g2) {
          v1 <- .x$Percent[.x$Group == g1]
          v2 <- .x$Percent[.x$Group == g2]
          if (length(v1) < 2 || length(v2) < 2 || dplyr::n_distinct(c(v1, v2)) < 2) return(NA_real_)
          stats::wilcox.test(v1, v2, exact = FALSE)$p.value
        }),
        y_max = y_max,
        y_range = y_range
      )
  }) %>%
  ungroup() %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH"), label = format_p(p.value))
write_csv_safe(wilcox_stats, file.path(result_dir, "Tcell_revision_subtype_proportion_wilcox_tests.csv"))

wilcox_plot <- wilcox_stats %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  group_by(T_celltype) %>%
  arrange(p.value, .by_group = TRUE) %>%
  mutate(
    y.position = y_max + row_number() * y_range * 0.16 + y_range * 0.06,
    tick.height = y_range * 0.035,
    label.y = y.position + y_range * 0.045
  ) %>%
  ungroup()

p_prop <- ggplot(prop_df, aes(Group, Percent, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.38, alpha = 0.82) +
  geom_point(
    position = position_jitter(width = 0.11, height = 0, seed = 20260816),
    shape = 21,
    colour = "black",
    stroke = 0.20,
    size = 1.35,
    alpha = 0.90
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.30
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = x, xend = x, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.30
  ) +
  geom_segment(
    data = wilcox_plot,
    aes(x = xend, xend = xend, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.30
  ) +
  geom_text(
    data = wilcox_plot,
    aes(x = (x + xend) / 2, y = label.y, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 2.65
  ) +
  facet_wrap(~T_celltype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.22))) +
  labs(x = NULL, y = "T cell subset frequency per sample (%)") +
  theme_pub(8.9) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7.2),
    axis.text.y = element_text(size = 7.8),
    strip.text = element_text(size = 8.2),
    legend.position = "top",
    legend.text = element_text(size = 8.2),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(4, 5, 6, 10)
  )

save_plot_safe(p_prop, "Tcell_subtype_proportion_compare_revision", width = 6.4, height = 3.2)
message("Saved T cell subtype proportion comparison boxplot.")

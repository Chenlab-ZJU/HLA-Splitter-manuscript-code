# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Required intermediate source-table / annotation generation.
# Original source: 08_revision_scRNA_celltype_composition_plots.R
# Figure/table scope: Fig2; FigS4; FigS5
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
  library(readr)
  library(gtools)
  library(patchwork)
})

env_project_dir <- Sys.getenv("CSF_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "08_revision_scRNA_celltype_composition_plots.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}
figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Read the saved object or source table used by this analysis.
scRNA <- readRDS(file.path(getOption("hla.data_root"), "csf", "CSF_scRNA.rds"))

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
celltype_levels <- c(
  "Naive T", "CD4+ T", "CD8+ T", "IFN+ T", "Tregs", "gd T", "Proliferating",
  "B", "Macrophages", "Monocytes", "cDC1", "cDC2", "pDC"
)
major_map <- c(
  "Naive T" = "T cells",
  "CD4+ T" = "T cells",
  "CD8+ T" = "T cells",
  "IFN+ T" = "T cells",
  "Tregs" = "T cells",
  "gd T" = "T cells",
  "Proliferating" = "Cycling",
  "B" = "B cells",
  "Macrophages" = "Myeloid",
  "Monocytes" = "Myeloid",
  "cDC1" = "DC",
  "cDC2" = "DC",
  "pDC" = "DC"
)
major_levels <- c("T cells", "B cells", "Myeloid", "DC", "Cycling")

group_colors <- c(
  "Metastases" = "#E42320",
  "Hemorrhage" = "#699ECA",
  "Encephalitis" = "#4DAF4A"
)
celltype_colors <- c(
  "Naive T" = "#E42320",
  "CD4+ T" = "#699ECA",
  "CD8+ T" = "#4DAF4A",
  "IFN+ T" = "#FF8C00",
  "Tregs" = "#F898CB",
  "gd T" = "#731A73",
  "Proliferating" = "#D65190",
  "B" = "#FFCB5B",
  "Macrophages" = "#E89B1E",
  "Monocytes" = "#0076B9",
  "cDC1" = "#3D505A",
  "cDC2" = "#0098B2",
  "pDC" = "#8A7355"
)

meta <- scRNA@meta.data %>%
  transmute(
    Sample = as.character(HLA_Splitter),
    Group = as.character(Group),
    Celltype = as.character(Celltype)
  ) %>%
  filter(
    !is.na(Sample), !Sample %in% c("Undefined", "Doublets", "Doublet"),
    !is.na(Group), Group %in% group_levels,
    !is.na(Celltype), Celltype %in% celltype_levels
  ) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Celltype = factor(Celltype, levels = celltype_levels),
    Major_celltype = factor(unname(major_map[as.character(Celltype)]), levels = major_levels)
  )

sample_group <- meta %>%
  distinct(Sample, Group)

sample_order <- sample_group %>%
  distinct(Sample, Group) %>%
  arrange(Group, mixedorder(Sample)) %>%
  pull(Sample)

celltype_comp <- meta %>%
  count(Sample, Celltype, name = "Count") %>%
  complete(Sample = sample_order, Celltype = celltype_levels, fill = list(Count = 0)) %>%
  left_join(sample_group, by = "Sample") %>%
  group_by(Sample) %>%
  mutate(
    Total_count = sum(Count),
    Percent = ifelse(Total_count > 0, Count / Total_count * 100, 0)
  ) %>%
  ungroup() %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    Group = factor(Group, levels = group_levels),
    Celltype = factor(Celltype, levels = rev(celltype_levels))
  )

major_comp <- meta %>%
  count(Sample, Major_celltype, name = "Count") %>%
  complete(Sample = sample_order, Major_celltype = major_levels, fill = list(Count = 0)) %>%
  left_join(sample_group, by = "Sample") %>%
  group_by(Sample) %>%
  mutate(
    Total_count = sum(Count),
    Percent = ifelse(Total_count > 0, Count / Total_count * 100, 0)
  ) %>%
  ungroup() %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    Group = factor(Group, levels = group_levels),
    Major_celltype = factor(Major_celltype, levels = major_levels)
  )

# Export the analysis source data for figure reproduction.
write_csv(celltype_comp, file.path(getOption("hla.data_root"), "csf", "scRNA_celltype_sample_percent_revision.csv"))
write_csv(major_comp, file.path(getOption("hla.data_root"), "csf", "scRNA_major_celltype_sample_percent_revision.csv"))

format_kw_p <- function(p) {
  ifelse(is.na(p), "KW p=NA", ifelse(p < 0.001, "KW p<0.001", paste0("KW p=", signif(p, 2))))
}

kw_test_by_type <- function(df, type_col) {
  df %>%
    group_by(.data[[type_col]]) %>%
    summarise(
      n_groups = n_distinct(Group),
      p_value = ifelse(
        n_groups >= 2 && dplyr::n_distinct(Percent) > 1,
        kruskal.test(Percent ~ Group)$p.value,
        NA_real_
      ),
      y.position = max(Percent, na.rm = TRUE) * 1.10 + 0.03,
      label = format_kw_p(p_value),
      .groups = "drop"
    )
}

major_kw <- kw_test_by_type(major_comp, "Major_celltype")
celltype_kw <- kw_test_by_type(celltype_comp, "Celltype")

write_csv(major_kw, file.path(result_dir, "scRNA_major_celltype_kw_test_revision.csv"))
write_csv(celltype_kw, file.path(result_dir, "scRNA_celltype_kw_test_revision.csv"))

format_wilcox_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", signif(p, 2))
  )
}

selected_celltypes <- c("CD4+ T", "B", "Macrophages")
selected_comp <- celltype_comp %>%
  filter(Celltype %in% selected_celltypes) %>%
  mutate(
    Celltype = factor(Celltype, levels = selected_celltypes),
    Group = factor(Group, levels = group_levels)
  )

selected_wilcox <- selected_comp %>%
  group_by(Celltype) %>%
  group_modify(function(.x, .y) {
    comparisons <- t(utils::combn(group_levels, 2))
    y_max <- max(.x$Percent, na.rm = TRUE)
    y_step <- max(y_max * 0.14, 0.75)
    tibble(
      group1 = comparisons[, 1],
      group2 = comparisons[, 2],
      x = match(comparisons[, 1], group_levels),
      xend = match(comparisons[, 2], group_levels),
      y.position = y_max + seq_len(nrow(comparisons)) * y_step,
      p.value = apply(comparisons, 1, function(pair) {
        values_1 <- .x$Percent[.x$Group == pair[1]]
        values_2 <- .x$Percent[.x$Group == pair[2]]
        if (length(values_1) < 2 || length(values_2) < 2) return(NA_real_)
        stats::wilcox.test(values_1, values_2, exact = FALSE)$p.value
      })
    ) %>%
      mutate(label = format_wilcox_p(p.value))
  }) %>%
  ungroup()

write_csv(selected_wilcox, file.path(result_dir, "scRNA_selected_celltype_wilcox_test_revision.csv"))

selected_wilcox_plot <- selected_wilcox %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  left_join(
    selected_comp %>%
      group_by(Celltype) %>%
      summarise(y_max = max(Percent, na.rm = TRUE), .groups = "drop"),
    by = "Celltype"
  ) %>%
  group_by(Celltype) %>%
  arrange(x, xend, .by_group = TRUE) %>%
  mutate(
    y_step = pmax(y_max * 0.26, 1.35),
    y.position = y_max + row_number() * y_step + y_step * 0.25,
    tick.height = pmax(y_step * 0.14, 0.18),
    label.y = y.position + pmax(y_step * 0.34, 0.38)
  ) %>%
  ungroup()

theme_pub <- function(base_size = 8.5) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = "grey40", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey45", linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

save_plot_safe <- function(plot, filename, width, height) {
  target_pdf <- file.path(figure_dir, paste0(filename, ".pdf"))
  target_png <- file.path(figure_dir, paste0(filename, ".png"))
  fallback_pdf <- file.path(figure_dir, paste0(filename, "_revision.pdf"))
  fallback_png <- file.path(figure_dir, paste0(filename, "_revision.png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  if (!file.copy(tmp_pdf, target_pdf, overwrite = TRUE)) {
    warning("Could not overwrite target PDF; writing fallback file instead: ", fallback_pdf)
    file.copy(tmp_pdf, fallback_pdf, overwrite = TRUE)
  }
  if (!file.copy(tmp_png, target_png, overwrite = TRUE)) {
    warning("Could not overwrite target PNG; writing fallback file instead: ", fallback_png)
    file.copy(tmp_png, fallback_png, overwrite = TRUE)
  }
}

p_heatmap <- ggplot(celltype_comp, aes(x = Sample, y = Celltype, fill = Percent)) +
  geom_tile(colour = "white", linewidth = 0.28) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#6F99C3", "#D7E4EF", "#F7F7F7", "#E79A94", "#B7192B"),
    name = "Cells (%)",
    limits = c(0, max(celltype_comp$Percent, na.rm = TRUE))
  ) +
  labs(x = "Samples grouped by diagnosis", y = "Celltype") +
  theme_pub(10.2) +
  theme(
    axis.title = element_text(size = 11.8, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9.0),
    axis.text.y = element_text(size = 10.0),
    strip.text = element_text(size = 10.8, face = "bold"),
    legend.title = element_text(size = 10.6, face = "bold"),
    legend.text = element_text(size = 9.6),
    panel.grid = element_blank(),
    legend.position = "right"
  )

plot_box <- function(df, type_col, type_levels, title, stat_df) {
  ggplot(df, aes(x = Group, y = Percent, fill = Group)) +
    geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.35, alpha = 0.78) +
    geom_point(
      position = position_jitter(width = 0.13, height = 0, seed = 20260813),
      size = 1.05,
      alpha = 0.78,
      shape = 21,
      colour = "black",
      stroke = 0.18
    ) +
    facet_wrap(stats::as.formula(paste("~", type_col)), scales = "free_y", ncol = 4) +
    geom_text(
      data = stat_df,
      aes(x = 2, y = y.position, label = label),
      inherit.aes = FALSE,
      size = 2.45,
      fontface = "bold",
      family = "Arial",
      colour = "black"
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(title = title, x = NULL, y = "Cells per sample (%)") +
    theme_pub(8.2) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 7.0),
      strip.text = element_text(size = 8.0),
      legend.position = "top",
      panel.grid.major.x = element_blank()
    )
}

p_major_box <- plot_box(major_comp, "Major_celltype", major_levels, "Major celltype composition", major_kw)
p_celltype_box <- plot_box(
  celltype_comp %>% mutate(Celltype = factor(Celltype, levels = celltype_levels)),
  "Celltype",
  celltype_levels,
  "Celltype composition",
  celltype_kw
)

p_selected_celltype_box <- ggplot(selected_comp, aes(x = Group, y = Percent, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.38, alpha = 0.82) +
  geom_point(
    position = position_jitter(width = 0.11, height = 0, seed = 20260813),
    size = 1.35,
    alpha = 0.90,
    shape = 21,
    colour = "black",
    stroke = 0.20
  ) +
  geom_segment(
    data = selected_wilcox_plot,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = selected_wilcox_plot,
    aes(x = x, xend = x, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = selected_wilcox_plot,
    aes(x = xend, xend = xend, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_text(
    data = selected_wilcox_plot,
    aes(x = (x + xend) / 2, y = label.y, label = label),
    inherit.aes = FALSE,
    size = 2.45,
    family = "Arial",
    fontface = "bold",
    colour = "black"
  ) +
  facet_wrap(~Celltype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_x_discrete(limits = group_levels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "Cells per sample (%)") +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7.8),
    strip.text = element_text(size = 8.2),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

save_plot_safe(p_heatmap, "scRNA_all_celltype_heatmap", width = 7.2, height = 4.6)
save_plot_safe(p_major_box, "scRNA_major_celltype_compare", width = 7.0, height = 3.8)
save_plot_safe(p_celltype_box, "scRNA_celltype_compare", width = 8.8, height = 6.4)
save_plot_safe(p_selected_celltype_box, "scRNA_selected_celltype_compare_wilcox", width = 5.2, height = 2.8)

message("Updated CSF celltype heatmap and boxplots in: ", figure_dir)

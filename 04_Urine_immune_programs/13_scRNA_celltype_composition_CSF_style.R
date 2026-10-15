# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final urine immune-cell composition, module and demultiplexing figures.
# Original source: 25_revision_scRNA_celltype_composition_CSF_style.R
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
  library(gtools)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "25_revision_scRNA_celltype_composition_CSF_style.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures", "scRNA", "Celltype_Composition")
result_dir <- file.path(project_dir, "Result", "scRNA", "Celltype_Composition")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
boxplot_group_levels <- c("NR", "TCMR", "ABMR", "BK")
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "Mix" = "#FF8C00",
  "BK" = "#F898CB"
)

major_levels <- c("CD4+ T", "CD8+ T", "Cycling T", "NK", "B cell", "Monocyte", "Macrophage", "Neutrophil", "cDC", "pDC")
celltype_levels <- c(
  "Naive CD4+ T cells",
  "Memory CD4+ T cells",
  "Naive CD8+ T cells",
  "CD8+ NKT-like cells",
  "ISG expressing immune cells",
  "Naive B cells",
  "Classical Monocytes",
  "Macrophages",
  "Myeloid Dendritic cells",
  "Plasmacytoid Dendritic cells",
  "Unknown"
)

read_csv_safe <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "KW p=NA",
    p < 0.001 ~ "KW p<0.001",
    p < 0.01 ~ paste0("KW p=", formatC(p, format = "f", digits = 3)),
    p < 0.1 ~ paste0("KW p=", formatC(p, format = "f", digits = 2)),
    TRUE ~ paste0("KW p=", formatC(p, format = "f", digits = 2))
  )
}

format_wilcox_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    p < 0.01 ~ formatC(p, format = "f", digits = 3),
    TRUE ~ formatC(p, format = "f", digits = 2)
  )
}

theme_pub <- function(base_size = 8.8) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      panel.grid.minor = element_blank(),
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

object_path <- file.path(getOption("hla.data_root"), "urine", "Urine_ALL.rds")
if (!file.exists(object_path)) {
  object_path <- file.path(project_dir, "Result", "Objects", "Urine_Allcells.rds")
}
if (!file.exists(object_path)) {
  stop("Urine scRNA object not found in Result/Objects.")
}

# Read the saved object or source table used by this analysis.
scRNA <- readRDS(object_path)
meta <- scRNA@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
    Group = as.character(Group),
    Major_celltype = as.character(Major_celltype),
    Celltype = as.character(celltype)
  ) %>%
  filter(
    !is.na(Sample), !Sample %in% c("Undefined", "Doublets", "Doublet", ""),
    !is.na(Group), Group %in% group_levels,
    !is.na(Major_celltype), !is.na(Celltype),
    Major_celltype %in% major_levels,
    Celltype %in% celltype_levels
  ) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Major_celltype = factor(Major_celltype, levels = major_levels),
    Celltype = factor(Celltype, levels = celltype_levels)
  )

sample_group <- meta %>%
  distinct(Sample, Group)
sample_order <- sample_group %>%
  arrange(Group, mixedorder(Sample)) %>%
  pull(Sample)

make_composition <- function(df, type_col, type_levels) {
  df %>%
    count(Sample, Group, Cell_label = .data[[type_col]], name = "Count") %>%
    complete(
      nesting(Sample, Group),
      Cell_label = factor(type_levels, levels = type_levels),
      fill = list(Count = 0)
    ) %>%
    group_by(Sample, Group) %>%
    mutate(
      Total_count = sum(Count),
      Percent = ifelse(Total_count > 0, Count / Total_count * 100, 0)
    ) %>%
    ungroup() %>%
    mutate(
      Sample = factor(Sample, levels = sample_order),
      Group = factor(Group, levels = group_levels),
      Cell_label = factor(as.character(Cell_label), levels = type_levels)
    )
}

major_comp <- make_composition(meta, "Major_celltype", major_levels)
celltype_comp <- make_composition(meta, "Celltype", celltype_levels)

# Export the analysis source data for figure reproduction.
write_csv_safe(major_comp, file.path(result_dir, "Urine_scRNA_major_celltype_sample_percent_revision.csv"))
write_csv_safe(celltype_comp, file.path(result_dir, "Urine_scRNA_celltype_sample_percent_revision.csv"))

kw_test_by_type <- function(df, label_col = "Cell_label") {
  df %>%
    group_by(.data[[label_col]]) %>%
    summarise(
      n_samples = n_distinct(Sample),
      n_groups = n_distinct(Group),
      p.value = ifelse(n_groups >= 2 && dplyr::n_distinct(Percent) > 1, kruskal.test(Percent ~ Group)$p.value, NA_real_),
      y_max = max(Percent, na.rm = TRUE),
      y_min = min(Percent, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      p.adjust = p.adjust(p.value, method = "BH"),
      label = format_p(p.value),
      y_range = pmax(y_max - y_min, 0.8),
      label.y = y_max + y_range * 0.16
    )
}

pairwise_wilcox_by_type <- function(df, compare_group_levels = group_levels) {
  comparisons <- as_tibble(t(utils::combn(compare_group_levels, 2)), .name_repair = "minimal") %>%
    setNames(c("group1", "group2"))
  df %>%
    group_by(Cell_label) %>%
    group_modify(function(.x, .y) {
      comparisons %>%
        mutate(
          p.value = purrr::map2_dbl(group1, group2, function(g1, g2) {
            v1 <- .x$Percent[.x$Group == g1]
            v2 <- .x$Percent[.x$Group == g2]
            if (length(v1) < 2 || length(v2) < 2 || dplyr::n_distinct(c(v1, v2)) < 2) return(NA_real_)
            stats::wilcox.test(v1, v2, exact = FALSE)$p.value
          })
        )
    }) %>%
    ungroup() %>%
    mutate(p.adjust = p.adjust(p.value, method = "BH"))
}

major_kw <- kw_test_by_type(major_comp)
celltype_kw <- kw_test_by_type(celltype_comp)
major_wilcox <- pairwise_wilcox_by_type(major_comp)
celltype_wilcox <- pairwise_wilcox_by_type(celltype_comp)
major_comp_no_mix <- major_comp %>%
  filter(as.character(Group) %in% boxplot_group_levels) %>%
  mutate(Group = factor(as.character(Group), levels = boxplot_group_levels))
major_wilcox_no_mix <- pairwise_wilcox_by_type(major_comp_no_mix, compare_group_levels = boxplot_group_levels)

write_csv_safe(major_kw, file.path(result_dir, "Urine_scRNA_major_celltype_KW_tests_revision.csv"))
write_csv_safe(celltype_kw, file.path(result_dir, "Urine_scRNA_celltype_KW_tests_revision.csv"))
write_csv_safe(major_wilcox, file.path(result_dir, "Urine_scRNA_major_celltype_pairwise_wilcox_revision.csv"))
write_csv_safe(celltype_wilcox, file.path(result_dir, "Urine_scRNA_celltype_pairwise_wilcox_revision.csv"))
write_csv_safe(major_comp_no_mix, file.path(result_dir, "Urine_scRNA_major_celltype_sample_percent_no_Mix_revision.csv"))
write_csv_safe(major_wilcox_no_mix, file.path(result_dir, "Urine_scRNA_major_celltype_pairwise_wilcox_no_Mix_revision.csv"))

heatmap_df <- celltype_comp %>%
  mutate(Cell_label = factor(as.character(Cell_label), levels = rev(celltype_levels)))
p_heatmap <- ggplot(heatmap_df, aes(Sample, Cell_label, fill = Percent)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#6F99C3", "#D7E4EF", "#F7F7F7", "#E79A94", "#B7192B"),
    name = "Cells (%)",
    limits = c(0, max(heatmap_df$Percent, na.rm = TRUE))
  ) +
  labs(x = "Samples grouped by diagnosis", y = NULL) +
  theme_pub(8.6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.8),
    axis.text.y = element_text(size = 7.6),
    axis.title = element_text(size = 8.6, face = "bold"),
    strip.text = element_text(size = 8.0),
    panel.grid = element_blank(),
    legend.title = element_text(size = 8.2, face = "bold"),
    legend.text = element_text(size = 7.2),
    plot.margin = margin(4, 8, 8, 8)
  )

plot_box <- function(
  comp_df,
  kw_df,
  filename,
  ncol,
  width,
  height,
  y_label,
  wilcox_df = NULL,
  p_cutoff = 0.1,
  plot_group_levels = group_levels,
  csf_style = FALSE,
  bar_style = FALSE
) {
  plot_df <- comp_df %>%
    filter(as.character(Group) %in% plot_group_levels) %>%
    mutate(
      Group = factor(as.character(Group), levels = plot_group_levels),
      Cell_label = factor(as.character(Cell_label), levels = unique(as.character(Cell_label)))
    )
  if (is.null(wilcox_df)) {
    kw_plot <- kw_df %>%
      mutate(Cell_label = factor(as.character(Cell_label), levels = levels(plot_df$Cell_label)))
  } else {
    kw_plot <- tibble(
      Cell_label = factor(character(), levels = levels(plot_df$Cell_label)),
      label.y = numeric(),
      label = character()
    )
  }
  y_stats <- plot_df %>%
    group_by(Cell_label) %>%
    summarise(
      y_max = max(Percent, na.rm = TRUE),
      y_min = min(Percent, na.rm = TRUE),
      y_range = pmax(y_max - y_min, 0.9),
      .groups = "drop"
    )
  wilcox_plot <- tibble(
    Cell_label = factor(character(), levels = levels(plot_df$Cell_label)),
    group1 = character(),
    group2 = character(),
    x = numeric(),
    xend = numeric(),
    y = numeric(),
    yend = numeric(),
    y_text = numeric(),
    tick_y = numeric(),
    label = character()
  )
  if (!is.null(wilcox_df)) {
    wilcox_plot <- wilcox_df %>%
      filter(!is.na(p.value), p.value < p_cutoff) %>%
      mutate(
        Cell_label = factor(as.character(Cell_label), levels = levels(plot_df$Cell_label)),
        group1 = factor(group1, levels = plot_group_levels),
        group2 = factor(group2, levels = plot_group_levels),
        x = as.numeric(group1),
        xend = as.numeric(group2),
        label = format_wilcox_p(p.value)
      ) %>%
      arrange(Cell_label, p.value, x, xend) %>%
      group_by(Cell_label) %>%
      mutate(stat_rank = row_number()) %>%
      ungroup() %>%
      left_join(y_stats, by = "Cell_label") %>%
      mutate(
        y = y_max + y_range * (0.18 + 0.24 * (stat_rank - 1)),
        yend = y,
        tick_y = y - y_range * 0.055,
        y_text = y + y_range * 0.095
      )
  }
  if (identical(filename, "Urine_scRNA_major_celltype_compare_revision")) {
    write_csv_safe(wilcox_plot, file.path(result_dir, "Urine_scRNA_major_celltype_pairwise_wilcox_no_Mix_p_lt_0.1_revision.csv"))
  }
  p <- ggplot(plot_df, aes(Group, Percent, fill = Group, colour = Group))
  if (bar_style) {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "col",
        width = 0.58,
        colour = "black",
        linewidth = 0.42,
        alpha = 0.42
      ) +
      stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        width = 0.22,
        linewidth = 0.42,
        colour = "black"
      )
  } else {
    p <- p +
      geom_boxplot(
        width = 0.58,
        outlier.shape = NA,
        colour = "black",
        linewidth = ifelse(csf_style, 0.42, 0.36),
        alpha = ifelse(csf_style, 0.34, 0.82)
      )
  }
  p <- p +
    geom_point(
      position = position_jitter(width = 0.10, height = 0, seed = ifelse(csf_style, 11, 20260816)),
      shape = 21,
      stroke = ifelse(csf_style, 0.52, 0.16),
      size = ifelse(csf_style, 2.10, 0.95),
      alpha = ifelse(csf_style, 0.92, 0.86)
    ) +
    geom_text(
      data = kw_plot,
      aes(x = 3, y = label.y, label = label),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "bold",
      size = 2.25
    ) +
    geom_segment(
      data = wilcox_plot,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.36,
      colour = "black"
    ) +
    geom_segment(
      data = wilcox_plot,
      aes(x = x, xend = x, y = tick_y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.36,
      colour = "black"
    ) +
    geom_segment(
      data = wilcox_plot,
      aes(x = xend, xend = xend, y = tick_y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.36,
      colour = "black"
    ) +
    geom_text(
      data = wilcox_plot,
      aes(x = (x + xend) / 2, y = y_text, label = label),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "bold",
      size = ifelse(bar_style, 3.8, 3.05)
    ) +
    facet_wrap(~Cell_label, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, ifelse(is.null(wilcox_df), 0.15, 0.34)))) +
    labs(x = NULL, y = y_label) +
    theme_pub(ifelse(csf_style, 10.2, 8.4)) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = ifelse(csf_style, 9.2, 6.2)),
      axis.text.y = element_text(size = ifelse(csf_style, 9.8, 7.2)),
      strip.text = element_text(size = ifelse(csf_style, 9.2, 7.8), face = "bold"),
      legend.position = ifelse(csf_style, "none", "top"),
      legend.text = element_text(size = 7.8),
      legend.title = element_text(size = 8.2, face = "bold"),
      panel.grid.major.x = element_blank(),
      axis.title.y = element_text(size = ifelse(csf_style, 11.8, 8.4), face = "bold", margin = margin(r = 5)),
      panel.spacing.x = unit(ifelse(csf_style, 0.75, 0.62), "lines"),
      panel.spacing.y = unit(ifelse(csf_style, 0.72, 0.62), "lines"),
      plot.margin = margin(4, 8, 8, 12)
    )
  save_plot_safe(p, filename, width = width, height = height)
  invisible(p)
}

save_plot_safe(p_heatmap, "Urine_scRNA_all_celltype_heatmap_revision", width = 8.6, height = 3.8)
plot_box(
  major_comp_no_mix,
  major_kw,
  "Urine_scRNA_major_celltype_compare_revision",
  ncol = 5,
  width = 11.0,
  height = 4.8,
  y_label = "Frequency (%)",
  wilcox_df = major_wilcox_no_mix,
  p_cutoff = 0.1,
  plot_group_levels = boxplot_group_levels,
  csf_style = TRUE,
  bar_style = TRUE
)
plot_box(
  celltype_comp,
  celltype_kw,
  "Urine_scRNA_celltype_compare_revision",
  ncol = 4,
  width = 8.6,
  height = 5.4,
  y_label = "Frequency (%)"
)

message("Saved urine scRNA celltype heatmap and composition boxplots.")

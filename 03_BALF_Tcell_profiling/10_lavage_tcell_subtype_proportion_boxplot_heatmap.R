# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Final BALF figures from saved objects and signatures.
# Original source: 18_revision_lavage_tcell_subtype_proportion_boxplot_heatmap.R
# Figure/table scope: Fig3b; Fig3c; Fig3d; Fig3e; Fig3f; Fig3g; Fig3h; Fig3i; FigS6a; FigS6b; FigS6c; FigS6d; FigS6e; FigS6f; FigS6g; FigS6h
# Data cohort: balf. Inputs use the compact English data layout.
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
  library(patchwork)
})

project_dir <- Sys.getenv("LAVAGE_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "balf")
}
project_dir <- gsub("\\\\", "/", project_dir)
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("CAP", "LC", "ILD")
group_colors <- c("CAP" = "#699ECA", "LC" = "#E42320", "ILD" = "#4DAF4A")
tcell_type_levels <- c(
  "Naive CD4", "Th1-like", "Tfh", "Th17", "Treg", "GZMK+ CD8",
  "ITGAE+ CD8", "NKT", "ISG+ T", "Cycling T", "gdT", "MAIT"
)
fallback_colors <- c(
  "#E42320", "#699ECA", "#4DAF4A", "#FF8C00", "#F898CB", "#D65190",
  "#731A73", "#FFCB5B", "#E89B1E", "#0076B9", "#3D505A", "#0098B2"
)
tcell_type_colors <- setNames(fallback_colors[seq_along(tcell_type_levels)], tcell_type_levels)
color_map_path <- file.path(getOption("hla.data_root"), "balf", "lavage_Tcell_UMAP_Cluster_Group_color_map.csv")
if (file.exists(color_map_path)) {
  color_map <- read.csv(color_map_path, stringsAsFactors = FALSE, check.names = FALSE)
  color_map <- color_map[color_map$category == "T_celltype", , drop = FALSE]
  if (nrow(color_map) > 0) {
    tcell_type_colors[color_map$label] <- color_map$color
  }
}

p_threshold <- 0.1

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
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

theme_box <- function(base_size = 9.2) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      strip.text = element_text(face = "bold", size = base_size - 0.2),
      strip.background = element_rect(fill = "grey94", colour = NA),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.position = "top",
      plot.margin = margin(4, 5, 5, 8)
    )
}

theme_heatmap <- function(base_size = 9.2) {
  theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      strip.text = element_text(face = "bold", size = base_size),
      strip.background = element_rect(fill = "grey94", colour = NA),
      panel.grid = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(4, 5, 5, 5)
    )
}

message("Loading lavage_Tcell.rds and calculating sample-level T-cell subtype proportions...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(file.path(getOption("hla.data_root"), "balf", "lavage_Tcell.rds"))
required_cols <- c("HLA_Splitter", "Group", "revised_celltype")
missing_cols <- setdiff(required_cols, colnames(tcell@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

meta <- tcell@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Sample = as.character(HLA_Splitter),
    Group = factor(as.character(Group), levels = group_levels),
    T_celltype = factor(as.character(revised_celltype), levels = tcell_type_levels)
  ) %>%
  filter(!is.na(Sample), !is.na(Group), !is.na(T_celltype))

present_types <- levels(droplevels(meta$T_celltype))
type_colors_present <- tcell_type_colors[present_types]

prop_df <- meta %>%
  count(Sample, Group, T_celltype, name = "Count") %>%
  complete(nesting(Sample, Group), T_celltype = present_types, fill = list(Count = 0)) %>%
  group_by(Sample, Group) %>%
  mutate(Total_T_cells = sum(Count), Percent = Count / Total_T_cells * 100) %>%
  ungroup() %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    T_celltype = factor(as.character(T_celltype), levels = present_types)
  ) %>%
  arrange(Group, Sample, T_celltype)

# Export the analysis source data for figure reproduction.
write.csv(prop_df, file.path(result_dir, "lavage_Tcell_subtype_sample_proportion_revision.csv"), row.names = FALSE)

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

wilcox_plot <- wilcox_stats %>%
  filter(!is.na(p.value), p.value < p_threshold) %>%
  group_by(T_celltype) %>%
  arrange(p.value, .by_group = TRUE) %>%
  mutate(
    y.position = y_max + row_number() * y_range * 0.16 + y_range * 0.05,
    tick.height = y_range * 0.035,
    label.y = y.position + y_range * 0.055
  ) %>%
  ungroup()

write.csv(wilcox_stats, file.path(result_dir, "lavage_Tcell_subtype_proportion_wilcox_tests_revision.csv"), row.names = FALSE)
write.csv(wilcox_plot, file.path(result_dir, "lavage_Tcell_subtype_proportion_wilcox_tests_p_lt_0.1_revision.csv"), row.names = FALSE)

plot_types <- if (nrow(wilcox_plot) > 0) {
  wilcox_plot %>%
    distinct(T_celltype) %>%
    pull(T_celltype) %>%
    as.character()
} else {
  present_types
}
plot_types <- intersect(present_types, plot_types)

box_df <- prop_df %>%
  filter(as.character(T_celltype) %in% plot_types) %>%
  mutate(T_celltype = factor(as.character(T_celltype), levels = plot_types))
wilcox_plot_box <- wilcox_plot %>%
  filter(as.character(T_celltype) %in% plot_types) %>%
  mutate(T_celltype = factor(as.character(T_celltype), levels = plot_types))

p_box <- ggplot(box_df, aes(Group, Percent, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.42, alpha = 0.34, colour = "black") +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 11),
    shape = 21,
    stroke = 0.52,
    size = 2.10,
    alpha = 0.92
  ) +
  geom_segment(
    data = wilcox_plot_box,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = wilcox_plot_box,
    aes(x = x, xend = x, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = wilcox_plot_box,
    aes(x = xend, xend = xend, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_text(
    data = wilcox_plot_box,
    aes(x = (x + xend) / 2, y = label.y, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    size = 3.0,
    vjust = -0.30,
    colour = "black"
  ) +
  facet_wrap(~T_celltype, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "T cell subset frequency per sample (%)") +
  theme_box(10.2) +
  theme(
    strip.background = element_rect(fill = "grey94", colour = NA),
    strip.text = element_text(size = 10.8, face = "bold"),
    axis.title.y = element_text(size = 11.8, face = "bold"),
    axis.text.x = element_text(size = 10.2, colour = "black"),
    axis.text.y = element_text(size = 9.8, colour = "black"),
    panel.spacing.x = unit(0.75, "lines"),
    legend.position = "none",
    plot.margin = margin(3, 4, 3, 4)
  )

heatmap_df <- prop_df %>%
  group_by(Group, T_celltype) %>%
  mutate(order_percent = mean(Percent, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Sample = factor(Sample, levels = unique(Sample[order(Group, Sample)])),
    T_celltype = factor(as.character(T_celltype), levels = rev(present_types))
  )

p_heatmap <- ggplot(heatmap_df, aes(Sample, T_celltype, fill = Percent)) +
  geom_tile(colour = "white", linewidth = 0.32) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#D7E4EF", "#F7F7F7", "#E79A94", "#B7192B"),
    values = scales::rescale(c(0, median(heatmap_df$Percent, na.rm = TRUE), quantile(heatmap_df$Percent, 0.85, na.rm = TRUE), max(heatmap_df$Percent, na.rm = TRUE))),
    name = "Cells (%)"
  ) +
  labs(x = "Samples grouped by diagnosis", y = "T cell subset") +
  theme_heatmap(9.2) +
  theme(
    axis.text.x = element_text(size = 7.2),
    axis.text.y = element_text(size = 8.2),
    legend.position = "right"
  )

save_plot_safe(p_box, "Tcell_subtype_proportion_compare_CAP_LC_ILD_revision", width = 5.0, height = 3.1)
save_plot_safe(p_heatmap, "Tcell_subtype_proportion_heatmap_CAP_LC_ILD_revision", width = 7.0, height = 4.8)

message("Saved lavage T-cell subtype proportion boxplot and heatmap.")

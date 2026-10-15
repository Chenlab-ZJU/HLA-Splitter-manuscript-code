# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Final urine T-cell, CXCR3 and expanded-clone figures.
# Original source: 31_revision_tcell_celltype_umap_TCMR_ABMR_facet.R
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
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "31_revision_tcell_celltype_umap_TCMR_ABMR_facet.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

object_file <- file.path(getOption("hla.data_root"), "urine", "Urine_Tcell.rds")
figure_dir <- file.path(project_dir, "Figures", "Tcell", "UMAP_Markers")
result_dir <- file.path(project_dir, "Result", "Tcell", "UMAP_Markers")
composition_figure_dir <- file.path(project_dir, "Figures", "Tcell", "Celltype_Composition")
composition_result_dir <- file.path(project_dir, "Result", "Tcell", "Celltype_Composition")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(composition_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(composition_result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
boxplot_group_levels <- c("NR", "TCMR", "ABMR", "BK")
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "Mix" = "#FF8C00",
  "BK" = "#F898CB"
)

tcell_levels <- c(
  "Naive CD4", "CD4 Tem", "CD4 Tex", "Treg",
  "Naive CD8", "CXCR4+ CD8", "GZMB+ CD8", "ISG+ CD8", "Cycling T"
)
tcell_colors <- c(
  "Naive CD4" = "#E42320",
  "CD4 Tem" = "#699ECA",
  "CD4 Tex" = "#8FB4DC",
  "Treg" = "#B85BA6",
  "Naive CD8" = "#4DAF4A",
  "CXCR4+ CD8" = "#FF8C00",
  "GZMB+ CD8" = "#D65190",
  "ISG+ CD8" = "#F898CB",
  "Cycling T" = "#731A73"
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
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

save_plot_to_dir <- function(plot, filename, width, height, out_dir) {
  target_pdf <- file.path(out_dir, paste0(filename, ".pdf"))
  target_png <- file.path(out_dir, paste0(filename, ".png"))
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

format_p <- function(p, prefix = "KW p") {
  ifelse(
    is.na(p), paste0(prefix, "=NA"),
    ifelse(p < 0.001, paste0(prefix, "<0.001"), paste0(prefix, "=", formatC(p, format = "f", digits = 2)))
  )
}

format_wilcox_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
  )
}

safe_kw <- function(df, value_col = "frequency_percent", group_col = "Group") {
  df <- df %>%
    filter(!is.na(.data[[value_col]]), !is.na(.data[[group_col]])) %>%
    mutate(.group_tmp = droplevels(factor(.data[[group_col]])))
  if (n_distinct(df$.group_tmp) < 2) return(NA_real_)
  if (length(unique(df[[value_col]])) < 2) return(NA_real_)
  suppressWarnings(kruskal.test(df[[value_col]] ~ df$.group_tmp)$p.value)
}

safe_wilcox_two_groups <- function(df, group_a, group_b, value_col = "frequency_percent", group_col = "Group") {
  df <- df %>%
    filter(.data[[group_col]] %in% c(group_a, group_b), !is.na(.data[[value_col]])) %>%
    mutate(.group_tmp = droplevels(factor(.data[[group_col]], levels = c(group_a, group_b))))
  if (n_distinct(df$.group_tmp) < 2) return(NA_real_)
  counts <- table(df$.group_tmp)
  if (any(counts[c(group_a, group_b)] < 1)) return(NA_real_)
  if (length(unique(df[[value_col]])) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(df[[value_col]] ~ df$.group_tmp, exact = FALSE)$p.value)
}

theme_umap <- function(base_size = 8.8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", size = base_size + 0.7, colour = "black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.36),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "grey93", colour = NA),
      strip.text = element_text(face = "bold", size = base_size + 0.3, colour = "black"),
      plot.margin = margin(4, 5, 4, 4)
    )
}

if (!file.exists(object_file)) {
  stop("Cannot find urine T-cell object: ", object_file)
}

message("Loading urine T-cell object...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(object_file)
if (!"umap" %in% Reductions(tcell)) {
  stop("UMAP reduction was not found in Urine_Tcell.rds.")
}

umap <- Embeddings(tcell, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell")
colnames(umap)[2:3] <- c("UMAP_1", "UMAP_2")

plot_df <- tcell@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Group = as.character(Group),
    Sample = as.character(Sample),
    T_celltype = as.character(revised_celltype),
    seurat_cluster = as.character(seurat_clusters)
  ) %>%
  left_join(umap, by = "cell") %>%
  filter(!is.na(UMAP_1), !is.na(UMAP_2), !is.na(T_celltype), T_celltype != "") %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    T_celltype = factor(T_celltype, levels = c(tcell_levels[tcell_levels %in% unique(T_celltype)], setdiff(sort(unique(T_celltype)), tcell_levels)))
  )

missing_cols <- setdiff(levels(plot_df$T_celltype), names(tcell_colors))
if (length(missing_cols) > 0) {
  extra_cols <- grDevices::hcl.colors(length(missing_cols), palette = "Dark 3")
  names(extra_cols) <- missing_cols
  tcell_colors <- c(tcell_colors, extra_cols)
}
tcell_colors <- tcell_colors[levels(plot_df$T_celltype)]

# Export the analysis source data for figure reproduction.
write_csv_safe(plot_df, file.path(result_dir, "TCell_celltype_UMAP_revision_source.csv"))

label_df <- plot_df %>%
  group_by(T_celltype) %>%
  summarise(
    n_cells = n(),
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

p_celltype <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = T_celltype)) +
  geom_point(size = 0.28, alpha = 0.80, stroke = 0) +
  scale_colour_manual(values = tcell_colors, drop = FALSE, name = "T-cell subtype") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 3.6, alpha = 1), ncol = 1)) +
  coord_equal() +
  theme_umap(8.8) +
  theme(legend.position = "right")

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_celltype <- p_celltype +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(x = UMAP_1, y = UMAP_2, label = T_celltype),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "bold",
      size = 2.25,
      colour = "black",
      min.segment.length = 0,
      segment.colour = "grey45",
      segment.size = 0.15,
      box.padding = 0.18,
      point.padding = 0.10,
      max.overlaps = Inf,
      seed = 20260826
    )
}

save_plot_safe(p_celltype, "TCell_celltype_UMAP_revision", width = 5.6, height = 4.6)

facet_groups <- c("NR", "TCMR", "ABMR", "BK")
panel_n <- plot_df %>%
  filter(Group %in% facet_groups) %>%
  count(Group, name = "n_cells") %>%
  mutate(
    Group = factor(as.character(Group), levels = facet_groups),
    Facet_group = paste0(as.character(Group), " (n=", n_cells, ")")
  )
facet_labels <- setNames(panel_n$Facet_group, as.character(panel_n$Group))

facet_df <- plot_df %>%
  filter(Group %in% facet_groups) %>%
  mutate(
    Group = factor(as.character(Group), levels = facet_groups),
    Facet_group = factor(facet_labels[as.character(Group)], levels = facet_labels)
  )

bg_df <- tidyr::crossing(
  plot_df %>% select(cell, UMAP_1, UMAP_2),
  Group = factor(facet_groups, levels = facet_groups)
) %>%
  mutate(Facet_group = factor(facet_labels[as.character(Group)], levels = facet_labels))

p_facet <- ggplot() +
  geom_point(
    data = bg_df,
    aes(UMAP_1, UMAP_2),
    colour = "grey78",
    size = 0.26,
    alpha = 0.58,
    stroke = 0
  ) +
  geom_point(
    data = facet_df,
    aes(UMAP_1, UMAP_2, colour = T_celltype),
    size = 0.58,
    alpha = 0.90,
    stroke = 0
  ) +
  facet_wrap(~Facet_group, nrow = 2) +
  scale_colour_manual(values = tcell_colors, drop = FALSE, name = "T-cell subtype") +
  guides(colour = guide_legend(override.aes = list(size = 5.2, alpha = 1), ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_umap(9.2) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8.8, colour = "black"),
    legend.title = element_text(face = "bold", size = 9.4),
    strip.text = element_text(face = "bold", size = 10.2, colour = "black"),
    panel.spacing.x = unit(0.85, "lines"),
    panel.spacing.y = unit(0.85, "lines")
  )

save_plot_safe(p_facet, "TCell_TCMR_ABMR_celltype_UMAP_facet_revision", width = 7.2, height = 6.4)

sample_group <- plot_df %>%
  filter(!is.na(Group), !is.na(Sample), Sample != "") %>%
  distinct(Sample, Group) %>%
  mutate(Group = factor(as.character(Group), levels = group_levels)) %>%
  arrange(Group, Sample)

sample_order <- sample_group$Sample

tcell_freq <- plot_df %>%
  filter(!is.na(Group), !is.na(Sample), Sample != "", !is.na(T_celltype)) %>%
  count(Group, Sample, T_celltype, name = "n_cells") %>%
  complete(
    nesting(Group, Sample),
    T_celltype = factor(levels(plot_df$T_celltype), levels = levels(plot_df$T_celltype)),
    fill = list(n_cells = 0)
  ) %>%
  group_by(Group, Sample) %>%
  mutate(
    total_t_cells = sum(n_cells),
    frequency = ifelse(total_t_cells > 0, n_cells / total_t_cells, NA_real_),
    frequency_percent = frequency * 100
  ) %>%
  ungroup() %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    Sample = factor(Sample, levels = sample_order),
    T_celltype = factor(as.character(T_celltype), levels = rev(levels(plot_df$T_celltype)))
  )

write_csv_safe(tcell_freq, file.path(composition_result_dir, "TCell_group_sample_subtype_frequency_heatmap_source.csv"))

tcell_freq_group_summary <- tcell_freq %>%
  group_by(T_celltype, Group) %>%
  summarise(
    n_samples = sum(!is.na(frequency_percent)),
    mean_frequency_percent = mean(frequency_percent, na.rm = TRUE),
    median_frequency_percent = median(frequency_percent, na.rm = TRUE),
    sd_frequency_percent = sd(frequency_percent, na.rm = TRUE),
    .groups = "drop"
  )

tcell_freq_kw <- tcell_freq %>%
  group_by(T_celltype) %>%
  summarise(p.value = safe_kw(pick(everything())), .groups = "drop") %>%
  mutate(
    p.adjust = p.adjust(p.value, method = "BH"),
    label = format_p(p.value, prefix = "KW p")
  )

tcell_freq_tcmr_wilcox <- tidyr::crossing(
  T_celltype = levels(tcell_freq$T_celltype),
  group1 = "TCMR",
  group2 = setdiff(group_levels, "TCMR")
) %>%
  rowwise() %>%
  mutate(
    p.value = safe_wilcox_two_groups(
      tcell_freq %>% filter(as.character(T_celltype) == as.character(.data$T_celltype)),
      group1,
      group2
    )
  ) %>%
  ungroup() %>%
  group_by(T_celltype) %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
  ungroup()

tcell_freq_pairwise_no_mix <- bind_rows(lapply(levels(tcell_freq$T_celltype), function(ct) {
  combos <- as.data.frame(t(utils::combn(boxplot_group_levels, 2)), stringsAsFactors = FALSE)
  colnames(combos) <- c("group1", "group2")
  combos$T_celltype <- ct
  combos$p.value <- vapply(seq_len(nrow(combos)), function(i) {
    safe_wilcox_two_groups(
      tcell_freq %>% filter(as.character(T_celltype) == ct),
      combos$group1[i],
      combos$group2[i]
    )
  }, numeric(1))
  combos
})) %>%
  group_by(T_celltype) %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
  ungroup()

freq_label_df <- tcell_freq %>%
  group_by(T_celltype) %>%
  summarise(
    y_min = min(frequency_percent, na.rm = TRUE),
    y_max = max(frequency_percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(tcell_freq_kw, by = "T_celltype") %>%
  mutate(
    y_span = pmax(y_max - y_min, 5),
    y = y_max + y_span * 0.12,
    x = 3
  )

write_csv_safe(tcell_freq_group_summary, file.path(composition_result_dir, "TCell_group_sample_subtype_frequency_group_summary.csv"))
write_csv_safe(tcell_freq_kw, file.path(composition_result_dir, "TCell_group_sample_subtype_frequency_KW_tests.csv"))
write_csv_safe(tcell_freq_tcmr_wilcox, file.path(composition_result_dir, "TCell_group_sample_subtype_frequency_TCMR_vs_other_wilcox_tests.csv"))
write_csv_safe(tcell_freq_pairwise_no_mix, file.path(composition_result_dir, "TCell_group_sample_subtype_frequency_pairwise_wilcox_no_Mix_tests.csv"))

tcell_freq_group_heat <- tcell_freq_group_summary %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    T_celltype = factor(as.character(T_celltype), levels = rev(levels(plot_df$T_celltype)))
  )
write_csv_safe(tcell_freq_group_heat, file.path(composition_result_dir, "TCell_group_mean_subtype_frequency_heatmap_source.csv"))

p_freq_heat <- ggplot(tcell_freq_group_heat, aes(Group, T_celltype, fill = mean_frequency_percent)) +
  geom_tile(colour = "white", linewidth = 0.28) +
  scale_fill_gradient2(
    low = "#3B76AF",
    mid = "white",
    high = "#B2182B",
    midpoint = median(tcell_freq_group_heat$mean_frequency_percent, na.rm = TRUE),
    limits = c(0, max(tcell_freq_group_heat$mean_frequency_percent, na.rm = TRUE)),
    name = "Mean frequency (%)"
  ) +
  labs(x = "Group", y = "T-cell subtype") +
  theme_classic(base_size = 8.6, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 7.8, colour = "black"),
    axis.text.y = element_text(size = 8.0, colour = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 8.3),
    legend.text = element_text(size = 7.8),
    plot.margin = margin(4, 5, 4, 4)
  )

save_plot_to_dir(
  p_freq_heat,
  "TCell_group_sample_subtype_frequency_heatmap_revision",
  width = 4.4,
  height = 3.3,
  out_dir = composition_figure_dir
)

freq_box_df <- tcell_freq %>%
  filter(as.character(Group) %in% boxplot_group_levels) %>%
  mutate(
    Group = factor(as.character(Group), levels = boxplot_group_levels),
    T_celltype = factor(as.character(T_celltype), levels = levels(plot_df$T_celltype))
  )

freq_y_stats <- freq_box_df %>%
  group_by(T_celltype) %>%
  summarise(
    y_max = max(frequency_percent, na.rm = TRUE),
    y_min = min(frequency_percent, na.rm = TRUE),
    y_range = pmax(y_max - y_min, 0.9),
    .groups = "drop"
  )

freq_wilcox_plot <- tcell_freq_pairwise_no_mix %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  mutate(
    T_celltype = factor(as.character(T_celltype), levels = levels(freq_box_df$T_celltype)),
    group1 = factor(group1, levels = boxplot_group_levels),
    group2 = factor(group2, levels = boxplot_group_levels),
    x = as.numeric(group1),
    xend = as.numeric(group2),
    label = format_wilcox_p(p.value)
  ) %>%
  arrange(T_celltype, p.value, x, xend) %>%
  group_by(T_celltype) %>%
  mutate(stat_rank = row_number()) %>%
  ungroup() %>%
  left_join(freq_y_stats, by = "T_celltype") %>%
  mutate(
    y = y_max + y_range * (0.18 + 0.24 * (stat_rank - 1)),
    yend = y,
    tick_y = y - y_range * 0.055,
    y_text = y + y_range * 0.095
  )
write_csv_safe(freq_wilcox_plot, file.path(composition_result_dir, "TCell_group_sample_subtype_frequency_pairwise_wilcox_no_Mix_p_lt_0.1_revision.csv"))

p_freq_box <- ggplot(freq_box_df, aes(Group, frequency_percent, fill = Group, colour = Group)) +
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
  ) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 20260816),
    shape = 21,
    stroke = 0.52,
    size = 2.10,
    alpha = 0.92
  ) +
  geom_segment(
    data = freq_wilcox_plot,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.36,
    colour = "black"
  ) +
  geom_segment(
    data = freq_wilcox_plot,
    aes(x = x, xend = x, y = tick_y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.36,
    colour = "black"
  ) +
  geom_segment(
    data = freq_wilcox_plot,
    aes(x = xend, xend = xend, y = tick_y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.36,
    colour = "black"
  ) +
  geom_text(
    data = freq_wilcox_plot,
    aes(x = (x + xend) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 3.8
  ) +
  facet_wrap(~T_celltype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, drop = FALSE, name = "Group") +
  scale_colour_manual(values = group_colors, drop = FALSE, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.34))) +
  labs(x = NULL, y = "Frequency (%)") +
  theme_classic(base_size = 10.2, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9.2, colour = "black"),
    axis.text.y = element_text(size = 9.8, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.45),
    axis.ticks = element_line(colour = "black", linewidth = 0.35),
    strip.background = element_rect(fill = "grey93", colour = NA),
    strip.text = element_text(face = "bold", size = 9.2, colour = "black"),
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 8.2),
    legend.text = element_text(size = 7.8),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    axis.title.y = element_text(size = 11.8, face = "bold", margin = margin(r = 5)),
    panel.spacing.x = unit(0.75, "lines"),
    panel.spacing.y = unit(0.72, "lines"),
    plot.margin = margin(4, 8, 8, 12)
  )

save_plot_to_dir(
  p_freq_box,
  "TCell_group_sample_subtype_frequency_boxplot_revision",
  width = 6.8,
  height = 7.4,
  out_dir = composition_figure_dir
)

message("Saved revised urine T-cell celltype UMAP, TCMR/ABMR faceted UMAP, and group-level subtype frequency heatmap/boxplot.")

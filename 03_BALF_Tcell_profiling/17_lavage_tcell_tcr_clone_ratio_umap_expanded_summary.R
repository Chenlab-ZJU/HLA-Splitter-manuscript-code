# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Final BALF figures from saved objects and signatures.
# Original source: 22_revision_lavage_tcell_tcr_clone_ratio_umap_expanded_summary.R
# Figure/table scope: Fig3b; Fig3c; Fig3d; Fig3e; Fig3f; Fig3g; Fig3h; Fig3i; FigS6a; FigS6b; FigS6c; FigS6d; FigS6e; FigS6f; FigS6g; FigS6h
# Data cohort: balf. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Expanded T cells: clone_ratio > 1; distinguish cell counts from clonotype counts.
# Keep within-sample TCR denominators and clonotype identifiers as defined below.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

project_dir <- Sys.getenv("LAVAGE_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "balf")
}
project_dir <- gsub("\\\\", "/", project_dir)
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures", "TCR")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("CAP", "LC", "ILD")
group_colors <- c("CAP" = "#699ECA", "LC" = "#E42320", "ILD" = "#4DAF4A")
clone_gradient <- c("#BDBDBD", "#FEE08B", "#F46D43", "#B2182B")

theme_pub <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey94", colour = NA),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      panel.grid = element_blank(),
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

message("Loading lavage_Tcell.rds...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(file.path(getOption("hla.data_root"), "balf", "lavage_Tcell.rds"))
if (!"umap" %in% names(tcell@reductions)) {
  stop("lavage_Tcell.rds does not contain a UMAP reduction.")
}

meta <- tcell@meta.data %>%
  rownames_to_column("cell")
required_cols <- c("Group", "HLA_Splitter", "clonalFrequency")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

umap <- as.data.frame(Embeddings(tcell, "umap")) %>%
  rownames_to_column("cell")
colnames(umap)[2:3] <- c("UMAP_1", "UMAP_2")

celltype_col <- intersect(c("revised_celltype", "Celltype", "celltype"), colnames(meta))[1]
if (is.na(celltype_col)) {
  meta$Celltype_for_plot <- NA_character_
} else {
  meta$Celltype_for_plot <- as.character(meta[[celltype_col]])
}
clone_id_col <- intersect(c("CTaa", "CTgene", "CTnt", "clonotype_id", "clonotype"), colnames(meta))[1]
if (is.na(clone_id_col)) {
  meta$Clone_ID_for_count <- NA_character_
} else {
  meta$Clone_ID_for_count <- as.character(meta[[clone_id_col]])
}

tcr_df <- meta %>%
  transmute(
    cell,
    Sample = as.character(HLA_Splitter),
    Group = factor(as.character(Group), levels = group_levels),
    clonalFrequency = suppressWarnings(as.numeric(clonalFrequency)),
    existing_clone_ratio = if ("clone_ratio" %in% colnames(meta)) suppressWarnings(as.numeric(clone_ratio)) else NA_real_,
    Celltype = Celltype_for_plot,
    Clone_ID = Clone_ID_for_count
  ) %>%
  filter(!is.na(Group), !is.na(Sample)) %>%
  group_by(Sample) %>%
  mutate(
    sample_baseline_clonal_frequency = mean(clonalFrequency[!is.na(clonalFrequency) & clonalFrequency > 0], na.rm = TRUE),
    clone_ratio_calculated = ifelse(
      is.finite(sample_baseline_clonal_frequency) & sample_baseline_clonal_frequency > 0,
      clonalFrequency / sample_baseline_clonal_frequency,
      0
    ),
    clone_ratio = ifelse(!is.na(existing_clone_ratio), existing_clone_ratio, clone_ratio_calculated),
    clone_ratio = ifelse(is.na(clone_ratio) | !is.finite(clone_ratio), 0, clone_ratio),
    expanded_clone = clone_ratio > 1
  ) %>%
  ungroup() %>%
  left_join(umap, by = "cell") %>%
  mutate(
    plot_ratio = pmin(clone_ratio, stats::quantile(clone_ratio[clone_ratio > 0], 0.995, na.rm = TRUE)),
    plot_size = ifelse(expanded_clone, rescale(pmin(clone_ratio, 8), to = c(0.55, 2.20), from = c(1, 8)), 0.18)
  )

# Export the analysis source data for figure reproduction.
write.csv(tcr_df, file.path(result_dir, "lavage_Tcell_TCR_clone_ratio_cell_level_revision.csv"), row.names = FALSE)

group_summary <- tcr_df %>%
  group_by(Group) %>%
  summarise(
    total_t_cells = n(),
    tcr_detected_cells = sum(clonalFrequency > 0, na.rm = TRUE),
    detected_tcr_clones = n_distinct(Clone_ID[clonalFrequency > 0 & !is.na(Clone_ID) & Clone_ID != ""]),
    expanded_tcr_clones = n_distinct(Clone_ID[expanded_clone & !is.na(Clone_ID) & Clone_ID != ""]),
    expanded_t_cells = sum(expanded_clone, na.rm = TRUE),
    expanded_fraction_all_t_cells = expanded_t_cells / total_t_cells,
    expanded_fraction_tcr_detected = expanded_t_cells / tcr_detected_cells,
    expanded_clone_fraction_detected_clones = expanded_tcr_clones / detected_tcr_clones,
    .groups = "drop"
  ) %>%
  mutate(Group = factor(as.character(Group), levels = group_levels))

sample_summary <- tcr_df %>%
  group_by(Sample, Group) %>%
  summarise(
    total_t_cells = n(),
    tcr_detected_cells = sum(clonalFrequency > 0, na.rm = TRUE),
    detected_tcr_clones = n_distinct(Clone_ID[clonalFrequency > 0 & !is.na(Clone_ID) & Clone_ID != ""]),
    expanded_tcr_clones = n_distinct(Clone_ID[expanded_clone & !is.na(Clone_ID) & Clone_ID != ""]),
    expanded_t_cells = sum(expanded_clone, na.rm = TRUE),
    expanded_fraction_all_t_cells = expanded_t_cells / total_t_cells,
    expanded_fraction_tcr_detected = ifelse(tcr_detected_cells > 0, expanded_t_cells / tcr_detected_cells, NA_real_),
    expanded_clone_fraction_detected_clones = ifelse(detected_tcr_clones > 0, expanded_tcr_clones / detected_tcr_clones, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(Group = factor(as.character(Group), levels = group_levels))

write.csv(group_summary, file.path(result_dir, "lavage_Tcell_TCR_expanded_clone_group_summary_revision.csv"), row.names = FALSE)
write.csv(sample_summary, file.path(result_dir, "lavage_Tcell_TCR_expanded_clone_sample_summary_revision.csv"), row.names = FALSE)

message("Drawing clone ratio UMAP...")
bg_df <- tcr_df %>% filter(!expanded_clone)
fg_df <- tcr_df %>% filter(expanded_clone) %>% arrange(plot_ratio)

p_umap <- ggplot() +
  geom_point(data = bg_df, aes(UMAP_1, UMAP_2), colour = "grey70", size = 0.16, alpha = 0.62, stroke = 0) +
  geom_point(data = fg_df, aes(UMAP_1, UMAP_2, colour = plot_ratio, size = plot_size), alpha = 0.93, stroke = 0) +
  scale_colour_gradientn(
    colours = clone_gradient,
    trans = "sqrt",
    name = "Clone ratio",
    guide = guide_colourbar(title.position = "top", barheight = unit(28, "mm"), barwidth = unit(4.5, "mm"))
  ) +
  scale_size_identity() +
  coord_fixed() +
  labs(
    title = "BALF T cells: TCR clone ratio",
    subtitle = "Expanded clone: clone ratio > 1",
    caption = "clone ratio = clonal frequency / mean non-zero clonal frequency within each sample",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_pub(9.8) +
  theme(
    plot.subtitle = element_text(size = 8.8, colour = "grey25", hjust = 0.5),
    plot.caption = element_text(size = 6.8, colour = "grey30", hjust = 0.5),
    legend.position = "right"
  )

p_expanded <- ggplot() +
  geom_point(data = bg_df, aes(UMAP_1, UMAP_2), colour = "grey70", size = 0.16, alpha = 0.62, stroke = 0) +
  geom_point(data = fg_df, aes(UMAP_1, UMAP_2, colour = Group), size = 0.56, alpha = 0.92, stroke = 0) +
  scale_colour_manual(values = group_colors, drop = FALSE, name = "Group") +
  coord_fixed() +
  labs(title = "Expanded T cells", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(
    legend.position = "right",
    legend.key.size = unit(5.5, "mm"),
    legend.text = element_text(size = 9.5),
    legend.title = element_text(size = 10, face = "bold")
  )

save_plot_safe(p_umap, "TCR_clone_ratio_UMAP_revision", width = 5.3, height = 4.8)
save_plot_safe(p_expanded, "TCR_expanded_clone_UMAP_revision", width = 5.3, height = 4.8)

group_summary_long <- group_summary %>%
  select(Group, total_t_cells, tcr_detected_cells, expanded_t_cells, expanded_fraction_all_t_cells) %>%
  mutate(label = paste0(expanded_t_cells, "/", total_t_cells, "\n", percent(expanded_fraction_all_t_cells, accuracy = 0.1)))

p_count <- ggplot(group_summary, aes(Group, expanded_t_cells, fill = Group)) +
  geom_col(width = 0.58, colour = "black", linewidth = 0.35, alpha = 0.82) +
  geom_text(aes(label = expanded_t_cells), vjust = -0.35, family = "Arial", fontface = "bold", size = 3.2) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Expanded T-cell number", x = NULL, y = "Expanded T cells") +
  theme_pub(9.8) +
  theme(legend.position = "none")

p_clone_count <- ggplot(group_summary, aes(Group, expanded_tcr_clones, fill = Group)) +
  geom_col(width = 0.58, colour = "black", linewidth = 0.35, alpha = 0.82) +
  geom_text(aes(label = expanded_tcr_clones), vjust = -0.35, family = "Arial", fontface = "bold", size = 3.2) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Expanded TCR clonotypes", x = NULL, y = "Expanded clonotypes") +
  theme_pub(9.8) +
  theme(legend.position = "none")

p_fraction <- ggplot(group_summary, aes(Group, expanded_fraction_all_t_cells, fill = Group)) +
  geom_col(width = 0.58, colour = "black", linewidth = 0.35, alpha = 0.82) +
  geom_text(aes(label = percent(expanded_fraction_all_t_cells, accuracy = 0.1)), vjust = -0.35, family = "Arial", fontface = "bold", size = 3.2) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Expanded T-cell fraction", x = NULL, y = "Expanded T cells / all T cells") +
  theme_pub(9.8) +
  theme(legend.position = "none")

p_clone_fraction <- ggplot(group_summary, aes(Group, expanded_clone_fraction_detected_clones, fill = Group)) +
  geom_col(width = 0.58, colour = "black", linewidth = 0.35, alpha = 0.82) +
  geom_text(aes(label = percent(expanded_clone_fraction_detected_clones, accuracy = 0.1)), vjust = -0.35, family = "Arial", fontface = "bold", size = 3.2) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Expanded clonotype fraction", x = NULL, y = "Expanded clonotypes / detected clonotypes") +
  theme_pub(9.8) +
  theme(legend.position = "none")

p_summary <- (p_count + p_fraction) / (p_clone_count + p_clone_fraction) + plot_layout(heights = c(1, 1))
save_plot_safe(p_summary, "TCR_expanded_clone_count_fraction_by_group_revision", width = 6.8, height = 6.0)

message("Saved TCR clone ratio UMAP and expanded clone summaries.")

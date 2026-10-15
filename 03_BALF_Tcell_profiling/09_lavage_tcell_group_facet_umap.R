# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Final BALF figures from saved objects and signatures.
# Original source: 17_revision_lavage_tcell_group_facet_umap.R
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
  tcell_colors_from_file <- color_map %>%
    filter(category == "T_celltype") %>%
    select(label, color)
  if (nrow(tcell_colors_from_file) > 0) {
    tcell_type_colors[tcell_colors_from_file$label] <- tcell_colors_from_file$color
  }
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

theme_umap <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = base_size + 1.2),
      strip.background = element_rect(fill = "grey94", colour = NA),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 0.3),
      legend.key = element_blank(),
      legend.key.height = unit(0.19, "in"),
      legend.key.width = unit(0.19, "in"),
      panel.spacing.x = unit(0.18, "in"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

message("Loading lavage_Tcell.rds...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(file.path(getOption("hla.data_root"), "balf", "lavage_Tcell.rds"))
if (!"umap" %in% Reductions(tcell)) {
  stop("lavage_Tcell.rds does not contain a UMAP reduction.")
}
required_cols <- c("revised_celltype", "Group")
missing_cols <- setdiff(required_cols, colnames(tcell@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

umap_df <- as.data.frame(Embeddings(tcell, "umap")[, 1:2])
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
plot_df <- bind_cols(
  umap_df,
  tcell@meta.data %>%
    transmute(
      Group = factor(as.character(Group), levels = group_levels),
      T_celltype = factor(as.character(revised_celltype), levels = tcell_type_levels)
    )
) %>%
  filter(!is.na(Group), !is.na(T_celltype))

type_colors_present <- tcell_type_colors[levels(droplevels(plot_df$T_celltype))]
facet_labels <- plot_df %>%
  count(Group, name = "n_cells") %>%
  mutate(label = paste0(Group, " (n=", format(n_cells, big.mark = ","), ")")) %>%
  select(Group, label)
facet_labeller <- setNames(facet_labels$label, facet_labels$Group)

# Export the analysis source data for figure reproduction.
write.csv(
  plot_df %>% count(Group, T_celltype, name = "n_cells") %>% arrange(Group, T_celltype),
  file.path(result_dir, "lavage_Tcell_group_facet_umap_cell_counts.csv"),
  row.names = FALSE
)

group_celltype_prop_df <- plot_df %>%
  count(Group, T_celltype, name = "n_cells") %>%
  group_by(Group) %>%
  mutate(
    total_cells = sum(n_cells),
    fraction = n_cells / total_cells,
    percent_label = scales::percent(fraction, accuracy = 0.1)
  ) %>%
  ungroup()

write.csv(
  group_celltype_prop_df,
  file.path(result_dir, "lavage_Tcell_group_facet_umap_group_celltype_proportion_revision.csv"),
  row.names = FALSE
)

background_df <- plot_df %>%
  select(UMAP_1, UMAP_2) %>%
  distinct()

p_facet <- ggplot() +
  geom_point(
    data = background_df,
    aes(UMAP_1, UMAP_2),
    colour = "grey82",
    size = 0.20,
    alpha = 0.36,
    stroke = 0
  ) +
  geom_point(
    data = plot_df,
    aes(UMAP_1, UMAP_2, colour = T_celltype),
    size = 0.32,
    alpha = 0.82,
    stroke = 0
  ) +
  facet_wrap(
    ~Group,
    nrow = 1,
    labeller = labeller(Group = facet_labeller)
  ) +
  coord_fixed(ratio = 1, expand = FALSE) +
  scale_colour_manual(values = type_colors_present, drop = TRUE, name = "T cell subset") +
  guides(colour = "none") +
  theme_umap(9.6) +
  theme(
    legend.position = "none"
  )

p_bar <- ggplot(group_celltype_prop_df, aes(Group, fraction, fill = T_celltype)) +
  geom_col(width = 0.50, colour = "white", linewidth = 0.22, alpha = 0.92) +
  scale_fill_manual(values = type_colors_present, drop = TRUE, name = "T cell subset") +
  scale_x_discrete(expand = expansion(mult = c(0.04, 0.04))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Cells (%)") +
  guides(fill = guide_legend(override.aes = list(alpha = 1), ncol = 1)) +
  theme_classic(base_size = 9.6, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    axis.title.y = element_text(face = "bold", size = 10.2),
    axis.text.x = element_text(size = 9.4, colour = "black"),
    axis.text.y = element_text(size = 8.8, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.35),
    axis.ticks = element_line(colour = "black", linewidth = 0.30),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9.0, face = "bold"),
    legend.text = element_text(size = 8.4),
    legend.key.size = unit(3.5, "mm"),
    legend.spacing.y = unit(0.6, "mm"),
    plot.margin = margin(30, 2, 2, 2)
  )

p_combined <- p_facet + p_bar + plot_layout(widths = c(5.0, 0.95))

save_plot_safe(p_combined, "TCell_UMAP_Group_Facet", width = 9.7, height = 3.55)
message("Saved TCell_UMAP_Group_Facet.pdf and .png")

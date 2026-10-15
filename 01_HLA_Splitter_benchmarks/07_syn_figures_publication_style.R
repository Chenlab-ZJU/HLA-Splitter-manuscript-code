# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Final benchmark / HLA-signal figures.
# Original source: 13_revision_syn_figures_publication_style.R
# Figure/table scope: Fig1e; Fig1f; Fig1g; Fig1h; FigS2a; FigS2b; FigS2c
# Data cohort: synthetic. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# Publication-style redraws for selected Syn_data figures.
# Outputs overwrite the requested PDF files and also save PNG previews.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggridges)
  library(scales)
})

project_dir <- file.path(getOption("hla.data_root"), "synthetic")
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

class_cols <- c(
  Correct = "#4DAF4A",
  Undefined = "#BDBDBD",
  Incorrect = "#E42320"
)

ridge_class_cols <- c(
  Correct = "#CEDE94",
  Undefined = "grey",
  Incorrect = "firebrick"
)

metric_cols <- c(
  HLA_count = "#7F7F7F",
  Accuracy = "#B2182B"
)

mixedsort_unique <- function(x) {
  x <- unique(as.character(x))
  if (requireNamespace("gtools", quietly = TRUE)) {
    gtools::mixedsort(x)
  } else {
    x[order(as.numeric(gsub("\\D+", "", x)), x)]
  }
}

theme_nc <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, color = "black", size = base_size + 1),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.40),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "grey96", color = "grey75", linewidth = 0.35),
      strip.text = element_text(face = "bold", color = "black"),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5.5, 7, 5.5, 5.5)
    )
}

save_nc_plot <- function(plot, filename, width, height, dpi = 600) {
  pdf_file <- file.path(figure_dir, paste0(filename, ".pdf"))
  png_file <- file.path(figure_dir, paste0(filename, ".png"))
  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white")
}

make_half_violin_density <- function(data, value_col, group_col, x_map, side = c("left", "right"),
                                     width = 0.34, n = 256) {
  side <- match.arg(side)
  split_data <- split(data[[value_col]], data[[group_col]])
  out <- lapply(names(split_data), function(g) {
    values <- split_data[[g]]
    values <- values[is.finite(values)]
    if (length(unique(values)) < 2) {
      return(NULL)
    }
    den <- density(values, n = n, na.rm = TRUE)
    x0 <- x_map[[g]]
    scaled <- den$y / max(den$y, na.rm = TRUE) * width
    x_density <- if (side == "left") x0 - scaled else x0 + scaled
    data.frame(
      cell_type = g,
      x = c(x0, x_density, x0),
      y = c(min(den$x), den$x, max(den$x)),
      side = side
    )
  })
  bind_rows(out)
}

# Read the saved object or source table used by this analysis.
demux <- read.csv(file.path(getOption("hla.data_root"), "synthetic", "demultiplex_result.csv"), row.names = 1, check.names = FALSE)
demux <- demux %>%
  mutate(class = factor(class, levels = c("Correct", "Undefined", "Incorrect"))) %>%
  filter(!is.na(class), !is.na(label), is.finite(HLA_count))

demux$label <- factor(demux$label, levels = mixedsort_unique(demux$label))

predict_fraction_df <- demux %>%
  count(label, class, name = "n") %>%
  group_by(label) %>%
  mutate(Fraction_pct = n / sum(n) * 100) %>%
  ungroup() %>%
  complete(label, class, fill = list(n = 0, Fraction_pct = 0)) %>%
  mutate(class = factor(class, levels = c("Correct", "Undefined", "Incorrect")))

predict_fraction_summary <- predict_fraction_df %>%
  group_by(class) %>%
  summarise(Mean_pct = mean(Fraction_pct), .groups = "drop")

p_predict_fraction <- ggplot(predict_fraction_df, aes(x = class, y = Fraction_pct)) +
  geom_boxplot(aes(fill = class), width = 0.58, outlier.shape = NA, linewidth = 0.48, alpha = 0.82) +
  geom_point(
    aes(fill = class),
    position = position_jitter(width = 0.10, height = 0),
    shape = 21,
    size = 1.8,
    stroke = 0.28,
    color = "black",
    alpha = 0.92
  ) +
  geom_text(
    data = predict_fraction_summary,
    aes(x = class, y = pmin(103, Mean_pct + 7), label = sprintf("%.1f%%", Mean_pct)),
    inherit.aes = FALSE,
    size = 2.7,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = class_cols, drop = FALSE) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 108), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Cells per sample (%)", title = "Prediction categories") +
  theme_nc(base_size = 8.8) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", angle = 25, hjust = 1)
  )

save_nc_plot(p_predict_fraction, "predict_fraction", width = 3.35, height = 3.45)
# Export the analysis source data for figure reproduction.
write.csv(predict_fraction_df, file.path(result_dir, "predict_fraction_publication_source.csv"), row.names = FALSE)

p_hla_ridge <- ggplot(demux, aes(x = HLA_count, y = class, fill = class)) +
  geom_density_ridges(
    scale = 1.05,
    rel_min_height = 0.01,
    alpha = 0.80,
    color = "white",
    linewidth = 0.35,
    quantile_lines = TRUE,
    quantiles = 2,
    vline_linetype = "dashed",
    vline_width = 0.28
  ) +
  scale_fill_manual(values = ridge_class_cols, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = expression(Log[10]~"(HLA counts + 1)"), y = NULL, title = "HLA abundance by prediction category") +
  theme_ridges(font_size = 8.8, grid = TRUE, center_axis_labels = TRUE) +
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0, size = 9.8, color = "black"),
    axis.title.x = element_text(face = "bold", size = 8.8, color = "black"),
    axis.text = element_text(size = 8.4, color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.45),
    axis.ticks.x = element_line(color = "black", linewidth = 0.40),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(5.5, 7, 5.5, 5.5)
  )

save_nc_plot(p_hla_ridge, "HLA_count_ridgeplot", width = 3.55, height = 2.75)

ct_stats <- demux %>%
  filter(!is.na(cell_type), cell_type != "") %>%
  group_by(cell_type) %>%
  summarise(
    n_cells = n(),
    correct_pct = mean(class == "Correct", na.rm = TRUE) * 100,
    median_hla = median(HLA_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 20) %>%
  arrange(desc(median_hla))

demux_ct <- demux %>%
  filter(cell_type %in% ct_stats$cell_type) %>%
  mutate(cell_type = factor(cell_type, levels = ct_stats$cell_type))

ct_stats <- ct_stats %>%
  mutate(cell_type = factor(cell_type, levels = levels(demux_ct$cell_type)))

set.seed(20260813)
demux_ct_points <- demux_ct %>%
  group_by(cell_type) %>%
  group_modify(~ {
    n_keep <- min(nrow(.x), 70)
    .x[sample.int(nrow(.x), n_keep), , drop = FALSE]
  }) %>%
  ungroup()

max_hla <- max(demux_ct$HLA_count, na.rm = TRUE)
scale_factor <- max_hla / 100

accuracy_sample <- demux_ct %>%
  group_by(cell_type, label) %>%
  summarise(
    n_cells = n(),
    accuracy_pct = mean(class == "Correct", na.rm = TRUE) * 100,
    accuracy_scaled = accuracy_pct * scale_factor,
    .groups = "drop"
  ) %>%
  filter(n_cells >= 5)

x_levels <- levels(demux_ct$cell_type)
x_map <- setNames(seq_along(x_levels), x_levels)

hla_half <- make_half_violin_density(
  data = mutate(demux_ct, cell_type_chr = as.character(cell_type)),
  value_col = "HLA_count",
  group_col = "cell_type_chr",
  x_map = x_map,
  side = "left",
  width = 0.32
)

acc_half <- make_half_violin_density(
  data = mutate(accuracy_sample, cell_type_chr = as.character(cell_type)),
  value_col = "accuracy_scaled",
  group_col = "cell_type_chr",
  x_map = x_map,
  side = "right",
  width = 0.32
)

half_medians <- bind_rows(
  demux_ct %>%
    group_by(cell_type) %>%
    summarise(y = median(HLA_count, na.rm = TRUE), .groups = "drop") %>%
    mutate(x = x_map[as.character(cell_type)], xmin = x - 0.30, xmax = x - 0.05, metric = "HLA count"),
  accuracy_sample %>%
    group_by(cell_type) %>%
    summarise(y = median(accuracy_scaled, na.rm = TRUE), .groups = "drop") %>%
    mutate(x = x_map[as.character(cell_type)], xmin = x + 0.05, xmax = x + 0.30, metric = "Accuracy")
)

demux_ct_points <- demux_ct_points %>%
  mutate(x = x_map[as.character(cell_type)] - runif(n(), 0.04, 0.29))

accuracy_points <- accuracy_sample %>%
  mutate(x = x_map[as.character(cell_type)] + runif(n(), 0.04, 0.29))

p_dual_axis <- ggplot() +
  geom_polygon(
    data = hla_half,
    aes(x = x, y = y, group = interaction(cell_type, side), fill = "HLA count"),
    color = "#4D4D4D",
    linewidth = 0.28,
    alpha = 0.78
  ) +
  geom_polygon(
    data = acc_half,
    aes(x = x, y = y, group = interaction(cell_type, side), fill = "Accuracy"),
    color = "#8F1424",
    linewidth = 0.28,
    alpha = 0.78
  ) +
  geom_point(
    data = demux_ct_points,
    aes(x = x, y = HLA_count),
    color = "#5F5F5F",
    alpha = 0.10,
    size = 0.22
  ) +
  geom_point(
    data = accuracy_points,
    aes(x = x, y = accuracy_scaled),
    color = metric_cols[["Accuracy"]],
    alpha = 0.82,
    size = 0.80
  ) +
  geom_segment(
    data = half_medians,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    linewidth = 0.45,
    color = "black"
  ) +
  scale_fill_manual(
    values = c("HLA count" = "#BDBDBD", "Accuracy" = metric_cols[["Accuracy"]]),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = seq_along(x_levels),
    labels = x_levels,
    expand = expansion(add = c(0.45, 0.45))
  ) +
  scale_y_continuous(
    name = expression(Log[10]~"(HLA counts + 1)"),
    expand = expansion(mult = c(0.02, 0.06)),
    sec.axis = sec_axis(~ . / scale_factor, name = "Correctly assigned cells (%)")
  ) +
  labs(x = "Cell type", title = "HLA abundance and demultiplexing accuracy by cell type") +
  theme_nc(base_size = 8.8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.title.y.right = element_text(color = metric_cols[["Accuracy"]], face = "bold"),
    axis.text.y.right = element_text(color = metric_cols[["Accuracy"]]),
    legend.position = c(0.82, 0.88),
    legend.background = element_rect(fill = alpha("white", 0.86), color = NA),
    legend.key.size = unit(7, "pt")
  )

save_nc_plot(p_dual_axis, "HLA_vs_Accuracy_DualAxis", width = 5.05, height = 3.65)
write.csv(accuracy_sample, file.path(result_dir, "HLA_vs_accuracy_celltype_sample_level_publication_source.csv"), row.names = FALSE)
write.csv(ct_stats, file.path(result_dir, "HLA_vs_accuracy_celltype_publication_source.csv"), row.names = FALSE)

message("Saved publication-style figures to: ", figure_dir)

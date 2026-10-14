# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Final BALF figures from saved objects and signatures.
# Original source: 23_revision_lavage_demultiplex_hla_tsne_predict_sample_number.R
# Figure/table scope: Fig3b; Fig3c; Fig3d; Fig3e; Fig3f; Fig3g; Fig3h; Fig3i; FigS6a; FigS6b; FigS6c; FigS6d; FigS6e; FigS6f; FigS6g; FigS6h
# Data cohort: 肺泡灌洗液. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(Rtsne)
  library(scales)
  library(grid)
})

project_dir <- Sys.getenv("LAVAGE_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "\u80ba\u6ce1\u704c\u6d17\u6db2")
}
project_dir <- gsub("\\\\", "/", project_dir)
data_dir <- file.path(project_dir, "Data")
result_dir <- file.path(project_dir, "Result")
demux_dir <- file.path(result_dir, "demultiplex")
figure_dir <- file.path(project_dir, "Figures", "Demultiplex")
dir.create(demux_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

batch_names <- c("Batch1", "Batch2", "Batch3", "Batch4")
group_levels <- c("CAP", "LC", "ILD", "Others")
group_colors <- c("CAP" = "#699ECA", "LC" = "#E42320", "ILD" = "#4DAF4A", "Others" = "#7A7A7A")
base_sample_colors <- c(
  "#E42320", "#699ECA", "#4DAF4A", "#FF8C00", "#F898CB", "#D65190", "#731A73",
  "#FFCB5B", "#E89B1E", "#0076B9", "#3D505A", "#0098B2", "#7E6148", "#B09C85",
  "#00A087", "#91D1C2", "#8491B4", "#DC0000", "#7CAE00", "#C77CFF", "#00BFC4",
  "#F8766D", "#619CFF", "#00BA38", "#B79F00", "#00B0F6", "#F564E3", "#A3A500",
  "#39B600", "#9590FF", "#00BF7D", "#D89000", "#FF62BC"
)

natural_sample_sort <- function(x) {
  x <- unique(as.character(x))
  idx <- suppressWarnings(as.integer(sub("^P", "", x)))
  x[order(is.na(idx), idx, x)]
}

save_plot_safe <- function(plot, path_no_ext, width, height) {
  tmp_pdf <- file.path(tempdir(), paste0(basename(path_no_ext), "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(basename(path_no_ext), "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
  file.copy(tmp_pdf, paste0(path_no_ext, ".pdf"), overwrite = TRUE)
  file.copy(tmp_png, paste0(path_no_ext, ".png"), overwrite = TRUE)
}

theme_pub <- function(base_size = 10.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.42),
      plot.title = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

read_hla_matrix <- function(batch, batch_id) {
  hla_path <- file.path(data_dir, "HLA_counts", batch)
  label_file <- file.path(hla_path, "sorted_genotypes.txt")
  barcode_file <- file.path(hla_path, "output.barcodes.txt")
  mtx_file <- file.path(hla_path, "HLA_type.mtx")
  if (!all(file.exists(label_file, barcode_file, mtx_file))) {
    stop("Missing HLA count files for ", batch)
  }
  labels <- read.delim(label_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
  barcodes <- read.delim(barcode_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
  mat <- as.matrix(Matrix::readMM(mtx_file))
  colnames(mat) <- labels
  rownames(mat) <- paste0(barcodes, "-", batch_id)
  mat
}

read_demux <- function(batch, batch_id) {
  path <- file.path(demux_dir, batch, "HLA_demultiplex.csv")
  if (!file.exists(path)) {
    stop("Missing demultiplex result: ", path)
  }
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  prediction_col <- if ("Adjusted_prediction" %in% colnames(x)) "Adjusted_prediction" else "Sample_prediction"
  x %>%
    transmute(
      Cell.barcode = paste0(.data[["Cell.barcode"]], "-", batch_id),
      Batch = batch,
      HLA_Splitter = as.character(.data[[prediction_col]]),
      Certainty = if ("Certainty" %in% colnames(x)) as.numeric(.data[["Certainty"]]) else NA_real_
    )
}

message("Reading HLA matrices and demultiplex calls...")
hla_mats <- lapply(seq_along(batch_names), function(i) read_hla_matrix(batch_names[i], i))
common_hla <- Reduce(intersect, lapply(hla_mats, colnames))
hla_all <- do.call(rbind, lapply(hla_mats, function(x) x[, common_hla, drop = FALSE]))
demux_all <- bind_rows(lapply(seq_along(batch_names), function(i) read_demux(batch_names[i], i)))

# Read the saved object or source table used by this analysis.
sample_info <- read.csv(file.path(data_dir, "SampleInfo.csv"), stringsAsFactors = FALSE, check.names = FALSE) %>%
  transmute(HLA_Splitter = as.character(.data[["Patient ID"]]), Group = as.character(Group))

plot_cells <- demux_all %>%
  filter(!is.na(HLA_Splitter), !HLA_Splitter %in% c("Undefined", "Doublets", "Undefiend")) %>%
  inner_join(sample_info, by = "HLA_Splitter") %>%
  filter(Group %in% group_levels)

common_cells <- intersect(plot_cells$Cell.barcode, rownames(hla_all))
plot_cells <- plot_cells %>% filter(Cell.barcode %in% common_cells)
hla_use <- hla_all[plot_cells$Cell.barcode, , drop = FALSE]
rownames(plot_cells) <- plot_cells$Cell.barcode
plot_cells <- plot_cells[rownames(hla_use), , drop = FALSE]

message(sprintf("Matched %d demultiplexed cells for plotting.", nrow(plot_cells)))
if (nrow(plot_cells) < 10) {
  stop("Too few matched cells for t-SNE.")
}

lib_size <- rowSums(hla_use)
keep_cells <- lib_size > 0
hla_use <- hla_use[keep_cells, , drop = FALSE]
plot_cells <- plot_cells[keep_cells, , drop = FALSE]

hla_norm <- sweep(hla_use, 1, rowSums(hla_use), "/") * 100
hla_norm <- log1p(hla_norm)
hla_norm[is.na(hla_norm)] <- 0
for (j in seq_len(ncol(hla_norm))) {
  pos <- hla_norm[, j] > 0
  if (sum(pos) > 1) {
    hla_norm[pos, j] <- as.numeric(scale(hla_norm[pos, j], center = FALSE, scale = TRUE))
  }
}

tsne_cache <- file.path(demux_dir, "HLA_TSNE_Adjusted_Prediction_revision_source.csv")
force_tsne <- tolower(Sys.getenv("FORCE_HLA_TSNE", unset = "false")) %in% c("true", "1", "yes")
if (file.exists(tsne_cache) && !force_tsne) {
  message("Using cached HLA t-SNE coordinates: ", tsne_cache)
  tsne_df <- read.csv(tsne_cache, stringsAsFactors = FALSE, check.names = FALSE) %>%
    as_tibble()
} else {
  message("Running HLA t-SNE...")
  set.seed(321)
  perplexity_val <- min(max(20, floor((nrow(hla_norm) - 1) / 3)), 50)
  tsne_out <- Rtsne(
    as.matrix(hla_norm),
    dims = 2,
    pca = TRUE,
    max_iter = 1000,
    theta = 0.4,
    perplexity = perplexity_val,
    verbose = FALSE,
    check_duplicates = FALSE
  )

  tsne_df <- as_tibble(tsne_out$Y, .name_repair = "minimal") %>%
    setNames(c("tSNE1", "tSNE2")) %>%
    bind_cols(plot_cells %>% rownames_to_column("Cell")) %>%
    mutate(HLA_count = log10(rowSums(hla_use) + 1))
# Export the analysis source data for figure reproduction.
  write.csv(tsne_df, tsne_cache, row.names = FALSE)
}

sample_levels <- natural_sample_sort(tsne_df$HLA_Splitter)
sample_colors <- setNames(rep(base_sample_colors, length.out = length(sample_levels)), sample_levels)
tsne_df <- tsne_df %>%
  mutate(
    HLA_Splitter = factor(as.character(HLA_Splitter), levels = sample_levels),
    Group = factor(as.character(Group), levels = group_levels)
  )

centroids <- tsne_df %>%
  group_by(HLA_Splitter) %>%
  summarise(tSNE1 = median(tSNE1), tSNE2 = median(tSNE2), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 20)

x_range <- range(tsne_df$tSNE1, na.rm = TRUE)
y_range <- range(tsne_df$tSNE2, na.rm = TRUE)
x_span <- diff(x_range)
y_span <- diff(y_range)
arrow_origin <- tibble(
  x = x_range[1] + x_span * 0.070,
  y = y_range[1] + y_span * 0.080
)
arrow_df <- bind_rows(
  arrow_origin %>% mutate(xend = x + x_span * 0.115, yend = y, label = "t-SNE 1", lx = x + x_span * 0.170, ly = y - y_span * 0.003),
  arrow_origin %>% mutate(xend = x, yend = y + y_span * 0.130, label = "t-SNE 2", lx = x, ly = y + y_span * 0.165)
)

p_tsne <- ggplot(tsne_df, aes(tSNE1, tSNE2, colour = HLA_Splitter)) +
  geom_point(size = 0.18, alpha = 0.66, stroke = 0) +
  ggrepel::geom_text_repel(
    data = centroids,
    aes(tSNE1, tSNE2, label = HLA_Splitter),
    inherit.aes = FALSE,
    family = "Arial",
    size = 2.4,
    colour = "black",
    box.padding = unit(0.16, "lines"),
    point.padding = unit(0.10, "lines"),
    segment.color = "grey62",
    segment.size = 0.16,
    max.overlaps = Inf
  ) +
  geom_segment(
    data = arrow_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.45,
    arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed")
  ) +
  geom_text(
    data = arrow_df,
    aes(x = lx, y = ly, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    size = 2.4,
    colour = "black",
    fontface = "bold"
  ) +
  scale_colour_manual(values = sample_colors, drop = FALSE) +
  labs(x = "t-SNE 1", y = "t-SNE 2") +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.margin = margin(3, 3, 3, 3)
  )

cell_count_df <- tsne_df %>%
  count(Group, HLA_Splitter, name = "Cellnum") %>%
  complete(Group = factor(group_levels, levels = group_levels), HLA_Splitter = factor(sample_levels, levels = sample_levels), fill = list(Cellnum = 0)) %>%
  filter(Cellnum > 0) %>%
  mutate(
    HLA_Splitter = factor(as.character(HLA_Splitter), levels = sample_levels),
    Group = factor(as.character(Group), levels = group_levels)
  )
write.csv(cell_count_df, file.path(demux_dir, "Predict_sample_number_revision_source.csv"), row.names = FALSE)

y_breaks <- c(10, 30, 100, 300, 1000, 3000)
p_count <- ggplot(cell_count_df, aes(HLA_Splitter, Cellnum, fill = Group)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.18, alpha = 0.82) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey35", linewidth = 0.38) +
  geom_text(aes(label = Cellnum), vjust = -0.30, size = 2.5, family = "Arial", colour = "black") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_log10(
    breaks = y_breaks,
    labels = comma,
    expand = expansion(mult = c(0.02, 0.16))
  ) +
  labs(x = "Predicted sample", y = "Number of cells", fill = "Group") +
  theme_pub(9.5) +
  theme(
    axis.title = element_text(size = 10.8, face = "bold"),
    axis.text.x = element_text(size = 7.0, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8.8),
    legend.position = "right",
    legend.title = element_text(size = 9.6, face = "bold"),
    legend.text = element_text(size = 9.0),
    legend.key.size = unit(0.36, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

save_plot_safe(p_tsne, file.path(demux_dir, "HLA_TSNE_Adjusted_Prediction"), width = 4.7, height = 4.2)
save_plot_safe(p_tsne, file.path(figure_dir, "HLA_TSNE_Adjusted_Prediction"), width = 4.7, height = 4.2)
save_plot_safe(p_count, file.path(demux_dir, "Predict_sample_number"), width = 7.4, height = 3.1)
save_plot_safe(p_count, file.path(figure_dir, "Predict_sample_number"), width = 7.4, height = 3.1)

message("Saved HLA_TSNE_Adjusted_Prediction and Predict_sample_number.")

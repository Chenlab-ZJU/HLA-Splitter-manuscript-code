# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final urine immune-cell composition, module and demultiplexing figures.
# Original source: 30_revision_demultiplex_HLA_tsne_predict_sample_number.R
# Figure/table scope: Fig4b; Fig4c; Fig4d; Fig4e; Fig4f; Fig4g; Fig4h; FigS7a; FigS7b; FigS7c; FigS7d; FigS7e
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(Rtsne)
  library(patchwork)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "30_revision_demultiplex_HLA_tsne_predict_sample_number.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

data_dir <- file.path(project_dir, "Data")
figure_dir <- file.path(project_dir, "Figures", "Demultiplex")
result_dir <- file.path(project_dir, "Result", "Demultiplex")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK", "Undefined")
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "Mix" = "#FF8C00",
  "BK" = "#F898CB",
  "Undefined" = "#9A9A9A"
)

urine_single_cell_colors <- c(
  "#E42320", "#699ECA", "#4DAF4A", "#FF8C00", "#F898CB", "#D65190",
  "#731A73", "#FFCB5B", "#E89B1E", "#0076B9", "#3D505A", "#0098B2",
  "#8A7355", "#CC5B45", "#008B8B", "#4E5689", "#F9B3AD", "#EDC66A",
  "#9FDAF7", "#C9A1CA", "#D9DEE7", "#C74546", "#EB7E60", "#8FB4DC"
)

make_urine_palette <- function(n) {
  if (n <= length(urine_single_cell_colors)) {
    urine_single_cell_colors[seq_len(n)]
  } else {
    grDevices::colorRampPalette(urine_single_cell_colors)(n)
  }
}

theme_pub <- function(base_size = 8.8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.36),
      axis.ticks = element_line(colour = "black", linewidth = 0.32),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 6, 4, 6)
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

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

sample_sort_key <- function(x) {
  prefix <- sub("_.*$", "", x)
  num <- suppressWarnings(as.integer(sub("^.*_", "", x)))
  prefix_rank <- dplyr::case_when(
    prefix == "NR" ~ 1L,
    prefix == "RE" ~ 2L,
    prefix == "BK" ~ 3L,
    x == "Undefined" ~ 9L,
    TRUE ~ 8L
  )
  prefix_rank * 1000 + dplyr::coalesce(num, 999L)
}

read_predict_number <- function() {
  predict_file <- file.path(getOption("hla.data_root"), "urine", "Predict_sample_number.csv")
  if (!file.exists(predict_file)) return(NULL)
  raw <- read.csv(predict_file, stringsAsFactors = FALSE, check.names = FALSE)
  raw <- raw[, nzchar(names(raw)), drop = FALSE]
  if (!all(c("Group", "predict_label", "Cellnum") %in% colnames(raw))) {
    stop("Predict_sample_number.csv does not contain Group, predict_label and Cellnum columns.")
  }
  raw %>%
    transmute(
      Group = as.character(Group),
      predict_label = as.character(predict_label),
      Cellnum = suppressWarnings(as.numeric(Cellnum))
    ) %>%
    filter(!is.na(Group), !is.na(predict_label), !is.na(Cellnum)) %>%
    mutate(
      Group = factor(Group, levels = group_levels),
      predict_label = factor(predict_label, levels = predict_label[order(sample_sort_key(predict_label))])
    )
}

plot_predict_sample_number <- function(df_plot) {
  p <- ggplot(df_plot, aes(predict_label, Cellnum, fill = Group)) +
    geom_col(width = 0.62, colour = "grey25", linewidth = 0.16, alpha = 0.82) +
    geom_text(
      aes(label = scales::comma(Cellnum)),
      vjust = -0.35,
      family = "Arial",
      fontface = "bold",
      size = 2.25,
      colour = "black"
    ) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "#234A7C", linewidth = 0.42) +
    annotate("text", x = 0.55, y = 58, label = "50 cells", hjust = 0, vjust = 0,
             family = "Arial", fontface = "bold", size = 2.5, colour = "#234A7C") +
    scale_fill_manual(values = group_colors, drop = FALSE, name = "Group") +
    scale_y_log10(
      breaks = c(10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000),
      labels = label_number(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.18))
    ) +
    labs(x = "Predicted sample", y = "Cell number") +
    theme_pub(8.6) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.6),
      axis.text.y = element_text(size = 7.4),
      axis.title = element_text(size = 9.0, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 7.4),
      legend.title = element_text(size = 7.8, face = "bold"),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25)
    )
  save_plot_safe(p, "Predict_sample_number", width = 9.2, height = 2.45)
}

read_hla_matrix_one <- function(base_dir, suffix) {
  genotype_file <- file.path(base_dir, "sorted_genotypes.txt")
  barcode_file <- file.path(base_dir, "filter_counts", "output.barcodes.txt")
  matrix_file <- file.path(base_dir, "filter_counts", "HLA_type.mtx")
  if (!all(file.exists(c(genotype_file, barcode_file, matrix_file)))) return(NULL)
  genes <- read.delim(genotype_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  barcodes <- read.delim(barcode_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  mat <- Matrix::readMM(matrix_file)
  if (ncol(mat) != length(genes)) stop("HLA matrix column number does not match sorted_genotypes.txt: ", base_dir)
  if (nrow(mat) != length(barcodes)) stop("HLA matrix row number does not match output.barcodes.txt: ", base_dir)
  barcodes <- paste0(substr(barcodes, 1, 16), suffix)
  mat <- as.matrix(mat)
  colnames(mat) <- genes
  rownames(mat) <- barcodes
  mat
}

combine_hla_matrices <- function(mats) {
  mats <- mats[!vapply(mats, is.null, logical(1))]
  all_genes <- unique(unlist(lapply(mats, colnames), use.names = FALSE))
  aligned <- lapply(mats, function(mat) {
    out <- matrix(0, nrow = nrow(mat), ncol = length(all_genes), dimnames = list(rownames(mat), all_genes))
    out[, colnames(mat)] <- mat
    out
  })
  do.call(rbind, aligned)
}

build_hla_tsne <- function() {
  object_path <- file.path(getOption("hla.data_root"), "urine", "Urine_ALL.rds")
  if (!file.exists(object_path)) object_path <- file.path(project_dir, "Result", "Objects", "Urine_Allcells.rds")
  if (!file.exists(object_path)) stop("Cannot find Urine_ALL.rds or Urine_Allcells.rds.")

  hla_mats <- list(
    read_hla_matrix_one(file.path(data_dir, "Group1", "HLA_counts"), "-B1"),
    read_hla_matrix_one(file.path(data_dir, "Group2", "HLA_counts", "Group1"), "-B2.1"),
    read_hla_matrix_one(file.path(data_dir, "Group2", "HLA_counts", "Group2"), "-B2.2"),
    read_hla_matrix_one(file.path(data_dir, "Group2", "HLA_counts", "Group3"), "-B2.3"),
    read_hla_matrix_one(file.path(data_dir, "Group3", "HLA_counts"), "-B3")
  )
  hla_matrix <- combine_hla_matrices(hla_mats)

  message("Loading urine object for HLA tSNE labels...")
  scRNA <- readRDS(object_path)
  meta <- scRNA@meta.data %>%
    as_tibble(rownames = "cell") %>%
    transmute(
      cell,
      predict = as.character(Sample),
      Group = as.character(Group)
    ) %>%
    filter(!is.na(predict), predict != "", !predict %in% c("Doublets", "Doublet"),
           !is.na(Group), Group != "Doublets")

  common_cells <- intersect(meta$cell, rownames(hla_matrix))
  if (length(common_cells) < 100) {
    stop("Too few cells overlap between HLA matrix and urine object: ", length(common_cells))
  }
  meta <- meta %>% filter(cell %in% common_cells)
  hla_matrix <- hla_matrix[meta$cell, , drop = FALSE]

  library_sizes <- rowSums(hla_matrix)
  keep <- library_sizes > 0
  hla_matrix <- hla_matrix[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  library_sizes <- library_sizes[keep]

  hla_norm <- sweep(hla_matrix, 1, library_sizes, "/") * 100
  hla_norm <- log1p(hla_norm)
  for (j in seq_len(ncol(hla_norm))) {
    positive <- hla_norm[, j] > 0
    if (sum(positive) > 1) {
      scaled <- scale(hla_norm[positive, j], center = FALSE, scale = TRUE)[, 1]
      hla_norm[positive, j] <- scaled
    }
  }
  hla_norm[!is.finite(hla_norm)] <- 0

  set.seed(321)
  tsne_out <- Rtsne(
    hla_norm,
    dims = 2,
    pca = TRUE,
    max_iter = 1000,
    theta = 0.4,
    perplexity = 20,
    verbose = FALSE,
    check_duplicates = FALSE
  )

  as_tibble(tsne_out$Y, .name_repair = "minimal") %>%
    setNames(c("tSNE1", "tSNE2")) %>%
    bind_cols(meta %>% dplyr::select(cell, predict, Group)) %>%
    mutate(
      Group = factor(Group, levels = group_levels),
      predict = factor(predict, levels = unique(predict[order(sample_sort_key(predict))])),
      HLA_count = log10(rowSums(hla_matrix) + 1)
    )
}

plot_hla_tsne <- function(tsne_df) {
  sample_levels <- levels(tsne_df$predict)
  sample_cols <- setNames(hue_pal(l = 58, c = 95)(length(sample_levels)), sample_levels)
  sample_cols["Undefined"] <- "#9A9A9A"

  centroid_df <- tsne_df %>%
    filter(predict != "Undefined") %>%
    group_by(predict) %>%
    summarise(
      n_cells = n(),
      tSNE1 = median(tSNE1, na.rm = TRUE),
      tSNE2 = median(tSNE2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_cells >= 50)

  p <- ggplot(tsne_df, aes(tSNE1, tSNE2, colour = predict)) +
    geom_point(size = 0.18, alpha = 0.58, stroke = 0) +
    scale_colour_manual(values = sample_cols, guide = "none") +
    labs(x = "tSNE 1", y = "tSNE 2") +
    theme_classic(base_size = 8.8, base_family = "Arial") +
    theme(
      axis.title = element_text(face = "bold", size = 9.3, colour = "black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.36),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9.4),
      plot.margin = margin(3, 3, 3, 3)
    )
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = centroid_df,
        aes(x = tSNE1, y = tSNE2, label = predict),
        inherit.aes = FALSE,
        family = "Arial",
        fontface = "bold",
        size = 1.85,
        colour = "black",
        min.segment.length = 0,
        segment.colour = "grey45",
        segment.size = 0.16,
        box.padding = 0.18,
        point.padding = 0.10,
        max.overlaps = Inf,
        seed = 20260826
      )
  } else {
    p <- p +
      geom_text(
        data = centroid_df,
        aes(x = tSNE1, y = tSNE2, label = predict),
        inherit.aes = FALSE,
        family = "Arial",
        fontface = "bold",
        size = 1.85,
        colour = "black",
        check_overlap = TRUE
      )
  }
  save_plot_safe(p, "HLA_TSNE_Adjusted_Prediction", width = 5.2, height = 4.4)
  invisible(p)
}

plot_hla_score_projection <- function(tsne_df) {
  score_df <- tsne_df %>%
    mutate(
      HLA_score = HLA_count,
      HLA_score_plot = pmin(HLA_score, quantile(HLA_score, 0.995, na.rm = TRUE))
    ) %>%
    arrange(HLA_score_plot)

# Export the analysis source data for figure reproduction.
  write_csv_safe(
    score_df %>% dplyr::select(cell, predict, Group, tSNE1, tSNE2, HLA_score),
    file.path(result_dir, "HLA_score_projection_on_HLA_TSNE_source.csv")
  )

  p <- ggplot(score_df, aes(tSNE1, tSNE2)) +
    geom_point(colour = "grey88", size = 0.16, alpha = 0.42, stroke = 0) +
    geom_point(aes(colour = HLA_score_plot), size = 0.22, alpha = 0.82, stroke = 0) +
    scale_colour_gradientn(
      colours = c("#D9D9D9", "#F6D95E", "#B2182B"),
      values = scales::rescale(c(
        min(score_df$HLA_score_plot, na.rm = TRUE),
        quantile(score_df$HLA_score_plot, 0.72, na.rm = TRUE),
        max(score_df$HLA_score_plot, na.rm = TRUE)
      )),
      name = "HLA score\nlog10(count + 1)",
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(0.28, "in"),
        barheight = unit(1.35, "in"),
        ticks.colour = "black",
        frame.colour = "black"
      )
    ) +
    labs(x = "tSNE 1", y = "tSNE 2") +
    theme_classic(base_size = 8.8, base_family = "Arial") +
    theme(
      axis.title = element_text(face = "bold", size = 9.3, colour = "black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.36),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9.4),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 7.7),
      legend.text = element_text(size = 7.0, colour = "black"),
      legend.key = element_blank(),
      plot.margin = margin(3, 5, 3, 3)
    )

  save_plot_safe(p, "HLA_score_projection_on_HLA_TSNE_Adjusted_Prediction", width = 5.7, height = 4.4)
  invisible(p)
}

plot_group_projection <- function(tsne_df) {
  group_df <- tsne_df %>%
    mutate(
      Group = factor(as.character(Group), levels = group_levels),
      draw_order = ifelse(Group == "Undefined", 0L, as.integer(Group))
    ) %>%
    arrange(draw_order)

  write_csv_safe(
    group_df %>% dplyr::select(cell, predict, Group, tSNE1, tSNE2, HLA_count),
    file.path(result_dir, "Group_projection_on_HLA_TSNE_source.csv")
  )

  p <- ggplot(group_df, aes(tSNE1, tSNE2)) +
    geom_point(
      data = group_df %>% filter(Group == "Undefined"),
      colour = "grey86",
      size = 0.18,
      alpha = 0.42,
      stroke = 0
    ) +
    geom_point(
      data = group_df %>% filter(Group != "Undefined"),
      aes(colour = Group),
      size = 0.22,
      alpha = 0.76,
      stroke = 0
    ) +
    scale_colour_manual(values = group_colors[names(group_colors) != "Undefined"], drop = FALSE, name = "Group") +
    labs(x = "tSNE 1", y = "tSNE 2") +
    guides(colour = guide_legend(override.aes = list(size = 4.2, alpha = 1))) +
    theme_classic(base_size = 8.8, base_family = "Arial") +
    theme(
      axis.title = element_text(face = "bold", size = 9.3, colour = "black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.36),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9.4),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 8.0),
      legend.text = element_text(size = 7.4, colour = "black"),
      legend.key = element_blank(),
      plot.margin = margin(3, 5, 3, 3)
    )

  save_plot_safe(p, "Group_projection_on_HLA_TSNE_Adjusted_Prediction", width = 5.6, height = 4.4)
  invisible(p)
}

add_umap_arrows <- function(p, df, x_col = "tSNE1", y_col = "tSNE2") {
  x_vals <- df[[x_col]]
  y_vals <- df[[y_col]]
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  y_min <- min(y_vals, na.rm = TRUE)
  y_max <- max(y_vals, na.rm = TRUE)
  x_rng <- x_max - x_min
  y_rng <- y_max - y_min
  x0 <- x_min + x_rng * 0.055
  y0 <- y_min + y_rng * 0.060
  arrow_len_x <- x_rng * 0.135
  arrow_len_y <- y_rng * 0.135

  p +
    annotate(
      "segment",
      x = x0,
      xend = x0 + arrow_len_x,
      y = y0,
      yend = y0,
      linewidth = 0.42,
      colour = "black",
      arrow = arrow(length = unit(0.08, "in"), type = "closed")
    ) +
    annotate(
      "segment",
      x = x0,
      xend = x0,
      y = y0,
      yend = y0 + arrow_len_y,
      linewidth = 0.42,
      colour = "black",
      arrow = arrow(length = unit(0.08, "in"), type = "closed")
    ) +
    annotate(
      "text",
      x = x0 + arrow_len_x * 0.68,
      y = y0 - y_rng * 0.030,
      label = "UMAP 1",
      family = "Arial",
      fontface = "bold",
      size = 2.55,
      colour = "black"
    ) +
    annotate(
      "text",
      x = x0 - x_rng * 0.030,
      y = y0 + arrow_len_y * 0.70,
      label = "UMAP 2",
      angle = 90,
      family = "Arial",
      fontface = "bold",
      size = 2.55,
      colour = "black"
    )
}

theme_embedding_no_axes <- function(base_size = 8.8) {
  theme_void(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9.4),
      legend.title = element_text(face = "bold", colour = "black"),
      legend.text = element_text(colour = "black"),
      legend.key = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

plot_hla_tsne_combined_panel <- function(tsne_df) {
  sample_levels <- unique(as.character(tsne_df$predict))
  sample_levels <- c(setdiff(sample_levels, "Undefined"), intersect(sample_levels, "Undefined"))
  sample_cols <- setNames(
    make_urine_palette(length(sample_levels)),
    sample_levels
  )
  if ("Undefined" %in% names(sample_cols)) sample_cols["Undefined"] <- "#9A9A9A"

  p <- ggplot(tsne_df, aes(tSNE1, tSNE2, colour = predict)) +
    geom_point(size = 0.18, alpha = 0.58, stroke = 0) +
    scale_colour_manual(values = sample_cols, drop = FALSE, name = "Predicted\nsample") +
    coord_equal() +
    labs(title = "Predicted sample") +
    guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1), ncol = 2)) +
    theme_embedding_no_axes(8.8) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 7.6, face = "bold"),
      legend.text = element_text(size = 6.4),
      legend.key.height = unit(0.20, "cm"),
      legend.key.width = unit(0.22, "cm")
    )

  add_umap_arrows(p, tsne_df)
}

plot_group_projection_combined_panel <- function(tsne_df) {
  keep_groups <- c("NR", "TCMR", "ABMR", "Mix", "BK")
  group_df <- tsne_df %>%
    filter(Group %in% keep_groups) %>%
    mutate(Group = factor(as.character(Group), levels = keep_groups)) %>%
    arrange(Group)

  p <- ggplot(group_df, aes(tSNE1, tSNE2, colour = Group)) +
    geom_point(size = 0.22, alpha = 0.76, stroke = 0) +
    scale_colour_manual(values = group_colors[keep_groups], drop = FALSE, name = "Group") +
    coord_equal() +
    labs(title = "Clinical group") +
    guides(colour = guide_legend(override.aes = list(size = 4.2, alpha = 1), ncol = 1)) +
    theme_embedding_no_axes(8.8) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 8.0),
      legend.text = element_text(size = 7.4),
      plot.margin = margin(4, 5, 4, 4)
    )

  add_umap_arrows(p, group_df)
}

plot_combined_hla_tsne <- function(tsne_df) {
  p_predict <- plot_hla_tsne_combined_panel(tsne_df)
  p_group <- plot_group_projection_combined_panel(tsne_df)

  p_combined <- (p_predict | p_group) +
    plot_layout(widths = c(1.10, 1.00), guides = "keep") +
    plot_annotation(tag_levels = "a") &
    theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 16, colour = "black"),
      plot.tag.position = c(0.015, 0.985)
    )

  save_plot_safe(
    p_combined,
    "Combined_HLA_TSNE_Adjusted_Prediction_Group_projection",
    width = 9.8,
    height = 4.4
  )
  invisible(p_combined)
}

predict_df <- read_predict_number()
if (!is.null(predict_df)) {
  write_csv_safe(predict_df, file.path(result_dir, "Predict_sample_number_publication_source.csv"))
  plot_predict_sample_number(predict_df)
}

tsne_source_file <- file.path(getOption("hla.data_root"), "urine", "HLA_TSNE_Adjusted_Prediction_source.csv")
if (file.exists(tsne_source_file)) {
  tsne_df <- read.csv(tsne_source_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
    mutate(
      Group = factor(Group, levels = group_levels),
      predict = factor(predict, levels = unique(predict[order(sample_sort_key(predict))]))
    )
} else {
  message("Reconstructing HLA tSNE coordinates...")
  tsne_df <- build_hla_tsne()
  write_csv_safe(tsne_df, tsne_source_file)
}
plot_hla_tsne(tsne_df)
plot_hla_score_projection(tsne_df)
plot_group_projection(tsne_df)
plot_combined_hla_tsne(tsne_df)

message("Saved revised HLA tSNE and predicted sample number plots.")

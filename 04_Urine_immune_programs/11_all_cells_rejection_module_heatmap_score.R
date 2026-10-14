#!/usr/bin/env Rscript
# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final urine immune-cell composition, module and demultiplexing figures.
# Original source: 40_revision_all_cells_rejection_module_heatmap_score.R
# Figure/table scope: Fig4b; Fig4c; Fig4d; Fig4e; Fig4f; Fig4g; Fig4h; FigS7a; FigS7b; FigS7c; FigS7d; FigS7e
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# Rejection-related gene-module heatmap and sample-level module scores
# across all urine scRNA-seq cells and clinical groups.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "40_revision_all_cells_rejection_module_heatmap_score.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

object_file <- file.path(project_dir, "Result", "Objects", "Urine_ALL.rds")
figure_dir <- file.path(project_dir, "Figures", "scRNA", "Group_DEG_Function")
result_dir <- file.path(project_dir, "Result", "scRNA", "Group_DEG_Function")
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

rejection_modules <- list(
  "Cytotoxicity" = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH", "GZMK", "CST7", "CTSW", "IFNG", "FASLG", "CCL5"),
  "IFN-CXCR3" = c("CXCR3", "CXCL9", "CXCL10", "CXCL11", "IFNG", "STAT1", "IRF1", "GBP1", "GBP2", "GBP5", "ISG15", "IFI6", "IFIT3", "MX1"),
  "MHC I" = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"),
  "MHC II" = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "CIITA"),
  "Fc/complement" = c("FCGR1A", "FCGR2A", "FCGR3A", "C1QA", "C1QB", "C1QC", "C3", "C5AR1", "ITGAM", "TYROBP"),
  "Myeloid inflammation" = c("LYZ", "S100A8", "S100A9", "IL1B", "TNF", "NLRP3", "CXCL8", "CCL2", "CCL3", "CCL4"),
  "B/plasma" = c("MS4A1", "CD79A", "CD79B", "BANK1", "JCHAIN", "MZB1", "XBP1", "SDC1", "IGHG1", "IGKC"),
  "Checkpoint/Treg" = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX", "FOXP3", "IL2RA", "IKZF2", "TNFRSF18"),
  "Cell cycle" = c("MKI67", "TOP2A", "TYMS", "STMN1", "HMGB2", "UBE2C", "CENPF")
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

get_data_matrix <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "data")
  )
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "p=NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
  )
}

safe_kw <- function(values, groups) {
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- droplevels(factor(groups[keep]))
  if (nlevels(groups) < 2 || dplyr::n_distinct(values) < 2) return(NA_real_)
  suppressWarnings(kruskal.test(values ~ groups)$p.value)
}

safe_wilcox <- function(values, groups, g1, g2) {
  keep <- is.finite(values) & groups %in% c(g1, g2)
  values <- values[keep]
  groups <- droplevels(factor(groups[keep], levels = c(g1, g2)))
  if (nlevels(groups) < 2 || any(table(groups) < 2) || dplyr::n_distinct(values) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(values ~ groups, exact = FALSE)$p.value)
}

theme_pub <- function(base_size = 8.8, style = c("classic", "bw")) {
  style <- match.arg(style)
  base_theme <- if (style == "bw") {
    theme_bw(base_size = base_size, base_family = "Arial")
  } else {
    theme_classic(base_size = base_size, base_family = "Arial")
  }
  base_theme +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.36),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 6, 4, 6)
    )
}

save_plot_safe <- function(plot, filename, width, height) {
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(filename, ".png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
  file.copy(tmp_pdf, pdf_path, overwrite = TRUE)
  file.copy(tmp_png, png_path, overwrite = TRUE)
  invisible(c(pdf = pdf_path, png = png_path))
}

save_plot_safe_with_png_pixels <- function(plot, filename, pdf_width, pdf_height, png_height_px, png_width_px = NULL, png_res = 96) {
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(filename, ".png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  if (is.null(png_width_px)) {
    png_width_px <- round(png_height_px * pdf_width / pdf_height)
  }
  ggsave(tmp_pdf, plot, width = pdf_width, height = pdf_height, units = "in", device = cairo_pdf)
  render_res <- 96
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = png_width_px, height = png_height_px, units = "px", dpi = render_res, bg = "white", device = ragg::agg_png)
  } else {
    grDevices::png(tmp_png, width = png_width_px, height = png_height_px, units = "px", res = render_res, bg = "white")
    print(plot)
    grDevices::dev.off()
  }
  if (png_res != render_res && requireNamespace("magick", quietly = TRUE)) {
    img <- magick::image_read(tmp_png)
    magick::image_write(img, path = tmp_png, format = "png", density = paste0(png_res, "x", png_res))
  }
  file.copy(tmp_pdf, pdf_path, overwrite = TRUE)
  file.copy(tmp_png, png_path, overwrite = TRUE)
  invisible(c(pdf = pdf_path, png = png_path))
}

read_group_deg <- function(project_dir) {
  deg_dir <- file.path(project_dir, "Result", "Signature_score_analysis", "Differential_Expression")
  files <- c(
    "NR" = file.path(deg_dir, "scRNA_NR_DEG.csv"),
    "TCMR" = file.path(deg_dir, "scRNA_TCMR_DEG.csv"),
    "ABMR" = file.path(deg_dir, "scRNA_ABMR_DEG.csv"),
    "BK" = file.path(deg_dir, "scRNA_BK_DEG.csv")
  )
  bind_rows(lapply(names(files), function(grp) {
    if (!file.exists(files[[grp]])) return(tibble())
    raw <- read.csv(files[[grp]], stringsAsFactors = FALSE, check.names = FALSE)
    raw <- raw[, nzchar(names(raw)), drop = FALSE]
    if (!"gene" %in% colnames(raw)) return(tibble())
    tibble(
      gene = as.character(raw$gene),
      DEG_group = if ("cluster" %in% colnames(raw)) as.character(raw$cluster) else grp,
      avg_log2FC = suppressWarnings(as.numeric(raw$avg_log2FC)),
      p_val_adj = suppressWarnings(as.numeric(raw$p_val_adj))
    )
  })) %>%
    filter(!is.na(gene), gene != "") %>%
    group_by(gene) %>%
    slice_min(order_by = p_val_adj, n = 1, with_ties = FALSE) %>%
    ungroup()
}

if (!file.exists(object_file)) stop("Cannot find all-cell object: ", object_file)

message("Loading all-cell urine object...")
# Read the saved object or source table used by this analysis.
obj <- readRDS(object_file)
DefaultAssay(obj) <- "RNA"

meta <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    Group = as.character(Group),
    Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
    Major_celltype = as.character(Major_celltype),
    revised_celltype = as.character(revised_celltype)
  ) %>%
  filter(!is.na(Group), Group %in% group_levels, !is.na(Sample), Sample != "")

obj <- obj[, meta$cell]
expr <- get_data_matrix(obj, assay = "RNA")
meta <- meta %>% mutate(Group = factor(Group, levels = group_levels))

module_gene_df <- tibble(
  module = rep(names(rejection_modules), lengths(rejection_modules)),
  gene = unlist(rejection_modules, use.names = FALSE)
) %>%
  distinct(module, gene) %>%
  filter(gene %in% rownames(expr)) %>%
  mutate(module = factor(module, levels = names(rejection_modules)))

deg_info <- read_group_deg(project_dir)
module_gene_df <- module_gene_df %>%
  left_join(deg_info, by = "gene") %>%
  group_by(module) %>%
  arrange(is.na(p_val_adj), p_val_adj, desc(abs(avg_log2FC)), gene, .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(module_gene = paste(module, gene, sep = "__"))

# Export the analysis source data for figure reproduction.
write_csv_safe(module_gene_df, file.path(result_dir, "All_cells_rejection_module_gene_list.csv"))

heatmap_genes <- unique(module_gene_df$gene)
message("Calculating group-level expression heatmap...")
group_avg <- AverageExpression(
  obj,
  assays = "RNA",
  features = heatmap_genes,
  group.by = "Group",
  slot = "data",
  verbose = FALSE
)$RNA %>%
  as.matrix()

group_avg <- group_avg[, intersect(group_levels, colnames(group_avg)), drop = FALSE]
group_z <- t(scale(t(group_avg)))
group_z[is.na(group_z)] <- 0
group_z <- pmax(pmin(group_z, 2.2), -2.2)

heat_df <- as.data.frame(group_z) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Group", values_to = "z_score") %>%
  left_join(module_gene_df, by = "gene", relationship = "many-to-many") %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    module = factor(as.character(module), levels = names(rejection_modules)),
    module_gene = factor(module_gene, levels = rev(module_gene_df$module_gene))
  )
write_csv_safe(heat_df, file.path(result_dir, "All_cells_rejection_module_group_heatmap_source.csv"))

p_heat <- ggplot(heat_df, aes(Group, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_gradient2(
    low = "#6F99C3",
    mid = "#F7F7F7",
    high = "#B7192B",
    midpoint = 0,
    limits = c(-2.2, 2.2),
    oob = scales::squish,
    name = "Z-score"
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub(9.0) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.5, face = "bold"),
    axis.text.y = element_text(size = 7.6),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_plot_safe(p_heat, "All_cells_rejection_module_group_DEG_heatmap_revision", width = 5.2, height = 8.2)

message("Calculating sample-level module scores...")
module_cell_scores <- lapply(names(rejection_modules), function(module_name) {
  genes <- intersect(rejection_modules[[module_name]], rownames(expr))
  if (length(genes) < 2) return(NULL)
  tibble(
    cell = colnames(expr),
    module = module_name,
    n_genes = length(genes),
    module_score = as.numeric(Matrix::colMeans(expr[genes, , drop = FALSE]))
  )
}) %>%
  bind_rows() %>%
  left_join(meta %>% select(cell, Sample, Group, Major_celltype), by = "cell") %>%
  filter(!is.na(Group), !is.na(Sample))

sample_scores <- module_cell_scores %>%
  group_by(Group, Sample, module) %>%
  summarise(
    n_cells = n(),
    n_genes = max(n_genes),
    module_score = mean(module_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    module = factor(module, levels = names(rejection_modules))
  )

module_kw <- sample_scores %>%
  group_by(module) %>%
  summarise(
    p_value = safe_kw(module_score, Group),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"), label = format_p(p_value))

pairwise_grid <- as_tibble(t(utils::combn(group_levels, 2)), .name_repair = "minimal") %>%
  setNames(c("group1", "group2"))

module_pairwise <- sample_scores %>%
  group_by(module) %>%
  group_modify(function(dat, key) {
    pairwise_grid %>%
      mutate(
        p_value = mapply(function(g1, g2) safe_wilcox(dat$module_score, as.character(dat$Group), g1, g2), group1, group2)
      )
  }) %>%
  ungroup() %>%
  group_by(module) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"), label = format_p(p_value)) %>%
  ungroup()

write_csv_safe(module_cell_scores, file.path(result_dir, "All_cells_rejection_module_cell_scores.csv"))
write_csv_safe(sample_scores, file.path(result_dir, "All_cells_rejection_module_sample_scores.csv"))
write_csv_safe(module_kw, file.path(result_dir, "All_cells_rejection_module_KW_tests.csv"))
write_csv_safe(module_pairwise, file.path(result_dir, "All_cells_rejection_module_pairwise_wilcox_tests.csv"))

message("Projecting rejection module scores onto UMAP...")
umap_reduction <- if ("umap" %in% names(obj@reductions)) "umap" else names(obj@reductions)[1]
umap_df <- Embeddings(obj, reduction = umap_reduction) %>%
  as.data.frame() %>%
  rownames_to_column("cell")
colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")

module_umap_df <- module_cell_scores %>%
  select(cell, module, module_score) %>%
  left_join(umap_df, by = "cell") %>%
  filter(is.finite(UMAP_1), is.finite(UMAP_2), is.finite(module_score)) %>%
  mutate(
    module = factor(as.character(module), levels = names(rejection_modules)),
    module_label = factor(as.character(module), levels = names(rejection_modules)),
    score_plot = pmin(
      module_score,
      ave(module_score, module, FUN = function(x) stats::quantile(x, 0.995, na.rm = TRUE))
    )
  ) %>%
  arrange(module, score_plot)

write_csv_safe(module_umap_df, file.path(result_dir, "All_cells_rejection_module_score_UMAP_source.csv"))

umap_x_breaks <- seq(
  floor(min(module_umap_df$UMAP_1, na.rm = TRUE) / 5) * 5,
  ceiling(max(module_umap_df$UMAP_1, na.rm = TRUE) / 5) * 5,
  by = 5
)
umap_y_breaks <- seq(
  floor(min(module_umap_df$UMAP_2, na.rm = TRUE) / 5) * 5,
  ceiling(max(module_umap_df$UMAP_2, na.rm = TRUE) / 5) * 5,
  by = 5
)

p_module_umap <- ggplot(module_umap_df, aes(UMAP_1, UMAP_2, colour = score_plot)) +
  geom_point(size = 0.08, alpha = 0.72, stroke = 0) +
  facet_wrap(~module_label, ncol = 3) +
  scale_colour_gradientn(
    colours = c("#D6D6D6", "#F6D85F", "#B7192B"),
    name = "Module\nscore"
  ) +
  coord_equal() +
  scale_x_continuous(breaks = umap_x_breaks) +
  scale_y_continuous(breaks = umap_y_breaks) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_pub(10.2, style = "classic") +
  theme(
    strip.background = element_rect(fill = "grey94", colour = NA),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.ticks = element_line(linewidth = 0.28, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.36),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.35),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0.18, "lines"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.36, "in"),
    legend.key.width = unit(0.12, "in"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(colour = guide_colourbar(
    title.position = "top",
    barheight = unit(1.20, "in"),
    barwidth = unit(0.14, "in"),
    ticks.linewidth = 0.25
  ))

save_plot_safe_with_png_pixels(
  p_module_umap,
  "All_cells_rejection_module_score_UMAP_revision",
  pdf_width = 10.4,
  pdf_height = 10.0,
  png_height_px = 1000,
  png_width_px = 1040,
  png_res = 300
)

sample_scores_boxplot <- sample_scores %>%
  filter(as.character(Group) %in% boxplot_group_levels) %>%
  mutate(
    Group = factor(as.character(Group), levels = boxplot_group_levels),
    module = factor(as.character(module), levels = names(rejection_modules))
  )

boxplot_pairwise_grid <- as_tibble(t(utils::combn(boxplot_group_levels, 2)), .name_repair = "minimal") %>%
  setNames(c("group1", "group2"))

module_pairwise_boxplot <- sample_scores_boxplot %>%
  group_by(module) %>%
  group_modify(function(dat, key) {
    boxplot_pairwise_grid %>%
      mutate(
        p_value = mapply(function(g1, g2) safe_wilcox(dat$module_score, as.character(dat$Group), g1, g2), group1, group2)
      )
  }) %>%
  ungroup() %>%
  group_by(module) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"), label = format_p(p_value)) %>%
  ungroup()

write_csv_safe(sample_scores_boxplot, file.path(result_dir, "All_cells_rejection_module_sample_scores_no_Mix.csv"))
write_csv_safe(module_pairwise_boxplot, file.path(result_dir, "All_cells_rejection_module_pairwise_wilcox_tests_no_Mix.csv"))

pairwise_plot <- module_pairwise_boxplot %>%
  filter(!is.na(p_value), p_value < 0.1) %>%
  mutate(
    module = factor(as.character(module), levels = names(rejection_modules)),
    group1 = factor(group1, levels = boxplot_group_levels),
    group2 = factor(group2, levels = boxplot_group_levels),
    x = as.numeric(group1),
    xend = as.numeric(group2),
    label = format_p(p_value)
  ) %>%
  left_join(
    sample_scores_boxplot %>%
      group_by(module) %>%
      summarise(
        y_min = min(module_score, na.rm = TRUE),
        y_max = max(module_score, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        y_range = pmax(y_max - y_min, 0.05),
        step_y = y_range * 0.27
      ),
    by = "module"
  ) %>%
  group_by(module) %>%
  arrange(p_value, .by_group = TRUE) %>%
  mutate(
    bracket_id = row_number(),
    y = y_max + step_y * bracket_id,
    y_text = y + step_y * 0.30,
    tick = step_y * 0.16
  ) %>%
  ungroup()

write_csv_safe(pairwise_plot, file.path(result_dir, "All_cells_rejection_module_pairwise_wilcox_tests_no_Mix_p_lt_0.1.csv"))

p_box <- ggplot(sample_scores_boxplot, aes(Group, module_score, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, colour = "black", linewidth = 0.42, alpha = 0.34) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 11),
    shape = 21,
    stroke = 0.52,
    size = 2.10,
    alpha = 0.92
  ) +
  geom_segment(
    data = pairwise_plot,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = pairwise_plot,
    aes(x = x, xend = x, y = y, yend = y - tick),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = pairwise_plot,
    aes(x = xend, xend = xend, y = y, yend = y - tick),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_text(
    data = pairwise_plot,
    aes(x = (x + xend) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 3.0,
    family = "Arial",
    fontface = "bold",
    colour = "black"
  ) +
  facet_wrap(~module, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, drop = FALSE, guide = "none") +
  scale_colour_manual(values = group_colors, drop = FALSE, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.24))) +
  labs(x = NULL, y = "Module score per sample") +
  theme_pub(10.2, style = "bw") +
  theme(
    strip.background = element_rect(fill = "grey94", colour = NA),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9.2),
    axis.text.y = element_text(size = 9.8),
    axis.title.y = element_text(size = 11.8, face = "bold"),
    strip.text = element_text(size = 9.2, face = "bold"),
    panel.spacing.x = unit(0.75, "lines"),
    panel.spacing.y = unit(0.72, "lines"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

save_plot_safe(p_box, "All_cells_rejection_module_sample_score_boxplots_revision", width = 7.8, height = 7.4)

message("Done. Outputs saved to: ", figure_dir)

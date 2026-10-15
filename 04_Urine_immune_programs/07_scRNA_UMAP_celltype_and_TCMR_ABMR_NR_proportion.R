# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final urine immune-cell composition, module and demultiplexing figures.
# Original source: 26_revision_scRNA_UMAP_celltype_and_TCMR_ABMR_NR_proportion.R
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
  library(patchwork)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "26_revision_scRNA_UMAP_celltype_and_TCMR_ABMR_NR_proportion.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_umap_dir <- file.path(project_dir, "Figures", "scRNA", "UMAP_Markers")
figure_comp_dir <- file.path(project_dir, "Figures", "scRNA", "Celltype_Composition")
result_comp_dir <- file.path(project_dir, "Result", "scRNA", "Celltype_Composition")
dir.create(figure_umap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_comp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_comp_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR")
group_colors <- c("NR" = "#E42320", "TCMR" = "#699ECA", "ABMR" = "#4DAF4A")

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

celltype_colors <- c(
  "Naive CD4+ T cells" = "#E42320",
  "Memory CD4+ T cells" = "#699ECA",
  "Naive CD8+ T cells" = "#4DAF4A",
  "CD8+ NKT-like cells" = "#FF8C00",
  "ISG expressing immune cells" = "#F898CB",
  "Naive B cells" = "#D65190",
  "Classical Monocytes" = "#731A73",
  "Macrophages" = "#FFCB5B",
  "Myeloid Dendritic cells" = "#E89B1E",
  "Plasmacytoid Dendritic cells" = "#0076B9",
  "Unknown" = "#8FB4DC"
)

revised_celltype_levels <- c(
  "Naive CD4+ T", "Naive CD8+ T", "CD4+ Tem", "CD8+ Tem", "CD8+ NKT", "NK",
  "Cycling T", "B", "C1QC+ Macro", "TREM2+ Macro", "MT1A+ Macro",
  "CD16+ Mono", "THBS1+ Mono", "NLRP3+ Mono", "Neutrophil", "cDC", "pDC"
)

revised_celltype_colors <- c(
  "Naive CD4+ T" = "#E42320",
  "Naive CD8+ T" = "#4DAF4A",
  "CD4+ Tem" = "#699ECA",
  "CD8+ Tem" = "#FF8C00",
  "CD8+ NKT" = "#F898CB",
  "NK" = "#D65190",
  "Cycling T" = "#731A73",
  "B" = "#FFCB5B",
  "C1QC+ Macro" = "#E89B1E",
  "TREM2+ Macro" = "#0076B9",
  "MT1A+ Macro" = "#3D505A",
  "CD16+ Mono" = "#0098B2",
  "THBS1+ Mono" = "#8A7355",
  "NLRP3+ Mono" = "#CC5B45",
  "Neutrophil" = "#008B8B",
  "cDC" = "#4E5689",
  "pDC" = "#C9A1CA"
)

major_celltype_levels <- c("CD4+ T", "CD8+ T", "Cycling T", "NK", "B cell", "Monocyte", "Macrophage", "Neutrophil", "cDC", "pDC")
major_celltype_colors <- c(
  "CD4+ T" = "#E42320",
  "CD8+ T" = "#699ECA",
  "Cycling T" = "#4DAF4A",
  "NK" = "#FF8C00",
  "B cell" = "#F898CB",
  "Monocyte" = "#D65190",
  "Macrophage" = "#731A73",
  "Neutrophil" = "#FFCB5B",
  "cDC" = "#E89B1E",
  "pDC" = "#0076B9"
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.margin = margin(4, 6, 4, 6)
    )
}

save_plot <- function(plot, dir, filename, width, height) {
  pdf_path <- file.path(dir, paste0(filename, ".pdf"))
  png_path <- file.path(dir, paste0(filename, ".png"))
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

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
  )
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
if (!"umap" %in% names(scRNA@reductions)) {
  stop("The urine scRNA object has no UMAP reduction. Please use the integrated object with existing UMAP coordinates.")
}

meta <- scRNA@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
    Group = as.character(Group),
    Major_celltype = as.character(Major_celltype),
    Celltype = as.character(celltype),
    Revised_celltype = as.character(revised_celltype)
  ) %>%
  filter(
    !is.na(Sample), !Sample %in% c("Undefined", "Doublets", "Doublet", ""),
    !is.na(Group), Group %in% c("NR", "TCMR", "ABMR", "Mix", "BK"),
    !is.na(Major_celltype), Major_celltype %in% major_celltype_levels,
    !is.na(Celltype), Celltype %in% celltype_levels,
    !is.na(Revised_celltype), Revised_celltype %in% revised_celltype_levels
  ) %>%
  mutate(
    Group = factor(Group, levels = c("NR", "TCMR", "ABMR", "Mix", "BK")),
    Major_celltype = factor(Major_celltype, levels = major_celltype_levels),
    Celltype = factor(Celltype, levels = celltype_levels),
    Revised_celltype = factor(Revised_celltype, levels = revised_celltype_levels)
  )

umap_df <- as.data.frame(Embeddings(scRNA, "umap")) %>%
  rownames_to_column("cell") %>%
  left_join(meta, by = "cell") %>%
  filter(!is.na(Celltype))
colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")

set.seed(20260826)
umap_df <- umap_df %>% sample_frac(1)

p_umap <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Major_celltype)) +
  geom_point(size = 0.22, alpha = 0.86, stroke = 0) +
  scale_colour_manual(values = major_celltype_colors, drop = FALSE, name = "Major celltype") +
  guides(colour = guide_legend(override.aes = list(size = 3.2, alpha = 1), ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_pub(8.6) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7.0),
    legend.title = element_text(size = 8.0, face = "bold"),
    axis.title = element_text(size = 8.6, face = "bold"),
    axis.text = element_text(size = 7.2),
    plot.margin = margin(8, 10, 4, 4)
  )

save_plot(p_umap, figure_umap_dir, "scRNA_UMAP_celltype", width = 6.2, height = 4.6)

sample_order <- meta %>%
  filter(Group %in% group_levels) %>%
  distinct(Sample, Group) %>%
  arrange(factor(Group, levels = group_levels), mixedorder(Sample)) %>%
  pull(Sample)

prop_df <- meta %>%
  filter(Group %in% group_levels) %>%
  count(Sample, Group, Major_celltype, name = "n_cells") %>%
  complete(
    nesting(Sample, Group),
    Major_celltype = factor(major_celltype_levels, levels = major_celltype_levels),
    fill = list(n_cells = 0)
  ) %>%
  group_by(Sample, Group) %>%
  mutate(
    total_cells = sum(n_cells),
    proportion = ifelse(total_cells > 0, n_cells / total_cells * 100, 0)
  ) %>%
  ungroup() %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    Sample = factor(as.character(Sample), levels = sample_order),
    Major_celltype = factor(as.character(Major_celltype), levels = major_celltype_levels)
  )

stat_df <- prop_df %>%
  group_by(Major_celltype) %>%
  summarise(
    p_TCMR_vs_NR = {
      v1 <- proportion[Group == "NR"]
      v2 <- proportion[Group == "TCMR"]
      if (length(v1) >= 2 && length(v2) >= 2 && dplyr::n_distinct(c(v1, v2)) > 1) wilcox.test(v2, v1, exact = FALSE)$p.value else NA_real_
    },
    p_ABMR_vs_NR = {
      v1 <- proportion[Group == "NR"]
      v2 <- proportion[Group == "ABMR"]
      if (length(v1) >= 2 && length(v2) >= 2 && dplyr::n_distinct(c(v1, v2)) > 1) wilcox.test(v2, v1, exact = FALSE)$p.value else NA_real_
    },
    median_NR = median(proportion[Group == "NR"], na.rm = TRUE),
    median_TCMR = median(proportion[Group == "TCMR"], na.rm = TRUE),
    median_ABMR = median(proportion[Group == "ABMR"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_TCMR_vs_NR = p.adjust(p_TCMR_vs_NR, method = "BH"),
    p_adj_ABMR_vs_NR = p.adjust(p_ABMR_vs_NR, method = "BH")
  )

# Export the analysis source data for figure reproduction.
write_csv_safe(prop_df, file.path(result_comp_dir, "cell_type_proportion_TCMR_ABMR_vs_NR_sample_percent.csv"))
write_csv_safe(stat_df, file.path(result_comp_dir, "cell_type_proportion_TCMR_ABMR_vs_NR_wilcox_tests.csv"))

bar_df <- prop_df %>%
  group_by(Group, Major_celltype) %>%
  summarise(
    mean_percent = mean(proportion, na.rm = TRUE),
    sem_percent = sd(proportion, na.rm = TRUE) / sqrt(dplyr::n()),
    n_samples = dplyr::n_distinct(Sample),
    .groups = "drop"
  )

label_df <- prop_df %>%
  group_by(Major_celltype) %>%
  summarise(y_max = max(proportion, na.rm = TRUE), .groups = "drop") %>%
  left_join(stat_df, by = "Major_celltype") %>%
  mutate(
    y1 = y_max + pmax(y_max * 0.14, 0.35),
    y2 = y_max + pmax(y_max * 0.28, 0.70),
    label_TCMR = ifelse(!is.na(p_TCMR_vs_NR) & p_TCMR_vs_NR < 0.1, paste0("vs NR\n", format_p(p_TCMR_vs_NR)), NA_character_),
    label_ABMR = ifelse(!is.na(p_ABMR_vs_NR) & p_ABMR_vs_NR < 0.1, paste0("vs NR\n", format_p(p_ABMR_vs_NR)), NA_character_)
  )

make_compare_plot <- function(target_group, target_color) {
  compare_colors <- setNames(c(group_colors[["NR"]], target_color), c("NR", target_group))

  df <- prop_df %>%
    filter(Group %in% c("NR", target_group)) %>%
    mutate(Group = factor(as.character(Group), levels = c("NR", target_group)))

  order_tbl <- df %>%
    filter(Group == target_group) %>%
    group_by(Major_celltype) %>%
    summarise(target_mean = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(target_mean))

  mean_df <- df %>%
    group_by(Major_celltype, Group) %>%
    summarise(mean_percent = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Major_celltype = factor(as.character(Major_celltype), levels = rev(as.character(order_tbl$Major_celltype))),
      Group = factor(as.character(Group), levels = c("NR", target_group))
    )

  stat_plot_df <- bind_rows(lapply(as.character(order_tbl$Major_celltype), function(ct) {
    sub <- df %>% filter(as.character(Major_celltype) == ct)
    tibble(
      Major_celltype = ct,
      p_value = if (n_distinct(sub$Group) == 2 && dplyr::n_distinct(sub$proportion) > 1) {
        suppressWarnings(wilcox.test(proportion ~ Group, data = sub, exact = FALSE)$p.value)
      } else {
        NA_real_
      }
    )
  })) %>%
    mutate(
      Major_celltype = factor(Major_celltype, levels = rev(as.character(order_tbl$Major_celltype))),
      label = case_when(
        is.na(p_value) ~ "p=NA",
        p_value < 0.001 ~ "p<0.001",
        TRUE ~ paste0("p=", signif(p_value, 2))
      )
    ) %>%
    left_join(
      mean_df %>%
        group_by(Major_celltype) %>%
        summarise(x = max(mean_percent, na.rm = TRUE), .groups = "drop"),
      by = "Major_celltype"
    ) %>%
    mutate(x = x + max(df$proportion, na.rm = TRUE) * 0.05)

  point_df <- df %>%
    mutate(
      Major_celltype = factor(as.character(Major_celltype), levels = rev(as.character(order_tbl$Major_celltype))),
      Group = factor(as.character(Group), levels = c("NR", target_group))
    )

  ggplot() +
    geom_col(
      data = mean_df,
      aes(x = mean_percent, y = Major_celltype, fill = Group),
      position = position_dodge(width = 0.68),
      width = 0.58,
      alpha = 0.33,
      color = NA
    ) +
    geom_point(
      data = point_df,
      aes(x = proportion, y = Major_celltype, color = Group),
      position = position_jitterdodge(jitter.width = 0, jitter.height = 0.10, dodge.width = 0.68),
      size = 1.45,
      alpha = 0.90,
      stroke = 0
    ) +
    geom_point(
      data = mean_df,
      aes(x = mean_percent, y = Major_celltype, color = Group),
      position = position_dodge(width = 0.68),
      size = 2.2,
      alpha = 1,
      stroke = 0
    ) +
    geom_text(
      data = stat_plot_df,
      aes(x = x, y = Major_celltype, label = label),
      hjust = 0,
      size = 2.35,
      color = "black",
      family = "Arial"
    ) +
    scale_fill_manual(values = compare_colors, drop = FALSE) +
    scale_color_manual(values = compare_colors, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.24))) +
    labs(
      title = paste0(target_group, " vs NR"),
      x = "Immune-cell proportion (%)",
      y = NULL
    ) +
    guides(
      fill = guide_legend(title = NULL, override.aes = list(alpha = 0.45)),
      color = guide_legend(title = NULL, override.aes = list(size = 2.8, alpha = 1))
    ) +
    theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      axis.text = element_text(color = "black", size = 7),
      axis.title.x = element_text(color = "black", size = 8),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.margin = margin(4, 18, 4, 4)
    )
}

p_tcmr <- make_compare_plot("TCMR", group_colors[["TCMR"]])
p_abmr <- make_compare_plot("ABMR", group_colors[["ABMR"]])
p_prop <- p_tcmr + p_abmr + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

save_plot(p_prop, figure_comp_dir, "cell_type_proportion_TCMR_ABMR_vs_NR_scatter_bars", width = 10.0, height = 4.6)

message("Saved urine scRNA celltype UMAP and TCMR/ABMR vs NR proportion scatter-bar plots.")

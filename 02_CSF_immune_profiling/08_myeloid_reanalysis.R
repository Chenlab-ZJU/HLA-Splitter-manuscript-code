# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Required intermediate source-table / annotation generation.
# Original source: 09_revision_myeloid_reanalysis.R
# Figure/table scope: Fig2; FigS4; FigS5
# Data cohort: csf. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(gtools)
  library(patchwork)
  library(ggrepel)
})

set.seed(20260813)

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  file.path(getwd(), "Code", "09_revision_myeloid_reanalysis.R")
}
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
group_colors <- c(
  "Metastases" = "#E42320",
  "Hemorrhage" = "#699ECA",
  "Encephalitis" = "#4DAF4A"
)
myeloid_colors <- c(
  "Mono" = "#0076B9",
  "HMOX1-Macro" = "#E89B1E",
  "MT1G-Macro" = "#699ECA",
  "PPARG-Macro" = "#4DAF4A",
  "Inflam-Macro" = "#D65190",
  "microglia-like" = "#8A7355",
  "IFN-Myeloid" = "#731A73",
  "Other" = "grey55"
)
signature_to_label <- c(
  "Mono_CD14" = "Mono",
  "Macro_FSCN1" = "HMOX1-Macro",
  "Macro_HMOX1" = "HMOX1-Macro",
  "Macro_IL1B" = "Inflam-Macro",
  "Macro_MT1G" = "MT1G-Macro",
  "Macro_PPARG" = "PPARG-Macro",
  "Macro_P2RY12" = "microglia-like",
  "Microglia_like" = "microglia-like"
)
manual_cluster_labels <- c(
  "0" = "microglia-like",
  "1" = "HMOX1-Macro",
  "2" = "MT1G-Macro",
  "3" = "Mono",
  "4" = "PPARG-Macro",
  "5" = "Inflam-Macro",
  "6" = "T-cell contamination",
  "7" = "IFN-Myeloid"
)

theme_pub <- function(base_size = 9.5) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      legend.title = element_text(face = "bold"),
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
}

keep_present <- function(object, genes) {
  intersect(unique(genes), rownames(object))
}

score_signature <- function(object, signatures, prefix) {
  for (nm in names(signatures)) {
    genes <- keep_present(object, signatures[[nm]])
    if (length(genes) < 2) {
      object[[paste0(prefix, nm)]] <- NA_real_
    } else {
      object <- AddModuleScore(object, features = list(genes), name = paste0(prefix, nm, "_"))
      score_col <- grep(paste0("^", prefix, nm, "_"), colnames(object@meta.data), value = TRUE)
      object[[paste0(prefix, nm)]] <- object@meta.data[[score_col[length(score_col)]]]
    }
  }
  object
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", signif(p, 2))
  )
}

# Read the saved object or source table used by this analysis.
scRNA <- readRDS(file.path(getOption("hla.data_root"), "csf", "CSF_scRNA.rds"))
myeloid_source <- c("Macrophages", "Monocytes", "Microglia-like")
myeloid_source <- intersect(myeloid_source, unique(as.character(scRNA$Celltype)))
Myeloid <- subset(scRNA, subset = Celltype %in% myeloid_source)
DefaultAssay(Myeloid) <- "RNA"

Myeloid <- NormalizeData(Myeloid, verbose = FALSE)
Myeloid <- FindVariableFeatures(Myeloid, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
Myeloid <- ScaleData(Myeloid, vars.to.regress = intersect(c("nCount_RNA", "percent.mt"), colnames(Myeloid@meta.data)), verbose = FALSE)
Myeloid <- RunPCA(Myeloid, npcs = 40, verbose = FALSE)

reduction_use <- "pca"
if (requireNamespace("harmony", quietly = TRUE) && "HLA_Splitter" %in% colnames(Myeloid@meta.data)) {
  harmony_obj <- tryCatch(
    {
      harmony::RunHarmony(Myeloid, group.by.vars = "HLA_Splitter", dims.use = 1:30, verbose = FALSE)
    },
    error = function(e) {
      message("Harmony failed; falling back to PCA. Reason: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(harmony_obj) && "harmony" %in% Reductions(harmony_obj)) {
    Myeloid <- harmony_obj
    reduction_use <- "harmony"
  }
}

Myeloid <- FindNeighbors(Myeloid, reduction = reduction_use, dims = 1:30, verbose = FALSE)
Myeloid <- FindClusters(Myeloid, resolution = 0.55, verbose = FALSE)
Myeloid <- RunUMAP(Myeloid, reduction = reduction_use, dims = 1:30, min.dist = 0.45, n.neighbors = 30, verbose = FALSE)
Myeloid$Group <- factor(as.character(Myeloid$Group), levels = group_levels)

myeloid_signatures <- list(
  Mono_CD14 = c("S100A8", "S100A9", "S100A12", "FCN1", "VCAN", "LST1", "LGALS2", "CDA"),
  Macro_FSCN1 = c("FSCN1", "APOC2", "CCL4", "CCR7", "LAMP3", "CD40", "CCL19"),
  Macro_HMOX1 = c("HMOX1", "CD163", "CLEC10A", "FCN1", "LGALS3", "CTSL", "FOLR2"),
  Macro_IL1B = c("IL1B", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CCL20", "NFKBIA"),
  Macro_MT1G = c("MT1G", "MT1E", "MT1F", "MT1H", "SELENOP", "FABP5", "SDC3"),
  Macro_PPARG = c("PPARG", "CDK14", "PLPP3", "MRC1", "MSR1", "APOE", "APOC1"),
  Macro_P2RY12 = c("P2RY12", "SPP1", "TREM2", "C1QA", "C1QB", "C1QC", "CX3CR1"),
  Microglia_like = c("P2RY12", "TMEM119", "CX3CR1", "GPR34", "SALL1", "C1QA", "C1QB")
)

Myeloid <- score_signature(Myeloid, myeloid_signatures, "Sig_")
signature_cols <- paste0("Sig_", names(myeloid_signatures))
cluster_scores <- Myeloid@meta.data %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(signature_cols), ~ mean(.x, na.rm = TRUE)), n_cells = n(), .groups = "drop")

cluster_annotation <- cluster_scores %>%
  rowwise() %>%
  mutate(
    Myeloid_celltype_revision = {
      values <- c_across(all_of(signature_cols))
      if (all(is.na(values))) "Other" else signature_to_label[[sub("^Sig_", "", signature_cols[which.max(values)])]]
    },
    best_score = max(c_across(all_of(signature_cols)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(seurat_clusters, Myeloid_celltype_revision, best_score, n_cells, everything())

cluster_annotation$Myeloid_celltype_revision <- dplyr::coalesce(
  unname(manual_cluster_labels[cluster_annotation$seurat_clusters]),
  cluster_annotation$Myeloid_celltype_revision
)

Myeloid$Myeloid_celltype_revision <- cluster_annotation$Myeloid_celltype_revision[
  match(as.character(Myeloid$seurat_clusters), cluster_annotation$seurat_clusters)
]
Myeloid <- subset(Myeloid, subset = Myeloid_celltype_revision != "T-cell contamination")
Myeloid$Myeloid_celltype_revision <- factor(Myeloid$Myeloid_celltype_revision, levels = names(myeloid_colors))

Idents(Myeloid) <- "seurat_clusters"
# Discover marker genes with the original assay and comparison settings.
markers <- FindAllMarkers(Myeloid, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.20, verbose = FALSE)
top10 <- markers %>%
  filter(!is.na(p_val_adj)) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>%
  ungroup()

# Export the analysis source data for figure reproduction.
write_csv(markers, file.path(result_dir, "Myeloid_revision_cluster_markers.csv"))
write_csv(top10, file.path(result_dir, "Myeloid_revision_cluster_top10_markers.csv"))
write_csv(cluster_annotation, file.path(result_dir, "Myeloid_revision_cluster_annotation.csv"))

m1_genes <- c(
  "CXCL9", "CXCL10", "CXCL11", "IDO1", "CCR7", "LAMP3", "CCL5", "CD40", "IRF8",
  "CD38", "CCL4", "APOL3", "HCK", "BCL2A1", "RSAD2", "CYP27B1", "KYNU",
  "SIGLEC1", "SLAMF1", "SLC2A6", "IFI44L", "SAMSN1", "LILRB2", "TLR2",
  "CCL8", "APOBEC3A", "AQP9", "C3AR1", "HLA-DQA1", "MNDA", "IL2RA",
  "SLC15A3", "CD80", "ADGRE1", "CHI3L1", "RASSF4", "AIM2", "MMP9",
  "TLR8", "SOCS1", "LST1", "NOD2", "BIRC3", "GPR183", "AIF1", "CD86"
)
m2_genes <- c(
  "MMP9", "MS4A6A", "CCL18", "AIF1", "NCF2", "ACP5", "CLEC4A", "CCL13",
  "LY86", "SLC15A3", "CLEC10A", "CCL23", "HLA-DQA1", "RNASE6", "ADAMDEC1",
  "HCK", "NPL", "TREM2", "IRF8", "SIGLEC1", "CLEC7A", "P2RY13", "CD86",
  "CD180", "CD68", "CD209", "C3AR1", "GPR183", "FPR3", "CD37", "LST1",
  "HRH1", "EGR2", "CHI3L1", "CFP", "C5AR1", "MNDA", "PTGER2", "CD300A",
  "PIK3IP1", "IGSF6", "PLA2G7", "FCGR2B", "SLAMF8", "QPCT", "ATP8B4"
)

# Calculate the specified module score using the retained gene set.
Myeloid <- AddModuleScore(Myeloid, features = list(keep_present(Myeloid, m1_genes)), name = "M1_score_revision_")
Myeloid <- AddModuleScore(Myeloid, features = list(keep_present(Myeloid, m2_genes)), name = "M2_score_revision_")
Myeloid$M1_score_revision <- Myeloid$M1_score_revision_1
Myeloid$M2_score_revision <- Myeloid$M2_score_revision_1
Myeloid$Polarization_delta_revision <- Myeloid$M1_score_revision - Myeloid$M2_score_revision
Myeloid$Polarization_state_revision <- ifelse(Myeloid$Polarization_delta_revision >= 0, "M1-like", "M2-like")

umap_mat <- Embeddings(Myeloid, "umap")
umap_df <- as.data.frame(umap_mat[, 1:2, drop = FALSE]) %>%
  setNames(c("UMAP_1", "UMAP_2")) %>%
  transmute(
    cell = rownames(.),
    UMAP_1 = UMAP_1,
    UMAP_2 = UMAP_2,
    Cluster = factor(Myeloid$seurat_clusters, levels = mixedsort(unique(as.character(Myeloid$seurat_clusters)))),
    Myeloid_celltype = Myeloid$Myeloid_celltype_revision,
    Group = Myeloid$Group,
    M1_score = Myeloid$M1_score_revision,
    M2_score = Myeloid$M2_score_revision,
    Polarization_delta = Myeloid$Polarization_delta_revision,
    Polarization_state = Myeloid$Polarization_state_revision,
    Sample = as.character(Myeloid$HLA_Splitter)
  )

score_df <- umap_df %>%
  select(cell, Sample, Group, Cluster, Myeloid_celltype, M1_score, M2_score, Polarization_delta, Polarization_state)
write_csv(score_df, file.path(result_dir, "Myeloid_revision_polarization_cell_scores.csv"))

sample_score <- score_df %>%
  group_by(Sample, Group) %>%
  summarise(
    n_cells = n(),
    mean_M1 = mean(M1_score, na.rm = TRUE),
    mean_M2 = mean(M2_score, na.rm = TRUE),
    mean_delta = mean(Polarization_delta, na.rm = TRUE),
    M1_like_fraction = mean(Polarization_state == "M1-like", na.rm = TRUE) * 100,
    .groups = "drop"
  )
write_csv(sample_score, file.path(result_dir, "Myeloid_revision_polarization_sample_summary.csv"))

pairwise_wilcox <- function(df, value_col) {
  comparisons <- t(utils::combn(group_levels, 2))
  tibble(
    metric = value_col,
    group1 = comparisons[, 1],
    group2 = comparisons[, 2],
    p.value = apply(comparisons, 1, function(pair) {
      x <- df[[value_col]][df$Group == pair[1]]
      y <- df[[value_col]][df$Group == pair[2]]
      if (length(x) < 2 || length(y) < 2) return(NA_real_)
      stats::wilcox.test(x, y, exact = FALSE)$p.value
    })
  ) %>%
    mutate(label = format_p(p.value))
}
wilcox_stats <- bind_rows(
  pairwise_wilcox(sample_score, "mean_M1"),
  pairwise_wilcox(sample_score, "mean_M2"),
  pairwise_wilcox(sample_score, "mean_delta"),
  pairwise_wilcox(sample_score, "M1_like_fraction")
)
write_csv(wilcox_stats, file.path(result_dir, "Myeloid_revision_polarization_wilcox_tests.csv"))

cluster_centers <- umap_df %>%
  group_by(Cluster) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
celltype_centers <- umap_df %>%
  group_by(Myeloid_celltype) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), n = n(), .groups = "drop") %>%
  filter(!is.na(Myeloid_celltype), n >= 30)
cluster_palette <- setNames(
  rep(unname(myeloid_colors), length.out = length(levels(umap_df$Cluster))),
  levels(umap_df$Cluster)
)

p_cluster <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Cluster)) +
  geom_point(size = 0.18, alpha = 0.82, stroke = 0) +
  geom_label_repel(
    data = cluster_centers,
    aes(label = Cluster),
    size = 3.0,
    family = "Arial",
    fontface = "bold",
    label.size = 0.15,
    label.padding = unit(0.10, "lines"),
    fill = scales::alpha("white", 0.82),
    colour = "black",
    seed = 20260813,
    max.overlaps = Inf
  ) +
  scale_colour_manual(values = cluster_palette) +
  labs(title = "Myeloid clusters", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(legend.position = "none", panel.grid = element_blank())

p_anno <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Myeloid_celltype)) +
  geom_point(size = 0.18, alpha = 0.82, stroke = 0) +
  scale_colour_manual(values = myeloid_colors, drop = TRUE, name = "Myeloid state") +
  guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1))) +
  labs(title = "Myeloid states", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.44, "cm"),
    legend.key.width = unit(0.44, "cm"),
    panel.grid = element_blank()
  )

p_group <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Group)) +
  geom_point(size = 0.18, alpha = 0.78, stroke = 0) +
  scale_colour_manual(values = group_colors, drop = FALSE, name = "Group") +
  guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1))) +
  labs(title = "Clinical groups", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.44, "cm"),
    legend.key.width = unit(0.44, "cm"),
    panel.grid = element_blank()
  )

save_plot_safe(p_cluster + p_anno + p_group + plot_layout(widths = c(1, 1.25, 1.05)), "Myeloid_UMAP_Cluster", width = 11.2, height = 3.6)

p_myeloid_publication <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Myeloid_celltype)) +
  geom_point(size = 0.18, alpha = 0.82, stroke = 0) +
  scale_colour_manual(values = myeloid_colors, drop = TRUE, name = "Myeloid state") +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 8.8, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(face = "bold", size = 10.2),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.36),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9.2),
    legend.text = element_text(size = 8.3),
    legend.key.height = unit(0.32, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2.7, alpha = 1)))
save_plot_safe(p_myeloid_publication, "Myeloid_UMAP_Cluster_publication", width = 3.85, height = 3.55)

polar_df <- score_df %>%
  sample_n(size = min(n(), 7000), replace = FALSE)
p_polar <- ggplot(polar_df, aes(M1_score, M2_score, colour = Group)) +
  geom_point(size = 0.52, alpha = 0.56, stroke = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.45, colour = "grey35") +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(x = "M1 module score", y = "M2 module score") +
  theme_pub(10.2) +
  theme(panel.grid = element_blank(), legend.position = "right")
save_plot_safe(p_polar, "Myeloid_Polarization", width = 4.4, height = 3.8)

plot_metric_box <- function(metric, y_label, filename) {
  df <- sample_score %>%
    mutate(Group = factor(as.character(Group), levels = group_levels))
  ymax <- max(df[[metric]], na.rm = TRUE)
  ymin <- min(df[[metric]], na.rm = TRUE)
  yrange <- max(ymax - ymin, 0.1)
  stat_df <- pairwise_wilcox(df, metric) %>%
    mutate(
      x = match(group1, group_levels),
      xend = match(group2, group_levels),
      y.position = ymax + row_number() * yrange * 0.16,
      label.y = y.position + yrange * 0.045
    )
  p <- ggplot(df, aes(Group, .data[[metric]], fill = Group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.80, linewidth = 0.42) +
    geom_point(
      position = position_jitter(width = 0.11, height = 0, seed = 20260813),
      size = 1.55,
      shape = 21,
      colour = "black",
      stroke = 0.22,
      alpha = 0.90
    ) +
    geom_segment(data = stat_df, aes(x = x, xend = xend, y = y.position, yend = y.position), inherit.aes = FALSE, linewidth = 0.32) +
    geom_segment(data = stat_df, aes(x = x, xend = x, y = y.position, yend = y.position - yrange * 0.035), inherit.aes = FALSE, linewidth = 0.32) +
    geom_segment(data = stat_df, aes(x = xend, xend = xend, y = y.position, yend = y.position - yrange * 0.035), inherit.aes = FALSE, linewidth = 0.32) +
    geom_text(data = stat_df, aes(x = (x + xend) / 2, y = label.y, label = label), inherit.aes = FALSE, size = 3.0, family = "Arial", fontface = "bold") +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(x = NULL, y = y_label) +
    theme_pub(10.2) +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    )
  save_plot_safe(p, filename, width = 3.1, height = 3.3)
}

plot_metric_box("mean_M1", "Mean M1 module score", "Myeloid_Polarization_M1")
plot_metric_box("mean_M2", "Mean M2 module score", "Myeloid_Polarization_M2")

saveRDS(Myeloid, file.path(getOption("hla.data_root"), "csf", "CSF_Myeloid_revision.rds"))
message("Updated myeloid re-analysis using reduction: ", reduction_use)

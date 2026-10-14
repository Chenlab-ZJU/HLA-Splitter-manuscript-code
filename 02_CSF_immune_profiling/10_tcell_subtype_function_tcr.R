# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Required intermediate source-table / annotation generation.
# Original source: 11_revision_tcell_subtype_function_tcr.R
# Figure/table scope: Fig2; FigS4; FigS5
# Data cohort: 脑脊液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.
# Expanded T cells: clone_ratio > 1; distinguish cell counts from clonotype counts.
# Keep within-sample TCR denominators and clonotype identifiers as defined below.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(gtools)
  library(patchwork)
  library(ggrepel)
  library(harmony)
  library(cowplot)
})

set.seed(20260813)

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "11_revision_tcell_subtype_function_tcr.R")
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
group_colors <- c("Metastases" = "#E42320", "Hemorrhage" = "#699ECA", "Encephalitis" = "#4DAF4A")
tcell_colors <- c(
  "Naive T" = "#E42320", "CD4+ Tcm" = "#699ECA", "Th17" = "#4DAF4A",
  "IFN+ T" = "#FF8C00", "CD4+ Tex" = "#F898CB", "Tregs" = "#731A73",
  "CD8+ Tmem" = "#FFCB5B", "CD8+ Trm" = "#E89B1E", "CD8+ Tact" = "#0076B9",
  "CD8+ Teff" = "#3D505A", "NK-like T" = "#0098B2", "gd T" = "#8A7355"
)
clone_colors <- c(
  "No TCR" = "grey82", "Un-expand" = "#D9DEE7", "Small (expand:1-2)" = "#9FDAF7",
  "Medium (expand:2-5)" = "#FFCB5B", "Large (expand:5-10)" = "#EB7E60",
  "Hyper (expand:>10)" = "#E42320"
)
tcell_signature_genes <- list(
  "Naive T" = c("CCR7", "LEF1", "TCF7", "IL7R", "SELL", "MAL"),
  "CD4+ Tcm" = c("IL7R", "CCR7", "LTB", "MAL", "CD40LG", "ANXA1"),
  "Th17" = c("CCR6", "AQP3", "KLRB1", "RORC", "IL7R", "S100A4"),
  "IFN+ T" = c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "RSAD2"),
  "CD4+ Tex" = c("PDCD1", "CTLA4", "ICOS", "MAF", "TOX", "TIGIT"),
  "Tregs" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT"),
  "CD8+ Tmem" = c("CD8A", "CD8B", "GZMK", "IL7R", "CD27", "CISH"),
  "CD8+ Trm" = c("ZNF683", "ITGAE", "CXCR6", "BHLHE40", "FGFBP2"),
  "CD8+ Tact" = c("CD8A", "GZMK", "HLA-DRA", "HLA-DQA1", "CCL5"),
  "CD8+ Teff" = c("GZMB", "GZMH", "PRF1", "NKG7", "GNLY", "KLRD1"),
  "NK-like T" = c("KLRD1", "NKG7", "GNLY", "XCL1", "KIR2DL4", "FCGR3A"),
  "gd T" = c("TRDC", "TRGC1", "TRGC2", "TRDV2", "TRGV9")
)

theme_pub <- function(base_size = 9.2) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      panel.grid.minor = element_blank(),
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

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", signif(p, 2))
  )
}

pairwise_wilcox <- function(df, value_col, facet_col = NULL) {
  if (is.null(facet_col)) df[[".facet"]] <- "all" else df[[".facet"]] <- as.character(df[[facet_col]])
  df %>%
    group_by(.facet) %>%
    group_modify(function(.x, .y) {
      comparisons <- t(utils::combn(group_levels, 2))
      y_max <- max(.x[[value_col]], na.rm = TRUE)
      y_min <- min(.x[[value_col]], na.rm = TRUE)
      y_range <- max(y_max - y_min, 0.08)
      tibble(
        group1 = comparisons[, 1],
        group2 = comparisons[, 2],
        x = match(comparisons[, 1], group_levels),
        xend = match(comparisons[, 2], group_levels),
        y.position = y_max + seq_len(nrow(comparisons)) * y_range * 0.18,
        tick.height = y_range * 0.035,
        label.y = y_max + seq_len(nrow(comparisons)) * y_range * 0.18 + y_range * 0.06,
        p.value = apply(comparisons, 1, function(pair) {
          a <- .x[[value_col]][.x$Group == pair[1]]
          b <- .x[[value_col]][.x$Group == pair[2]]
          if (length(a) < 2 || length(b) < 2) return(NA_real_)
          stats::wilcox.test(a, b, exact = FALSE)$p.value
        })
      ) %>% mutate(label = format_p(p.value))
    }) %>%
    ungroup()
}

score_signature <- function(object, signatures, prefix) {
  for (nm in names(signatures)) {
    genes <- intersect(signatures[[nm]], rownames(object))
    if (length(genes) < 2) {
      object[[paste0(prefix, nm)]] <- NA_real_
    } else {
      object <- AddModuleScore(object, features = list(genes), name = paste0(prefix, make.names(nm), "_"))
      score_col <- grep(paste0("^", prefix, make.names(nm), "_"), colnames(object@meta.data), value = TRUE)
      object[[paste0(prefix, nm)]] <- object@meta.data[[score_col[length(score_col)]]]
    }
  }
  object
}

# Read the saved object or source table used by this analysis.
previous_tcell <- readRDS(file.path(result_dir, "CSF_Tcell.rds"))
previous_annotation <- previous_tcell@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(cell, previous_T_celltype = as.character(T_celltype))

scRNA <- readRDS(file.path(result_dir, "CSF_scRNA.rds"))
Tcell <- subset(scRNA, subset = Celltype %in% c("Naive T", "CD4+ T", "CD8+ T", "IFN+ T", "Tregs", "Proliferating", "gd T"))
Tcell@meta.data <- Tcell@meta.data %>%
  as_tibble(rownames = "cell") %>%
  left_join(previous_annotation, by = "cell") %>%
  tibble::column_to_rownames("cell") %>%
  as.data.frame()
DefaultAssay(Tcell) <- "RNA"
Tcell <- NormalizeData(Tcell, verbose = FALSE)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
Tcell <- ScaleData(Tcell, vars.to.regress = intersect(c("nCount_RNA", "percent.mt"), colnames(Tcell@meta.data)), verbose = FALSE)
Tcell <- RunPCA(Tcell, npcs = 40, verbose = FALSE)
reduction_use <- "pca"
if ("orig.ident" %in% colnames(Tcell@meta.data)) {
  harmony_obj <- tryCatch(
    harmony::RunHarmony(Tcell, group.by.vars = "orig.ident", dims.use = 1:30, verbose = FALSE),
    error = function(e) {
      message("Harmony failed; falling back to PCA. Reason: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(harmony_obj) && "harmony" %in% Reductions(harmony_obj)) {
    Tcell <- harmony_obj
    reduction_use <- "harmony"
  }
}
Tcell <- FindNeighbors(Tcell, reduction = reduction_use, dims = 1:30, verbose = FALSE)
Tcell <- FindClusters(Tcell, resolution = 0.8, verbose = FALSE)
Tcell <- RunUMAP(Tcell, reduction = reduction_use, dims = 1:30, min.dist = 0.45, n.neighbors = 30, verbose = FALSE)

Tcell$T_celltype <- as.character(Tcell$previous_T_celltype)
Tcell$T_celltype[is.na(Tcell$T_celltype) | Tcell$T_celltype %in% c("T_cycling", "Contamination", "Doublet")] <- NA_character_
Tcell$Group <- factor(as.character(Tcell$Group), levels = group_levels)
Tcell$T_celltype <- factor(as.character(Tcell$T_celltype), levels = names(tcell_colors))
Tcell <- subset(Tcell, subset = !is.na(T_celltype))

cluster_annotation <- Tcell@meta.data %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  count(seurat_clusters, T_celltype, name = "n") %>%
  group_by(seurat_clusters) %>%
  mutate(cluster_cells = sum(n), fraction = n / cluster_cells) %>%
  arrange(seurat_clusters, desc(n), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()
# Export the analysis source data for figure reproduction.
write_csv(cluster_annotation, file.path(result_dir, "Tcell_revision_harmony_cluster_annotation.csv"))
Idents(Tcell) <- "seurat_clusters"
# Discover marker genes with the original assay and comparison settings.
tcell_markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.20, verbose = FALSE)
tcell_top10 <- tcell_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>%
  ungroup()
write_csv(tcell_markers, file.path(result_dir, "Tcell_revision_harmony_cluster_markers.csv"))
write_csv(tcell_top10, file.path(result_dir, "Tcell_revision_harmony_cluster_top10_markers.csv"))

tcr_path <- file.path(result_dir, "Table", "TCR_information_all.csv")
if (file.exists(tcr_path)) {
  tcr_raw <- read_csv(tcr_path, show_col_types = FALSE, name_repair = "unique_quiet")
  barcode_cols <- grep("^barcode", colnames(tcr_raw), value = TRUE)
  tcr_cell <- tcr_raw %>%
    transmute(
      cell = as.character(.data[[barcode_cols[1]]]),
      Sample = as.character(sample),
      clone_id = as.character(raw_clonotype_id),
      CTaa = if ("cdr3" %in% colnames(tcr_raw)) as.character(cdr3) else NA_character_,
      chain = as.character(chain),
      v_gene = as.character(v_gene),
      clonalFrequency = as.numeric(clonalFrequency),
      clone_ratio = as.numeric(clone_ratio),
      cloneSize = as.character(cloneSize)
    ) %>%
    group_by(cell, Sample) %>%
    summarise(
      clone_id = paste(unique(na.omit(clone_id)), collapse = ";"),
      CTaa = paste(unique(na.omit(CTaa)), collapse = ";"),
      TCR_chain = paste(sort(unique(na.omit(chain))), collapse = "+"),
      TRBV = paste(sort(unique(v_gene[grepl("^TRBV", v_gene)])), collapse = ";"),
      clonalFrequency = max(clonalFrequency, na.rm = TRUE),
      clone_ratio = max(clone_ratio, na.rm = TRUE),
      cloneSize = cloneSize[which.max(replace_na(clone_ratio, -Inf))],
      .groups = "drop"
    ) %>%
    mutate(
      clonalFrequency = ifelse(is.infinite(clonalFrequency), 0, clonalFrequency),
      clone_ratio = ifelse(is.infinite(clone_ratio), 0, clone_ratio),
      cloneSize = ifelse(is.na(cloneSize) | cloneSize == "", "No TCR", cloneSize)
    )
  meta_tcr <- Tcell@meta.data %>%
    as_tibble(rownames = "cell") %>%
    transmute(cell, Sample = as.character(HLA_Splitter)) %>%
    left_join(tcr_cell, by = c("cell", "Sample"))
  Tcell$clone_ratio_revision <- replace_na(meta_tcr$clone_ratio, 0)
  Tcell$cloneSize_revision <- replace_na(meta_tcr$cloneSize, "No TCR")
  Tcell$clonalFrequency_revision <- replace_na(meta_tcr$clonalFrequency, 0)
  Tcell$TRBV_revision <- replace_na(meta_tcr$TRBV, "")
  Tcell$TCR_chain_revision <- replace_na(meta_tcr$TCR_chain, "None")
}
Tcell$cloneSize_revision <- factor(as.character(Tcell$cloneSize_revision), levels = names(clone_colors))

umap_mat <- Embeddings(Tcell, "umap")[, 1:2, drop = FALSE]
umap_df <- as.data.frame(umap_mat) %>%
  setNames(c("UMAP_1", "UMAP_2")) %>%
  mutate(
    cell = rownames(.),
    Cluster = factor(Tcell$seurat_clusters, levels = mixedsort(unique(as.character(Tcell$seurat_clusters)))),
    Group = Tcell$Group,
    T_celltype = Tcell$T_celltype,
    clone_ratio = Tcell$clone_ratio_revision,
    cloneSize = Tcell$cloneSize_revision
  )
write_csv(umap_df, file.path(result_dir, "Tcell_revision_umap_metadata.csv"))

umap_xlim <- range(umap_df$UMAP_1, finite = TRUE)
umap_ylim <- range(umap_df$UMAP_2, finite = TRUE)
coord_tcell_umap <- coord_fixed(xlim = umap_xlim, ylim = umap_ylim, expand = FALSE)
cluster_centers <- umap_df %>%
  group_by(Cluster) %>%
  summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE), UMAP_2 = median(UMAP_2, na.rm = TRUE), .groups = "drop")

p_tcell_cluster <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Cluster)) +
  geom_point(size = 0.18, alpha = 0.82, stroke = 0) +
  geom_label_repel(
    data = cluster_centers,
    aes(x = UMAP_1, y = UMAP_2, label = Cluster),
    inherit.aes = FALSE,
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
  scale_colour_manual(
    values = setNames(rep(unname(tcell_colors), length.out = length(levels(umap_df$Cluster))), levels(umap_df$Cluster))
  ) +
  coord_tcell_umap +
  labs(title = "T cell clusters", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

p_tcell_type <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = T_celltype)) +
  geom_point(size = 0.18, alpha = 0.82, stroke = 0) +
  scale_colour_manual(values = tcell_colors, drop = TRUE, name = "T cell subset") +
  guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1))) +
  coord_tcell_umap +
  labs(title = "T cell subsets", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.44, "cm"),
    legend.key.width = unit(0.44, "cm"),
    panel.grid = element_blank()
  )
p_tcell_group <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Group)) +
  geom_point(size = 0.18, alpha = 0.78, stroke = 0) +
  scale_colour_manual(values = group_colors, drop = FALSE, name = "Group") +
  guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1))) +
  coord_tcell_umap +
  labs(title = "Clinical groups", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.44, "cm"),
    legend.key.width = unit(0.44, "cm"),
    panel.grid = element_blank()
  )
clone_umap <- umap_df %>% arrange(clone_ratio)
p_tcr_umap <- ggplot(clone_umap, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = clone_ratio, size = pmin(clone_ratio, 12)), alpha = 0.84, stroke = 0) +
  scale_colour_gradientn(colours = c("grey82", "#F9B3AD", "#EB7E60", "#E42320"), name = "Expansion\nratio") +
  scale_size_continuous(range = c(0.10, 0.85), guide = "none") +
  coord_tcell_umap +
  labs(title = "TCR expansion", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.2) +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
tcell_subset_legend <- cowplot::get_legend(p_tcell_type)
tcell_group_legend <- cowplot::get_legend(p_tcell_group)
tcell_umap_panels <- (
  p_tcell_cluster +
    p_tcell_type + theme(legend.position = "none") +
    p_tcell_group + theme(legend.position = "none") +
    plot_layout(nrow = 1, widths = c(1, 1, 1))
)
tcell_umap_legends <- cowplot::plot_grid(
  tcell_subset_legend,
  tcell_group_legend,
  ncol = 1,
  rel_heights = c(1.7, 0.8),
  align = "v"
)
tcell_umap_combined <- cowplot::plot_grid(
  tcell_umap_panels,
  tcell_umap_legends,
  nrow = 1,
  rel_widths = c(3.0, 0.85),
  align = "h",
  axis = "tb"
)
save_plot_safe(tcell_umap_combined, "CSF_UMAP_Tcell_celltype_revision", width = 11.2, height = 3.6)
save_plot_safe(p_tcr_umap, "Tcell_TCR_expansion_UMAP_revision", width = 4.0, height = 3.6)

tcell_dotplot_markers <- list(
  "Naive T" = c("CCR7", "LEF1", "IL7R"),
  "CD4+ Tcm" = c("MAL", "LTB", "CD40LG"),
  "Th17" = c("CCR6", "AQP3", "KLRB1"),
  "IFN+ T" = c("ISG15", "IFIT3", "MX1"),
  "CD4+ Tex" = c("PDCD1", "CTLA4", "ICOS"),
  "Tregs" = c("FOXP3", "IL2RA", "TIGIT"),
  "CD8+ Tmem" = c("GZMK", "CD27", "CISH"),
  "CD8+ Trm" = c("ZNF683", "FGFBP2", "BHLHE40"),
  "CD8+ Tact" = c("CD8A", "HLA-DRA", "HLA-DQA1"),
  "CD8+ Teff" = c("GZMB", "PRF1", "NKG7"),
  "NK-like T" = c("KLRD1", "GNLY", "XCL1"),
  "gd T" = c("TRDC", "TRGC1", "TRGC2")
)
dotplot_marker_df <- tibble(
  T_celltype = rep(names(tcell_dotplot_markers), lengths(tcell_dotplot_markers)),
  gene = unlist(tcell_dotplot_markers, use.names = FALSE)
) %>%
  filter(T_celltype %in% levels(droplevels(Tcell$T_celltype)), gene %in% rownames(Tcell)) %>%
  distinct(T_celltype, gene)
write_csv(dotplot_marker_df, file.path(result_dir, "Tcell_revision_marker_gene_dotplot_markers.csv"))

dotplot_features <- unique(dotplot_marker_df$gene)
p_tcell_marker_dot <- DotPlot(
  Tcell,
  features = dotplot_features,
  group.by = "T_celltype",
  assay = "RNA",
  dot.scale = 5.8
) +
  scale_colour_gradientn(colours = c("#6F99C3", "#F7F7F7", "#B7192B"), name = "Average\nexpression") +
  scale_size(range = c(0.6, 4.8), name = "Percent\nexpressed") +
  labs(x = NULL, y = NULL) +
  theme_pub(9.8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.4),
    axis.text.y = element_text(size = 8.8),
    legend.title = element_text(size = 9.2, face = "bold"),
    legend.text = element_text(size = 8.4),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.position = "right"
  )
save_plot_safe(p_tcell_marker_dot, "Tcell_marker_gene_dotplot_revision", width = 8.4, height = 4.7)

meta <- Tcell@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Sample = as.character(HLA_Splitter),
    Group = factor(as.character(Group), levels = group_levels),
    T_celltype = factor(as.character(T_celltype), levels = names(tcell_colors)),
    clone_ratio = as.numeric(clone_ratio_revision),
    cloneSize = factor(as.character(cloneSize_revision), levels = names(clone_colors))
  )

prop_df <- meta %>%
  count(Sample, Group, T_celltype, name = "Count") %>%
  group_by(Sample, Group) %>%
  mutate(Total = sum(Count), Percent = Count / Total * 100) %>%
  ungroup() %>%
  tidyr::complete(
    tidyr::nesting(Sample, Group),
    T_celltype,
    fill = list(Count = 0, Percent = 0)
  ) %>%
  filter(!is.na(T_celltype))
write_csv(prop_df, file.path(result_dir, "Tcell_revision_subtype_sample_percent.csv"))

sample_order <- meta %>%
  distinct(Sample, Group) %>%
  arrange(Group, mixedorder(Sample)) %>%
  pull(Sample)

tcell_prop_heat_df <- prop_df %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    Group = factor(as.character(Group), levels = group_levels),
    T_celltype = factor(as.character(T_celltype), levels = rev(names(tcell_colors)))
  )
write_csv(tcell_prop_heat_df, file.path(result_dir, "Tcell_revision_subtype_proportion_heatmap_source.csv"))

p_tcell_prop_heat <- ggplot(tcell_prop_heat_df, aes(Sample, T_celltype, fill = Percent)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#D7E4EF", "#F7F7F7", "#E79A94", "#B7192B"),
    name = "Cells (%)",
    limits = c(0, max(tcell_prop_heat_df$Percent, na.rm = TRUE))
  ) +
  labs(x = "Samples grouped by diagnosis", y = "T cell subset") +
  theme_pub(9.6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.2),
    axis.text.y = element_text(size = 8.6),
    strip.text = element_text(size = 9.2, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )
save_plot_safe(p_tcell_prop_heat, "Tcell_celltype_proportion_heatmap_revision", width = 6.2, height = 4.2)

prop_kw <- prop_df %>%
  group_by(T_celltype) %>%
  summarise(p.value = kruskal.test(Percent ~ Group)$p.value, .groups = "drop") %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH"), label = format_p(p.value))
write_csv(prop_kw, file.path(result_dir, "Tcell_revision_subtype_KW_tests.csv"))

p_prop <- ggplot(prop_df, aes(Group, Percent, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35, alpha = 0.82) +
  geom_point(position = position_jitter(width = 0.11, height = 0, seed = 20260813), shape = 21, colour = "black", stroke = 0.18, size = 1.15, alpha = 0.86) +
  facet_wrap(~T_celltype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(x = NULL, y = "Cells per sample (%)") +
  theme_pub(8.7) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.9), strip.text = element_text(size = 7.5), legend.position = "top", panel.grid.major.x = element_blank())
save_plot_safe(p_prop, "Tcell_subtype_proportion_compare_revision", width = 7.2, height = 5.4)

tcell_modules <- list(
  "Cytotoxic" = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX"),
  "IFN response" = c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "RSAD2", "OAS1"),
  "Activation" = c("IL2RA", "CD69", "ICOS", "TNFRSF9", "HLA-DRA", "HLA-DQA1"),
  "Treg" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT"),
  "Memory/naive" = c("CCR7", "IL7R", "LEF1", "TCF7", "SELL", "MAL")
)
module_score_cell <- lapply(names(tcell_modules), function(module_name) {
  genes <- intersect(tcell_modules[[module_name]], rownames(Tcell))
  if (length(genes) < 2) return(NULL)
  mat <- FetchData(Tcell, vars = genes)
  tibble(cell = rownames(mat), module = module_name, module_score = rowMeans(as.matrix(mat), na.rm = TRUE))
}) %>% bind_rows() %>% left_join(meta %>% dplyr::select(cell, Sample, Group), by = "cell")
module_score_sample <- module_score_cell %>%
  group_by(Sample, Group, module) %>%
  summarise(n_cells = n(), module_score = mean(module_score, na.rm = TRUE), .groups = "drop") %>%
  mutate(module = factor(module, levels = names(tcell_modules)))
write_csv(module_score_sample, file.path(result_dir, "Tcell_revision_functional_module_sample_scores.csv"))
module_kw <- module_score_sample %>%
  group_by(module) %>% summarise(p.value = kruskal.test(module_score ~ Group)$p.value, .groups = "drop") %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH"), label = format_p(p.value))
module_wilcox <- pairwise_wilcox(module_score_sample, "module_score", "module")
write_csv(module_kw, file.path(result_dir, "Tcell_revision_functional_module_KW_tests.csv"))
write_csv(module_wilcox, file.path(result_dir, "Tcell_revision_functional_module_wilcox_tests.csv"))
p_module <- ggplot(module_score_sample, aes(Group, module_score, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.36, alpha = 0.82) +
  geom_point(position = position_jitter(width = 0.11, height = 0, seed = 20260813), shape = 21, colour = "black", stroke = 0.18, size = 1.25, alpha = 0.88) +
  facet_wrap(~module, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(x = NULL, y = "Module score per sample") +
  theme_pub(8.8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7.0), strip.text = element_text(size = 7.8), legend.position = "top", panel.grid.major.x = element_blank())
save_plot_safe(p_module, "Tcell_functional_module_score_compare_revision", width = 6.1, height = 4.3)

gene_expr <- FetchData(Tcell, vars = intersect(c("GZMB", "PRF1", "NKG7", "GNLY", "GZMK", "PDCD1", "LAG3", "CXCR3", "IFNG", "ISG15"), rownames(Tcell))) %>%
  as_tibble(rownames = "cell") %>%
  left_join(meta %>% dplyr::select(cell, Sample, Group), by = "cell") %>%
  pivot_longer(cols = -c(cell, Sample, Group), names_to = "gene", values_to = "expression") %>%
  group_by(Sample, Group, gene) %>%
  summarise(mean_expression = mean(expression, na.rm = TRUE), pct_expressing = mean(expression > 0) * 100, .groups = "drop")
write_csv(gene_expr, file.path(result_dir, "Tcell_revision_function_gene_sample_expression.csv"))
p_gene <- ggplot(gene_expr, aes(Group, mean_expression, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35, alpha = 0.82) +
  geom_point(position = position_jitter(width = 0.11, height = 0, seed = 20260813), shape = 21, colour = "black", stroke = 0.18, size = 1.05, alpha = 0.86) +
  facet_wrap(~gene, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(x = NULL, y = "Mean expression per sample") +
  theme_pub(8.5) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.6), strip.text = element_text(size = 7.4), legend.position = "top", panel.grid.major.x = element_blank())
save_plot_safe(p_gene, "Tcell_function_gene_expression_compare_revision", width = 7.6, height = 4.2)

tcell_marker_modules <- list(
  "Naive/memory" = c("CCR7", "IL7R", "LEF1", "TCF7", "SELL", "MAL", "LTB"),
  "CD4 activation/Treg" = c("CD40LG", "CCR6", "AQP3", "FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "Cytotoxic" = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK"),
  "Exhaustion" = c("PDCD1", "LAG3", "HAVCR2", "TOX", "ICOS", "MAF"),
  "IFN/chemokine" = c("IFNG", "CXCR3", "CCL3", "CCL4", "ISG15", "IFI6", "IFIT3", "MX1"),
  "Antigen presentation" = c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DPA1", "CD74")
)
tcell_marker_df <- tibble(
  module = rep(names(tcell_marker_modules), lengths(tcell_marker_modules)),
  gene = unlist(tcell_marker_modules, use.names = FALSE)
) %>%
  filter(gene %in% rownames(Tcell)) %>%
  distinct(module, gene) %>%
  group_by(module) %>%
  mutate(module_gene = paste(module, gene, sep = "__")) %>%
  ungroup()
write_csv(tcell_marker_df, file.path(result_dir, "Tcell_revision_heatmap_gene_modules.csv"))

heatmap_genes <- unique(tcell_marker_df$gene)
celltype_avg <- AverageExpression(Tcell, assays = "RNA", features = heatmap_genes, group.by = "T_celltype", slot = "data")$RNA %>%
  as.matrix()
celltype_z <- t(scale(t(celltype_avg)))
celltype_z[is.na(celltype_z)] <- 0
celltype_heat_df <- as.data.frame(celltype_z) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "T_celltype", values_to = "z_score") %>%
  left_join(tcell_marker_df, by = "gene", relationship = "many-to-many") %>%
  mutate(
    module = factor(module, levels = names(tcell_marker_modules)),
    module_gene = factor(module_gene, levels = rev(tcell_marker_df$module_gene)),
    T_celltype = factor(T_celltype, levels = names(tcell_colors))
  )
write_csv(celltype_heat_df, file.path(result_dir, "Tcell_revision_celltype_heatmap_source.csv"))

p_celltype_heat <- ggplot(celltype_heat_df, aes(T_celltype, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_gradient2(low = "#6F99C3", mid = "#F7F7F7", high = "#B7192B", midpoint = 0, name = "Z-score") +
  labs(x = NULL, y = NULL) +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.3),
    axis.text.y = element_text(size = 7.8),
    strip.text.y = element_text(size = 7.8, angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey93", colour = NA),
    panel.grid = element_blank(),
    legend.position = "right"
  )
save_plot_safe(p_celltype_heat, "Tcell_celltype_heatmap_revision", width = 6.5, height = 7.0)

group_avg <- AverageExpression(Tcell, assays = "RNA", features = heatmap_genes, group.by = "Group", slot = "data")$RNA %>%
  as.matrix()
group_z <- t(scale(t(group_avg)))
group_z[is.na(group_z)] <- 0
group_heat_df <- as.data.frame(group_z) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = all_of(group_levels), names_to = "Group", values_to = "z_score") %>%
  left_join(tcell_marker_df, by = "gene", relationship = "many-to-many") %>%
  mutate(
    module = factor(module, levels = names(tcell_marker_modules)),
    module_gene = factor(module_gene, levels = rev(tcell_marker_df$module_gene)),
    Group = factor(Group, levels = group_levels)
  )
write_csv(group_heat_df, file.path(result_dir, "Tcell_revision_group_DEG_key_gene_heatmap_source.csv"))

p_group_heat <- ggplot(group_heat_df, aes(Group, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_gradient2(low = "#6F99C3", mid = "#F7F7F7", high = "#B7192B", midpoint = 0, name = "Z-score") +
  labs(x = NULL, y = NULL) +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.2),
    axis.text.y = element_text(size = 7.8),
    strip.text.y = element_text(size = 7.8, angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey93", colour = NA),
    panel.grid = element_blank(),
    legend.position = "right"
  )
save_plot_safe(p_group_heat, "Tcell_group_DEG_key_gene_heatmap_revision", width = 4.8, height = 7.0)

tcr_sample <- meta %>%
  group_by(Sample, Group) %>%
  summarise(
    n_tcell = n(),
    n_tcr = sum(clone_ratio > 0, na.rm = TRUE),
    expanded_fraction = mean(clone_ratio > 1, na.rm = TRUE) * 100,
    mean_clone_ratio = mean(clone_ratio, na.rm = TRUE),
    max_clone_ratio = max(clone_ratio, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(tcr_sample, file.path(result_dir, "Tcell_revision_TCR_expansion_sample_summary.csv"))

tcr_wilcox <- bind_rows(
  pairwise_wilcox(tcr_sample, "expanded_fraction") %>% mutate(metric = "expanded_fraction"),
  pairwise_wilcox(tcr_sample, "mean_clone_ratio") %>% mutate(metric = "mean_clone_ratio"),
  pairwise_wilcox(tcr_sample, "max_clone_ratio") %>% mutate(metric = "max_clone_ratio")
)
write_csv(tcr_wilcox, file.path(result_dir, "Tcell_revision_TCR_expansion_wilcox_tests.csv"))

tcr_long <- tcr_sample %>%
  select(Sample, Group, expanded_fraction, mean_clone_ratio, max_clone_ratio) %>%
  pivot_longer(cols = c(expanded_fraction, mean_clone_ratio, max_clone_ratio), names_to = "metric", values_to = "value")
p_tcr_box <- ggplot(tcr_long, aes(Group, value, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.36, alpha = 0.82) +
  geom_point(position = position_jitter(width = 0.11, height = 0, seed = 20260813), shape = 21, colour = "black", stroke = 0.18, size = 1.35, alpha = 0.88) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(x = NULL, y = "TCR expansion metric") +
  theme_pub(8.8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7.0), strip.text = element_text(size = 7.8), legend.position = "top", panel.grid.major.x = element_blank())
save_plot_safe(p_tcr_box, "Tcell_TCR_expansion_group_compare_revision", width = 5.8, height = 2.8)

tcr_dist <- meta %>%
  count(Sample, Group, T_celltype, cloneSize, name = "n") %>%
  group_by(Sample, Group, T_celltype) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()
write_csv(tcr_dist, file.path(result_dir, "Tcell_revision_TCR_cloneSize_by_subtype_sample.csv"))
p_tcr_dist <- tcr_dist %>%
  group_by(Group, T_celltype, cloneSize) %>%
  summarise(percent = mean(percent, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(T_celltype, percent, fill = cloneSize)) +
  geom_col(width = 0.74, colour = "white", linewidth = 0.15) +
  facet_grid(Group ~ ., scales = "free_y") +
  scale_fill_manual(values = clone_colors, drop = FALSE, name = "Clone size") +
  labs(x = NULL, y = "Mean cells (%)") +
  theme_pub(8.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.8), strip.text.y = element_text(angle = 0), legend.position = "right", panel.grid.major.x = element_blank())
save_plot_safe(p_tcr_dist, "Tcell_TCR_clone_expansion_distribution_revision", width = 6.8, height = 5.2)

saveRDS(Tcell, file.path(result_dir, "CSF_Tcell_revision_TCR_annotated.rds"))
message("Updated CSF T cell subtype, function and TCR revision analysis.")

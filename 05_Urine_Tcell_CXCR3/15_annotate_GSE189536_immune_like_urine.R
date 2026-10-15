#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: External Harmony annotation and final CXCR3/module validation.
# Original source: annotate_GSE189536_immune_like_urine.R
# Figure/table scope: FigS10a; FigS10b; FigS10c; FigS10d; FigS10e; FigS10f
# Data cohort: gse189536. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# Harmonize GSE189536 biopsy immune-cell annotations to the urine scRNA-seq
# major-celltype and T-cell/myeloid subtype naming style.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

env_data_dir <- Sys.getenv("GSE189536_DATA_DIR", unset = "")
if (nzchar(env_data_dir)) {
  data_dir <- gsub("\\\\", "/", env_data_dir)
  data_dir <- sub("/$", "", data_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_dir <- if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  } else {
    file.path(getOption("hla.data_root"), "gse189536")
  }
  if (!dir.exists(script_dir)) script_dir <- file.path(getOption("hla.data_root"), "gse189536")
  data_dir <- normalizePath(script_dir, winslash = "/", mustWork = TRUE)
}
manual_dir <- file.path(data_dir, "Harmony_immune_cells_manual_annotation_results")
out_dir <- file.path(manual_dir, "Urine_style_annotation")
fig_dir <- file.path(out_dir, "Figures")
tab_dir <- file.path(out_dir, "Tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(getOption("hla.data_root"), "gse189536", "Immune_filtered.rds")

major_levels <- c("CD4+ T", "CD8+ T", "Cycling T", "NK", "B cell", "Monocyte", "Macrophage", "Neutrophil", "cDC", "pDC")
urine_revised_levels <- c(
  "Naive CD4+ T", "Naive CD8+ T", "CD4+ Tem", "CD8+ Tem", "CD8+ NKT",
  "NK", "Cycling T", "B", "C1QC+ Macro", "TREM2+ Macro", "MT1A+ Macro",
  "CD16+ Mono", "THBS1+ Mono", "NLRP3+ Mono", "Neutrophil", "cDC", "pDC"
)
tcell_levels <- c("Naive CD4", "CD4 Tem", "CD4 Tex", "CD8 Tex-like", "Treg", "Naive CD8", "CXCR4+ CD8", "GZMB+ CD8", "ISG+ CD8", "Cycling T")
group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")

major_colors <- c(
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
revised_colors <- c(
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
tcell_colors <- c(
  "Naive CD4" = "#E42320",
  "CD4 Tem" = "#699ECA",
  "CD4 Tex" = "#D95F02",
  "CD8 Tex-like" = "#B15928",
  "Treg" = "#7570B3",
  "Naive CD8" = "#4DAF4A",
  "CXCR4+ CD8" = "#1F78B4",
  "GZMB+ CD8" = "#FF8C00",
  "ISG+ CD8" = "#F898CB",
  "Cycling T" = "#731A73"
)
group_colors <- c("NR" = "#E42320", "TCMR" = "#699ECA", "ABMR" = "#4DAF4A", "Mix" = "#FF8C00", "BK" = "#F898CB")

marker_sets <- list(
  "Naive CD4" = c("CD4", "CCR7", "SELL", "IL7R", "TCF7", "LEF1"),
  "CD4 Tem" = c("CD4", "IL7R", "CD40LG", "ANXA1", "S100A4"),
  "CD4 Tex" = c("CD4", "PDCD1", "CTLA4", "TIGIT", "LAG3", "HAVCR2", "TOX"),
  "Treg" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT"),
  "Naive CD8" = c("CD8A", "CD8B", "CCR7", "SELL", "TCF7", "LEF1"),
  "CXCR4+ CD8" = c("CD8A", "CD8B", "CXCR4"),
  "GZMB+ CD8" = c("CD8A", "CD8B", "GZMB", "PRF1", "GNLY", "NKG7"),
  "ISG+ CD8" = c("CD8A", "CD8B", "IFIT1", "IFIT3", "ISG15", "MX1"),
  "Cycling T" = c("MKI67", "TOP2A", "STMN1"),
  "CD16+ Mono" = c("FCGR3A", "MS4A7", "LST1", "AIF1"),
  "THBS1+ Mono" = c("THBS1", "VCAN", "S100A8", "S100A9", "FCN1"),
  "NLRP3+ Mono" = c("NLRP3", "IL1B", "S100A8", "S100A9", "CXCL8"),
  "C1QC+ Macro" = c("C1QA", "C1QB", "C1QC", "APOE", "MRC1"),
  "TREM2+ Macro" = c("TREM2", "SPP1", "APOE", "TYROBP", "C1QA"),
  "MT1A+ Macro" = c("MT1A", "MT2A", "MT1E", "MT1X", "C1QA")
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

get_data_matrix <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "data")
  )
}

zscore_safe <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

score_marker_set <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  if (length(genes) == 0) return(rep(NA_real_, ncol(expr)))
  Matrix::colMeans(expr[genes, , drop = FALSE])
}

pick_max_label <- function(score_df, labels) {
  mat <- as.matrix(score_df[, labels, drop = FALSE])
  mat[!is.finite(mat)] <- -Inf
  labels[max.col(mat, ties.method = "first")]
}

theme_umap <- function(base_size = 9.0) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(
    file.path(fig_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    device = function(filename, width, height, units, res, bg, ...) {
      grDevices::png(filename, width = width, height = height, units = units, res = res, bg = bg, type = "cairo")
    }
  )
}

message("Loading existing GSE189536 immune object...")
obj <- tryCatch(
  readRDS(input_rds),
  error = function(e) stop("Cannot read input RDS: ", input_rds, "\n", conditionMessage(e))
)
DefaultAssay(obj) <- "RNA"
expr <- get_data_matrix(obj, assay = "RNA")

meta <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    Group = case_when(
      sample_type %in% group_levels ~ as.character(sample_type),
      diagnosis == "No rejection" ~ "NR",
      diagnosis == "Mixed" ~ "Mix",
      TRUE ~ as.character(diagnosis)
    ),
    Group = ifelse(Group %in% group_levels, Group, NA_character_)
  )

available_marker_df <- tibble(
  module = rep(names(marker_sets), lengths(marker_sets)),
  gene = unlist(marker_sets, use.names = FALSE)
) %>%
  mutate(found = gene %in% rownames(expr)) %>%
  arrange(module, desc(found), gene)
# Export the analysis source data for figure reproduction.
write_csv_safe(available_marker_df, file.path(tab_dir, "GSE189536_urine_style_marker_genes_used.csv"))

message("Scoring urine-style marker modules...")
score_df <- as.data.frame(lapply(marker_sets, function(genes) score_marker_set(expr, genes)), check.names = FALSE)
score_df[] <- lapply(score_df, zscore_safe)
score_df$cell <- colnames(expr)

meta2 <- meta %>%
  left_join(score_df, by = "cell")

t_like <- c("CD4 T", "CD8 T", "Treg", "gamma-delta T")
tcell_allowed <- list(
  "CD4 T" = c("Naive CD4", "CD4 Tem", "CD4 Tex", "Cycling T"),
  "Treg" = c("Treg", "CD4 Tex", "CD4 Tem", "Cycling T"),
  "CD8 T" = c("Naive CD8", "CXCR4+ CD8", "GZMB+ CD8", "ISG+ CD8", "Cycling T"),
  "gamma-delta T" = c("GZMB+ CD8", "ISG+ CD8", "CXCR4+ CD8", "Cycling T")
)

meta2$Urine_Tcell_revised_celltype <- NA_character_
for (raw_type in names(tcell_allowed)) {
  idx <- which(meta2$manual_cell_type == raw_type)
  if (length(idx) > 0) {
    meta2$Urine_Tcell_revised_celltype[idx] <- pick_max_label(meta2[idx, , drop = FALSE], tcell_allowed[[raw_type]])
  }
}
meta2$Urine_Tcell_revised_celltype[meta2$manual_cell_type == "Treg"] <- "Treg"
meta2$Urine_Tcell_revised_celltype[meta2$manual_cell_type == "gamma-delta T" & is.na(meta2$Urine_Tcell_revised_celltype)] <- "GZMB+ CD8"

lineage_genes <- intersect(c("CD4", "CD8A", "CD8B"), rownames(obj))
if (length(lineage_genes) > 0) {
  lineage_df <- FetchData(obj, vars = lineage_genes) %>%
    rownames_to_column("cell")
  meta2 <- meta2 %>%
    left_join(lineage_df, by = "cell") %>%
    mutate(
      CD4_expr = if ("CD4" %in% colnames(.)) as.numeric(CD4) else 0,
      CD8A_expr = if ("CD8A" %in% colnames(.)) as.numeric(CD8A) else 0,
      CD8B_expr = if ("CD8B" %in% colnames(.)) as.numeric(CD8B) else 0,
      Urine_Tcell_revised_celltype = if_else(
        Urine_Tcell_revised_celltype == "CD4 Tex" &
          (CD8A_expr + CD8B_expr) > CD4_expr &
          (CD8A_expr + CD8B_expr) > 0,
        "CD8 Tex-like",
        as.character(Urine_Tcell_revised_celltype)
      )
    ) %>%
    select(-any_of(c("CD4", "CD8A", "CD8B", "CD4_expr", "CD8A_expr", "CD8B_expr")))
}

myeloid_allowed <- c("CD16+ Mono", "THBS1+ Mono", "NLRP3+ Mono", "C1QC+ Macro", "TREM2+ Macro", "MT1A+ Macro")
myeloid_idx <- which(
  (meta2$manual_cell_type == "Myeloid" | meta2$manual_immune_cell_type %in% c("Myeloid")) &
    !meta2$manual_immune_cell_type %in% c("mDC", "pDC", "Neutrophils")
)
myeloid_label <- character(length(myeloid_idx))
if (length(myeloid_idx) > 0) {
  myeloid_label <- pick_max_label(meta2[myeloid_idx, , drop = FALSE], myeloid_allowed)
}

meta2$Urine_revised_celltype <- case_when(
  meta2$manual_cell_type == "B cell" ~ "B",
  meta2$manual_cell_type == "Plasma" ~ "B",
  meta2$manual_cell_type == "NK" ~ "NK",
  meta2$manual_cell_type == "Neutrophil" ~ "Neutrophil",
  meta2$manual_immune_cell_type == "mDC" ~ "cDC",
  meta2$manual_immune_cell_type == "pDC" ~ "pDC",
  meta2$manual_cell_type == "gamma-delta T" ~ "CD8+ NKT",
  meta2$manual_cell_type %in% t_like & meta2$Urine_Tcell_revised_celltype == "Naive CD4" ~ "Naive CD4+ T",
  meta2$manual_cell_type %in% t_like & meta2$Urine_Tcell_revised_celltype %in% c("CD4 Tem", "CD4 Tex", "Treg") ~ "CD4+ Tem",
  meta2$manual_cell_type %in% t_like & meta2$Urine_Tcell_revised_celltype == "Naive CD8" ~ "Naive CD8+ T",
  meta2$manual_cell_type %in% t_like & meta2$Urine_Tcell_revised_celltype %in% c("CXCR4+ CD8", "GZMB+ CD8", "ISG+ CD8", "CD8 Tex-like") ~ "CD8+ Tem",
  meta2$manual_cell_type %in% t_like & meta2$Urine_Tcell_revised_celltype == "Cycling T" ~ "Cycling T",
  TRUE ~ NA_character_
)
meta2$Urine_revised_celltype[myeloid_idx] <- myeloid_label

is_unassigned_isg <- meta2$manual_cell_type == "ISG+ immune" | is.na(meta2$Urine_revised_celltype)
if (any(is_unassigned_isg)) {
  broad_scores <- meta2[is_unassigned_isg, , drop = FALSE] %>%
    transmute(
      `CD8+ Tem` = pmax(`GZMB+ CD8`, `ISG+ CD8`, na.rm = TRUE),
      NK = `GZMB+ CD8`,
      `NLRP3+ Mono` = `NLRP3+ Mono`,
      `C1QC+ Macro` = `C1QC+ Macro`,
      B = 0
    )
  meta2$Urine_revised_celltype[is_unassigned_isg] <- pick_max_label(broad_scores, colnames(broad_scores))
}

meta2$Urine_major_celltype <- case_when(
  meta2$Urine_revised_celltype %in% c("Naive CD4+ T", "CD4+ Tem") ~ "CD4+ T",
  meta2$Urine_revised_celltype %in% c("Naive CD8+ T", "CD8+ Tem", "CD8+ NKT") ~ "CD8+ T",
  meta2$Urine_revised_celltype == "Cycling T" ~ "Cycling T",
  meta2$Urine_revised_celltype == "NK" ~ "NK",
  meta2$Urine_revised_celltype == "B" ~ "B cell",
  meta2$Urine_revised_celltype %in% c("CD16+ Mono", "THBS1+ Mono", "NLRP3+ Mono") ~ "Monocyte",
  meta2$Urine_revised_celltype %in% c("C1QC+ Macro", "TREM2+ Macro", "MT1A+ Macro") ~ "Macrophage",
  meta2$Urine_revised_celltype == "Neutrophil" ~ "Neutrophil",
  meta2$Urine_revised_celltype == "cDC" ~ "cDC",
  meta2$Urine_revised_celltype == "pDC" ~ "pDC",
  TRUE ~ NA_character_
)

meta2 <- meta2 %>%
  mutate(
    Urine_major_celltype = factor(Urine_major_celltype, levels = major_levels),
    Urine_revised_celltype = factor(Urine_revised_celltype, levels = urine_revised_levels),
    Urine_Tcell_revised_celltype = factor(Urine_Tcell_revised_celltype, levels = tcell_levels),
    Group = factor(Group, levels = group_levels)
  )

obj <- AddMetaData(obj, setNames(meta2$Group, meta2$cell), col.name = "Group")
obj <- AddMetaData(obj, setNames(meta2$Urine_major_celltype, meta2$cell), col.name = "Urine_major_celltype")
obj <- AddMetaData(obj, setNames(meta2$Urine_revised_celltype, meta2$cell), col.name = "Urine_revised_celltype")
obj <- AddMetaData(obj, setNames(meta2$Urine_Tcell_revised_celltype, meta2$cell), col.name = "Urine_Tcell_revised_celltype")

cell_metadata <- meta2 %>%
  select(
    cell, biopsy, diagnosis, sample_type, Group, seurat_clusters,
    cell_type, immune_cell_type, manual_cell_type,
    Urine_major_celltype, Urine_revised_celltype, Urine_Tcell_revised_celltype
  )
write_csv_safe(cell_metadata, file.path(tab_dir, "GSE189536_urine_style_annotation_cell_metadata.csv"))

cluster_summary <- cell_metadata %>%
  group_by(seurat_clusters, manual_cell_type) %>%
  summarise(
    cells = n(),
    Urine_major_celltype = names(sort(table(Urine_major_celltype), decreasing = TRUE))[1],
    Urine_revised_celltype = names(sort(table(Urine_revised_celltype), decreasing = TRUE))[1],
    Urine_Tcell_revised_celltype = ifelse(all(is.na(Urine_Tcell_revised_celltype)), NA_character_, names(sort(table(Urine_Tcell_revised_celltype), decreasing = TRUE))[1]),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(as.character(seurat_clusters)))
write_csv_safe(cluster_summary, file.path(tab_dir, "GSE189536_urine_style_annotation_cluster_summary.csv"))

composition_by_group <- cell_metadata %>%
  filter(!is.na(Group), !is.na(Urine_major_celltype)) %>%
  count(Group, biopsy, Urine_major_celltype, name = "cells") %>%
  group_by(Group, biopsy) %>%
  mutate(freq = cells / sum(cells) * 100) %>%
  ungroup()
write_csv_safe(composition_by_group, file.path(tab_dir, "GSE189536_urine_style_major_celltype_composition_by_biopsy.csv"))

saveRDS(obj, file.path(getOption("hla.data_root"), "gse189536", "Immune_annotated.rds"))

if (!"umap" %in% names(obj@reductions)) stop("UMAP reduction is missing in the input object.")
umap_df <- Embeddings(obj, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  left_join(cell_metadata, by = "cell")
colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")
write_csv_safe(umap_df, file.path(tab_dir, "GSE189536_urine_style_annotation_umap_metadata.csv"))

p_major <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Urine_major_celltype)) +
  geom_point(size = 0.22, alpha = 0.82, stroke = 0) +
  coord_equal() +
  scale_colour_manual(values = major_colors, na.value = "grey82", drop = TRUE) +
  labs(title = "Major celltype", colour = "Major celltype", x = "UMAP 1", y = "UMAP 2") +
  theme_umap(9.0) +
  guides(colour = guide_legend(override.aes = list(size = 3.2, alpha = 1)))

p_revised <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Urine_revised_celltype)) +
  geom_point(size = 0.22, alpha = 0.82, stroke = 0) +
  coord_equal() +
  scale_colour_manual(values = revised_colors, na.value = "grey82", drop = TRUE) +
  labs(title = "Urine-style revised celltype", colour = "Revised celltype", x = "UMAP 1", y = "UMAP 2") +
  theme_umap(9.0) +
  guides(colour = guide_legend(override.aes = list(size = 3.0, alpha = 1), ncol = 1))

tcell_umap <- umap_df %>%
  mutate(Tcell_plot = ifelse(is.na(Urine_Tcell_revised_celltype), "Non-T immune", as.character(Urine_Tcell_revised_celltype)))
tcell_cols <- c(tcell_colors, "Non-T immune" = "#D9D9D9")
p_tcell <- ggplot(tcell_umap, aes(UMAP_1, UMAP_2)) +
  geom_point(data = tcell_umap %>% filter(Tcell_plot == "Non-T immune"), colour = "#D9D9D9", size = 0.16, alpha = 0.38, stroke = 0) +
  geom_point(data = tcell_umap %>% filter(Tcell_plot != "Non-T immune"), aes(colour = Tcell_plot), size = 0.24, alpha = 0.88, stroke = 0) +
  coord_equal() +
  scale_colour_manual(values = tcell_cols, drop = TRUE) +
  labs(title = "T-cell subtype", colour = "T-cell subtype", x = "UMAP 1", y = "UMAP 2") +
  theme_umap(9.0) +
  guides(colour = guide_legend(override.aes = list(size = 3.0, alpha = 1)))

p_overview <- p_major + p_revised + p_tcell + plot_layout(ncol = 3, widths = c(1, 1.08, 1))
save_plot(p_overview, "GSE189536_urine_style_annotation_UMAP_overview", width = 13.2, height = 4.4)

message("Reclustering standalone T cells with Harmony-corrected embeddings...")
tcell_cells <- cell_metadata %>%
  filter(!is.na(Urine_Tcell_revised_celltype)) %>%
  pull(cell)
if (length(tcell_cells) < 50) {
  warning("Too few annotated T cells for standalone reclustering.")
} else if (!"harmony" %in% names(obj@reductions)) {
  warning("Harmony reduction is missing; standalone T-cell UMAP was not regenerated.")
} else {
  tobj <- subset(obj, cells = tcell_cells)
  DefaultAssay(tobj) <- "RNA"
  dims_use <- seq_len(min(20, ncol(Embeddings(tobj, "harmony"))))
  set.seed(20260827)
  tobj <- FindNeighbors(tobj, reduction = "harmony", dims = dims_use, verbose = FALSE)
  tobj <- FindClusters(tobj, resolution = 0.65, verbose = FALSE)
  tobj <- RunUMAP(
    tobj,
    reduction = "harmony",
    dims = dims_use,
    reduction.name = "tcell_harmony_umap",
    reduction.key = "tHUMAP_",
    verbose = FALSE
  )
  saveRDS(tobj, file.path(getOption("hla.data_root"), "gse189536", "Tcell_harmony.rds"))

  t_umap <- Embeddings(tobj, "tcell_harmony_umap") %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    left_join(
      tobj@meta.data %>%
        rownames_to_column("cell") %>%
        transmute(
          cell,
          seurat_clusters = as.character(seurat_clusters),
          Tcell_subtype = factor(as.character(Urine_Tcell_revised_celltype), levels = tcell_levels),
          Group = factor(as.character(Group), levels = group_levels),
          biopsy = as.character(biopsy)
        ),
      by = "cell"
    )
  colnames(t_umap)[2:3] <- c("UMAP_1", "UMAP_2")
  write_csv_safe(t_umap, file.path(tab_dir, "GSE189536_Tcell_standalone_recluster_UMAP_metadata.csv"))
  write_csv_safe(
    t_umap %>% count(seurat_clusters, Tcell_subtype, Group, name = "cells"),
    file.path(tab_dir, "GSE189536_Tcell_harmony_cluster_celltype_group_counts.csv")
  )

  cluster_centers <- t_umap %>%
    group_by(seurat_clusters) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), cells = n(), .groups = "drop")

  cluster_levels <- sort(unique(t_umap$seurat_clusters))
  cluster_cols <- setNames(hue_pal(l = 58, c = 90)(length(cluster_levels)), cluster_levels)
  p_t_cluster <- ggplot(t_umap, aes(UMAP_1, UMAP_2, colour = seurat_clusters)) +
    geom_point(size = 0.30, alpha = 0.82, stroke = 0) +
    geom_text(
      data = cluster_centers,
      aes(x = UMAP_1, y = UMAP_2, label = seurat_clusters),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "bold",
      size = 3.1,
      colour = "black"
    ) +
    coord_equal() +
    scale_colour_manual(values = cluster_cols, guide = "none") +
    labs(title = "T-cell clusters", subtitle = "Harmony neighbors/clusters/UMAP", x = "UMAP 1", y = "UMAP 2") +
    theme_umap(9.0)

  p_t_subtype <- ggplot(t_umap, aes(UMAP_1, UMAP_2, colour = Tcell_subtype)) +
    geom_point(size = 0.32, alpha = 0.86, stroke = 0) +
    coord_equal() +
    scale_colour_manual(values = tcell_colors, na.value = "grey82", drop = TRUE) +
    labs(title = "T-cell subtype", colour = "T-cell subtype", x = "UMAP 1", y = "UMAP 2") +
    theme_umap(9.0) +
    guides(colour = guide_legend(override.aes = list(size = 3.2, alpha = 1)))

  p_t_group <- ggplot(t_umap, aes(UMAP_1, UMAP_2, colour = Group)) +
    geom_point(size = 0.28, alpha = 0.78, stroke = 0) +
    coord_equal() +
    scale_colour_manual(values = group_colors, na.value = "grey82", drop = TRUE) +
    labs(title = "Clinical group", colour = "Group", x = "UMAP 1", y = "UMAP 2") +
    theme_umap(9.0) +
    guides(colour = guide_legend(override.aes = list(size = 3.2, alpha = 1)))

  p_tcell_standalone <- p_t_cluster + p_t_subtype + p_t_group + plot_layout(ncol = 3, widths = c(1, 1.1, 1))
  save_plot(p_tcell_standalone, "GSE189536_urine_style_Tcell_standalone_recluster_UMAP", width = 12.8, height = 4.2)
}

dot_markers <- c(
  "CD3D", "CD4", "CD8A", "CCR7", "IL7R", "GZMK", "GZMB", "NKG7", "FOXP3", "MKI67",
  "MS4A1", "MZB1", "LYZ", "FCGR3A", "THBS1", "NLRP3", "C1QA", "TREM2", "S100A8", "CLEC9A", "LILRA4"
)
dot_markers <- intersect(dot_markers, rownames(obj))
p_dot <- DotPlot(obj, features = dot_markers, group.by = "Urine_revised_celltype", cols = c("#E6E6E6", "#B7192B"), dot.scale = 5.0) +
  RotatedAxis() +
  labs(x = NULL, y = NULL, colour = "Average\nexpression", size = "Percent\nexpressed") +
  theme_classic(base_size = 9.0, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.0, colour = "black"),
    axis.text.y = element_text(size = 8.0, colour = "black"),
    legend.title = element_text(face = "bold", size = 8.0),
    legend.text = element_text(size = 7.2)
  )
save_plot(p_dot, "GSE189536_urine_style_annotation_marker_dotplot", width = 8.8, height = 4.8)

message("Saved urine-style annotated object: ", file.path(getOption("hla.data_root"), "gse189536", "Immune_annotated.rds"))
message("Saved tables and figures to: ", out_dir)

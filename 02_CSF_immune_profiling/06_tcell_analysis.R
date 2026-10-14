# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Upstream cohort processing and reusable helper definitions.
# Original source: 05_tcell_analysis.R
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

# ==============================================================================
# File: 05_tcell_analysis.R
# Description: T 细胞亚群降维、可视化与热图分析
# ==============================================================================
rm(list = ls())
source(file.path(getOption("hla.code_root"), "02_CSF_immune_profiling/01_setup_and_helpers.R"))

# Read the saved object or source table used by this analysis.
scRNA <- readRDS("./Result/CSF_scRNA.rds")

# 1. Tcell 提取与降维 ----------------------------------------------------------
Tcell <- subset(scRNA, subset = Celltype %in% c("Naive T","CD4+ T", "CD8+ T", "IFN+ T", "Tregs", "Proliferating", "gd T"))
Tcell <- NormalizeData(Tcell) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%
  RunHarmony(group.by.vars = "orig.ident") %>% 
  RunUMAP(dims = 1:20, reduction = "harmony") %>%
  FindNeighbors(dims = 1:20, reduction = "harmony") %>% 
  FindClusters(resolution = 0.8, reduction = "harmony")

# UMAP 绘图
export_top_markers(Tcell, "seurat_clusters", "Tcell")

# 假设手动注释已完成
df_anno <- read.csv("./Result/Tcell_makers_summary.csv")
Tcell$T_celltype <- df_anno$Celltype[match(Tcell$seurat_clusters, df_anno$seurat_clusters)]
unique(Tcell$T_celltype) %>% mixedsort() %>% print()
Tcell = Tcell[,!Tcell$T_celltype %in% c("Contamination","T_cycling")]
Tcell$T_celltype <- factor(Tcell$T_celltype,
                           levels = c("Naive T", "CD4+ Tcm", "Th17", "IFN+ T", "CD4+ Tex", "Tregs",
                                      "CD8+ Tmem", "CD8+ Trm","CD8+ Tact", "CD8+ Teff",  "NK-like T", "gd T"))

# 绘制全局UMAP三联图 --------------------------------------------------------
p1 <- DimPlot(Tcell, group.by = "seurat_clusters", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p2 <- DimPlot(Tcell, group.by = "T_celltype", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p3 <- DimPlot(Tcell, group.by = "Group", label = TRUE, alpha = 0.5) + scale_color_manual(values = GROUP_COLORS)
ggsave("./Figures/CSF_UMAP_Tcell_celltype.pdf", p1 + p2 + p3, height = 5, width = 18)

# Tcell Celltype 基因绘图
Tcell_mk <- c("CCR7", "LEF1", "IL7R", "MAL", "CCR6", "AQP3",
              "IFIT3", "ISG15", "CTLA4", "PDCD1", "FOXP3", "IL2RA",
              "CD27", "CISH", "ZNF683", "FGFBP2", "GZMK",
              "GZMH", "GZMB", "KLRD1", "NKG7", "GNLY","TRDC")
p_dot_Tcell <- DotPlot(Tcell, group.by = "T_celltype", features = Tcell_mk) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("navy","grey","firebrick3"))
ggsave("./Figures/Tcell_Dotplot_clusters.pdf", plot = p_dot_Tcell, height = 4, width = 9)

# 细胞亚群比例图 -------------------------------------------------
comparison <- list(c("Metastases","Hemorrhage"), c("Metastases","Encephalitis"), c("Hemorrhage","Encephalitis"))
p_cell_ratio <- plot_cell_ratio_comparison(cell_subset = Tcell,
                                           sample_col = "HLA_Splitter", cluster_col =  "T_celltype", group_col = "Group",
                                           facet = T, comparison_groups = comparison,
                                           write_to_csv = "./Result/Tcell_celltype_ratio", plot_colors = GLOBAL_COLORS)
ggsave("./Figures/Tcell_celltype_compare.pdf", p_cell_ratio, width = 10, height = 10)

# Tcell Celltype 热图
p_list <- plot_celltype_ratio(seurat_obj = Tcell,
                              sample_col = "HLA_Splitter", group_col = "Group", cluster_col = "T_celltype",
                              comparison_groups = comparison, group_colors = GROUP_COLORS, cluster_colors = GLOBAL_COLORS)

ggsave("./Figures/Tcell_celltype_stacked_bar.pdf", p_list$stacked_bar, width = 5, height = 4)
ggsave("./Figures/Tcell_celltype_heatmap.pdf", p_list$heatmap, width = 6, height = 6)


# 2. Tcell 不同疾病分组中的差异基因 ------------------------------------------
Idents(Tcell) <- "Group"
# Discover marker genes with the original assay and comparison settings.
CELLDEG <- FindMarkers(Tcell, ident.1 ="Metastases", ident.2 = "Encephalitis", logfc.threshold = 0.1, min.pct = 0.1)
CELLDEG <- CELLDEG[!is.na(CELLDEG$p_val),]
CELLDEG$gene <- rownames(CELLDEG)
# Export the analysis source data for figure reproduction.
write.csv(CELLDEG, "./Result/Tcell_Group_marker_Metastases-Encephalitis.csv")

# 火山图 & GO 富集
p_volcano <- plot_volcano(diffexp = CELLDEG,
                          top_n_labels = 10,
                          pos_group_name = "Metastases", 
                          neg_group_name = "Encephalitis")

ggsave("./Figures/Tcell_Group_Volcano_plot_Metastases-Encephalitis.pdf", p_volcano, height = 6, width = 8)

p_enrich <- Enrichment_analysis(diffexp = CELLDEG, 
                                dbs = c("GO_Biological_Process_2023"), 
                                pos_group_name = "Metastases", 
                                neg_group_name = "Hemorrhage",
                                p_val_adj_threshold = 0.05, 
                                log2fc_threshold = 1, 
                                top_n_pathways = 6, 
                                plot_colors = GLOBAL_COLORS, 
                                output_csv_path = "./Result/Tcell_Group_enriched_GO.csv")
ggsave("./Figures/Tcell_Group_enriched_GO_Metastases-Encephalitis.pdf", p_enrich, height = 6, width = 4)

# 3. T细胞 特征基因可视化 ------------------------------------------------------
# T细胞特征基因小提琴图
selected_gene_vln <- c("GZMK","GZMA","GZMH","GZMB","IFITM1","IFI6","IRF7","ISG15")
plots_violins <- VlnPlot(Tcell, cols = GLOBAL_COLORS[1:3], pt.size = 0, group.by = "Group",
                         features = selected_gene_vln, ncol = 3, combine = FALSE)

for(i in seq_along(plots_violins)) {
  plots_violins[[i]] <- plots_violins[[i]] + theme_classic() + 
    geom_boxplot(width = 0.1, fill = "white", outliers = FALSE) +
    theme(axis.text.x = element_text(size = 10, color="black", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 10, color="black"),
          axis.title.y= element_text(size=12, color="black"),
          axis.title.x = element_blank(), legend.position='none') +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    stat_compare_means(aes(label = ..p.signif..), comparisons = comparison, step.increase = 0.2, vjust = -0.05)
}
P <- CombinePlots(plots_violins, ncol = 4)
ggsave("./Figures/Tcell_markergenes_violinplot.pdf", P, height = 5, width = 8)

# 3. Heatmap (Bulk & Single-cell) ----------------------------------------------
tcell_selected_gene <- c("GZMK","CISH","PTGER4","TIGIT","CXCR3","CXCR6")
Tcell_avg <- AverageExpression(Tcell, group.by = "Group", assays = "RNA", return.seurat = TRUE)
gene_matrix_t <- Tcell_avg@assays$RNA$scale.data %>% as.matrix()

plot_heatmap_t_bulk <- pheatmap(gene_matrix_t[tcell_selected_gene, ], cluster_rows = TRUE, cluster_cols = FALSE, angle_col = 45, cellwidth = 40, cellheight = 15)
ggsave("./Figures/Tcell_Heatmap_bulk.pdf", plot_heatmap_t_bulk, height = 10, width = 10)

saveRDS(Tcell, "./Result/CSF_Tcell.rds")

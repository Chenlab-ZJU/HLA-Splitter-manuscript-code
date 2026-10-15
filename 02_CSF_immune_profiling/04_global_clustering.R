# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Upstream cohort processing and reusable helper definitions.
# Original source: 03_global_clustering.R
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

# ==============================================================================
# File: 03_global_clustering.R
# Description: 全局降维聚类、注释与全局可视化
# ==============================================================================
rm(list = ls())
source(file.path(getOption("hla.code_root"), "02_CSF_immune_profiling/01_setup_and_helpers.R"))

# Read the saved object or source table used by this analysis.
scRNA_Intergrate <- readRDS("Result/CSF_scRNA_Intergrate.rds")
#scRNA = readRDS("./Result/CSF_scRNA.rds")

# 1. 数据筛选与重新聚类 --------------------------------------------------------
scRNA <- subset(scRNA_Intergrate, subset = !HLA_Splitter %in% c("Doublets", "Undefined"))
SampleInfo <- read.csv(file.path(getOption("hla.data_root"), "csf", "SampleInfo.csv"))
Sample_levels <- c(mixedsort(SampleInfo$Patient.ID))
scRNA$Group = factor(scRNA$Group, setdiff(scRNA$Group,c("Undefined","Doublet")))
# 
# valid_samples <- Sample_levels[table(scRNA$HLA_Splitter) > 100]
# scRNA <- subset(scRNA, subset = HLA_Splitter %in% valid_samples)

scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%
  RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.8)

# 2. 全局 Marker 提取与可视化 --------------------------------------------------
export_top_markers(scRNA, "seurat_clusters", "scRNA")

# 注意：运行到此处请暂停，确保 ./Result/scRNA_makers_summary.csv 已填写 Celltype
df <- read.csv("./Result/scRNA_makers_summary.csv")
Celltype <- df$Celltype
names(Celltype) <- as.character(df$seurat_clusters)
scRNA <- SetIdent(scRNA, value = "seurat_clusters")
scRNA <- RenameIdents(scRNA, Celltype)
scRNA$Celltype <- Idents(scRNA)
unique(scRNA$Celltype) %>% print()

# 为了演示，假定 Celltype 已经映射好
scRNA = scRNA[,scRNA$Celltype != "Neuronal-like"]
scRNA$Celltype <- factor(scRNA$Celltype, levels = c("Naive T","CD4+ Tact", "IFN+ T", "Tregs", "CD8+ CTL","CD8+ Trm", "CD8+ Tmem",
                                                    "gd T", "B cells", "NK cells","Microglia-like", "Macrophages",
                                                    "Monocytes", "cDC1", "cDC2", "pDC","Proliferating"))

# 3. 绘制全局UMAP三联图 --------------------------------------------------------
p1 <- DimPlot(scRNA, group.by = "seurat_clusters", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p2 <- DimPlot(scRNA, group.by = "Celltype", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p3 <- DimPlot(scRNA, group.by = "Group", label = TRUE, alpha = 0.5) + scale_color_manual(values = GROUP_COLORS)
ggsave("./Figures/CSF_UMAP_celltype.pdf", p1 + p2 + p3, height = 5, width = 18)

p_F = FeaturePlot(scRNA, features = c("C1QC", "FCGBP", "TREM2", "SPP1"),ncol=2,cols = c("grey","firebrick"))
ggsave("./Figures/CSF_FeaturePlot.pdf", p_F, height = 8, width = 8)

# UMAP 分面图 
combined_plot <- ggarrange(plotlist = plot_grouped_umap(scRNA, "Celltype", "Group", GLOBAL_COLORS), common.legend = TRUE, legend = "right")
ggsave("./Figures/scRNA_all_celltype_facet.pdf", combined_plot, height = 6, width = 8)

# Celltype 热图
p_list <- plot_celltype_ratio(seurat_obj = scRNA,
                              sample_col = "HLA_Splitter", group_col = "Group", cluster_col = "Celltype",
                              comparison_groups = comparison, group_colors = GROUP_COLORS, cluster_colors = GLOBAL_COLORS)

ggsave("./Figures/scRNA_celltype_stacked_bar.pdf", p_list$stacked_bar, width = 5, height = 4)
ggsave("./Figures/scRNA_celltype_heatmap.pdf", p_list$heatmap, width = 5, height = 5)

# 全局 Marker 基因绘图
global_mk <- c("TCF7","IL7R", "KLRB1", "CCR6", "IFIT3","ISG15",
               "CTLA4","IL2RA","CD8A", "PRF1", "ZNF683", "FGFBP2","CD27","CISH",
               "XCL1","TRDC","MS4A1","BANK1", "KLRD1", "GNLY", 
               "P2RY12","TREM2", "FCGR3A", "CD163","S100A8","CD14",
               "XCR1", "CLEC9A","CD1C","FCER1A","LILRA4","CLEC4C","MKI67","TYMS")
p_dot <- DotPlot(scRNA, group.by = "Celltype", features = global_mk) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x=NULL, y=NULL) +
  scale_color_gradientn(colours = c("navy","white","firebrick3"))
ggsave("./Figures/CSF_Dotplot_clusters.pdf", plot = p_dot, height = 4, width = 9)

# 细胞亚群比例图 -------------------------------------------------
p_cell_ratio <- plot_cell_ratio_comparison(cell_subset = scRNA,
                                           sample_col = "HLA_Splitter", 
                                           cluster_col =  "Celltype",
                                           group_col = "Group",
                                           facet = T, comparison_groups = comparison,
                                           stat_method = "t.test",
                                           write_to_csv = "./Result/scRNA_celltype_ratio", 
                                           plot_colors = GLOBAL_COLORS)
ggsave("./Figures/scRNA_celltype_compare.pdf", p_cell_ratio, width = 10, height = 10)

# 2. 疾病两两分组中的差异基因 ------------------------------------------
Idents(scRNA) <- "Group"
# Discover marker genes with the original assay and comparison settings.
CELLDEG <- FindMarkers(scRNA, ident.1 ="Metastases", ident.2 = "Encephalitis", logfc.threshold = 0.1, min.pct = 0.1)
CELLDEG <- CELLDEG[!is.na(CELLDEG$p_val),]
CELLDEG$gene <- rownames(CELLDEG)
# Export the analysis source data for figure reproduction.
write.csv(CELLDEG, "./Result/scRNA_Group_marker_Metastases-Encephalitis.csv")

# 火山图 & GO 富集
p_volcano <- plot_volcano(diffexp = CELLDEG,
                          top_n_labels = 10,
                          pos_group_name = "Metastases", 
                          neg_group_name = "Encephalitis")

ggsave("./Figures/scRNA_Group_Volcano_plot_Metastases-Encephalitis.pdf", p_volcano, height = 6, width = 8)

# 3. 不同疾病相较于平均值的差异基因 ------------------------------------------
Idents(scRNA) <- "Group"

CELLDEG <- FindAllMarkers(scRNA, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
CELLDEG <- CELLDEG[!is.na(CELLDEG$p_val), ]
write.csv(CELLDEG, "./Result/scRNA_Group_All_Markers.csv", row.names = FALSE)

# Pseudobulk分析，输出 Lollipop 图 + 分组专属 CSV + 分组专属 Bar&Dot 图
source(file.path(getOption("hla.code_root"), "00_Shared/plot_top_pseudobulk.R"))

plot_top_pseudobulk(scRNA, CELLDEG, group_by = "Group",sample_col = "HLA_Splitter",
                    top_n_lollipop = 20, top_n_bardot = 10, my_colors = GROUP_COLORS)

# 保存                                           
saveRDS(scRNA, file.path(getOption("hla.data_root"), "csf", "CSF_scRNA.rds"))

# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Upstream cohort processing and reusable helper definitions.
# Original source: 04_myeloid_analysis.R
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
# File: 04_myeloid_analysis.R
# Description: 巨噬细胞/髓系细胞亚群深入分析、M1/M2极化与 GSVA
# ==============================================================================
rm(list = ls())
source(file.path(getOption("hla.code_root"), "02_CSF_immune_profiling/01_setup_and_helpers.R"))

# Read the saved object or source table used by this analysis.
scRNA <- readRDS(file.path(getOption("hla.data_root"), "csf", "CSF_scRNA.rds"))

# 1. Myeloid 数据提取与降维 ----------------------------------------------------
Myeloid <- subset(scRNA, subset = Celltype %in% c("Macrophages", "Monocytes"))
Myeloid <- NormalizeData(Myeloid) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%
  RunUMAP(dims = 1:30, min.dist = 0.6) %>% FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.6)

export_top_markers(Myeloid, "seurat_clusters", "Myeloid")

# 假设手动注释已完成
df_myeloid_anno <- read.csv("./Result/Myeloid_makers_summary.csv")
Myeloid$Myeloid_celltype <- df_myeloid_anno$Celltype[match(Myeloid$seurat_clusters, df_myeloid_anno$seurat_clusters)]
unique(Myeloid$Myeloid_celltype) %>% mixedsort() %>% print()
Myeloid = Myeloid[,Myeloid$Myeloid_celltype != "Contamination"]
Myeloid$Myeloid_celltype <- factor(Myeloid$Myeloid_celltype,
                                   levels = c("Mono_CD14", "Macro_FSCN1","Macro_HMOX1", "Macro_IL1B","Macro_MT1G", "Macro_PPARG","Macro_P2RY12"))

# 绘制全局UMAP三联图 --------------------------------------------------------
p1 <- DimPlot(Myeloid, group.by = "seurat_clusters", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p2 <- DimPlot(Myeloid, group.by = "Myeloid_celltype", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p3 <- DimPlot(Myeloid, group.by = "Group", label = TRUE, alpha = 0.5) + scale_color_manual(values = GROUP_COLORS)
ggsave("./Figures/CSF_UMAP_Macrophages_celltype.pdf", p1 + p2 + p3, height = 5, width = 18)

# Myeloid Celltype 基因绘图
myeloid_mk <- c("S100A8", "S100A9", "FSCN1", "APOC2", "HMOX1", "CLEC10A", 
                "IL1B", "CXCL5", "MT1G", "SELENOP", "PPARG", "CDK14",
                "P2RY12", "SPP1")
p_dot_myeloid <- DotPlot(Myeloid, group.by = "Myeloid_celltype", features = myeloid_mk) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("navy","grey","firebrick3"))
ggsave("./Figures/Myeloid_Dotplot_clusters.pdf", plot = p_dot_myeloid, height = 4, width = 9)

# Myeloid Celltype 热图
comparison <- list(c("Metastases","Hemorrhage"), c("Metastases","Encephalitis"), c("Hemorrhage","Encephalitis"))
p_list <- plot_celltype_ratio(seurat_obj = Myeloid,
                              sample_col = "HLA_Splitter", group_col = "Group", cluster_col = "Myeloid_celltype",
                              comparison_groups = comparison, group_colors = GROUP_COLORS, cluster_colors = GLOBAL_COLORS)

ggsave("./Figures/Myeloid_celltype_stacked_bar.pdf", p_list$stacked_bar, width = 5, height = 4)
ggsave("./Figures/Myeloid_celltype_heatmap.pdf", p_list$heatmap, width = 6, height = 6)

# 2. Myeloid 不同疾病分组中的差异基因 ------------------------------------------
Idents(Myeloid) <- "Group"
# Discover marker genes with the original assay and comparison settings.
CELLDEG <- FindMarkers(Myeloid, ident.1 ="Metastases", ident.2 = "Hemorrhage", logfc.threshold = 0.1, min.pct = 0.1)
CELLDEG <- CELLDEG[!is.na(CELLDEG$p_val),]
CELLDEG$gene <- rownames(CELLDEG)
# Export the analysis source data for figure reproduction.
write.csv(CELLDEG, "./Result/Myeloid_Group_marker_Metastases-Hemorrhage.csv")

# 火山图 & GO 富集
p_volcano <- plot_volcano(diffexp = CELLDEG,
                          top_n_labels = 10,
                          pos_group_name = "Metastases", 
                          neg_group_name = "Hemorrhage")
                          
ggsave("./Figures/Myeloid_Group_Volcano_plot.pdf", p_volcano, height = 6, width = 8)

p_enrich <- Enrichment_analysis(diffexp = CELLDEG, 
                                dbs = c("GO_Biological_Process_2023"), 
                                pos_group_name = "Metastases", 
                                neg_group_name = "Hemorrhage",
                                p_val_adj_threshold = 0.05, 
                                log2fc_threshold = 1, 
                                top_n_pathways = 6, 
                                plot_colors = GLOBAL_COLORS, 
                                output_csv_path = "./Result/Myeloid_Group_enriched_GO.csv")
ggsave("./Figures/Myeloid_Group_enriched_GO.pdf", p_enrich, height = 6, width = 4)

# Myeloid Marker 基因绘图
myeloid_mk <- c("CD163", "MRC1", "MSR1", "LGALS3", "TREM2", "SPP1", "APOE", "APOC1", "CD68")

VlnPlot_myeloid <- VlnPlot(Myeloid, 
                           group.by = "Group", 
                           features = myeloid_mk, 
                           stack = TRUE,
                           pt.size = 0,
                           flip = TRUE) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(), 
    legend.position = "none" 
  )

ggsave("./Figures/Myeloid_VlnPlot_groups.pdf", plot = VlnPlot_myeloid, height = 6, width = 8)

myeloid_mk <- c(CD163, MRC1, MSR1, LGALS3,"TREM2", "SPP1","APOE","APOC1", "CD68")
VlnPlot_myeloid <- VlnPlot(Myeloid, group.by = "Group", features = myeloid_mk, stack = T) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("./Figures/Myeloid_VlnPlot_groups.pdf", plot = p_dot_myeloid, height = 6, width = 6)


# 2. M1/M2 极化分析 ------------------------------------------------------------
m1_genes <- c("CXCL9","CCL19","CXCL10","IDO1","EBI3","TNFAIP6","CCR7","LAMP3","CCL5","CD40","IRF8","CD38","CCL4","APOL3","CXCL11","ADAMDEC1","HCK","BCL2A1","RSAD2","CYP27B1","KYNU","SIGLEC1","SLAMF1","SLC2A6","IF144L","SAMSN1","LILRB2","TLR2","CCL8","APOBEC3A","AQP9","C3AR1","HLA-DQA1","MNDA","IL2RA","SLC15A3","CD80","ADGRE1","CHI3L1","RASSF4","IL4R","AIM2","PLA1A","MMP9","TLR8","HESX1","SOCS1","LST1","NOD2","IL7R","FLVCR2","BIRC3","PTGIR","CHI3L2","PLA2G7","GPR183","AIF1","CD4","HLA-DOB","MSC","ADGRE2","CCR5","CLEC2D","INFRSF4","CD86","SPIB","LAG3","DHX58","CLIC2","ACHE","CCL18","LY86","CCL14","TNIP3","Clorf54","MMP25","ACP5","APOBEC3G")
m2_genes <- c("MMP9","MS4A6A","CCL18","AIF1","NCF2","CD4","ACP5","CLEC4A","CCL13","LY86","SLC15A3","CLEC10A","CCL23","HLA-DQA1","RNASE6","ADAMDEC1","HCK","NPL","TREM2","IRF8","SIGLEC1","CCL4","SAMSN1","CLEC7A","TLR2","P2RY13","CCL8","CD86","CD180","CD68","CD209","C3AR1","GPR183","CCL14","DPEP2","FPR3","CD37","LST1","CLIC2","HRH1","EGR2","CHI3L1","Clorf54","CFP","C5AR1","MNDA","PTGER2","CD300A","PIK3IP1","IGSF6","RENBP","PLA2G7","FCGR2B","SLAMF8","RASGRP3","QPCT","LILRB2","ATP8B4","PTGIR","FRMD4A","TLR8","KYNU","CHST15")
# Calculate the specified module score using the retained gene set.
Myeloid <- AddModuleScore(Myeloid, features = list(m1_genes), name = "M1_score")
Myeloid <- AddModuleScore(Myeloid, features = list(m2_genes), name = "M2_score")

plot_data_polar <- data.frame(M1_score = Myeloid$M1_score1, M2_score = Myeloid$M2_score1, Group = Myeloid$Group)
p_polar <- ggplot(plot_data_polar, aes(x = M1_score, y = M2_score, color = Group)) +
  geom_point(size = 1.5, alpha = 0.7) + scale_color_manual(values = GROUP_COLORS) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 1) + theme_bw()
p_margin <- ggMarginal(p_polar, type = "boxplot", margins = "both", groupColour = TRUE, groupFill = T, alpha = 0.5)
ggsave("./Figures/Myeloid_Polarization.pdf", p_margin, width = 8, height = 6)

# 3. GSVA 免疫通路分析 ---------------------------------------------------------
Myeloid_avg <- AverageExpression(Myeloid, group.by = "Group", assays = "RNA", return.seurat = TRUE)
gene_matrix <- Myeloid_avg@assays$RNA$scale.data %>% as.matrix()
gmtFile <- read.gmt(file.path(file.path(getOption("hla.data_root"), "references"), "KEGG_2021_Human.txt"))
genesets <- lapply(unique(gmtFile$term), function(x) gmtFile$gene[gmtFile$term == x])
names(genesets) <- unique(gmtFile$term)

gsva_es <- gsva(gsvaParam(gene_matrix, genesets, minSize = 1))
select_pathway <- c("ECM-receptor interaction", "Antigen processing and presentation", "Allograft rejection", 
                    "Endocytosis", "IL-17 signaling pathway", "Fc gamma R-mediated phagocytosis", "Complement and coagulation cascades")
plot_gsva_heatmap <- pheatmap(gsva_es[select_pathway, ], cluster_rows = TRUE, cluster_cols = TRUE, cellwidth = 20, cellheight = 20)
ggsave("./Figures/Myeloid_GSVA_heatmap.pdf", plot_gsva_heatmap, height = 10, width = 10)

saveRDS(Myeloid, "./Result/CSF_Myeloid.rds")

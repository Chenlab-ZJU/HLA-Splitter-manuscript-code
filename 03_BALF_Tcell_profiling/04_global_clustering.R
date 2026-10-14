# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: 03_global_clustering.R
# Figure/table scope: Fig3; FigS6; Supplementary Table 13
# Data cohort: 肺泡灌洗液. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# File: 02_global_clustering.R
# Description: 全局降维聚类、细胞群注释及 T 细胞的再次分群
# ==============================================================================
rm(list = ls())
source(file.path(getOption("hla.code_root"), "03_BALF_Tcell_profiling/01_setup_and_helpers.R"))

# Read the saved object or source table used by this analysis.
scRNA_Intergrate = readRDS("Result/lavage_scRNA_Intergrate.rds")
#Tcell = readRDS("./Result/lavage_Tcell.rds")

# 1. 数据筛选与重新聚类 --------------------------------------------------------
scRNA <- subset(scRNA_Intergrate, subset = !HLA_Splitter %in% c("Doublets", "Undefined"))
SampleInfo <- read.csv("./Data/SampleInfo.csv")
Sample_levels <- c(mixedsort(SampleInfo$Patient.ID))
scRNA$Group = factor(scRNA$Group, setdiff(scRNA$Group,c("Undefined","Doublet")))
# 过滤细胞数小于50的样本
valid_samples <- Sample_levels[table(scRNA$HLA_Splitter) > 50]
scRNA <- subset(scRNA, subset = HLA_Splitter %in% valid_samples)

# 1. 全局细胞聚类与注释 --------------------------------------------------------
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
Tcell = scRNA[,scRNA$Celltype != "Contament"]
Tcell$Celltype <- factor(Tcell$Celltype, levels = c("Naive T", "Quiescent T", "CD4+ Tact", "Th17-like", "Tph", "IFN+ T", "Tregs", 
                                                    "CD8+ Tem", "CD8+ Trm", "CD8+ Teff", "Cycling T", "MAIT cells"))

# 重新聚类
Tcell = Tcell[,Tcell$Group != "Others"]
Tcell$Group = factor(Tcell$Group, levels = c("CAP","LC","ILD"))
                       
Tcell <- RunUMAP(Tcell,dims = 1:30,  min.dist = 0.1)
#FeaturePlot(Tcell, features = c("CD4","CD8A","KLRD1","ISG15"))

# 3. 绘制全局UMAP三联图 --------------------------------------------------------
p1 <- DimPlot(Tcell, group.by = "seurat_clusters", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p2 <- DimPlot(Tcell, group.by = "Celltype", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p3 <- DimPlot(Tcell, group.by = "Group", label = TRUE, alpha = 0.5) + scale_color_manual(values = GROUP_COLORS)
ggsave("./Figures/Tcell_UMAP_celltype.pdf", p1 + p2 + p3, height = 5, width = 18)

# 全局 Marker 基因绘图
global_mk <- c("CCR7", "LEF1", "BRD1", "GRK5", "CD40LG", "HLA-DR", 
               "CCR6", "KLRB1", "CXCL13", "TNFRSF18", "IFIT3", "ISG15", "FOXP3", "IL2RA",
               "GZMK", "EOMES", "ZNF683","ITGAE",  "GZMH", "PRF1", 
               "PDCD1", "MKI67","TYMS", "SLC4A10", "TRAV1-2")
p_dot <- DotPlot(Tcell, group.by = "Celltype", features = global_mk) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x=NULL, y=NULL) +
  scale_color_gradientn(colours = c("navy","white","firebrick3"))
ggsave("./Figures/Tcell_Dotplot_clusters.pdf", plot = p_dot, height = 4, width = 9)

# UMAP 分面图
combined_plot <- ggarrange(plotlist = plot_grouped_umap(Tcell, "Group"), ncol = 3,common.legend = TRUE, legend = "right")

ggsave("./Figures/Tcell_all_celltype_facet.pdf", combined_plot, height = 5, width = 12)

# Celltype 热图
p_list <- plot_celltype_ratio(seurat_obj = Tcell,
                              sample_col = "HLA_Splitter", group_col = "Group", cluster_col = "Celltype",
                              comparison_groups = comparison, group_colors = GROUP_COLORS, cluster_colors = GLOBAL_COLORS)

ggsave("./Figures/Tcell_celltype_stacked_bar.pdf", p_list$stacked_bar, width = 5, height = 4)
ggsave("./Figures/Tcell_celltype_heatmap.pdf", p_list$heatmap, width = 5, height = 5)


# 细胞亚群比例图 -------------------------------------------------
comparison <- list(c("CAP","LC"), c("CAP","ILD"), c("ILD","LC"))
p_cell_ratio <- plot_cell_ratio_comparison(cell_subset = Tcell,
                                           sample_col = "HLA_Splitter", cluster_col =  "Celltype", group_col = "Group",
                                           facet = T, comparison_groups = comparison,
                                           write_to_csv = "./Result/Tcell_celltype_ratio", plot_colors = GLOBAL_COLORS)
ggsave("./Figures/Tcell_celltype_compare.pdf", p_cell_ratio, width = 10, height = 10)

#保存
# Export the analysis source data for figure reproduction.
saveRDS(Tcell, "./Result/lavage_Tcell.rds")

# ==============================================================================
# Description: 细胞组成比例对比、组间差异表达及富集分析
# ==============================================================================
Tcell = readRDS("./Result/lavage_Tcell.rds")
Tcell = Tcell[,Tcell$Group!= "Others"]

comparison <- list(c("CAP","LC"), c("CAP","ILD"), c("ILD","LC"))

# ==============================================================================
# 1. Pseudobulk分析，输出 Lollipop 图 + 分组专属 CSV + 分组专属 Bar&Dot 图
# ==============================================================================
source(file.path(getOption("hla.code_root"), "00_Shared/plot_top_pseudobulk.R"))

out_dir = "./Result/Marker_Analysis"

plot_top_pseudobulk(Tcell, 
                    group_by = "Group",
                    sample_col = "HLA_Splitter",
                    top_n_lollipop = 10, 
                    top_n_bardot = 10, 
                    out_dir = out_dir,
                    my_colors = GROUP_COLORS)

unique_groups = unique(Tcell$Group)

Group_DEG = data.frame()
for (grp in unique_groups) {
  top_gene_tmp = read.csv(paste0(out_dir,"/Marker_Stats_", grp, ".csv"))
  Group_DEG = rbind(Group_DEG, top_gene_tmp)
}
write.csv(Group_DEG, paste0(out_dir, "/Tcell_Group_DEG.csv"), row.names = FALSE)

# ==============================================================================
# 3. Heatmap (Pseudobulk) - 聚合所有分组的 Top 差异基因绘制全局热图
# ==============================================================================
message("Generating global pseudobulk heatmap...")

# 自动提取每个分组 log2FC 最高的基因
top_heatmap_genes <- Group_DEG %>%
  group_by(Group_Col) %>%
  dplyr::slice_max(order_by = PB_Log2FC, n = 5, with_ties = FALSE) %>%
  pull(gene) %>%
  unique()

Tcell_avg <- AverageExpression(Tcell, group.by = "Group", assays = "RNA", return.seurat = TRUE)
Tcell_avg <- ScaleData(Tcell_avg, features = rownames(Tcell_avg))
gene_matrix_t <- Tcell_avg@assays$RNA$scale.data %>% as.matrix()

# 取交集，防止某些基因在 Pseudobulk 的 scale.data 中丢失
genes_to_plot <- intersect(top_heatmap_genes, rownames(gene_matrix_t))

plot_heatmap_t_bulk <- pheatmap(gene_matrix_t[genes_to_plot, ], 
                                cluster_rows = TRUE, cluster_cols = FALSE, 
                                angle_col = 45, cellwidth = 40, cellheight = 12) # 动态调整 cellheight 防止图太长

ggsave("./Figures/Tcell_Heatmap_bulk_AllGroups.pdf", plot_heatmap_t_bulk, height = 12, width = 8)

# ==============================================================================
# 4. 循环遍历每个分组：绘制专属小提琴图和富集分析
# ==============================================================================
unique_groups <- as.character(unique(Group_DEG$Group_Col))

for (grp in unique_groups) {
  message("\n>> 正在执行亚群精细分析: ", grp)
  
  # 画小提琴图
  plot_top_violins(
    seu_obj = Tcell,
    deg_data = Group_DEG,
    target_group = grp,
    out_dir = out_dir,
    comparisons_list = comparison,
    plot_colors = GROUP_COLORS,
    group_col = "Group_Col",
    lfc_col = "PB_Log2FC",
    lfc_threshold = 0.1,
    top_n = 8
  )
  
  # 跑 GO 富集并画图
  run_go_enrichment(
    deg_data = Group_DEG,
    target_group = grp,
    out_dir = out_dir,
    group_col = "Group_Col",
    lfc_col = "PB_Log2FC",
    lfc_threshold = 0.1,
    top_n_pathways = 8
  )
}

# ==============================================================================
# 5. 综合性pesudobulk差异分析
# ==============================================================================
comparison <- list(c("CAP","LC"), c("CAP","ILD"), c("ILD","LC"))

source(file.path(getOption("hla.code_root"), "03_BALF_Tcell_profiling/07_run_sc_pseudobulk_functional_analysis.R"))

#Tcell$Major_type = "Tcell"

res <- run_sc_pseudobulk_functional_analysis(
  seu = Tcell,
  group_col = "Group",
  sample_col = "HLA_Splitter",
  celltype_col = "Major_type",
  default_celltype = "Tcell",
  comparisons = comparison,
  organism = "human",
  gene_id_type = "SYMBOL",
  assay = "RNA",
  layer = "counts",
  min_cells_per_sample = 20,
  min_samples_per_group = 2,
  min_count = 20,
  min_samples_expr = 3,
  padj_cutoff = 0.05,
  lfc_cutoff = 0.25,
  output_dir = "./Result/Tcell_pseudobulk_DEG_GSEA"
)

# 1. 画基于 pvalue 的热图
heatmap_plot <- plot_lrt_global_heatmap(
  lrt_table = res$overall_LRT,
  pb_counts = res$pseudobulk$counts,
  pb_meta = res$pseudobulk$metadata,
  target_celltype = "Tcell",
  filter_by = "pvalue",    # 明确指定使用原始 pvalue
  p_cutoff = 0.01,         # 将阈值放宽到 0.05
  k_clusters = 4
)

# 2. 画基于 pvalue 的趋势聚类图
trend_res <- plot_lrt_cluster_trends(
  lrt_table = res$overall_LRT,
  pb_counts = res$pseudobulk$counts,
  pb_meta = res$pseudobulk$metadata,
  target_celltype = "Tcell",
  filter_by = "pvalue",    # 明确指定使用原始 pvalue
  p_cutoff = 0.01,         # 将阈值放宽到 0.05
  k_clusters = 4
)

# 查看图形
print(trend_res$plot)

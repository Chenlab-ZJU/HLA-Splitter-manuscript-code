# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Upstream HLA typing, resolution and PBMC demultiplexing.
# Original source: 30PBMC.R
# Figure/table scope: FigS1; FigS2f; Supplementary Tables 1 and 9
# Data cohort: 30PBMC. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

setwd(file.path(getOption("hla.data_root"), "30PBMC"))
rm(list = ls())

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(reshape2)
library(R.utils)
library(pheatmap)
library(ggridges)
library(stringr)
library(ggsci)
library(gtools)

dir.create("Figures", showWarnings = FALSE)
dir.create("Result", showWarnings = FALSE)
dir.create("Result/demultiplex", showWarnings = FALSE)

# 定义全局颜色盘 (使用 ggsci 的 igv 配色，样本多时比较适用)
color <- pal_igv("default")(51)

filename = "30PBMC"

# ==========================================
# 1. 读取数据与基础质控 (QC)
# ==========================================
Counts <- Read10X(data.dir = paste0("./Data/GEM/filtered_feature_bc_matrix/"))
scRNA_all <- CreateSeuratObject(Counts, project = filename, assay = "RNA")

scRNA_all[["percent.mt"]] <- PercentageFeatureSet(scRNA_all, pattern = "^MT-")

p1 = VlnPlot(scRNA_all, features=c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=2, pt.size=0)
p2 = FeatureScatter(scRNA_all, feature1="nCount_RNA", feature2="nFeature_RNA")
ggsave(paste0("./Figures/QC.pdf"), p1+p2, dpi=300, height=10, width = 12)

# 过滤低质量细胞
scRNA_all <- subset(scRNA_all, subset = nFeature_RNA > 200 & nFeature_RNA < 100000 & percent.mt < 20)

# ==========================================
# 2. 聚类与去双细胞 (Scrublet)
# ==========================================
scRNA_all <- NormalizeData(scRNA_all) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA()
scRNA_all <- RunTSNE(scRNA_all, dims = 1:30)
scRNA_all <- RunUMAP(scRNA_all, dims = 1:30)
scRNA_all <- FindNeighbors(scRNA_all, dims = 1:30) %>% FindClusters(resolution=0.8)

# install.packages("remotes")
# remotes::install_github("r-lib/devtools")
# devtools::install_github("furlan-lab/scrubletR")
library(scrubletR)
scRNA_all <- scrublet(scRNA_all)
p1 = FeaturePlot(scRNA_all, features = "doublet_scores", reduction = "tsne")
ggsave("./Figures/scrublet_doublet_score.pdf", p1, height=5, width = 6)

scRNA_all$doublets <- scRNA_all$doublet_scores > 0.2   
p2 = DimPlot(scRNA_all, group.by = "doublets", cols=c("goldenrod", "navy"), reduction = "tsne")
ggsave("./Figures/scrublet_doublet.pdf", p2, height=5, width = 6)

# ==========================================
# 3. 提取单细胞并重新聚类
# ==========================================
scRNA = scRNA_all[, !scRNA_all$doublets]
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA()
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA <- FindNeighbors(scRNA, dims = 1:30) %>% FindClusters(resolution=0.3)

# ==========================================
# 4. 细胞亚群差异基因分析与导出
# ==========================================
scRNA <- SetIdent(scRNA, value = "seurat_clusters")
# Discover marker genes with the original assay and comparison settings.
cluster.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Export the analysis source data for figure reproduction.
write.csv(cluster.markers, "./Result/scRNA_markers.csv", row.names = FALSE)

df_summary <- data.frame(seurat_clusters = levels(scRNA$seurat_clusters))
# 提取每个 cluster 的 Top 20 基因
df_top_genes <- cluster.markers %>%
  filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>% 
  summarise(top_genes = paste0(gene, collapse = ", ")) %>% 
  ungroup() %>% 
  rename(seurat_clusters = cluster)

final_df <- left_join(df_summary, df_top_genes, by = "seurat_clusters")
final_df$Celltype = c("")
write.csv(final_df, "./Result/scRNA_celltype.csv", row.names = FALSE)

message("请手动打开 ./Result/scRNA_markers_summary.csv 文件，添加 Celltype 列后保存！")

# ==========================================
# 5. 基于高表达基因对细胞亚群进行修正与重命名
# ==========================================
# Read the saved object or source table used by this analysis.
df = read.csv("./Result/scRNA_celltype.csv")

# 构建映射字典 (Named Vector)
Celltype = df$Celltype
names(Celltype) = as.character(df$seurat_clusters)

# 将数字 cluster 替换为细胞类型
scRNA <- SetIdent(scRNA, value = "seurat_clusters")
scRNA <- RenameIdents(scRNA, Celltype)
scRNA$Celltype <- Idents(scRNA)

# # 规范化因子的顺序，方便画图时图例按这个顺序排列
# scRNA$Celltype <- factor(scRNA$Celltype, levels = c("CD4+ T", "CD8+ T", "B cell", "NK",
#                                           "cDC", "pDC", "Monocyte", "Macrophage",
#                                           "Neutrophil", "Cycling cell", "Low-quality"))

# ==========================================
# 6. 绘制 UMAP 拼图
# ==========================================
my_colors <- pal_igv("default")(length(levels(scRNA$Celltype)))

# 图1：原始数字 Cluster
scRNA <- SetIdent(scRNA, value = "seurat_clusters")
p1 <- DimPlot(scRNA, reduction = "umap", label = TRUE, alpha = 0.5, label.size = 3) +
  theme(legend.position = "right", legend.box = "vertical", legend.key.size = unit(3, "pt")) +
  scale_color_discrete() 

# 图2：注释后的 Celltype
scRNA <- SetIdent(scRNA, value = "Celltype")
p2 <- DimPlot(scRNA, reduction = "umap", label = TRUE, alpha = 0.5, label.size = 3) +
  theme(legend.position = "right", legend.box = "vertical", legend.key.size = unit(3, "pt")) +
  scale_color_manual(values = my_colors) # 使用我们定义的调色板

# 拼图并保存
combined_plot <- p1 + p2
ggsave("./Figures/scRNA_UMAP_celltype.pdf", combined_plot, height = 5, width = 12)


# ==========================================
# 7. 基于 HLA-Splitter 的样本拆分
# ==========================================
# 读取HLA counts信息
sampleInfo = read.csv("./Data/HLA_list/HLA-genotypes_ONT.csv", header=TRUE, row.names = 1)
namelist = colnames(sampleInfo)

# 构建 HLA list
HLA_list = list()
for (i in 1:length(namelist)){ HLA_list[[i]] = unique(sampleInfo[,i]) }
names(HLA_list) = namelist

# 读取矩阵
HLA_label = read.delim("./Data/HLA_counts/sorted_genotypes.txt", header = F)
barcodes = read.delim(paste0("./Data/HLA_counts/output.barcodes.txt"), header=F)
HLA_matrix = readMM(paste0("./Data/HLA_counts/HLA_type.mtx")) %>% as.matrix() %>% as.data.frame()
colnames(HLA_matrix) = HLA_label$V1
rownames(HLA_matrix) = paste0(barcodes$V1, "-1")

# 加载函数并进行拆分
source(file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks/01_HLA_demultiplex_v2.R"))
sample_predict <- HLA_demultiplex(HLA_list = HLA_list, 
                                  HLA_matrix = HLA_matrix, 
                                  graph_corr = TRUE, 
                                  outdir = "Result/demultiplex")

# ==========================================
# 8. 将拆分结果映射回 Seurat 对象并绘图
# ==========================================
scRNA$HLA_splitter = sample_predict[match(names(sample_predict), colnames(scRNA))] # 根据细胞barcode精准匹配
scRNA$HLA_splitter[is.na(scRNA$HLA_splitter)] = "Undefined"

# 绘制 UMAP: Cluster vs Celltype
scRNA = SetIdent(scRNA, value = "seurat_clusters")
p1 = DimPlot(scRNA, reduction = "umap", label=T) +
  scale_color_discrete(type = c(color)) +
  theme(legend.position = "right", legend.box = "vertical", legend.key.size = unit(3,"pt"))

scRNA = SetIdent(scRNA, value = "HLA_splitter")
p2 = DimPlot(scRNA, reduction = "umap", label=F) +
  scale_color_discrete(type = c(color)) +
  theme(legend.position = "right", legend.box = "vertical", legend.key.size = unit(3,"pt"))

# 拼接并保存
ggsave(paste0("./Figures/UMAP_Cluster_Splitter.pdf"), p1 + p2, height = 5, width = 12)

# ==========================================
# 9. 细胞比例直方图
# ==========================================
df_ratio = data.frame(Sample = scRNA$HLA_splitter, Celltype = scRNA$Celltype) %>% 
  table() %>% prop.table(margin = 1) %>% data.frame() # 使用 margin = 1 按样本计算比例
df_ratio$Freq = df_ratio$Freq * 100

p3 = ggplot(df_ratio, aes(x = Sample, y = Freq, fill = Celltype)) +
  geom_col(position = "fill", alpha=0.7) +
  scale_fill_manual(values = color) +
  theme_classic() +
  labs(y = "Frequency of cell types (%)", x = "") +
  theme(
    text = element_text(size = 15),
    axis.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
  )
ggsave(paste0("./Figures/Celltype_ratio.pdf"), p3, height = 5, width = 10)

# ==========================================
# 10. 保存文件
# ==========================================
saveRDS(scRNA_all, "./Code/30PBMC_all.rds")
saveRDS(scRNA, "./Code/30PBMC.rds")

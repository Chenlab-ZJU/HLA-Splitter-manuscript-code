# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Upstream synthetic-data preparation and benchmark computation.
# Original source: 00_scRNA_analysis.R
# Figure/table scope: Fig1b; Fig1c; Fig1d; FigS2d; FigS2e
# Data cohort: synthetic. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# 0. 环境设置与依赖加载
# ==============================================================================
setwd(file.path(getOption("hla.data_root"), "synthetic"))
rm(list = ls())

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(harmony)
library(reshape2)
library(pheatmap)
library(sctransform)
library(RColorBrewer)
library(grid)
library(stringr)
library(ggsci)
library(scales)
library(scrubletR)

# 自动创建所需的所有输出目录
dir.create("Figures/QC", recursive = TRUE, showWarnings = FALSE)
dir.create("Result", showWarnings = FALSE)
dir.create("demultiplex", showWarnings = FALSE)

# ==============================================================================
# 1. 数据读取、去双细胞与基础质控 (Per-sample)
# ==============================================================================
filename <- "Syn_data"
# Read the saved object or source table used by this analysis.
sampleInfo <- read.csv("./Data/HLA_list/HLA-genotypes-2field.csv", header=TRUE, row.names = 1)
namelist <- colnames(sampleInfo)

scRNAlist <- list()

for (i in seq_along(namelist)) {
  sample_name <- namelist[i]
  counts <- Read10X(data.dir = paste0("./Data/GEX/", sample_name, "/filtered_feature_bc_matrix"))
  obj <- CreateSeuratObject(counts, project = sample_name)
  obj <- RenameCells(obj, add.cell.id = sample_name)
  # 计算线粒体比例
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  # 基础过滤：去除低基因表达、高线粒体、以及双细胞
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20 & doublets == FALSE)
  scRNAlist[[sample_name]] <- obj
}

#去除不同样本中cellbarcode相同的cell
barcodes_all = c()
for (i in 1:length(namelist)){
  barcodes_all = c(barcodes_all,colnames(scRNAlist[[i]]))
}
barcodes_duplicate = barcodes_all[duplicated(barcodes_all)] %>% unique()
for (i in 1:length(namelist)){
  scRNAlist[[i]]=scRNAlist[[i]][,!colnames(scRNAlist[[i]]) %in% barcodes_duplicate]
}

# 保存初步质控后的单样本列表
# Export the analysis source data for figure reproduction.
saveRDS(scRNAlist, "./Result/scRNAlist.RDS")
#scRNAlist <- readRDS("./Result/scRNAlist.RDS")

# ==============================================================================
# 2. 数据合并与 QC 可视化
# ==============================================================================
scRNA_merge <- merge(scRNAlist[[1]], y = scRNAlist[-1], project = "Synthetic_Data")

# Seurat v5 推荐在整合或标准化前执行 JoinLayers
scRNA_merge <- JoinLayers(scRNA_merge)
meta <- scRNA_merge@meta.data
color <- pal_igv("default")(length(namelist) + 1)

# 准备细胞数量统计框
CellNum <- as.data.frame(table(meta$orig.ident))
colnames(CellNum) <- c("Sample", "nCount_Cell")

# 画图：细胞数
p1 <- ggplot(CellNum, aes(x = Sample, y = nCount_Cell)) +
  geom_col(aes(fill = Sample), alpha=0.7) +
  geom_text(aes(label = nCount_Cell), position = position_stack(), vjust = 1.5) +
  theme_classic() +
  scale_fill_discrete(type = c(color)) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)) # 样本名太长时倾斜

# 画图：UMI 数 (nCount_RNA)
p2 <- ggplot(meta, aes(x = orig.ident, y = nCount_RNA)) +
  geom_violin(aes(fill = orig.ident), alpha=0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_classic() +
  ylim(c(0, 10000)) +
  scale_fill_discrete(type = c(color)) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# 画图：基因数 (nFeature_RNA)
p3 <- ggplot(meta, aes(x = orig.ident, y = nFeature_RNA)) +
  geom_violin(aes(fill = orig.ident), alpha=0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_classic() +
  ylim(c(0, 5000)) +
  scale_fill_discrete(type = c(color)) +
  theme(legend.position = 'none',
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./Figures/QC/CellNum.pdf", p1, height=10, width=20)
ggsave("./Figures/QC/UMI.pdf", p2, height=8, width=20)
ggsave("./Figures/QC/Gene.pdf", p3, height=8, width=20)


# ==============================================================================
# 3. 标准化、降维与批次矫正 (Harmony)
# ==============================================================================
scRNA_merge <- NormalizeData(scRNA_merge) %>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE)

# 使用 Harmony 进行批次矫正
scRNA_merge <- RunHarmony(scRNA_merge, group.by.vars = "orig.ident")

# 去除Doublets
scRNA_merge <- scrublet(scRNA_merge)
scRNA_merge$doublets <- scRNA_merge$doublet_scores > 0.2

# ==============================================================================
# 4. 聚类与 UMAP
# ==============================================================================
scRNA_merge <- RunUMAP(scRNA_merge, reduction = "harmony", dims = 1:30, min.dist = 0.3)
scRNA_merge <- FindNeighbors(scRNA_merge, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)

# ==============================================================================
# 5. 全局 Marker 提取与可视化
# ==============================================================================
#remotes::install_github('satijalab/seurat-wrappers')
export_top_markers <- function(seurat_obj, ident_col, out_prefix) {
  Idents(seurat_obj) <- ident_col
  markers <- SeuratWrappers::RunPrestoAll(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers, paste0("./Result/", out_prefix, "_makers.csv"), row.names = FALSE)
  
  df_summary <- seurat_obj@meta.data %>% select(all_of(ident_col)) %>% distinct() %>% arrange(across(all_of(ident_col)))
  df_top_genes <- markers %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>%
    summarise(top_genes = paste0(gene, collapse = ", ")) %>% rename(!!ident_col := cluster)
  
  final_df <- left_join(df_summary, df_top_genes, by = ident_col)
  write.csv(final_df, paste0("./Result/", out_prefix, "_makers_summary.csv"), row.names = FALSE)
}

export_top_markers(scRNA_merge, "seurat_clusters", "scRNA")

# 注意：运行到此处请暂停，确保 ./Result/scRNA_makers_summary.csv 已填写 Celltype
df <- read.csv("./Result/scRNA_makers_summary.csv")
Celltype <- df$Celltype
names(Celltype) <- as.character(df$seurat_clusters)
scRNA_merge <- SetIdent(scRNA_merge, value = "seurat_clusters")
scRNA_merge <- RenameIdents(scRNA_merge, Celltype)
scRNA_merge$Celltype <- Idents(scRNA_merge)
unique(scRNA_merge$Celltype) %>% print()

# 为了演示，假定 Celltype 已经映射好
scRNA_merge = scRNA_merge[,scRNA_merge$Celltype != "Erythrocytes"]
scRNA_merge$Celltype <- factor(scRNA_merge$Celltype, levels = c( "T cells", "B cells", "NK cells", "Monocytes", "MAIT", "DCs", "Platelets"))


# 3. 绘制全局UMAP三联图 --------------------------------------------------------
p1 <- DimPlot(scRNA_merge, group.by = "seurat_clusters", label = TRUE, alpha = 0.5) + scale_color_manual(values = color )
p2 <- DimPlot(scRNA_merge, group.by = "Celltype", label = TRUE, alpha = 0.5) + scale_color_manual(values = color )
ggsave("./Figures/scRNA_merge_UMAP_celltype.pdf", p1 + p2, height = 5, width = 12)

# 全局 Marker 基因绘图
global_mk <- c("CD3E", "CD3D", "MS4A1","BANK1",
               "NKG7", "KLRF1", "CD14", "S100A8",
               "SLC4A10", "TRAV1-2", "CD1C", "FCER1A",
               "PPBP", "GP9")
p_dot <- DotPlot(scRNA_merge, group.by = "Celltype", features = global_mk) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x=NULL, y=NULL) +
  scale_color_gradientn(colours = c("navy","white","firebrick3"))
ggsave("./Figures/scRNA_merge_Dotplot_clusters.pdf", plot = p_dot, height = 4, width = 8)

df_celltype = scRNA_merge@meta.data
write.csv(df_celltype, "./Result/scRNA_celltype.csv")
# ==============================================================================
# 6. 保存最终 Seurat 对象
# ==============================================================================
# 统一保存到 Result 文件夹
saveRDS(scRNA_merge, "./Result/scRNA_merge.RDS")
#scRNA_merge = readRDS("./Result/scRNA_merge.RDS")
print("Pipeline execution completed successfully!")

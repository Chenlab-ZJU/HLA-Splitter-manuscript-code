# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Experimental sorted-cell mixture validation.
# Original source: 4mix_PBMC.R
# Figure/table scope: FigS3a; FigS3b; FigS3c; FigS3d
# Data cohort: pbmc4. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

setwd(file.path(getOption("hla.data_root"), "pbmc4"))
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
library(ggsci)

#save.image(file = "./4mix_PBMC.RData")
#load(file = "./4mix_PBMC.RData")

# 定义全局颜色盘 (使用 ggsci 的 igv 配色，样本多时比较适用)
color <- pal_igv("default")(51)

filename = "4mix_PBMC"

## 构建seurat对象 ####
Counts <- Read10X(data.dir = paste0("./GEM/filtered_cell_gene_matrix/"))
colnames(Counts)=gsub('[_]','',colnames(Counts))
scRNA_all <- CreateSeuratObject(Counts,project = filename)
scRNA_all[["percent.mt"]] <- PercentageFeatureSet(scRNA_all, pattern = "^MT-")

## 质检数据 ####
p1 = VlnPlot(scRNA_all, features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3,pt.size=0)
p2 = FeatureScatter(scRNA_all, feature1="nCount_RNA", feature2="nFeature_RNA")
scRNA_all <- subset(scRNA_all, subset= nFeature_RNA >200 & nFeature_RNA < 30000 & percent.mt <20)
ggsave(paste0("./QC/VlnPlot.pdf"),p1,dpi=300)
ggsave(paste0("./QC/FeatureScatter.pdf"),p2,dpi=300)

## Normalize Data ####
scRNA_all <- NormalizeData(scRNA_all) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA()
ElbowPlot(scRNA_all,ndims = 50, reduction = "pca")
VizDimLoadings(scRNA_all, dims = 1:2, reduction = "pca")
scRNA_all <- RunUMAP(scRNA_all, dims = 1:30)
scRNA_all <- RunTSNE(scRNA_all, dims = 1:30)
scRNA_all <- FindNeighbors(scRNA_all, dims = 1:30) %>% FindClusters(resolution=0.3)

## remove doublets #### 
#devtools::install_github("furlan-lab/scrubletR")
library(scrubletR)

scRNA_all <- scrublet(scRNA_all)
p1 = FeaturePlot(scRNA_all, features = "doublet_scores",reduction = "tsne")
ggsave("./demultiplex/scrublet_doublet_score.pdf",p1,height=5,width = 6)

hist(scRNA_all$doublet_scores)
scRNA_all$doublets<-scRNA_all$doublet_scores > 0.1   #(You pick this)
p2 = DimPlot(scRNA_all, group.by = "doublets", cols=c("goldenrod", "navy"),reduction = "tsne")
ggsave("./demultiplex/scrublet_doublet.pdf",p2,height=5,width = 6)

## TSNE plot ####
scRNA = scRNA_all[,!scRNA_all$doublets]
scRNA = scRNA[,!scRNA$seurat_clusters %in% c("4")] #去除不确定的细胞类群
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA()
scRNA <- RunTSNE(scRNA, dims = 1:30)
scRNA <- FindNeighbors(scRNA, dims = 1:30) %>% FindClusters(resolution=0.3)

scRNA = SetIdent(scRNA,value = "seurat_clusters")
p1 = DimPlot(scRNA, reduction = "tsne", label=T)+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(3,"pt"))
scRNA = SetIdent(scRNA,value = "celltype")
p2 = DimPlot(scRNA, reduction = "tsne", label=T)+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(3,"pt"))
ggsave(paste0("./Figures/TSNE_cluster.pdf"),p1+p2,dpi=300,height = 5,width = 12)

## 寻找差异基因 ####
scRNA = SetIdent(scRNA,value = "seurat_clusters")
# Discover marker genes with the original assay and comparison settings.
scRNA.markers <- FindAllMarkers(scRNA, only.pos= TRUE, min.pct= 0.25,
                                logfc.threshold= 0.25)
# Export the analysis source data for figure reproduction.
write.csv(scRNA.markers,paste0("./Rawdata/Marker.csv"))

## heatmap ####
scRNA.markers %>%
  group_by(cluster) %>%
  top_n(n= 10, wt=avg_log2FC) ->top10
subobj <- ScaleData(object = scRNA, features = top10$gene)
subobj <- subset(subobj, downsample = 300)
set.seed(42)
options(repr.plot.width = 13, repr.plot.height=13)
plot_heatmap = DoHeatmap(subobj, features = top10$gene,
                         disp.min=-2.5, disp.max=2.5)+
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))

ggsave(paste0("./Figures/heatmap.pdf"),plot_heatmap,height = 20,width = 15)

## 细胞亚群差异基因绘图 ####
classmk <-c("CD3D","KLRF1","LYZ","MS4A1")
scRNA = SetIdent(scRNA,value = "seurat_clusters")
p1 = FeaturePlot(scRNA, features="CD3D",reduction = "tsne",label = F)+
  scale_color_continuous(type = "viridis")
p2 = FeaturePlot(scRNA, features="KLRF1",reduction = "tsne",label = F)+
  scale_color_continuous(type = "viridis")
p3 = FeaturePlot(scRNA, features="LYZ",reduction = "tsne",label = F)+
  scale_color_continuous(type = "viridis")
p4 = FeaturePlot(scRNA, features="MS4A1",reduction = "tsne",label = F)+
  scale_color_continuous(type = "viridis")
p = ggarrange(p1,p2,p3,p4,common.legend = T)
ggsave(paste0("./Figures/marker_clusters.pdf"),p,height = 10,width = 9)

plot_violin = VlnPlot(scRNA, features=classmk,
                      pt.size=0,stack = T,same.y.lims = T, flip = T)+NoLegend()
ggsave(paste0("./Figures/violin_clusters.pdf"),plot=plot_violin,height = 6,width = 5)

## 细胞亚群分群 ####
scRNA = SetIdent(scRNA,value = "seurat_clusters")
celltype_major <- c("Monocyte","NK","B_cell","T_cell","Monocyte","NK",
                    "T_cell","T_cell","B_cell","T_cell","Monocyte","Monocyte","Monocyte")

names(celltype_major) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, celltype_major)
scRNA$celltype_major = Idents(scRNA)
scRNA$celltype_major = factor(scRNA$celltype_major,levels=c("T_cell","NK","Monocyte","B_cell","Doublet"))

scRNA = SetIdent(scRNA,value = "seurat_clusters")
p1 = DimPlot(scRNA, reduction = "tsne", label=T)+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(3,"pt"))
scRNA = SetIdent(scRNA,value = "celltype_major")
p2 = DimPlot(scRNA, reduction = "tsne", label=T)+
  scale_color_discrete(type = c(color,"grey","black"))+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(3,"pt"))
scRNA = SetIdent(scRNA,value = "celltype")
p3 = DimPlot(scRNA, reduction = "tsne", label=T)+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(3,"pt"))
ggsave(paste0("./Figures/TSNE_cluster.pdf"),p1+p2+p3 ,height=5,width=18,dpi=300)

#分群后violin_plot
scRNA = SetIdent(scRNA,value = "celltype_major")
plot_violin = VlnPlot(scRNA, features=classmk,
                      pt.size=0,stack = T,same.y.lims = T, flip = T)+NoLegend()+
  scale_fill_discrete(type = c(color,"grey","black"))+
  theme(
    legend.position = "none", # Remove the legend
    axis.title.y = element_text(size = 15, face = "bold"), # X-axis title (expression)
    axis.title.x = element_blank(), # No Y-axis title if feature names are on the y-axis
    axis.text.x = element_text(size = 12), # Adjust x-axis tick label size
    axis.text.y = element_text(size = 12), # Adjust y-axis tick label size
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5) # Center and bold plot title
  ) +
  labs(y = "Normalized Hash Tag expression")
ggsave(paste0("./Figures/violin_clusters.pdf"),plot=plot_violin,height = 5,width = 5)

#保存数据
saveRDS(scRNA,"./Rawdata/4mix_PBMC.rds")

#基于HLA-Splitter样本拆分######
# 读取HLAcounts信息
# Read the saved object or source table used by this analysis.
sampleInfo = read.csv(file.path(getOption("hla.data_root"), "pbmc4", "HLA-genotypes.csv"),header=TRUE,row.names = 1)
namelist=colnames(sampleInfo)
length(namelist)
# HLA list
HLA_list=list()
for (i in 1:length(namelist)){HLA_list[[i]] = unique(sampleInfo[,i])}
names(HLA_list) = namelist

HLA_label = read.delim("./Data/HLA_counts/sorted_genotypes.txt",header = F)
barcodes=read.delim(paste0("./Data/HLA_counts/output.barcodes.txt"),header=F)
HLA_matrix = readMM(paste0("./Data/HLA_counts/HLA_type.mtx"))  %>% as.matrix() %>% as.data.frame()
colnames(HLA_matrix) = HLA_label$V1
rownames(HLA_matrix) = paste0(barcodes$V1)

#细胞样本来源标记
outdir = "./Result/demultiplex/"
scRNA = readRDS(file.path(getOption("hla.data_root"), "pbmc4", "4mix_PBMC.rds"))
sample_label = scRNA$celltype_major
barcodes = intersect(barcodes$V1,names(sample_label))
sample_label = sample_label[barcodes]
HLA_matrix = HLA_matrix[barcodes,]

# Sample demultiplexing
source(file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks/01_HLA_demultiplex_v2.R"))
sample_predict <- HLA_demultiplex(HLA_list = HLA_list, 
                                  HLA_matrix = HLA_matrix, 
                                  graph_corr = TRUE, 
                                  outdir = outdir)
table(sample_predict)

#绘制TSNE图####
colnames(tsne_result) <- gsub("adjust_predict", "Predict_label", colnames(tsne_result))
tsne_result$True_label <- sample_label

p1 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = True_label)) +
  geom_point(size = 0.2, alpha = 0.8) + # 稍微增加透明度和大小，质感更好
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 4))) +  
  scale_color_manual(values = c(color[1:length(namelist)]))+
  labs(title = "True Labels", color = "Sample")

p2 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = Predict_label)) +
  geom_point(size = 0.2, alpha = 0.8) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_color_manual(values = c(color[1:length(namelist)], "grey"))+
  labs(title = "Predicted Labels", color = "Sample")

combined_plot <- p1 + p2
print(combined_plot)

ggsave(paste0(outdir, "/Combined_TSNE_plot.pdf"), plot = combined_plot, height = 5, width = 11)

#比较两种方法的一致程度####
scRNA$HLA_Splitter = sample_predict
scRNA$HLA_Splitter[is.na(scRNA$HLA_Splitter)] = "Undefined"
df = data.frame(Celltype = scRNA$celltype_major,HLA_demutliplex = scRNA$HLA_Splitter)
df = table(df) %>% prop.table(1) %>% as.data.frame()
df = df[!is.na(df$Freq),]

p1 <- ggplot(df, aes(x = Celltype, y = Freq, fill = HLA_demutliplex)) +
  geom_col(alpha = 0.8, position = "stack", width = 0.8) + # Adjusted alpha, added width for bar appearance
  geom_text(aes(label = paste0(round(Freq * 100, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "black", # Set text color to black for better contrast
            size = 4) + 
  ylab("Cell Frequency (%)") + xlab("") + 
  theme_classic() +
  theme(
    text = element_text(size = 14), # Base font size
    axis.title.y = element_text(size = 15, face = "bold"), # Y-axis title bold and larger
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), # Rotate x-axis labels for readability
    axis.text.y = element_text(size = 12), # Y-axis tick label size
    legend.title = element_text(size = 14, face = "bold"), # Legend title bold
    legend.text = element_text(size = 12), # Legend text size
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5) # Optional: Add and center a plot title
  ) +
  scale_fill_manual(values = c(color[1:length(namelist)], "grey"))
p1
ggsave(paste0(outdir,"/label_consistency.pdf"),p1,height=6,width = 6)

df = data.frame(Celltype = scRNA$celltype_major,HLA_demutliplex = scRNA$HLA_Splitter)
df = table(df) %>% prop.table(1)*100
write.csv(df,paste0(outdir,"/demultiplex_res.csv"))

#绘制score热图######
score_df <- data.frame(score,
                       sample_predict = sample_predict,
                       sample_label = sample_label)

score_df <- score_df[order(score_df$sample_predict), ]

mat <- score_df %>%
  select(all_of(paste0(namelist,".score"))) %>% # Use select(all_of(namelist)) for safer column selection
  t() # Transpose to get features (namelist) as rows and samples as columns

annotation_col_df <- data.frame(
  sample_label = score_df$sample_predict # Using 'sample_predict' for annotation
)
rownames(annotation_col_df) <- colnames(mat) # Assign sample IDs as row names

desired_levels <- c(namelist, "Undefined")
annotation_col_df$sample_label <- factor(annotation_col_df$sample_label, levels = desired_levels)

ann_colors_list <- list(
  sample_label = c(color[1:length(namelist)], "grey")
)
names(ann_colors_list$sample_label) <- desired_levels # Assign names to colors to map to factor levels

p1 <- pheatmap(
  mat,
  cluster_rows = FALSE,      # Don't cluster rows (features/namelist items)
  cluster_cols = FALSE,      # Don't cluster columns (samples), as they are pre-ordered
  scale = "column",          # Scale values within each column (sample)
  rect_gp = gpar(col = 'white', lwd = 0.5), # Slightly thinner white lines between cells
  row_names_side = "left",   # Place row names on the left
  show_colnames = FALSE,     # Don't show column names (sample IDs) at the top
  annotation_col = annotation_col_df, # Add the column annotation
  annotation_colors = ann_colors_list, # Apply custom colors for annotation
  main = "Score Heatmap by Predicted Sample" 
)
ggsave(paste0(outdir,"/score_heatmap.pdf"),p1,dpi=300,height=5,width =7)

## 细胞亚群HLA-score绘图 ####
library(ggpubr)
score = apply(score,1,scale) %>% t()
colnames(score) = c(paste0(namelist,".score"))
score[is.na(score)]=0
tsne_result = cbind(tsne_result,score)

p1 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=Donor1.score))+
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")

p2 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=Donor2.score)) +
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")

p3 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=Donor3.score)) +
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")

p4 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=Donor4.score)) +
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")
p = ggarrange(p1,p2,p3,p4,common.legend = T)
ggsave(paste0(outdir,"/HLA-score_clusters.pdf"),p,height = 10,width = 10)

## 细胞亚群特征基因绘图 ####
library(ggpubr)
#打分正则化
marker = FetchData(scRNA, c("CD3D","KLRF1","LYZ","MS4A1"))
marker = marker[rownames(tsne_result),]
tsne_result = cbind(tsne_result,marker)

p1 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=CD3D))+
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")

p2 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=KLRF1)) +
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")

p3 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=LYZ)) +
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")

p4 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=MS4A1)) +
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")
p = ggarrange(p1,p2,p3,p4,common.legend = T)
ggsave(paste0(outdir,"/Genemarker_clusters.pdf"),p,height = 10,width = 10)


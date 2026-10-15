# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Experimental cell-hashing validation.
# Original source: 10X-hash.R
# Figure/table scope: FigS3e; FigS3f; FigS3g; FigS3h
# Data cohort: hash10x. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

setwd(file.path(getOption("hla.data_root"), "hash10x"))
rm(list = ls())
#remotes::install_github("satijalab/seurat", "seurat5", quiet = T)
#remotes::install_github("mojaveazure/seurat-object", "seurat5")

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
library(ggpubr)

# 定义全局颜色盘 (使用 ggsci 的 igv 配色，样本多时比较适用)
color <- pal_igv("default")(51)

filename = "10X-hash"
 
#读取umi和hash数据
Counts <- Read10X(data.dir = paste0("./Data/GEM/filtered_feature_bc_matrix/"))
scRNA_all.umis <- Counts$`Gene Expression`
scRNA_all.htos <- Counts$`Antibody Capture`
joint.bcs <- intersect(colnames(scRNA_all.umis), colnames(scRNA_all.htos))
scRNA_all.umis <- scRNA_all.umis[, joint.bcs]
scRNA_all.htos <- as.matrix(scRNA_all.htos[, joint.bcs])
scRNA_all <- CreateSeuratObject(Counts$`Gene Expression`,project = filename,assay = "RNA")

#add HTO assay
scRNA_all[["HTO"]] <- CreateAssayObject(scRNA_all.htos)


#计算线粒体比例
scRNA_all[["percent.mt"]] <- PercentageFeatureSet(scRNA_all, pattern = "^MT-")

#绘制质控图
p1 = VlnPlot(scRNA_all, features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=2,pt.size=0)
p2 = FeatureScatter(scRNA_all, feature1="nCount_RNA", feature2="nFeature_RNA")
ggsave(paste0("./Figures/QC.pdf"),p1+p2,dpi=300,height=10,width = 12)

#过滤低质量细胞
scRNA_all <- subset(scRNA_all, subset = nFeature_RNA > 200 & nFeature_RNA < 30000 & percent.mt <20)

#normalize data
scRNA_all <- NormalizeData(scRNA_all) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA()
scRNA_all <- NormalizeData(scRNA_all, assay = "HTO", normalization.method = "CLR")

#Demultiplex cells based on HTO enrichment
scRNA_all <- MULTIseqDemux(scRNA_all, assay = "HTO")

#Visualize demultiplexing results
scRNA_all$HTO_classification = scRNA_all$MULTI_ID
table(scRNA_all$HTO_classification)
Idents(scRNA_all) <- "hash.ID"
p1 = RidgePlot(scRNA_all, assay = "HTO", features = rownames(scRNA_all[["HTO"]])[1:4], ncol = 2)
p2 = FeatureScatter(scRNA_all, feature1 = "Hash-tag1", feature2 = "Hash-tag2")+
  FeatureScatter(scRNA_all, feature1 = "Hash-tag1", feature2 = "Hash-tag3")+
  FeatureScatter(scRNA_all, feature1 = "Hash-tag2", feature2 = "Hash-tag3")+
  FeatureScatter(scRNA_all, feature1 = "Hash-tag3", feature2 = "Hash-tag4")
Idents(scRNA_all) <- "HTO_classification"
p3 = VlnPlot(scRNA_all, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
ggsave(paste0("./demultiplex/HTO_RidgePlot.pdf"),p1,dpi=300,height=10,width = 12)
ggsave(paste0("./demultiplex/HTO_Scatter.pdf"),p2,dpi=300,height=10,width = 12)
ggsave(paste0("./demultiplex/HTO_VlnPlot.pdf"),p3,dpi=300,height=10,width = 12)

#单细胞转录组数据处理
DefaultAssay(scRNA_all) <- "RNA"
scRNA_all <- NormalizeData(scRNA_all) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA()
scRNA_all <- RunTSNE(scRNA_all, reduction = "pca", dims = 1:30)
scRNA_all <- FindNeighbors(scRNA_all,  dims = 1:30) %>% FindClusters(resolution=0.3)

## 去除Doublets (doublet_scores > 0.2)细胞 #### 
#devtools::install_github("furlan-lab/scrubletR")
library(scrubletR)
scRNA_all <- scrublet(scRNA_all)
p1 = FeaturePlot(scRNA_all, features = "doublet_scores",reduction = "tsne")
ggsave("./Figures/scrublet_doublet_score.pdf",p1,height=5,width = 6)

hist(scRNA_all$doublet_scores)
scRNA_all$doublets<-scRNA_all$doublet_scores > 0.2   #(You pick this)
p2 = DimPlot(scRNA_all, group.by = "doublets", cols=c("goldenrod", "navy"),reduction = "tsne")
ggsave("./Figures/scrublet_doublet.pdf",p2,height=5,width = 6)

## TSNE plot ####
scRNA = scRNA_all[,!scRNA_all$doublets]
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA()
scRNA <- RunTSNE(scRNA, dims = 1:30)
scRNA <- FindNeighbors(scRNA, dims = 1:30) %>% FindClusters(resolution=0.3)

scRNA = SetIdent(scRNA,value = "seurat_clusters")
p1 = DimPlot(scRNA, reduction = "tsne", label=T)+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(3,"pt"))
scRNA = SetIdent(scRNA,value = "celltype")
p2 = DimPlot(scRNA, reduction = "tsne", label=T)+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(3,"pt"))

ggsave(paste0("./Figures/TSNE_cluster_celltype.pdf"),p1+p2,dpi=300,height = 5,width = 15)

#TSNE图标注样本
scRNA = SetIdent(scRNA,value = "HTO_classification")
scRNA$HTO_classification = factor(scRNA$HTO_classification,
                                  levels = c("Hash-tag1", "Hash-tag2", "Hash-tag3", "Hash-tag4", "Doublet", "Negative"))
p3 = DimPlot(scRNA, reduction = "tsne", label=F, alpha = 0.8)+
  theme_void()+
  theme(legend.position = "right",legend.box = "vertical",legend.key.size = unit(3,"pt"))+
  guides(color=guide_legend(override.aes =list(size = 5)))+
  scale_color_discrete(type = c(color,"grey","black"))

ggsave(paste0("./Figures/TSNE_cluster_sample.pdf"),p3,dpi=300,height = 5,width = 6)

#HTO violin plot
DefaultAssay(scRNA) <- "HTO"
plot_violin = VlnPlot(scRNA, features=c("Hash-tag1", "Hash-tag2", "Hash-tag3", "Hash-tag4"),
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
ggsave(paste0("./Figures/HTO_violin_clusters.pdf"),plot=plot_violin,height = 5,width = 5)


#保存文件
# Export the analysis source data for figure reproduction.
saveRDS(scRNA_all,"./Rawdata/10X_hash_all_cells.rds")
saveRDS(scRNA,"./Rawdata/10X_hash_cells.rds")

#基于HLA-Splitter样本拆分######
# Read the saved object or source table used by this analysis.
sampleInfo = read.csv(file.path(getOption("hla.data_root"), "hash10x", "HLA-genotypes.csv"),header=TRUE,row.names = 1)
namelist=colnames(sampleInfo)
length(namelist)
HLA_list=list()
for (i in 1:length(namelist)){HLA_list[[i]] = unique(sampleInfo[,i])}
names(HLA_list) = namelist
HLA_label = read.delim("./Data/HLA_counts/sorted_genotypes.txt",header = F)
barcodes=read.delim(paste0("./Data/HLA_counts/output.barcodes.txt"),header=F)
HLA_matrix = readMM(paste0("./Data/HLA_counts/HLA_type.mtx"))  %>% as.matrix() %>% as.data.frame()
colnames(HLA_matrix) = HLA_label$V1
rownames(HLA_matrix) = paste0(barcodes$V1,"-1")

#细胞样本来源标记
outdir = "./Result/demultiplex/"
scRNA = readRDS(file.path(getOption("hla.data_root"), "hash10x", "10X_hash_cells.rds"))
sample_label=scRNA$HTO_classification
barcodes = intersect(paste0(barcodes$V1,"-1"),names(sample_label))
sample_label = sample_label[barcodes]
HLA_matrix = HLA_matrix[barcodes,]
length(barcodes)/ncol(scRNA)

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
  scale_color_manual(values = c(color[1:length(namelist)], "grey","black"))+
  labs(title = "True Labels", color = "Sample")

p2 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = Predict_label)) +
  geom_point(size = 0.2, alpha = 0.8) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_color_manual(values = c(color[1:length(namelist)], "grey","black"))+
  labs(title = "Predicted Labels", color = "Sample")

combined_plot <- p1 + p2
print(combined_plot)

ggsave(paste0(outdir, "/Combined_TSNE_plot.pdf"), plot = combined_plot, height = 5, width = 11)


#比较两种方法的一致程度####
scRNA$HLA_Splitter = sample_predict
scRNA$HLA_Splitter[is.na(scRNA$HLA_Splitter)] = "Undefined"
df = data.frame(MULTI_demutliplex = scRNA$HTO_classification,HLA_demutliplex = scRNA$HLA_Splitter)
df = table(df) %>% prop.table(1) %>%as.data.frame()
df$MULTI_demutliplex = factor(df$MULTI_demutliplex, levels=c("Hash-tag1", "Hash-tag2", "Hash-tag3", "Hash-tag4", "Doublet", "Negative"))
df$HLA_demutliplex = factor(df$HLA_demutliplex, levels=c(namelist,"Undefined"))

df = df[!df$MULTI_demutliplex %in% c("Doublet","Negative"),]

p1 = ggplot(df, aes(x = MULTI_demutliplex,y = Freq, fill = HLA_demutliplex))+
  geom_col(alpha = 0.7)+
  geom_text(aes(label = paste0(round(Freq*100,1),"%")),position = position_stack(vjust = 0.5))+
  ylab("Cell frequency (%)")+xlab("")+
  theme_classic()+
  scale_fill_discrete(type = c(color[1:length(namelist)],"grey"))
p1
ggsave(paste0(outdir,"/label_consistency.pdf"),p1,height=6,width = 6)

df = data.frame(MULTI_demutliplex = scRNA$HTO_classification,HLA_demutliplex = scRNA$HLA_Splitter)
df = table(df) %>% prop.table(1)*100
write.csv(outdir,"/demultiplex_res.csv")

#绘制score流形图######
library(ggridges) 
score_mat = data.frame(score,sample_predict=sample_predict,sample_label=sample_label)
p1 = ggplot(score_mat,aes(x=Donor1.score,y=sample_predict,fill=sample_predict))+
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none")
p2 = ggplot(score_mat,aes(x=Donor2.score,y=sample_predict,fill=sample_predict))+
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none")
p3 = ggplot(score_mat,aes(x=Donor3.score,y=sample_predict,fill=sample_predict))+
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none")
p4 = ggplot(score_mat,aes(x=Donor4.score,y=sample_predict,fill=sample_predict))+
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none")
p1+p2+p3+p4
ggsave(paste0(outdir,"/score_ridgeplot.pdf"),p1+p2+p3+p4,width = 8,height = 5)

#绘制score单细胞heatmap
library(grid)
score_mat = data.frame(score,sample_predict=sample_predict,sample_label=sample_label)
score_mat = score_mat[order(score_mat$sample_predict),]
#score_mat$sample_predict = factor(score_mat$sample_predict,levels=c(namelist,"Undefined","Doublet"))
mat = score_mat[,c(1:length(namelist))] %>% t()
annotation = data.frame(sample_label = score_mat$sample_predict)
rownames(annotation)=colnames(mat)
annotation$sample_label = factor(annotation$sample_label,levels = c(namelist,"Undefined","Doublet"))
ann_colors = c(color[1:length(namelist)],"black","grey")
names(ann_colors) = c(namelist,"Undefined","Doublet")
ann_colors = list(sample_label = ann_colors)
p1 = pheatmap(mat,cluster_rows = F,cluster_cols = F,scale="column",
              rect_gp = gpar(col='white',lwd=2),
              row_names_side = "left",show_colnames=F,annotation_col=annotation,
              annotation_colors = ann_colors)
ggsave(paste0(outdir,"/score_heatmap.pdf"),p1,dpi=300,height=5,width =8)

##绘制HTO umi在HLA TSNE上的分布###
HTO_score = scRNA@assays$HTO@data %>% t()
HTO_score = HTO_score[rownames(tsne_result),]
colnames(HTO_score) = c("Hash.tag1","Hash.tag2","Hash.tag3","Hash.tag4")
tsne_result = cbind(tsne_result,HTO_score)

p1 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=Hash.tag1))+
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")

p2 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=Hash.tag2)) +
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")

p3 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=Hash.tag3)) +
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")

p4 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=Hash.tag4)) +
  geom_point(size=0.1)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(type = "viridis")
p = ggarrange(p1,p2,p3,p4,common.legend = T)
ggsave(paste0(outdir,"/HTO_score_TSNE.pdf"),p ,dpi=300,height=10,width = 10)

## 细胞亚群HLA-score绘图 ####
#打分正则化
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

# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: lavage_find_signature.R
# Figure/table scope: Fig3; FigS6; Supplementary Table 13
# Data cohort: 肺泡灌洗液. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

## 不同分组差异基因分析(RE vs NR) #########
cell_subset = Tcell[,Tcell$Group %in% c("Infection","Tumor","Interlung")]
#cell_subset = Myeloid

table(cell_subset$HLA_Splitter)
#寻找差异基因
cell_subset$is.Interlung= ifelse(cell_subset$Group == "Interlung","Interlung","Others")
Idents(cell_subset) <- "is.Interlung"
# Discover marker genes with the original assay and comparison settings.
CELLDEG <- FindMarkers(cell_subset,ident.1 = "Interlung",ident.2 = "Others",logfc.threshold = 0.5, min.pct = 0.1, test.use = "wilcox")
CELLDEG = CELLDEG[CELLDEG$avg_log2FC>1 & CELLDEG$p_val_adj<0.05,]
# Export the analysis source data for figure reproduction.
write.csv(CELLDEG,"./Rawdata/Tcell_Interlung_VS_Others_DEG.csv")

## 细胞信号通路分析（AUCell)
library(AUCell)
library(clusterProfiler)
library(GSEABase)

#cell_subset = scRNA
cell_subset = Tcell
#cell_subset = Myeloid

counts <- GetAssayData(cell_subset,layer = "data")
cells_rankings <- AUCell_buildRankings(counts)

# 定义排异特征基因
# Read the saved object or source table used by this analysis.
geneSets <- read.csv("./MSigDB/Lung_disease_DEGs.csv")
geneSets = as.list(geneSets)
geneSets = lapply(geneSets,function(x)x[nchar(x)>0])

# 计算AUC score
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

# 可视化(RE vs. NR)
names(geneSets)
cell_subset$Infection_Score = as.numeric(getAUC(cells_AUC)["Infection", ])
cell_subset$Tumor_Score = as.numeric(getAUC(cells_AUC)["Cancer", ])
cell_subset$Interlung_Score = as.numeric(getAUC(cells_AUC)["Interstitial", ]) 

hist(cell_subset$Infection_Score)

#绘制各样本signature score
Sample_Info = read.csv("./样本分组/样本信息.csv")
df = FetchData(cell_subset, c("Group","HLA_Splitter","celltype",
                              "Infection_Score","Tumor_Score","Interlung_Score"))
df_tmp = aggregate(df[,4:6], list(df$Group,df$HLA_Splitter),mean)
df_tmp[3:5] = apply(df_tmp[3:5],2,scale)
colnames(df_tmp) = c("Group","SampleName","Infection_Score","Tumor_Score","Interlung_Score")

#Heatmap
df_tmp = melt(df_tmp,id.vars = c("Group","SampleName"),
              measure.vars = c("Infection_Score","Tumor_Score","Interlung_Score"),
              variable.name = "Signatrue",value.name = "Score")

p = ggplot(df_tmp, aes(x = Signatrue, y = SampleName, fill = Score))+
  geom_raster()+theme_classic()+
  scale_fill_viridis_c()
p
ggsave("./Figures/scRNA_Heatmap_Signature_Score.pdf",p,height = 6,width = 4.5)

#箱型图
comparison = list(c("Infection","Interlung"),c("Infection","Tumor"),c("Interlung","Tumor"))
p1 = ggplot(df_tmp,aes(x=Group,y=Score))+ facet_grid(.~Signatrue)+
  geom_boxplot(fill="white",outliers = F)+
  geom_point(aes(color=Group),size=3,position = position_jitter(0.2),alpha = 0.7)+
  stat_compare_means(aes(label = ..p.signif..),comparisons = comparison)+
  ylab("Signature score of Rejection")+
  scale_color_discrete(type = c(color))+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_classic()+theme(text = element_text (size =15),
                        axis.title.x = element_blank(),
                        plot.title = element_text(hjust = 0.5))
p1
ggsave("./Figures/Tcell_Signature_sample.pdf",p1,height = 4,width = 8)

#计算AUC曲线####
#install.packages("pROC")
library(pROC)
Roc_All = data.frame()
AUC = c()
for (group_tmp in c("Infection","Tumor","Interlung")){
  tmp = df_tmp[df_tmp$Signatrue == paste0(group_tmp,"_Score"),]
  tmp = data.frame(True_label = ifelse(tmp$Group == group_tmp, 1,0),
                   score = tmp$Score)
  Roc = plot.roc(tmp$True_label, tmp$score, percent = TRUE, main="Smoothing")
  AUC = c(AUC,paste0(group_tmp,"_Score:",round(Roc$auc,2),"%"))
  tmp = data.frame(sensitivity = Roc$sensitivities, 
                   specificity = Roc$specificities,
                   Sample = paste0(group_tmp,"_Score"),
                   AUC = paste0(group_tmp,"_Score:",round(Roc$auc,2),"%"))
  Roc_All = rbind(Roc_All,tmp)
}

p = ggplot(Roc_All,aes(x=specificity,y=sensitivity,color=Sample))+
  geom_line(size = 1)+
  annotate("text", x=20, y=5, label= AUC[1],size=4)+
  annotate("text", x=20, y=15, label= AUC[2],size=4)+
  annotate("text", x=20, y=25, label= AUC[3],size=4)+
  theme_classic()+
  ylab("Sensitivity (%)")+xlab("1-Specificity (%)")+
  guides(color=guide_legend(override.aes =list(size = 5)))+
  scale_color_discrete(type = c(color))
p
ggsave(paste0("./Figures/Tcell_Signature_score_Plot.pdf"),p,dpi=300,height=4,width = 5)

#绘制signature score在UMAP上分布
p1 = plot_density(cell_subset, c("Infection_Score","Tumor_Score","Interlung_Score"),pal = "inferno",reduction = "umap")
p1
FeaturePlot(cell_subset, c("Infection_Score"),reduction = "umap")
FeaturePlot(cell_subset, c("Tumor_Score"),reduction = "umap")
FeaturePlot(cell_subset, c("Interlung_Score"),reduction = "umap")
ggsave("./Figures/Tcell_Signature_Featureplot.pdf",p1,height = 3.5,width = 12)

#绘制signature score在celltypes上分布
p1 = VlnPlot(cell_subset, group.by = "celltype",features = c("Infection_Score"), alpha =0.7,
             pt.size=0, sort=T)+NoLegend()+
  coord_flip()+
  geom_boxplot(fill = "white",position = position_dodge(width = 0.8),width = 0.3)+
  ylab("Signature score of BK")+xlab(NULL)+
  scale_color_discrete(type = color)+
  scale_fill_discrete(type = color)
p1
ggsave("./Figures/Tcell_Infection_Score_vlnplot.pdf",p1,height = 6,width = 4)

#Heatmap plot (celltype)
df = FetchData(cell_subset, c("Group","celltype","Infection_Score","Tumor_Score","Interlung_Score"))
df_tmp = aggregate(df[,c(3:5)], list(df$Group,df$celltype),mean)
colnames(df_tmp) = c("Group","celltype","Infection_Score","Tumor_Score","Interlung_Score")
Cell_precent = table(df[,c(1,2)]) %>% prop.table(1) %>% data.frame()
Cell_precent = Cell_precent[!is.na(Cell_precent$Freq) & Cell_precent$Freq!=0,]
df_tmp$Cell_precent = Cell_precent$Freq*100
p1 = ggplot(df_tmp, aes(x = Group, y = celltype, color = Infection_Score, size = Cell_precent))+
  geom_point()+theme_bw()+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.5))+
  scale_color_viridis_c()
p2 = ggplot(df_tmp, aes(x = Group, y = celltype, color = Tumor_Score, size = Cell_precent))+
  geom_point()+theme_bw()+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.5))+
  scale_color_viridis_c()
p3 = ggplot(df_tmp, aes(x = Group, y = celltype, color = Interlung_Score, size = Cell_precent))+
  geom_point()+theme_bw()+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.5))+
  scale_color_viridis_c()
p1+p2+p3
ggsave(paste0("./Figures/Tcell_celltype_Score.pdf"),p1+p2+p3,height = 5,width = 14)

#scatter plot (Sample)
library(ggrepel)
df = FetchData(cell_subset, c("Group","HLA_Splitter","Infection_Score","Tumor_Score","Interlung_Score"))
df_tmp = aggregate(df[,3:5], list(df$Group,df$HLA_Splitter),mean)
colnames(df_tmp) = c("Group","Sample","Infection_Score","Tumor_Score","Interlung_Score")
df_tmp[,3:5] = apply(df_tmp[,3:5],2,scale)
p = ggplot(df_tmp,aes(Infection_Score,Tumor_Score,color = Group,label = Sample))+
  geom_point(size = 5, alpha = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  geom_text_repel(size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE,max.overlaps =30)+
  scale_color_discrete(type = c(color))+
  theme_classic()
p  
ggsave(paste0("./Figures/TCell_signature_scatter.pdf"),p,height = 4,width = 5.5)

## 绘制基因表达heatmap
library(ComplexHeatmap)
library(circlize)

cell_subset = Tcell[,Tcell$Group %in% c("Infection","Interlung","Tumor")]
gene_matrix <- AverageExpression(cell_subset,group.by = c("Group"),assays = "RNA", return.seurat = T)
gene_matrix <- gene_matrix@assays$RNA$scale.data %>% as.matrix()

#selected_gene <- read.csv("./MSigDB/Lung_disease_DEGs.csv")
#selected_gene = melt(selected_gene,measure.vars = c("Infection","Interstitial","Cancer"),variable.names("Signature","Score"))

selected_gene <- read.csv("./Rawdata/Tcell_AllGroup_DEG_wilcox.csv")
selected_gene = intersect(selected_gene$X,rownames(gene_matrix))
gene_matrix_tmp = gene_matrix[selected_gene,]

annotated_gene = c("CMTM8","PELI2","HSPA1B","KLRC4","CD8B",
                   "SELL","IFI44L","CMPK2","XCL2","LY96",
                   "GIMAP8","RTP4","HLA-DPB1","HLA-DPA1",
                   "LGALS1","HLA-DQB1","HLA-DRA","GZMB","PRF1",
                   "IFNG","FCRL6","CCL4L2","CCL3","IFIT3",
                   "HSPA6","NMB","PASK","FCGR3A","TUBB2A","GZMB","KIR3DL1",
                   "IL10","CX3CR1","LILRB1","CLU","CD80","HFE","FOXP3")

gene_pos = c()
for (i in 1:length(annotated_gene)){
  tmp = which(rownames(gene_matrix_tmp) == annotated_gene[i])
  names(tmp) = rep(annotated_gene[i],length(tmp))
  gene_pos = append(gene_pos,tmp[])
}
row_anno = rowAnnotation(gene = anno_mark(at=gene_pos,labels=names(gene_pos)))
col_fun = colorRamp2(c(-2,0,2),c("navy","white","firebrick3"))

pdf(paste0("./Figures/Tcell_Heatmap_markergene.pdf"), width = 5, height = 8)
Heatmap(gene_matrix_tmp ,cluster_rows = F,cluster_columns = F, 
        show_column_names = T, show_row_names = F,
        right_annotation = row_anno,
        column_title = NULL,row_title = NULL,
        heatmap_legend_param = list(title="Expression",title_position="leftcenter-rot"),
        col=col_fun)
dev.off()

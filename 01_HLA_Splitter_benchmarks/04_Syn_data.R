# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Upstream synthetic-data preparation and benchmark computation.
# Original source: 01_Syn_data.R
# Figure/table scope: Fig1b; Fig1c; Fig1d; FigS2d; FigS2e
# Data cohort: synthetic. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

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

dir.create("Figures")

## 读取40例PBMC样本，构建seurat对象 ####
# Read the saved object or source table used by this analysis.
scRNAlist = readRDS("./Result/scRNAlist.RDS")

## 自定义随机样本生成，随机选择x个细胞 ####
cell_num_select = c(10,25,50,100,250,500,1000,2500,5000)
for (cell_num_tmp in cell_num_select){
  for (j in 1:3){
    df = data.frame()
    for (i in 1:length(namelist)){
      barcode = colnames(scRNAlist[[i]])
      barcode_select = sample(barcode,min(length(barcode),cell_num_tmp))
      df_tmp = data.frame(barcode_select,sample = namelist[i])
      df = rbind(df,df_tmp)
    }
    df = df[!duplicated(df$barcode_select),]
# Export the analysis source data for figure reproduction.
    write.csv(df, paste0("./Data/Synthetic data/Barcodes/",cell_num_tmp,"cells-",j,".csv"))
    write.table(df$barcode_select, paste0("./Data/Synthetic data/Barcodes/",cell_num_tmp,"cells-",j,".txt"),row.names = F, col.names = F, quote = F)
  }
}

#从txt文件重新生成csv文件
barcode_sample_list = data.frame(barcode = colnames(scRNA_merge), sample = scRNA_merge$orig.ident)
table(barcode_sample_list$sample)
for (cell_num_tmp in cell_num_select){
  for (j in 1:3){
    barcode_select = read.delim(paste0("./Rawdata/",cell_num_tmp,"cells-",j,".txt"),header = F)[,1]
    df = barcode_sample_list[barcode_select,]
    colnames(df) = c("barcode_select","sample")
    write.csv(df, paste0("./Rawdata/",cell_num_tmp,"cells-",j,".csv"))
  }
}

# HLA matrix读取######
sampleInfo = read.csv(file.path(getOption("hla.data_root"), "synthetic", "HLA-genotypes.csv"),header=TRUE,row.names = 1)
namelist=colnames(sampleInfo)
length(namelist)
# HLA list
HLA_list=list()
for (i in 1:length(namelist)){HLA_list[[i]] = unique(sampleInfo[,i])}
names(HLA_list) = namelist
#HLA_matrix_All
HLA_matrix_All=matrix()
for (i in 1:length(namelist)){
  HLA_label = read.delim(paste0("./Data/HLA_counts/",namelist[i],"/HLA_genotype.tsv"),header = F)
  barcodes=read.delim(paste0("./Data/HLA_counts/",namelist[i],"/barcodes.tsv"),header=F)
  HLA_matrix = readMM(paste0("./Data/HLA_counts/",namelist[i],"/HLA_count_matrix.mtx"))  %>% as.matrix() %>% as.data.frame()
  colnames(HLA_matrix) = HLA_label$V1
  barcodes = paste0(barcodes$V1,"-1")
  rownames(HLA_matrix) = paste0(barcodes)
  if(nrow(HLA_matrix_All)>1){
    HLA_matrix_All = rbind(HLA_matrix_All,HLA_matrix)}else{
      HLA_matrix_All = HLA_matrix}
}

# 批量样本测试(cell number 50-2500)######
cell_num_select = c(10,25,50,100,250,500,1000,2500)
performance_list <- data.frame()

for (cell_num_tmp in cell_num_select){
  for (j in 1:3){
    #细胞样本来源标记
    data_tmp = read.csv(paste0("./Data/Synthetic data/Barcodes/",cell_num_tmp,"cells-",j,".csv"))
    sample_label = data_tmp$sample
    barcodes = data_tmp$barcode_select
    names(sample_label) = barcodes
    
    HLA_matrix = HLA_matrix_All[barcodes,]
    HLA_matrix[is.na(HLA_matrix)] = 0
    rownames(HLA_matrix) = barcodes
    nrow(HLA_matrix)
    
    #新建文件
    outdir = paste0("./demultiplex/",cell_num_tmp,"cells-",j)
    dir.create(outdir)
    
    #HLA-Splitter demultiplex
    source(file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks/01_HLA_demultiplex_v2.R"))
    sample_predict <- HLA_demultiplex(HLA_list = HLA_list, 
                                      HLA_matrix = HLA_matrix, 
                                      graph_corr = TRUE, 
                                      outdir = outdir)
    table(sample_predict)
    
    #evaluation
    sample_predict = read.csv(paste0(outdir,"/HLA_demultiplex.csv"))
    sample_predict = sample_predict$Sample_prediction
    
    source(file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks/02_Demultiplex_evaluation.R"))
    demultiplex_evaluation(sample_label,sample_predict,outdir)
    
    # 计算正确比对的比例
    sample_predict[is.na(sample_predict)] <- "Undefined"
    total_cells <- length(sample_label)
    count_undefined <- sum(sample_predict == "Undefined")
    count_correct <- sum(sample_predict == sample_label & sample_predict != "Undefined")
    count_incorrect <- total_cells - count_correct - count_undefined 
    
    # 将结果整理为单行 Data Frame 并存入列表
    current_res <- data.frame(
      Cell_Number = cell_num_tmp,
      Replicate = j,
      Total_Cells = total_cells,
      Correct_Ratio = count_correct / total_cells,
      Incorrect_Ratio = count_incorrect / total_cells,
      Undefined_Ratio = count_undefined / total_cells
    )
    performance_list <- rbind(performance_list,current_res)
    
  }
}
write.csv(performance_list, file.path(getOption("hla.data_root"), "synthetic", "Demultiplex_performance.csv"), row.names = FALSE)

#绘制HLA-Splitter拆分结果统计图(cell number 50-2500)
library(tidyr)
df <- read.csv(file.path(getOption("hla.data_root"), "synthetic", "Demultiplex_performance.csv"))
df_long <- df %>%
  select(Cell_Number, Replicate, Correct_Ratio, Incorrect_Ratio, Undefined_Ratio) %>%
  pivot_longer(
    cols = ends_with("_Ratio"),
    names_to = "Label",
    values_to = "Percent"
  ) %>%
  mutate(
    Percent = Percent * 100, # 将比例(0-1)转换为百分比(0-100)
    Label = gsub("_Ratio", "", Label), # 去掉后缀，留下 Correct, Incorrect, Undefined
    Cell_Number = factor(Cell_Number)
  )

df_sum <- df_long %>%
  group_by(Cell_Number, Label) %>%
  summarise(
    Mean_Percent = mean(Percent, na.rm = TRUE),
    SD_Percent = sd(Percent, na.rm = TRUE),
    .groups = "drop"
  )

df_sum$Label <- factor(df_sum$Label, levels = c("Undefined", "Incorrect", "Correct"))

df_line <- df_sum %>%
  filter(Label == "Correct") %>%
  mutate(
    ymin_SD = Mean_Percent - SD_Percent,
    ymax_SD = Mean_Percent + SD_Percent
  )

p <- ggplot(df_sum, aes(x = Cell_Number, y = Mean_Percent)) + 
  geom_col(aes(fill = Label), position = "stack", alpha = 0.8, linewidth = 0.3) + 
  geom_line(data = df_line, aes(y = Mean_Percent, group = 1), color = "black", linewidth = 0.8) + 
  geom_point(data = df_line, aes(y = Mean_Percent, group = 1), color = "black", size = 2) + 
  geom_errorbar(data = df_line, aes(ymin = ymin_SD, ymax = ymax_SD),
                width = 0.2, color = "black", linewidth = 0.8) + 
  geom_text(data = df_line, aes(y = ymax_SD, label = paste0(round(Mean_Percent, 1), "%")),
            vjust = +2, color = "black", size = 4) + 
  labs(
    title = "Demultiplexing accuracy across samples",
    x = "Multiplexed cells per sample",
    y = "Cell percentage (%)",
    fill = "Prediction"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 15), 
    axis.title = element_text(size = 16, face = "bold"), 
    axis.text = element_text(size = 14, color = "black"), 
    axis.text.x = element_text(angle = 30, vjust = 0.5), 
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5) 
  ) +
  scale_fill_manual(values = c("Correct" = "#CEDE94", "Incorrect" = "firebrick", "Undefined" = "grey80"))

print(p)
ggsave("./Figures/Cellnum_Benchmark_Plots.pdf", plot = p, width = 6.5, height = 5, dpi = 300)

df_sum = aggregate(df_long$Percent,list(df_long$Cell_Number,df_long$Label),mean)
colnames(df_sum) = c("Cell_Number", "Label", "Percent(%)")
SD_Value = aggregate(df_long$Percent, list(df_long$Cell_Number,df_long$Label),sd)
df_sum$St.d = SD_Value$x
write.csv(df_sum,"./Result/Demutliplex performance summary.csv")

#对40样本，500细胞混样结果进行细致刻画#######
data_tmp = read.csv("./Data/Synthetic data/Barcodes/500cells-1.csv")
sample_label = data_tmp$sample
barcodes = data_tmp$barcode_select
names(sample_label) = barcodes

HLA_matrix = HLA_matrix_All[barcodes,]
HLA_matrix[is.na(HLA_matrix)] = 0
rownames(HLA_matrix) = barcodes
nrow(HLA_matrix)

outdir = "./Result/500cells-1/"
dir.create(outdir)

# Sample demultiplexing
source(file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks/01_HLA_demultiplex_v2.R"))
sample_predict <- HLA_demultiplex(HLA_list = HLA_list, 
                                  HLA_matrix = HLA_matrix, 
                                  graph_corr = T, 
                                  outdir = outdir)
table(sample_predict)

#label proportion
df = data.frame(label = sample_label,predict = sample_predict)
df$class = "Correct"
df$class[as.character(df$label)!=as.character(df$predict)] = "Incorrect"
df$class[df$predict=="Undefined"] = "Undefined"

df$class = factor(df$class, levels = c("Correct","Undefined","Incorrect"))
df$label = factor(df$label, levels = c(namelist))
df_tmp = data.frame(prop.table(table(df$class))*100)
df_tmp$label = paste0(df_tmp$Var1,": ",round(df_tmp$Freq,2),"%")

df$HLA_count = log(rowSums(HLA_matrix)+1,10)

p = ggplot(df_tmp,aes(x= "ratio",y=Freq,fill=label))+
  geom_col(position="stack",alpha=0.5)+coord_polar(theta = 'y')+
  theme_void()+theme(text = element_text (size =15))+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_fill_discrete(type = c("#CEDE94","firebrick","grey"))
p
ggsave(paste0("./Figures/label_proportion.pdf"),p,height=5,width = 8)

# HLA-tsne result of true label and predict label
tsne_result$True_label = sample_label
all_samples <- unique(c(as.character(tsne_result$True_label), as.character(tsne_result$sample_predict)))

valid_samples <- setdiff(all_samples, "Undefined")
color_mapping <- setNames(color[1:length(valid_samples)], valid_samples)
color_mapping["Undefined"] <- "grey80" 

# 绘制 TSNE 图并应用手动颜色映射
p1 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = True_label)) +
  geom_point(size = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  # 使用 scale_color_manual 强制应用我们的命名颜色字典
  scale_color_manual(values = color_mapping, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "True Labels")

p2 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = sample_predict)) +
  geom_point(size = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  scale_color_manual(values = color_mapping, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "Predicted Labels")

combined_plot <- p1 + p2
print(combined_plot)
ggsave("./Figures/HLA_TSNE_Prediction.pdf", plot = combined_plot, dpi = 300, height = 5, width = 13)

#计算回收样本比例
df_tmp = data.frame(prop.table(table(df[,c(1,3)]),1)*100)
df_tmp$label = factor(df_tmp$label,levels = namelist)
df_tmp = df_tmp[!is.na(df_tmp$label),]
p = ggplot(df_tmp ,aes(x = class,y=Freq))+
  geom_boxplot(fill="white",color="black", size = 1, outliers = F)+
  geom_point(aes(color=label), position = position_jitter(width = 0.3),size=3, alpha = 0.8)+
  theme_bw()+ylab("Percentage in each sample(%)")+
  geom_text(aes(x=1,y=50,label = round(mean(df_tmp$Freq[df_tmp$class == "Correct"]),2)))+
  geom_text(aes(x=2,y=50,label = round(mean(df_tmp$Freq[df_tmp$class == "Undefined"]),2)))+
  geom_text(aes(x=3,y=50,label = round(mean(df_tmp$Freq[df_tmp$class == "Incorrect"]),2)))+
  scale_color_manual(values = color_mapping, drop = FALSE) +
  theme_classic()+theme(text = element_text (size =15),
                        plot.title = element_text(hjust = 0.5),
                        axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75))
p
ggsave(paste0("./Figures/predict_fraction.pdf"),p,dpi=300,height=5,width=5)

# HLA预测结果在TSNE图中的分布
p1 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=HLA_count)) +
  geom_point(size=0.3)+theme_classic()+
  theme(text = element_text (size = 15))+
  scale_color_continuous(low = "#304F80", high = "#FDFF00")

p2 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=class)) +
  geom_point(size=0.3)+theme_classic()+
  theme(text = element_text (size = 15))+
  guides(color=guide_legend(override.aes =list(size = 5)))+
  scale_color_discrete(type = c("#CEDE94","grey","firebrick"))
p1+p2
ggsave(paste0("./Figures/HLA_TSNE_counts.pdf"),p1+p2,dpi=300,height=5,width = 13)

#绘制score热图
library(gtools)
score_df <- data.frame(score,
                       sample_predict = sample_predict,
                       sample_label = sample_label)
score_df <- score_df[mixedorder(score_df$sample_label), ]

mat <- score_df %>%
  select(all_of(paste0(namelist,".score"))) %>%
  t() 

annotation_col_df <- data.frame(sample_label = score_df$sample_label)
rownames(annotation_col_df) <- colnames(mat) 

desired_levels <- c(namelist, "Undefined")
annotation_col_df$sample_label <- factor(annotation_col_df$sample_label, levels = desired_levels)

ann_colors_list <- list(
  sample_label = c(color[1:length(namelist)], "grey")
)
names(ann_colors_list$sample_label) <- desired_levels

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
ggsave(paste0("./Figures/score_heatmap.pdf"),p1,dpi=300,height=5,width =8)


# 不同预测结果细胞的HLA丰度流形图
library(ggridges)
df$HLA_count = log(rowSums(HLA_matrix)+1,10)
tsne_result$HLA_count = df$HLA_count
df = df[is.finite(df$HLA_count),]

p3 <- ggplot(df, aes(x = HLA_count, y = class, fill = class)) +
  geom_density_ridges(alpha = 0.7, color = "white", scale = 0.9) + # 
  scale_fill_manual(values = c("#CEDE94", "grey", "firebrick")) +
  labs(
    x = "Log10 normalized HLA count", # Renamed x-axis title for clarity
    y = "Sample Class",               # Renamed y-axis title for clarity
    title = "Density Distribution of HLA Count by Sample Class" # Added a descriptive plot title
  ) +
  theme_ridges(font_size = 12, grid = TRUE) + # Apply theme_ridges first, set base font size, enable grid
  theme(
    legend.position = "none",        # Remove legend if 'class' is self-explanatory or colors are fixed
    text = element_text(size = 15),  # Base text size for all plot elements
    axis.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5), # Make axis titles bolder and slightly larger
    axis.text = element_text(size = 14), # Axis tick label size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )
p3
ggsave(paste0("./Figures/HLA_count_ridgeplot.pdf"),p3,dpi=300,height=5,width=6)
aggregate(df$HLA_count,list(df$class),mean)

# 检验HLA表达量低的细胞和HLA表达量高的细胞拆分准确率差异
library(tidyverse)
library(ggpubr)
df$HLA_group <- ifelse(df$HLA_count > 1, "High", "Low")

sample_stats <- df %>%
  group_by(HLA_group, label) %>%
  summarise(
    Correct = sum(class == "Correct") / n() * 100,
    Undefined = sum(class == "Undefined") / n() * 100,
    Incorrect = sum(class == "Incorrect") / n() * 100,
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(Correct, Undefined, Incorrect), 
               names_to = "Status", 
               values_to = "Percentage")

p_stat <- ggplot(sample_stats, aes(x = HLA_group, y = Percentage, fill = HLA_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 2, shape = 21, color = "black", alpha = 0.8) +
  facet_wrap(~Status, scales = "free_y") + 
  scale_fill_manual(values = c("Low" = "#CEDE94", "High" = "firebrick")) +
  labs(
    x = "HLA expression levels (HLA counts >10)",
    y = "Cell proportion per sample (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 14),
    strip.background = element_rect(fill = "white"), # 分面标签背景
    strip.text = element_text(face = "bold", size = 12)
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.3, label.y = 120)

print(p_stat)
ggsave(paste0("./Figures/HLA_Consistency_Comparison_Facets.pdf"), p_stat, width = 6, height = 5)

# 比较不同细胞类型中的HLA-counts和拆分准确性
df_Celltype <- read.csv("./Result/scRNA_celltype.csv", row.names = 1)

df$cell_type <- df_Celltype[match(rownames(df), rownames(df_Celltype)), "Celltype"]
df = df[!is.na(df$cell_type),]
ct_stats <- df %>%
  group_by(cell_type) %>%
  summarise(
    correct_pct = sum(class == "Correct", na.rm = TRUE) / n() * 100,
    median_hla = median(HLA_count, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(median_hla)) # 按 HLA 中位数降序排序

df$cell_type <- factor(df$cell_type, levels = ct_stats$cell_type)
max_hla <- max(df$HLA_count, na.rm = TRUE)
scale_factor <- (max_hla * 1.05) / 100 

p_combined <- ggplot(df, aes(x = cell_type)) +
  geom_boxplot(aes(y = HLA_count, fill = cell_type), alpha = 0.4, outlier.size = 0.5) +
  geom_line(data = ct_stats, aes(y = correct_pct * scale_factor, group = 1), 
            color = "firebrick", linewidth = 1) + 
  geom_point(data = ct_stats, aes(y = correct_pct * scale_factor), 
             color = "firebrick", size = 3) +
  scale_y_continuous(
    name = "Log10 HLA Counts",
    sec.axis = sec_axis(~ . / scale_factor, name = "Correct labeled cells (%)")
  ) +
  labs(x = "Cell types") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "firebrick", face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

print(p_combined)
ggsave("./Figures/HLA_vs_Accuracy_DualAxis.pdf", p_combined, width = 6, height = 4)

write.csv(df,file.path(getOption("hla.data_root"), "synthetic", "demultiplex_result.csv"))

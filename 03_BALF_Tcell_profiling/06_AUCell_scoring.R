# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: 05_AUCell_scoring.R
# Figure/table scope: Fig3; FigS6; Supplementary Table 13
# Data cohort: balf. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# File: 05_AUCell_scoring.R
# Description: 基于数据驱动的差异基因提取、AUCell细胞亚群定位与疾病评分模型构建
# ==============================================================================
rm(list = ls())
source(file.path(getOption("hla.code_root"), "03_BALF_Tcell_profiling/01_setup_and_helpers.R"))
# Read the saved object or source table used by this analysis.
Tcell <- readRDS(file.path(getOption("hla.data_root"), "balf", "lavage_Tcell.rds"))
Sample_Info <- read.csv(file.path(getOption("hla.data_root"), "balf", "SampleInfo.csv"))

# 确保只保留目标分组并去除无用因子水平
cell_subset <- subset(Tcell, subset = Group %in% c("LC", "ILD", "CAP"))
cell_subset$HLA_Splitter %>% table()
cell_subset$Group <- factor(cell_subset$Group, levels = c("CAP", "LC", "ILD"))
comparison <- list(c("CAP","LC"), c("CAP","ILD"), c("ILD","LC"))

# ==============================================================================
# 1. 寻找各疾病分组的关键特征基因 (Data-driven Signatures)
# ==============================================================================
Idents(cell_subset) <- "Group"
message("正在计算各分组的差异表达基因...")

# AUCell 需要的 list 格式
#geneSets <- read.csv("./Result/MSigDB/Lung_disease_DEGs.csv")
#geneSets = as.list(geneSets)
#geneSets = lapply(geneSets,function(x)x[nchar(x)>0])

# geneSets <- split(filtered_markers$gene, filtered_markers$cluster)

Group_DEG = read.csv('./Result/Tcell_Group_DEG.csv')
filtered_markers = Group_DEG %>%
  filter(PB_Log2FC>0.1)  %>%
  group_by(Group_Col)
geneSets <- split(filtered_markers$gene, filtered_markers$Group_Col)


# ==============================================================================
# 2. 运行 AUCell 疾病评分模型
# ==============================================================================
message("正在构建 AUCell 排序并计算评分...")
names(geneSets) <- paste0(names(geneSets), "_Signature")
counts <- GetAssayData(cell_subset, layer = "data")
cells_rankings <- AUCell_buildRankings(counts, plotStats = FALSE, splitByBlocks = TRUE)

# 计算 AUC (aucMaxRank 通常设为基因总数的 5%-10%)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.1)

# 将计算结果注入到 Seurat 的 meta.data 中
signature_names <- names(geneSets)
for (sig in signature_names) {
  cell_subset[[sig]] <- as.numeric(getAUC(cells_AUC)[sig, ])
}

# ==============================================================================
# 3. 寻找特征性细胞亚群 (Cell-level Evaluation)
# ==============================================================================
# 3.1 UMAP 分布图
p_umap_sig <- FeaturePlot(cell_subset, features = signature_names, ncol = 3, pt.size = 0.5, reduction = "umap") & 
  scale_color_gradientn(colors = c("lightgrey", "firebrick3")) & theme_void()
ggsave("./Figures/Tcell_Signature_UMAP.pdf", p_umap_sig, height = 4, width = 12)

# 3.2 哪种细胞亚型具有最高的疾病特征(Violin Plot)
p_vln_sig <- VlnPlot(cell_subset, group.by = "Celltype", features = signature_names, stack = TRUE, pt.size = 0, flip = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"), legend.position = "none")
ggsave("./Figures/Tcell_Signature_Celltype_Vln.pdf", p_vln_sig, height = 6, width = 8)

# 3.3 DotPlot 综合气泡图 (展示细胞亚型 vs 疾病评分)
p_dot_sig <- DotPlot(cell_subset, group.by = "Celltype", features = signature_names) +
  coord_flip() +
  scale_color_gradientn(colors = c("navy", "white", "firebrick3")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(x = "Disease Signature", y = "Cell Type")
ggsave("./Figures/Tcell_Signature_Celltype_Dotplot.pdf", p_dot_sig, height = 4, width = 7)


# ==============================================================================
# 4. 建立样本级疾病评估模型与 ROC 检验 (Sample-level Evaluation)
# ==============================================================================
# 提取样本均值
df_sample <- FetchData(cell_subset, c("Group", "HLA_Splitter", signature_names))
df_agg <- aggregate(df_sample[, signature_names], by = list(Group = df_sample$Group, SampleName = df_sample$HLA_Splitter), FUN = mean)

# 4.1 样本级评分 Boxplot
df_melt <- reshape2::melt(df_agg, id.vars = c("Group", "SampleName"), variable.name = "Signature", value.name = "Score")
p_box_samp <- ggplot(df_melt, aes(x = Group, y = Score, fill = Group)) +
  facet_wrap(~ Signature, scales = "free_y") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  stat_compare_means(comparisons = comparison, label = "p.signif") +
  my_theme()+
  scale_fill_manual(values = GLOBAL_COLORS) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Sample-level Disease Signature Scores", x = "Actual Group", y = "Mean AUCell Score")
p_box_samp
ggsave("./Figures/Tcell_Signature_Sample_Boxplot.pdf", p_box_samp, height = 5, width = 6)

# scatter plot (Sample)
df_agg_scale = df_agg
df_agg[,3:5] = apply(df_agg[,3:5],2,scale)
p_scatter_samp = ggplot(df_agg,aes(CAP_Signature, LC_Signature,color = Group,label = SampleName))+
  geom_point(size = 3, alpha = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  geom_text_repel(size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE,max.overlaps =30)+
  scale_color_manual(values = GLOBAL_COLORS) +
  my_theme()+ coord_fixed()+
  xlim(c(-2,3))+ylim(c(-2,3))
p_scatter_samp  
ggsave(paste0("./Figures/TCell_signature_scatter.pdf"),p_scatter_samp,height = 4,width = 5.5)

# 4.2 ROC 曲线评估评分模型的诊断能力
Roc_All <- data.frame()
AUC_labels <- c()

# 循环三种疾病状态计算 ROC
for (grp in levels(cell_subset$Group)) {
  sig_name <- paste0(grp, "_Signature")
  # 真实标签：是否属于当前疾病组
  true_labels <- ifelse(df_agg$Group == grp, 1, 0)
  scores <- df_agg[[sig_name]]
  
  # 如果该组存在样本，则计算 ROC
  if (sum(true_labels) > 0) {
    Roc <- plot.roc(true_labels, scores, percent = TRUE, quiet = TRUE)
    AUC_labels <- c(AUC_labels, paste0(grp, " Model AUC: ", round(Roc$auc, 2), "%"))
    
    tmp_roc_df <- data.frame(
      sensitivity = Roc$sensitivities,
      specificity = Roc$specificities,
      Model = sig_name
    )
    Roc_All <- rbind(Roc_All, tmp_roc_df)
  }
}
Roc_All <- Roc_All %>% 
  arrange(Model, sensitivity)

# 绘制多分类 ROC 曲线
p_roc <- ggplot(Roc_All, aes(x = specificity, y = sensitivity, color = Model)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 100, linetype = "dashed", color = "grey50") +
  theme_classic(base_size = 14) +
  scale_x_reverse() +
  labs(title = "Diagnostic Performance of AUCell Signatures",
       x = "Specificity (%)", y = "Sensitivity (%)") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = c(0.5, 0.25), 
        legend.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 在图上添加 AUC 文本注释
y_pos <- c(35, 30, 25)
for(i in seq_along(AUC_labels)) {
  p_roc <- p_roc + ggplot2::annotate("text", x = 30, y = y_pos[i], label = AUC_labels[i], size = 4.5, fontface = "bold")
}
p_roc
ggsave("./Figures/Tcell_Signature_ROC_Plot.pdf", p_roc, dpi = 300, height = 5, width = 5)

# 保存包含评分的 Seurat 对象
# Export the analysis source data for figure reproduction.
saveRDS(cell_subset, "./Result/lavage_Tcell_Scored.rds")

# 保存每个样本的疾病评分
write.csv(df_agg, "./Result/Disease_signature_score.csv")

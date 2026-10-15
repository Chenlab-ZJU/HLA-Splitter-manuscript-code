# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Upstream HLA typing, resolution and PBMC demultiplexing.
# Original source: Nanopore.R
# Figure/table scope: FigS1; FigS2f; Supplementary Tables 1 and 9
# Data cohort: pbmc30. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

setwd(file.path(getOption("hla.data_root"), "pbmc30"))
rm(list = ls())

# 确保加载所有必需的包
library(tidyverse)
library(ggplot2)
library(gtools)    
library(reshape2)  
library(ggpubr) 
library(tidyr)

# ==========================================
# 1. 导入数据并处理 (Part 1: 1-field & 2-field 分级评估)
# ==========================================
# Read the saved object or source table used by this analysis.
df <- read.csv(file.path(getOption("hla.data_root"), "pbmc30", "HLA_typing.csv"))

# 辅助函数 1：提取 2-field 精度 (例如 A*11:01:01 -> A*11:01)
truncate_to_2field <- function(allele) {
  if(is.na(allele) || allele == "") return(NA)
  match <- str_match(allele, "^([A-Z0-9]+\\*[0-9]+:[0-9]+)")
  if(!is.na(match[1,2])) return(match[1,2]) else return(allele)
}

# 辅助函数 2：提取 1-field 精度 (例如 A*11:01:01 -> A*11)
truncate_to_1field <- function(allele) {
  if(is.na(allele) || allele == "") return(NA)
  match <- str_match(allele, "^([A-Z0-9]+\\*[0-9]+)")
  if(!is.na(match[1,2])) return(match[1,2]) else return(allele)
}

# 转换为长数据，并同时保留 1-field 和 2-field 的解析结果
df_long <- df %>%
  pivot_longer(cols = -c(Sample.ID, Method), 
               names_to = c("Gene", "Allele_Num"), 
               names_sep = "_", 
               values_to = "Allele_Raw") %>%
  rowwise() %>%
  mutate(
    Allele_2Field = truncate_to_2field(Allele_Raw),
    Allele_1Field = truncate_to_1field(Allele_Raw)
  ) %>%
  ungroup()

df_compare <- df_long %>%
  # 按样本、方法、基因分组排序，确保 Allele 1 和 2 正确对应
  group_by(Sample.ID, Method, Gene) %>%
  arrange(Allele_2Field, .by_group = TRUE) %>%
  mutate(Sorted_Allele_Idx = paste0("Allele_", row_number())) %>%
  ungroup() %>%
  
  # 将 SBT 和 ONT 的 1-field 与 2-field 结果并列排开
  pivot_wider(
    id_cols = c(Sample.ID, Gene, Sorted_Allele_Idx),
    names_from = Method, 
    values_from = c(Allele_2Field, Allele_1Field)
  ) %>%
  
  # 【核心修改】临床级 HLA 不一致性评估逻辑
  mutate(
    Concordance = case_when(
      is.na(Allele_2Field_SBT) | is.na(`Allele_2Field_ONT-Rapid`) | 
        grepl("NA", Allele_2Field_SBT) | grepl("NA", `Allele_2Field_ONT-Rapid`) ~ "Missing",
      
      # 1. 如果连 1-field 都没对上 (如 A*11 vs A*24) -> 严重错误 (红色)
      Allele_1Field_SBT != `Allele_1Field_ONT-Rapid` ~ "Mismatch (1-field)", 
      
      # 2. 1-field 对上了，但 2-field 没对上 (如 A*11:01 vs A*11:02) -> 精度瑕疵 (黄色)
      Allele_2Field_SBT != `Allele_2Field_ONT-Rapid` ~ "Mismatch (2-field)", 
      
      # 3. 2-field 完全一致 -> 完美匹配 (绿色)
      TRUE ~ "Match (≥2-field)"                                               
    ),
    # 提取 2-field 用于格子内文本展示
    Display_Allele = coalesce(`Allele_2Field_ONT-Rapid`, Allele_2Field_SBT)
  )

# 控制因子顺序，美化标签
df_compare$Gene <- factor(df_compare$Gene, levels = c("A", "B", "C", "DRB1", "DQB1"))
df_compare$Sample.ID = factor(df_compare$Sample.ID, levels = mixedsort(unique(df_compare$Sample.ID), decreasing = F))
df_compare$Sorted_Allele_Idx <- factor(df_compare$Sorted_Allele_Idx, 
                                       levels = c("Allele_1", "Allele_2"), 
                                       labels = c("Allele 1", "Allele 2"))

# 定义红、黄、绿颜色映射
my_colors <- c("Match (≥2-field)" = "#749B58FF",      # 绿色
               "Mismatch (2-field)" = "#E69F00",      # 黄色 (深琥珀)
               "Mismatch (1-field)" = "#CE3D32FF",    # 红色
               "Missing" = "grey80")                  # 灰色

# 绘制热图 (p1)
p1 <- ggplot(df_compare, aes(x = Sorted_Allele_Idx, y = Sample.ID, fill = Concordance)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Display_Allele), color = "white", size = 3.5, fontface = "bold") +
  
  # 分面按基因分块
  facet_grid(~ Gene, scales = "free_x", space = "free_x") +
  
  scale_fill_manual(values = my_colors) +
  labs(
    x = NULL, 
    y = "Sample ID"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.2, "lines") 
  )

# # ==========================================
# # 2. 导入数据并处理 (Part 2: Coverage)
# # ==========================================
# data = read.csv("./nanopore/coverage.csv")
# df_cov = melt(data, id.vars=c("SampleID"),
#               measure.vars = c("A_1","B_1","C_1","DRB1_1","DPB1_1","DQB1_1"),
#               variable.name = "HLA_locus", value.name = "locus_coverage")
# df_cov$SampleID = factor(df_cov$SampleID, levels = mixedsort(unique(df_cov$SampleID)))
# 
# locus_colors <- c("#CE3D32FF", "#F0E685FF","#749B58FF","#5050FFFF", "#466983FF", "#FF9E4A")
# 
# # 绘制箱线图 (p2)
# p2 <- ggplot() + 
#   facet_grid(HLA_locus ~ ., scales = "free_y") +
#   geom_boxplot(data = df_cov, aes(x=SampleID, y=locus_coverage, color=HLA_locus),
#                position = position_dodge2(width = 1), fill="white", alpha = 0.3, outlier.shape = NA) +
#   geom_point(data = df_cov, aes(x=SampleID, y=locus_coverage, color=HLA_locus),
#              position = position_dodge2(width = 0.75), size = 0.5, alpha = 0.5) +
#   ylab("Average locus coverage") +
#   xlab("Sample ID") +
#   scale_color_manual(values = locus_colors) + 
#   theme_bw() + 
#   theme(
#     text = element_text(size = 14),
#     axis.title = element_text(face = "bold"),
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
#     legend.position = "none",
#     strip.text.y = element_text(size = 12, face = "bold", angle = 270),
#     strip.background = element_rect(fill = "grey90")
#   )

# # ==========================================
# # 3. 对比纳米孔typing和二代测序HLAtyping对HLA-Splitter拆分结果的影响
# # ==========================================

HLA_res_ONT = read.csv("./Result/HLA_Splitter_ONT/HLA_demultiplex.csv")
HLA_res_SBT = read.csv("./Result/HLA_Splitter_SBT/HLA_demultiplex.csv")

# 绘制并排柱状图
ont_counts <- HLA_res_ONT %>%
  filter(!is.na(Sample_prediction)) %>%
  group_by(Sample_prediction) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  mutate(Method = "ONT-Rapid")

sbt_counts <- HLA_res_SBT %>%
  filter(!is.na(Sample_prediction)) %>%
  group_by(Sample_prediction) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  mutate(Method = "SBT")

df_counts <- bind_rows(ont_counts, sbt_counts) %>%
  rename(Sample_ID = Sample_prediction)

df_counts$Sample_ID <- factor(df_counts$Sample_ID, levels = mixedsort(unique(df_counts$Sample_ID)))

df_counts$Method <- factor(df_counts$Method, levels = c("SBT", "ONT-Rapid"))

method_colors <- c("SBT" = "#CE3D32FF", "ONT-Rapid" = "#466983FF")

p_bar <- ggplot(df_counts, aes(x = Sample_ID, y = Cell_Count, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7, alpha = 0.7) +
  # geom_text(aes(label = Cell_Count), 
  #           position = position_dodge(width = 0.8), 
  #           vjust = -0.5, size = 3.5, fontface = "bold", color = "black") +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Assigned Sample ID",
    y = "Number of Cells",
    fill = "Method"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "top",          # 将图例放在顶部节省横向空间
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
print(p_bar)

# 绘制heatmap热图
HLA_res_df = data.frame("HLA_res_ONT" = HLA_res_ONT$Sample_prediction,
                        "HLA_res_SBT" = HLA_res_SBT$Sample_prediction)

heatmap_data <- HLA_res_df %>%
  filter(!is.na(HLA_res_SBT) & !is.na(HLA_res_ONT)) %>%
  group_by(HLA_res_SBT, HLA_res_ONT) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(HLA_res_SBT) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

all_samples <- unique(c(as.character(heatmap_data$HLA_res_SBT), as.character(heatmap_data$HLA_res_ONT)))
sorted_samples <- mixedsort(all_samples)

heatmap_data$HLA_res_SBT <- factor(heatmap_data$HLA_res_SBT, levels = sorted_samples)
heatmap_data$HLA_res_ONT <- factor(heatmap_data$HLA_res_ONT, levels = rev(sorted_samples))


p_heatmap <- ggplot(heatmap_data, aes(x = HLA_res_SBT, y = HLA_res_ONT, fill = Proportion)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # geom_text(data = filter(heatmap_data, Proportion > 0),
  #           aes(label = scales::percent(Proportion, accuracy = 0.1)), 
  #           size = 3.5, color = "black") +
  scale_fill_gradient(low = "white", high = "#CE3D32FF", 
                      labels = scales::percent_format(accuracy = 1),
                      limits = c(0, 1)) +
  labs(
    x = "SBT Assigned Samples (Reference)",
    y = "ONT Assigned Samples",
    fill = "Concordance"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),     # 去除网格线
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )
print(p_heatmap)

# ==========================================
# 3. 图像合并 (Combined Plot)
# ==========================================
message("Combining plots...")

# 上下合并
final_plot <- ggarrange(p1, p_bar, p_heatmap,
                        ncol = 1, nrow = 3, 
                        labels = c("A", "B" , "C"), 
                        heights = c(1.5, 1, 1.5))

# 建议宽度设定为 14
ggsave("./Figures/Combined_HLA_Concordance_Coverage.pdf", final_plot, width = 12, height = 15)

# ==============================================================================
# 4. 绘制不同拆分结果的TSNE图
# ==============================================================================
# 加载函数并进行拆分
source(file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks/01_HLA_demultiplex_v2.R"))
sample_predict <- HLA_demultiplex(HLA_list = HLA_list, 
                                  HLA_matrix = HLA_matrix, 
                                  graph_corr = TRUE, 
                                  outdir = "Result/demultiplex")

# HLA-tsne result of true label and predict label
HLA_res_ONT = read.csv("./Result/HLA_Splitter_ONT/HLA_demultiplex.csv")
HLA_res_SBT = read.csv("./Result/HLA_Splitter_SBT/HLA_demultiplex.csv")
HLA_res_ONT$Cell.barcode = paste0(HLA_res_ONT$Cell.barcode, "-1")
HLA_res_SBT$Cell.barcode = paste0(HLA_res_SBT$Cell.barcode, "-1")

tsne_result$HLA_res_ONT = HLA_res_ONT$Sample_prediction[match(HLA_res_ONT$Cell.barcode, rownames(tsne_result))]
tsne_result$HLA_res_SBT = HLA_res_SBT$Sample_prediction[match(HLA_res_SBT$Cell.barcode, rownames(tsne_result))]

all_samples <- unique(c(as.character(tsne_result$HLA_res_ONT), as.character(tsne_result$HLA_res_SBT)))

valid_samples <- setdiff(all_samples, "Undefined")
color_mapping <- setNames(color[1:length(valid_samples)], valid_samples)
color_mapping["Undefined"] <- "grey80" 

# 绘制 TSNE 图并应用手动颜色映射
p1 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = HLA_res_SBT)) +
  geom_point(size = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  # 使用 scale_color_manual 强制应用我们的命名颜色字典
  scale_color_manual(values = color_mapping, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "HLA_res_SBT")

p2 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = HLA_res_ONT)) +
  geom_point(size = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  scale_color_manual(values = color_mapping, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = "HLA_res_ONT")

combined_plot <- p1 + p2
print(combined_plot)
ggsave("./Figures/HLA_TSNE_ONT_SBT.pdf", plot = combined_plot, dpi = 300, height = 5, width = 14)

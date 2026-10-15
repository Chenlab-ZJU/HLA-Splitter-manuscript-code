# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Upstream urine processing and sample-level signature analysis.
# Original source: 02_demultiplexing.R
# Figure/table scope: Fig4; FigS7; FigS8
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# File: 02_demultiplexing_and_doublets.R
# Description: 整合 HLA-Splitter 拆分结果，去除 Doublets
# ==============================================================================
source(file.path(getOption("hla.code_root"), "04_Urine_immune_programs/01_setup_and_helpers.R"))
# Read the saved object or source table used by this analysis.
scRNA_Intergrate <- readRDS(find_object_file("Urine_scRNA_Intergrate.rds", "Urine_ALL_Intergrate.rds"))

# 读取HLA样本拆分结果
samplenames <- ANALYSIS_BATCHES

list_1_field <- list()
list_2_field <- list()
counter_1 <- 1
counter_2 <- 1

for (batch in samplenames) {
  demux_files <- find_hla_demultiplex_files(batch)
  if (nrow(demux_files) == 0) {
    warning("No HLA demultiplex result found for batch: ", batch)
    next
  }

  for (row_id in seq_len(nrow(demux_files))) {
    tmp <- read.csv(demux_files$file[row_id])
    run_index <- demux_files$run_index[row_id]
    tmp$Cell.barcode <- paste0(str_sub(tmp$Cell.barcode, 1, 16), "-", batch, "-", run_index)

    if (demux_files$method[row_id] == "2-field") {
      list_2_field[[counter_2]] <- tmp
      counter_2 <- counter_2 + 1
    } else {
      list_1_field[[counter_1]] <- tmp
      counter_1 <- counter_1 + 1
    }
  }
}

if (length(list_1_field) == 0 && length(list_2_field) > 0) list_1_field <- list_2_field
if (length(list_2_field) == 0 && length(list_1_field) > 0) list_2_field <- list_1_field
if (length(list_1_field) == 0 || length(list_2_field) == 0) {
  stop("No usable HLA demultiplex CSV files were found under Data/Group*/HLA_counts.")
}
# 组合成完整的数据框
demultiplex_res_1_field <- bind_rows(list_1_field)
demultiplex_res_2_field <- bind_rows(list_2_field)
rownames(demultiplex_res_1_field) = demultiplex_res_1_field$Cell.barcode
rownames(demultiplex_res_2_field) = demultiplex_res_2_field$Cell.barcode

message("已成功合并 1-field 数据：", nrow(demultiplex_res_1_field), " 行")
message("已成功合并 2-field 数据：", nrow(demultiplex_res_2_field), " 行")

# ==============================================================================
# 2. 比较 1-field 和 2-field 精度的拆分结果
# ==============================================================================
comparison_df <- inner_join(
  demultiplex_res_1_field %>% select(Cell.barcode, prediction_1field = Sample_prediction),
  demultiplex_res_2_field %>% select(Cell.barcode, prediction_2field = Sample_prediction),
  by = "Cell.barcode"
)

# 计算两者的一致性
comparison_df <- comparison_df %>%
  mutate(is_consistent = (prediction_1field == prediction_2field))

# 统计总体一致率
consistency_summary <- comparison_df %>%
  summarise(
    Total_Cells = n(),
    Consistent_Cells = sum(is_consistent, na.rm = TRUE),
    Consistency_Rate = Consistent_Cells / Total_Cells
  )

print("=== 1-field 与 2-field 拆分结果一致性统计 ===")
print(consistency_summary)

# ==============================================================================
# 3. 绘制1-field 和 2-field 精度的拆分结果的柱状图
# ==============================================================================
counts_1f <- comparison_df %>%
  dplyr::count(prediction_1field, name = "Count") %>%
  mutate(Method = "1-field", Sample = prediction_1field) %>%
  select(Sample, Method, Count)

counts_2f <- comparison_df %>%
  dplyr::count(prediction_2field, name = "Count") %>%
  mutate(Method = "2-field", Sample = prediction_2field) %>%
  select(Sample, Method, Count)

df_counts <- bind_rows(counts_1f, counts_2f)

df_consistency <- comparison_df %>%
  filter(!grepl("Doublet|Negative", prediction_1field, ignore.case = TRUE)) %>%
  group_by(Sample = prediction_1field) %>%
  summarise(
    Total = n(),
    Consistent = sum(is_consistent, na.rm = TRUE),
    Consistency_Rate = Consistent / Total
  )

# 为了让 X 轴样本排序一致，统一设定 Factor levels
sample_levels <- mixedsort(unique(df_counts$Sample))
df_counts$Sample <- factor(df_counts$Sample, levels = sample_levels)
df_consistency$Sample <- factor(df_consistency$Sample, levels = sample_levels)

p_counts <- ggplot(df_counts, aes(x = Sample, y = Count, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", alpha = 0.85) +
  scale_fill_manual(values = c("1-field" = "#4DAF4A", "2-field" = "#699ECA")) +
  labs(
    title = "Cell Counts per Sample (1-field vs 2-field)",
    y = "Number of Cells"
  )+ my_theme()

p_rate <- ggplot(df_consistency, aes(x = Sample, y = Consistency_Rate)) +
  geom_col(fill = "#E42320", color = "black", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = percent(Consistency_Rate, accuracy = 0.1)), 
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.1)) + # Y轴留点空间给标签
  labs(
    title = "Demultiplexing Consistency Rate per Sample",
    x = "Predicted Sample (Baseline: 1-field)",
    y = "Consistency Rate"
  ) + my_theme()

combined_plot <- p_counts / p_rate + 
  plot_layout(heights = c(1, 1))

print(combined_plot)

ggsave(fig_file("Demultiplex", 1, "Counts_and_Consistency"), plot = combined_plot, width = 14, height = 6)


# ==============================================================================
# 4. 添加分组信息与保存
# ==============================================================================
demultiplex_result = demultiplex_res_1_field
rownames(demultiplex_result) <- demultiplex_result$Cell.barcode

scRNA_Intergrate$HLA_Splitter <- demultiplex_result[colnames(scRNA_Intergrate), "Sample_prediction"]
scRNA_Intergrate$HLA_Splitter[is.na(scRNA_Intergrate$HLA_Splitter)] <- "Undefined"

SampleInfo <- read.csv(data_path("metadata", "Sample_Info.csv"))

scRNA_Intergrate$Group <- SampleInfo$Group[match(scRNA_Intergrate$HLA_Splitter, SampleInfo$Patient.ID)]
scRNA_Intergrate$Group[is.na(scRNA_Intergrate$Group)] <- 'Undefined'
scRNA_Intergrate$Group <- SampleInfo$Group[match(scRNA_Intergrate$HLA_Splitter, SampleInfo$Patient.ID)]
scRNA_Intergrate$Sample.ID <- SampleInfo$Sample.ID[match(scRNA_Intergrate$HLA_Splitter, SampleInfo$Patient.ID)]
scRNA_Intergrate$Sample.ID[is.na(scRNA_Intergrate$Sample.ID)] <- 'Undefined'

Group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
scRNA = scRNA_Intergrate[,scRNA_Intergrate$Group %in% Group_levels]
scRNA$Group <- factor(scRNA$Group, levels = Group_levels)
scRNA$HLA_Splitter = factor(scRNA$HLA_Splitter, levels = mixedsort(unique(scRNA$HLA_Splitter)))
scRNA$Sample.ID = factor(scRNA$Sample.ID, levels = mixedsort(unique(scRNA$Sample.ID)))

# 保存结果
df <- data.frame(HLA_demutliplex = scRNA$Sample.ID)
df <- table(df)
# Export the analysis source data for figure reproduction.
write.csv(df, paste0("./Result/Demultiplex/demultiplex_res.csv"))

saveRDS(scRNA_Intergrate, "./Result/Objects/Urine_scRNA_Intergrate.rds")
saveRDS(scRNA, "./Result/Objects/Urine_scRNA.rds")

table(scRNA_Intergrate$orig.ident)
table(scRNA$orig.ident)
# ==============================================================================
# 5. 整合 HLA-Splitter 信息并绘制 tSNE
# ==============================================================================
library(tibble)
# 读取全局 HLA 列表
hla_genotype_file <- data_path("HLA_list", "HLA-genotypes_1field.csv")
if (!file.exists(hla_genotype_file)) {
  hla_genotype_file <- data_path("metadata", "HLA-genotypes.csv")
}
sampleInfo <- read.csv(hla_genotype_file, header = TRUE, row.names = 1)
HLA_list <- lapply(sampleInfo, unique)

# 收集并合并所有样本和批次的矩阵
HLA_matrix_list <- list()

for (samplename in samplenames) {
  matrix_sets <- find_hla_matrix_sets(samplename)
  if (nrow(matrix_sets) == 0) {
    message("No HLA matrix files found for: ", samplename)
    next
  }

  for (row_id in seq_len(nrow(matrix_sets))) {
    label_file <- matrix_sets$label_file[row_id]
    barcode_file <- matrix_sets$barcode_file[row_id]
    mtx_file <- matrix_sets$mtx_file[row_id]
    run_index <- matrix_sets$run_index[row_id]

    if (file.exists(label_file) && file.exists(barcode_file) && file.exists(mtx_file)) {
      HLA_label <- read.delim(label_file, header = FALSE)
      barcodes <- read.delim(barcode_file, header = FALSE)
      HLA_matrix <- as.data.frame(as.matrix(readMM(mtx_file)))
      colnames(HLA_matrix) <- HLA_label$V1
      rownames(HLA_matrix) <- paste0(str_sub(barcodes$V1, 1, 16), "-", samplename, "-", run_index)
      HLA_matrix_list[[paste0(samplename, "_", run_index)]] <- HLA_matrix
    } else {
      warning("Incomplete HLA matrix set skipped: ", dirname(mtx_file))
    }
  }
}
if (length(HLA_matrix_list) > 0) {
  list_with_rownames <- lapply(HLA_matrix_list, function(df) {
    as.data.frame(df) %>% rownames_to_column(var = "Barcode")
  })
  HLA_matrix_All <- bind_rows(list_with_rownames)
  HLA_matrix_All <- column_to_rownames(HLA_matrix_All, var = "Barcode")
  HLA_matrix_All[is.na(HLA_matrix_All)] <- 0
  
  message(sprintf("Successfully combined %d matrices. Total unique HLA alleles: %d", 
                  length(HLA_matrix_list), ncol(HLA_matrix_All)))
} else {
  stop("Error: No HLA matrices were found and loaded.")
}

# Filter Seurat object and match barcodes
message("Filtering Seurat object and matching barcodes...")
# Subset to remove Doublets (assuming 'Group' is in meta.data)
common_barcodes <- intersect(rownames(HLA_matrix_All), colnames(scRNA))
HLA_matrix <- HLA_matrix_All[common_barcodes, ]

message(sprintf("Matched %d barcodes.", length(common_barcodes)))

# Normalize HLA expression
message("Normalizing and scaling HLA counts...")
library_sizes <- rowSums(HLA_matrix)
HLA_matrix_norm <- HLA_matrix / library_sizes * 1e2
HLA_matrix_norm <- log(HLA_matrix_norm + 1)
HLA_matrix_norm[is.na(HLA_matrix_norm)] <- 0

# Scale HLA expression (positive values only per column)
for (i in seq_len(ncol(HLA_matrix_norm))) {
  exp_vals <- HLA_matrix_norm[, i]
  positive_idx <- exp_vals > 0
  if (sum(positive_idx) > 1) {
    exp_positive <- scale(exp_vals[positive_idx], center = FALSE, scale = TRUE)[, 1]
    HLA_matrix_norm[names(exp_positive), i] <- exp_positive
  }
}

# Prepare TSNE inputs
message("Running Rtsne...")
score_mat <- data.frame(HLA_Splitter_1field = demultiplex_res_1_field[common_barcodes,"Sample_prediction"],
                        HLA_Splitter_2field = demultiplex_res_2_field[common_barcodes,"Sample_prediction"])
namelist = mixedsort(unique(scRNA$HLA_Splitter))
score_mat$HLA_Splitter_1field <- factor(score_mat$HLA_Splitter_1field, levels = c(namelist))
score_mat$HLA_Splitter_2field <- factor(score_mat$HLA_Splitter_2field, levels = c(namelist))

score_mat$HLA_count <- log10(rowSums(HLA_matrix) + 1)

# Run TSNE
perplexity_val = min(max(10, floor((nrow(HLA_matrix_norm) - 1) / 3)), 50)
tsne_out = Rtsne(as.matrix(HLA_matrix_norm), dims = 2, pca = TRUE,
                 max_iter = 1000, theta = 0.4, perplexity = perplexity_val, 
                 verbose = FALSE, check_duplicates = FALSE)

tsne_result <- as.data.frame(tsne_out$Y)
colnames(tsne_result) <- c("tSNE1", "tSNE2")
tsne_result <- cbind(tsne_result, score_mat)


tsne_result$HLA_Splitter_1field_Sample <- SampleInfo$Sample.ID[match(tsne_result$HLA_Splitter_1field, SampleInfo$Patient.ID)]
tsne_result$HLA_Splitter_2field_Sample <- SampleInfo$Sample.ID[match(tsne_result$HLA_Splitter_2field, SampleInfo$Patient.ID)]

# TSNE plot
message("Generating plot...")
namelist = mixedsort(unique(scRNA$Sample.ID))
expanded_colors <- c(colorRampPalette(pal_igv("default")(51))(length(namelist)-1), "grey")
plot_colors <- setNames(expanded_colors, c(namelist))

# 分别计算 1-field 和 2-field 的中心点 (Centroids)
df_centroids_1 <- aggregate(cbind(tSNE1, tSNE2) ~ HLA_Splitter_1field, data = tsne_result, FUN = mean)
df_centroids_2 <- aggregate(cbind(tSNE1, tSNE2) ~ HLA_Splitter_2field, data = tsne_result, FUN = mean)

p1 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = HLA_Splitter_1field_Sample)) +
  geom_point(size = 0.1, alpha = 0.5) +
  #geom_text(data = df_centroids_1, aes(label = HLA_Splitter_1field), size = 4, color = "black", position = position_nudge(y = 2)) +
  theme_void() + 
  ggtitle("Demultiplex result using 1_field HLA genotypes") +
  theme(text = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Sample label")+
  scale_color_manual(values = plot_colors) # 传入命名向量，ggplot 会自动根据名字匹配颜色

p2 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = HLA_Splitter_2field_Sample)) +
  geom_point(size = 0.1, alpha = 0.5) +
  #geom_text(data = df_centroids_2, aes(label = HLA_Splitter_2field), size = 4, color = "black", position = position_nudge(y = 2)) +
  theme_void() +
  theme(text = element_text(size = 15)) +
  ggtitle("Demultiplex result using 2_field HLA genotypes") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = plot_colors)

combined_plot <- ggarrange(p1, p2, common.legend = T, legend = "right")
print(combined_plot) 

ggsave(fig_file("Demultiplex", 2, "HLA_TSNE_1field_2field"), plot = combined_plot, width = 14, height = 6)


# ==============================================================================
# 6. 统计拆分样本的细胞数量
# ==============================================================================
df_plot <- data.frame(predict_label = scRNA$Sample.ID, Group = scRNA$Group) %>%
  dplyr::count(Group, predict_label, name = "Cellnum")
group_list = levels(scRNA_filter$Group)
thresh = 50

p3 <- ggplot(data = df_plot, aes(x = predict_label, y = Cellnum, fill = Group)) +
  geom_col(alpha = 0.8, width = 0.5) + 
  geom_text(aes(label = Cellnum), vjust = -0.5, size = 4, color = "black") + 
  geom_hline(yintercept = thresh, linetype = "dashed", color = "darkblue", linewidth = 0.8) + # Added dashed line at y=100
  theme_classic() + scale_y_log10() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), # Rotate x-axis labels for readability
    axis.text.y = element_text(size = 12)) +
  scale_fill_manual(values = GROUP_COLORS) + 
  ylab("Number of Cells") + # Clearer y-axis label
  xlab("Predicted Sample Label")

# Save plot
ggsave(fig_file("Demultiplex", 3, "Sample_Prediction_Cell_Count"), plot = p3, width = 14, height = 6)

# ==============================================================================
# 7. 统计不同分组中总细胞量和免疫细胞量的比例
# ==============================================================================
SampleInfo <- read.csv(data_path("metadata", "Sample_Info.csv"))

df <- data.frame(
  Sample = SampleInfo$Sample.ID,
  Group = SampleInfo$Group,
  Rejection = SampleInfo$Rejection,
  TotalCell = SampleInfo$TotalCell * 1e6, # Convert million cells to absolute count
  PTD = SampleInfo$PTD
)

# 过滤掉无法分配的样本
df <- df %>% filter(Sample != "/")
df <- df %>%
  left_join(df_plot %>% select(predict_label, Cellnum),
            by = c("Sample" = "predict_label")) %>%
  dplyr::rename(ImmuneCell = Cellnum)

df$ImmuneCell[is.na(df$ImmuneCell)] <- 0

df <- df %>% mutate(CellRatio = ImmuneCell / TotalCell)
df <- df[mixedorder(df$Sample),]
df$Group <- factor(df$Group, levels = c("NR", "TCMR", "ABMR", "Mix", "BK"))

p_scatter <- ggplot(df, aes(x = TotalCell, y = ImmuneCell, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = Group), alpha = 0.15) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,5000000)+ ylim(0,1000) +
  labs(
    title = "Correlation between Total Cells and Immune Cells",
    x = "Total Cells in Urine (Absolute Count)",
    y = "Immune Cells (scRNA-seq Count)",
    color = "Clinical Group",
    fill = "Clinical Group"
  ) +
  scale_color_manual(values = GROUP_COLORS) +
  scale_fill_manual(values = GROUP_COLORS) +
  my_theme()

print(p_scatter)
ggsave(fig_file("Demultiplex", 4, "Immune_vs_Total_Cell_LinearFit"), plot = p_scatter, dpi = 300, height = 5, width = 7)


#保存
saveRDS(scRNA_Intergrate, "./Result/Objects/Urine_scRNA_Intergrate.rds")

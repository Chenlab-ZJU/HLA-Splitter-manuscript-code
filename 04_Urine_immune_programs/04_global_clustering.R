# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Upstream urine processing and sample-level signature analysis.
# Original source: 03_global_clustering.R
# Figure/table scope: Fig4; FigS7; FigS8
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# File: 03_global_clustering.R
# Description: 全局降维聚类, 细胞群注释及 T 细胞的再次分群
# ==============================================================================
source(file.path(getOption("hla.code_root"), "04_Urine_immune_programs/01_setup_and_helpers.R"))
# Read the saved object or source table used by this analysis.
scRNA <- readRDS(find_object_file("Urine_scRNA.rds", "Urine_ALL.rds"))

# ==============================================================================
# 1. 数据筛选与重新聚类
# ==============================================================================
scRNA <- subset(scRNA, subset = !HLA_Splitter %in% c("Doublets", "Undefined"))

SampleInfo <- read.csv(data_path("metadata", "Sample_Info.csv"))
Sample_levels <- c(mixedsort(unique(scRNA$HLA_Splitter)))

# 剔除未定义和双细胞组别
scRNA$Group <- factor(scRNA$Group, levels = setdiff(unique(scRNA$Group), c("Undefined", "Doublet")))

# 过滤细胞数小于 50 的样本
valid_samples <- Sample_levels[table(scRNA$HLA_Splitter) > 50]
scRNA <- subset(scRNA, subset = HLA_Splitter %in% valid_samples)


# ==============================================================================
# 2. 全局细胞聚类与注释
# ==============================================================================
scRNA <- scRNA %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = "orig.ident") %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, min.dist = 0.3) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)

# ==============================================================================
# 3. 全局 Marker 提取与可视化
# ==============================================================================
export_top_markers(scRNA, "seurat_clusters", "scRNA_all")

# 注意：运行到此处请暂停，确保 ./Result/scRNA/Markers/scRNA_all_makers_summary.csv 已填写 Celltype
df <- read.csv("./Result/scRNA/Markers/scRNA_all_makers_summary.csv")

# 大类注释 (Major_celltype)
Major_celltype <- df$Major_celltype
names(Major_celltype) <- as.character(df$seurat_clusters)

scRNA <- SetIdent(scRNA, value = "seurat_clusters")
scRNA <- RenameIdents(scRNA, Major_celltype)
scRNA$Major_celltype <- Idents(scRNA)

message("当前的大类细胞包括：")
print(unique(scRNA$Major_celltype))

# 剔除杂细胞并排序 (注: Contament 拼写保持原样，若为 Contaminant 可自行修改)
scRNA <- subset(scRNA, Major_celltype != "Contament")

# 转换为因子并固定顺序
scRNA$Major_celltype <- factor(
  scRNA$Major_celltype, 
  levels = c("T/NK", "B cell", "Myeloid", "DC")
)

scRNA = RunUMAP(scRNA, reduction = "harmony", dims = 1:30, min.dist = 0.3)

# 整合细分亚群 (Celltype)
scRNA$Celltype <- as.character(scRNA$Major_celltype)

# 假设 Tcell 和 Myeloid 对象已在环境中（或由前面代码生成）
scRNA@meta.data[colnames(Tcell), "Celltype"] <- as.character(Tcell$Celltype)
scRNA@meta.data[colnames(Myeloid), "Celltype"] <- as.character(Myeloid$Celltype)
scRNA <- subset(scRNA, !Celltype %in% c("NA","T cell","Myeloid"))

message("整合细分亚群后的细胞包括：")
print(unique(scRNA$Celltype))

# 转换为因子并固定顺序
scRNA$Celltype <- factor(
  scRNA$Celltype, 
  levels = c("Naive T", "CD4+ Trm","CD4+ Tph", "Tregs", "Th17", "CD8+ Trm", "CD8+ Tact", "ISG+ T","Proliferating T",
             "B cell","NK-like",
             "Monocyte", "Macro_C1QA", "Macro_CXCL9", "Macro_CCL3",
             "Macro_TREM2","Macro_TXNRD1", "Macro_NLRP3", 
             "Macro_MT", "DC")
)

# ==============================================================================
# 4. 绘制全局可视化图表 (UMAP、FeaturePlot、DotPlot等)
# ==============================================================================
# 全局 UMAP 三联图
p1 <- DimPlot(scRNA, group.by = "Major_celltype", label = T, alpha = 0.5) + 
  scale_color_manual(values = GLOBAL_COLORS)+ my_theme()
p2 <- DimPlot(scRNA, group.by = "Celltype", label = FALSE, alpha = 0.5) + 
  scale_color_manual(values = GLOBAL_COLORS)+ my_theme()
p3 <- DimPlot(scRNA, group.by = "Group", label = FALSE, alpha = 0.5) + 
  scale_color_manual(values = GROUP_COLORS)+ my_theme()

ggsave(fig_file("scRNA/UMAP_Markers", 1, "UMAP_Major_Subtype_Group"), p1 + p2 + p3, height = 5, width = 18)

# 关键基因 FeaturePlot
p_F <- FeaturePlot(scRNA, features = c("CXCR3", "CXCR4", "GZMB", "CD27"), ncol = 2, cols = c("grey", "firebrick"))
ggsave(fig_file("scRNA/UMAP_Markers", 2, "FeaturePlot_CXCR3_CXCR4_GZMB_CD27"), p_F, height = 8, width = 8)

# 全局 Marker 基因气泡图
global_mk <- c("CCR7", "LEF1", "CD27", "IGFL2",
               "KLRB1","IL21", "CD40LG", "FOXP3", "IL2RA",
               "CCR6", "IL23R", "GZMK", "GZMB" ,
               "CXCR4", "IFNG","IFIT3", "ISG15","MKI67", "TYMS",
               "MS4A1", "KLRD1",
               "S100A8", "S100A9", "C1QA", "C1QB", "CXCL9", "CXCL10", 
               "CCL3", "TREM2", "APOE","GNG11", "TXNRD1", 
               "IL1B", "NLRP3", "MT1A", "MT1G", "CD1C", "FCER1A")

p_dot <- DotPlot(scRNA, group.by = "Celltype", features = global_mk) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = NULL, y = NULL) +
  scale_color_gradientn(colours = c("navy", "white", "firebrick3"))

ggsave(fig_file("scRNA/UMAP_Markers", 3, "DotPlot_Celltype_Marker_Genes"), plot = p_dot, height = 5, width = 12)

# UMAP 分面图
combined_plot <- ggarrange(
  plotlist = plot_grouped_umap(scRNA, "Group"), 
  ncol = 5, common.legend = TRUE, legend = "right"
)
ggsave(fig_file("scRNA/UMAP_Markers", 4, "Grouped_UMAP_by_Rejection_Group"), combined_plot, height = 5, width = 20)

# 不同 Group 的细胞类型比例柱状图
p4 <- ggplot(scRNA@meta.data, aes(x = Group, fill = Major_celltype)) +
  geom_bar(position = "fill", width = 0.8, color = "black", linewidth = 0, alpha = 0.5) + # position="fill" 自动计算百分比
  scale_fill_manual(values = GLOBAL_COLORS) +
  scale_y_continuous(labels = scales::percent) + # 将 y 轴转换为百分比格式
  labs(x = "Group", y = "Cell Proportion", fill = "Major Celltype") +
  my_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 倾斜 x 轴标签防止文字重叠
    panel.grid = element_blank(), # 移除网格线
    axis.line = element_line(color = "black")
  )

ggsave(fig_file("scRNA/Celltype_Composition", 1, "Major_Celltype_Proportion_by_Group"), p4, height = 5, width = 4)


# 细胞构成比例与热图
comparison <- list(c("NR", "TCMR"), c("NR", "ABMR"), c("NR", "Mix"))

p_list <- plot_celltype_ratio(
  seurat_obj = scRNA,
  sample_col = "Sample.ID", 
  group_col = "Group", 
  cluster_col = "Celltype",
  comparison_groups = comparison, 
  group_colors = GROUP_COLORS, 
  cluster_colors = GLOBAL_COLORS
)

ggsave(fig_file("scRNA/Celltype_Composition", 2, "Celltype_Stacked_Bar_by_Sample"), p_list$stacked_bar, width = 10, height = 3)
ggsave(fig_file("scRNA/Celltype_Composition", 3, "Celltype_Proportion_Heatmap"), p_list$heatmap, width = 5, height = 5)

# 细胞亚群比例图
p_cell_ratio <- plot_cell_ratio_comparison(
  cell_subset = scRNA,
  sample_col = "HLA_Splitter", 
  cluster_col = "Celltype", 
  group_col = "Group",
  facet = TRUE, 
  comparison_groups = comparison, 
  write_to_csv = "./Result/scRNA_celltype_ratio", 
  plot_colors = GLOBAL_COLORS
)

ggsave(fig_file("scRNA/Celltype_Composition", 4, "Celltype_Proportion_Group_Comparison"), p_cell_ratio, width = 8, height = 8)
# Export the analysis source data for figure reproduction.
saveRDS(scRNA, "./Result/Objects/Urine_scRNA.rds")

# ==============================================================================
# 5. Pseudobulk 分析: 输出 Lollipop 图 + 分组专属 CSV + Bar&Dot 图
# ==============================================================================
source_if_exists(file.path(HELPER_DIR, "plot_top_pseudobulk.R"), required = TRUE)

# 设定新的输出路径，避免覆盖之前的疾病分组结果
out_dir <- "./Result/scRNA/Group_Markers"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cell_subset = scRNA

plot_top_pseudobulk(
  seurat_obj = cell_subset, 
  group_by = "Group",  
  sample_col = "HLA_Splitter",
  top_n_lollipop = 10, 
  top_n_bardot = 10, 
  my_colors = GROUP_COLORS, 
  out_dir = out_dir
)

# 合并各亚群的 Top DEG
Group_DEG <- data.frame()
unique_groups <- as.character(unique(cell_subset$Group))

for (grp in unique_groups) {
  csv_path <- paste0(out_dir, "/Marker_Stats_", grp, ".csv")
  if (file.exists(csv_path)) { 
    top_gene_tmp <- read.csv(csv_path)
    Group_DEG <- rbind(Group_DEG, top_gene_tmp)
  }
}
write.csv(Group_DEG, paste0(out_dir, "/Top_DEGs.csv"), row.names = FALSE)

# ==============================================================================
# 6. 绘制全局 Pseudobulk 热图
# ==============================================================================
message("Generating global pseudobulk heatmap...")

# 自动提取每个细胞类型 log2FC 最高的基因
top_heatmap_genes <- Group_DEG %>%
  dplyr::group_by(Group_Col) %>%
  filter(PB_Pval<0.05)  %>%
  dplyr::slice_max(order_by = PB_Log2FC, n = 8, with_ties = FALSE) %>%
  dplyr::pull(gene) %>%
  unique()

top_heatmap_genes <- c("GNLY", "NKG7", "CCL5", "GZMA", "CTSW", 
                       "HLA-DRB5", "HLA-DQA2", "APOE", "APOC1", 
                       "MT1H", "MT1G", "CXCL10", "CXCL9", "CXCL8", 
                       "CCL2", "CCL7", "CCL3L1", "CCL4", "S100A8", 
                       "S100A9", "S100A11", "S100A12", "IL1B", "IL1R2", 
                       "IFITM2", "IFITM3", "SERPINB2", "THBS1", "F3", 
                       "CXCL3", "CCL3", "GZMH", "AREG", "HSPA6", "HSPD1", 
                       "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA8", "MS4A1", 
                       "CXCR4", "HLA-DRA", "HLA-DPA1", "HLA-DQA1")

top_heatmap_genes <- c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1",
                       "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "CD74",
                       "CXCL9", "CXCL10", "CCL2", "CCL3", "CCL4", "CCL5", "CCL7",
                       "NKG7", "GNLY", "GZMA", "GZMB", "CST7", "CTSW", "IFNG",
                       "CD3E", "CD69", "LAG3", "CTLA4",
                       "MS4A1", "CD79A", "IGHM", "BANK1",
                       "S100A8", "S100A9", "S100A12", "IL1B", "CXCL8", "HMOX1",
                       "PLAUR", "TYMP", "BCL2A1")

cell_subset_avg <- AverageExpression(cell_subset, group.by = "Group", assays = "RNA", return.seurat = TRUE)
cell_subset_avg <- ScaleData(cell_subset_avg, features = rownames(cell_subset_avg))
gene_matrix_sc <- cell_subset_avg@assays$RNA$scale.data %>% as.matrix()

genes_to_plot <- intersect(top_heatmap_genes, rownames(gene_matrix_sc))

library(pheatmap)
plot_heatmap_sc_bulk <- pheatmap(
  gene_matrix_sc[genes_to_plot, ], 
  cluster_rows = T, 
  cluster_cols = FALSE, 
  angle_col = 45, 
  cellwidth = 40, 
  cellheight = 12 
)

ggsave(fig_file("scRNA/Pseudobulk", 1, "Group_DEG_Heatmap"), plot_heatmap_sc_bulk, height = 20, width = 8)

# ==============================================================================
# 7. 循环遍历每个细胞类型：绘制专属小提琴图和富集分析
# ==============================================================================
valid_groups <- as.character(unique(Group_DEG$Group_Col))

for (grp in valid_groups) {
  message("\n>> 正在执行亚群精细分析: ", grp)
  
  # 绘制小提琴图
  plot_top_violins(
    seu_obj = cell_subset,             
    deg_data = Group_DEG,
    target_group = grp,
    out_dir = out_dir,
    comparisons_list = comparison,   
    plot_colors = GROUP_COLORS,
    group_col = "Group_Col",
    lfc_col = "PB_Log2FC",
    vln_group_by = "Group",
    lfc_threshold = 0.1,
    top_n = 8
  )
  
  # 运行 GO 富集并绘图
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
# 8. 核心代表性基因的单细胞表达与 Pseudo-bulk 样本级验证 (组合图优化版)
# ==============================================================================
target_genes <- c("CXCL9", "CXCL10","CCL5", "CCL7", 
                  "HLA-DRA", "HLA-DPB1", "CD74",
                  "NKG7", "GNLY", "GZMB", "GZMA",
                  "CD3E")


sample_col <- "HLA_Splitter"
group_col <- "Group"

# 1. 检查基因是否存在
valid_genes <- intersect(target_genes, rownames(cell_subset))
if (length(valid_genes) == 0) stop("所有目标基因均未在表达矩阵中找到，请检查基因名。")

# 2. 提取数据并清理 NA
expr_data <- FetchData(cell_subset, vars = c(valid_genes, group_col, sample_col))
expr_data <- expr_data[!is.na(expr_data[[group_col]]), ]

# ------------------------------------------------------------------------------
# 8.1 绘制组合单细胞小提琴图 (Single-cell Level Combined)
# ------------------------------------------------------------------------------
sc_data_long <- expr_data %>%
  pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = valid_genes)) 

p_vln_combined <- ggplot(sc_data_long, aes(x = .data[[group_col]], y = Expression)) +
  geom_violin(aes(fill = .data[[group_col]], color = .data[[group_col]]), 
              alpha = 0.2, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", color = "black", 
               outlier.shape = NA, linewidth = 0.3) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = GROUP_COLORS) +
  scale_color_manual(values = GROUP_COLORS) +
  theme_classic(base_size = 14) +
  stat_compare_means(
    comparisons = comparison, 
    method = "wilcox.test", 
    label = "p.signif", 
    step.increase = 0.1, 
    tip.length = 0.02, 
    bracket.size = 0.3
  ) +
  labs(title = "Single-cell Normalized Expression", x = NULL, y = "Expression Level") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    strip.text = element_text(face = "bold.italic", size = 14), 
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  )

ggsave(fig_file("scRNA/Pseudobulk", 2, "Core_Genes_Single_Cell_Violin"), p_vln_combined, width = 14, height = 12)

# ------------------------------------------------------------------------------
# 8.2 绘制组合 Pseudo-bulk 样本级柱状图 (Sample Level Combined)
# ------------------------------------------------------------------------------
pb_data_long <- expr_data %>%
  group_by(!!sym(group_col), !!sym(sample_col)) %>%
  summarise(across(all_of(valid_genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Mean_Expression") %>%
  mutate(Gene = factor(Gene, levels = valid_genes))

p_bar_combined <- ggplot(pb_data_long, aes(x = .data[[group_col]], y = Mean_Expression)) +
  # 柱状图：计算均值 (修复了旧版 size 参数警告，使用 linewidth)
  stat_summary(fun = mean, geom = "bar", 
               aes(fill = .data[[group_col]], color = .data[[group_col]]),
               alpha = 0.4, width = 0.8, linewidth = 0.6) +
  
  # 误差棒：直接使用匿名函数避免作用域找不到变量报错，修复 size 警告
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
      data.frame(y = m, ymin = m, ymax = m + se)
    }, 
    geom = "errorbar", 
    aes(color = .data[[group_col]]), 
    width = 0.2, linewidth = 0.7
  ) +
  
  # 散点抖动展示具体样本点
  geom_jitter(aes(color = .data[[group_col]]), width = 0.18, size = 1.5, alpha = 0.8) +
  
  # 分面及排版样式
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = GROUP_COLORS) +
  scale_color_manual(values = GROUP_COLORS) +
  theme_classic(base_size = 14) +
  stat_compare_means(
    comparisons = comparison, 
    method = "wilcox.test", 
    label = "p.signif", 
    step.increase = 0.12, 
    tip.length = 0.02, 
    bracket.size = 0.3
  ) +
  labs(title = "Pseudo-bulk Mean Expression per Sample", x = NULL, y = "Mean Expression") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    strip.text = element_text(face = "bold.italic", size = 14), 
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  )

ggsave(fig_file("scRNA/Pseudobulk", 3, "Core_Genes_Pseudobulk_Barplot"), p_bar_combined, width = 8, height = 6)


# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Upstream urine processing and sample-level signature analysis.
# Original source: 09_Pseudobulk_analysis.R
# Figure/table scope: Fig4; FigS7; FigS8
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# File: 04_group_marker_analysis.R
# Description: 灵活交互式的 Pseudobulk 分析、差异基因提取、富集与组合图绘制
# ==============================================================================
source(file.path(getOption("hla.code_root"), "04_Urine_immune_programs/01_setup_and_helpers.R"))
source_if_exists(file.path(HELPER_DIR, "plot_top_pseudobulk.R"), required = TRUE)

if (!exists("scRNA")) scRNA <- readRDS(find_object_file("Urine_scRNA.rds", "Urine_ALL.rds"))
if (!exists("Tcell")) Tcell <- readRDS(find_object_file("Urine_Tcell.rds"))
if (!exists("Myeloid")) Myeloid <- readRDS(find_object_file("Urine_Myeloid.rds"))

# ==============================================================================
# 【配置区】更换分析对象时，仅需修改这里的参数，下方代码完全无需改动
# ==============================================================================
current_obj   <- scRNA 
analysis_name <- "scRNA"     # [输出] 文件夹及命名前缀 (例如换成 "Tcell_Subset")

group_col     <- "Group"            # 分组所在的 meta.data 列名
sample_col    <- "HLA_Splitter"     # 样本所在的 meta.data 列名

my_colors     <- GROUP_COLORS       # 当前图表使用的调色盘
my_comparisons <- list(c("NR","TCMR"), c("NR","ABMR"), c("NR","Mix")) # 统计比较组合
# ==============================================================================
# 1. Pseudobulk分析: Lollipop 图 + 分组 CSV + Bar&Dot 图
# ==============================================================================
# 动态生成专属输出路径，避免相互覆盖
out_dir <- result_path(analysis_name, "Group_Markers")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_top_pseudobulk(
  seurat_obj = current_obj, 
  group_by = group_col,  
  sample_col = sample_col,
  top_n_lollipop = 10, 
  top_n_bardot = 10, 
  my_colors = my_colors, 
  out_dir = out_dir
)

# 合并各亚群的 Top DEG
Group_DEG <- data.frame()
unique_groups <- as.character(unique(current_obj@meta.data[[group_col]]))

for (grp in unique_groups) {
  csv_path <- file.path(out_dir, paste0("Marker_Stats_", grp, ".csv"))
  if (file.exists(csv_path)) { 
    top_gene_tmp <- read.csv(csv_path)
    Group_DEG <- rbind(Group_DEG, top_gene_tmp)
  }
}
# Export the analysis source data for figure reproduction.
write.csv(Group_DEG, file.path(out_dir, paste0(analysis_name, "_Top_DEGs.csv")), row.names = FALSE)

# ==============================================================================
# 2. 绘制全局 Pseudobulk 热图
# ==============================================================================
message("Generating global pseudobulk heatmap for: ", analysis_name)

# 自动提取每个细胞类型 log2FC 最高的基因
top_heatmap_genes <- Group_DEG %>%
  dplyr::group_by(Group_Col) %>%
  dplyr::slice_max(order_by = PB_Log2FC, n = 5, with_ties = FALSE) %>%
  dplyr::pull(gene) %>%
  unique()


# 自定义基因(scRNA)
top_heatmap_genes <- c(
  "GNLY", "NKG7", "MT1M", "MT1H", "MT1G",
  "GZMB", "CXCL9", "CXCL10", "ICAM1",
  "S100A8", "C1QA", "STAT1", "IFITM1", "HMOX1","IFNG","CD3E", "CD7",
  "STAT2", "IRF7", "PDCD1", "CXCL8", "IFIT2", "CXCL2", "CXCL3", "ICAM1"
)

# 自定义基因(Tcell)
top_heatmap_genes <- c(
  "CD4", "CD8A","CXCR3", "CD27", "STAT1",
  "TNF", "IFNG", "IL21", "CCL4", 
  "GZMK", "GZMB", "NKG7",
  "HLA-A","HLA-B" ,"HLA-C","HLA-DRA","HLA-DPA","HLA-DQA",
  "HAVCR2","TIGIT",  "PDCD1", "CTLA4"
)

# 自定义基因(Myeloid)
top_heatmap_genes <- c("C1QC", "NLRP3", "TREM2", "SPP1", "APOC1", "APOE",
                  "HLA-A","HLA-B" ,"HLA-C","HLA-DRA","HLA-DPA","HLA-DQA",
                  "CXCL9", "CXCL10", "IFNG", "TNF","CCL3",  "CCL5")


seu_avg <- AverageExpression(current_obj, group.by = group_col, assays = "RNA", return.seurat = TRUE)
seu_avg <- ScaleData(seu_avg, features = rownames(seu_avg))
gene_matrix_sc <- seu_avg@assays$RNA$scale.data %>% as.matrix()

genes_to_plot <- intersect(top_heatmap_genes, rownames(gene_matrix_sc))

library(pheatmap)
plot_heatmap_sc_bulk <- pheatmap(
  gene_matrix_sc[genes_to_plot, ], 
  cluster_rows = TRUE, 
  cluster_cols = FALSE, 
  angle_col = 45, 
  cellwidth = 40, 
  cellheight = 12 
)

ggsave(fig_file(file.path(analysis_name, "Pseudobulk"), 1, "Group_DEG_Heatmap"), plot_heatmap_sc_bulk, height = 12, width = 8)

# ==============================================================================
# 3. 循环遍历每个细胞类型：专属小提琴图和 GO 富集分析
# ==============================================================================
valid_groups <- as.character(unique(Group_DEG$Group_Col))

for (grp in valid_groups) {
  message("\n>> 正在执行亚群精细分析: ", grp)
  
  # 绘制小提琴图
  plot_top_violins(
    seu_obj = current_obj,             
    deg_data = Group_DEG,
    target_group = grp,
    out_dir = out_dir,
    comparisons_list = my_comparisons,   
    plot_colors = my_colors,
    group_col = "Group_Col",
    lfc_col = "PB_Log2FC",
    vln_group_by = group_col,
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
# 4. 核心代表性基因的组合图 (Single-cell & Pseudo-bulk)
# ==============================================================================
# 核心靶基因 (scRNA)
target_genes <- c(
  # ABMR (内皮损伤、NK细胞、趋化)
  "NKG7", "MT1H","CXCL9", "CXCL10",
  "GNLY", "GZMB", "ICAM1",
  # TCMR (巨噬细胞/T细胞浸润、组织损伤)
  "S100A8", "C1QA", "STAT1", "IFITM1", "HMOX1","IFNG","CD3E", "CD7",
  # BK肾病 (抗病毒、耗竭、慢性炎症)
  "STAT2", "IRF7", "PDCD1", "CXCL8", "IFIT2", "CXCL2", "CXCL3", "ICAM1"
)

# 核心靶基因 (Tcell)
target_genes <- c(
  "CXCR3", "CD27", "STAT1", "TNF",
  "IL21", "HAVCR2", "PDCD1", "CTLA4"
)

# 核心靶基因 (Myeloid)
target_genes <- c("TREM2", "SPP1", "APOE",
                  "CXCL9", "CXCL10", "IFNG")

valid_genes <- intersect(target_genes, rownames(current_obj))
if (length(valid_genes) == 0) stop("所有目标基因均未在表达矩阵中找到，请检查基因名或输入对象。")

expr_data <- FetchData(current_obj, vars = c(valid_genes, group_col, sample_col))
expr_data <- expr_data[!is.na(expr_data[[group_col]]), ]

# ------------------------------------------------------------------------------
# 4.1 单细胞小提琴组合图 (Single-cell Level Combined)
# ------------------------------------------------------------------------------
sc_data_long <- expr_data %>%
  pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = valid_genes)) 

p_vln_combined <- plot_violin_comparison(seurat_obj = current_obj, 
                                         group_col = "Group", 
                                         comparison_groups = comparison,
                                         features_to_plot = valid_genes, 
                                         plot_colors = GLOBAL_COLORS)

ggsave(fig_file(file.path(analysis_name, "Pseudobulk"), 2, "Core_Genes_Single_Cell_Violin"), p_vln_combined, width = 12, height = 6)

# ------------------------------------------------------------------------------
# 4.2 Pseudo-bulk 柱状组合图 (Sample Level Combined)
# ------------------------------------------------------------------------------
pb_data_long <- expr_data %>%
  group_by(!!sym(group_col), !!sym(sample_col)) %>%
  summarise(across(all_of(valid_genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Mean_Expression") %>%
  mutate(Gene = factor(Gene, levels = valid_genes))

p_bar_combined <- ggplot(pb_data_long, aes(x = .data[[group_col]], y = Mean_Expression)) +
  stat_summary(
    fun = mean, geom = "bar", 
    aes(fill = .data[[group_col]], color = .data[[group_col]]),
    alpha = 0.5, width = 0.65, linewidth = 0.6
  ) +
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      data.frame(y = m, ymin = m, ymax = m + sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    },
    geom = "errorbar", aes(color = .data[[group_col]]), width = 0.2, linewidth = 0.7
  ) +
  geom_jitter(aes(color = .data[[group_col]]), width = 0.18, size = 1.5, alpha = 0.8) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_color_manual(values = my_colors, aesthetics = c("color", "fill")) +
  theme_classic(base_size = 14) +
  stat_compare_means(
    comparisons = my_comparisons, method = "wilcox.test", 
    label = "p.signif", step.increase = 0.12, tip.length = 0.02, bracket.size = 0.3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
  labs(
    title = paste0(analysis_name, " - Pseudo-bulk Mean Expression"), 
    x = NULL, 
    y = "Mean Expression"
  ) +
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

ggsave(fig_file(file.path(analysis_name, "Pseudobulk"), 3, "Core_Genes_Pseudobulk_Barplot"), p_bar_combined, width = 8, height = 12)

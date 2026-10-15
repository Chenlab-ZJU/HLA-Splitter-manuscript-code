# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis.
# Original source: 05_Tcell_subclustering.R
# Figure/table scope: Fig5; FigS9; Supplementary Table 18; Methods CXCR3-axis analysis
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# File: 04_celltype_subclustering.R
# Description: Myeloid和Tcell亚群重新聚类、Marker发掘、M1/M2极化评分
# ==============================================================================
source(file.path(getOption("hla.code_root"), "04_Urine_immune_programs/01_setup_and_helpers.R"))
#scRNA <- readRDS("./Result/Objects/Urine_scRNA.rds")
# Read the saved object or source table used by this analysis.
Tcell = readRDS(file.path(getOption("hla.data_root"), "urine", "Urine_Tcell.rds"))

# ======================= B. T细胞 (Tcell) 分析 =======================
Tcell = scRNA[,scRNA$Major_celltype %in% c("T cell", "NK-like", "Cycling")]

Tcell <- Tcell %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 800) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = "orig.ident") %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, min.dist = 0.3) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)

# 全局 Marker 提取与可视化 --------------------------------------------------
export_top_markers(Tcell, "seurat_clusters", "Tcell")

# 注意：运行到此处请暂停，确保 ./Result/Tcell/Markers/Tcell_makers_summary.csv 已填写 Celltype
df <- read.csv("./Result/Tcell/Markers/Tcell_makers_summary.csv")
Celltype <- df$Celltype
names(Celltype) <- as.character(df$seurat_clusters)
Tcell <- SetIdent(Tcell, value = "seurat_clusters")
Tcell <- RenameIdents(Tcell, Celltype)
Tcell$Celltype <- Idents(Tcell)
unique(Tcell$Celltype) %>% print()

#除去杂细胞
Tcell = Tcell[,!Tcell$Celltype %in% c("Contament", "T_stress", "T_MT")]
Tcell$Celltype <- factor(Tcell$Celltype, levels = c("NK","CD4_Tcm","T_Quiescent","CD4_Tph","Treg","T_Activated",
                                                    "CD8_Teff","CD8_Trm","T_ISG", "T_cycling"))
Tcell = RunUMAP(Tcell, reduction = "harmony", dims = 1:30, min.dist = 0.3)

# 4. 绘制全局UMAP三联图 --------------------------------------------------------
p1 <- DimPlot(Tcell, group.by = "seurat_clusters", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p2 <- DimPlot(Tcell, group.by = "Celltype", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p3 <- DimPlot(Tcell, group.by = "Group", alpha = 0.5) + scale_color_manual(values = GROUP_COLORS)
ggsave(fig_file("Tcell/UMAP_Markers", 1, "UMAP_Cluster_Celltype_Group"), p1 + p2 + p3, height = 5, width = 16)

p_F = FeaturePlot(Tcell, features = c("CXCR3", "CXCR4", "IFNG", "FCGR3A"),ncol=2,cols = c("grey","firebrick"))
ggsave(fig_file("Tcell/UMAP_Markers", 2, "FeaturePlot_CXCR3_CXCR4_IFNG_FCGR3A"), p_F, height = 8, width = 8)

FeaturePlot(Tcell, features = c("CD3D", "KLRD1"),ncol=2,cols = c("grey","firebrick"))
# 全局 Marker 基因绘图
global_mk <- c("CCR7", "LEF1", "CD27", "IGFL2",
               "KLRB1","IL21", "CD40LG", "FOXP3", "IL2RA",
               "CCR6", "IL23R", "GZMK", "GZMB" ,
               "CXCR4", "IFNG","IFIT3", "ISG15","MKI67", "TYMS","KLRD1","FCGR3A")
p_dot <- DotPlot(Tcell, group.by = "Celltype", features = global_mk) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x=NULL, y=NULL) +
  scale_color_gradientn(colours = c("navy","white","firebrick3"))
ggsave(fig_file("Tcell/UMAP_Markers", 3, "DotPlot_Celltype_Marker_Genes"), plot = p_dot, height = 4, width = 9)

# UMAP 分面图 & 细胞比例直方图
combined_plot <- ggarrange(plotlist = plot_grouped_umap(Tcell, "Celltype", "Group", GLOBAL_COLORS), 
                           common.legend = TRUE, legend = "right")
ggsave(fig_file("Tcell/UMAP_Markers", 4, "Grouped_UMAP_by_Rejection_Group"), combined_plot, height = 6, width = 8)

comparison <- list(c("NR","TCMR"),c("NR","ABMR"),c("NR","Mix"))
p_list = plot_celltype_ratio(seurat_obj = Tcell,
                             sample_col = "Sample.ID",group_col = "Group",
                             cluster_col = "Celltype",
                             comparison_groups = comparison,
                             group_colors = GROUP_COLORS, 
                             cluster_colors = GLOBAL_COLORS)

ggsave(fig_file("Tcell/Celltype_Composition", 1, "Celltype_Stacked_Bar_by_Sample"), p_list$stacked_bar, width = 8, height = 4)
ggsave(fig_file("Tcell/Celltype_Composition", 2, "Celltype_Proportion_Heatmap"), p_list$heatmap, width = 5, height = 5)

# 细胞亚群比例图 -------------------------------------------------
comparison <- list(c("NR","TCMR"),c("NR","ABMR"),c("NR","Mix"))
p_cell_ratio <- plot_cell_ratio_comparison(cell_subset = Tcell,
                                           sample_col = "HLA_Splitter", cluster_col =  "Celltype", group_col = "Group",
                                           facet = T, comparison_groups = comparison, 
                                           write_to_csv = "./Result/Tcell_celltype_ratio", plot_colors = GLOBAL_COLORS)
ggsave(fig_file("Tcell/Celltype_Composition", 3, "Celltype_Proportion_Group_Comparison"), p_cell_ratio, width = 10, height = 10)
# Export the analysis source data for figure reproduction.
saveRDS(Tcell, file.path(getOption("hla.data_root"), "urine", "Urine_Tcell.rds"))


# ==============================================================================
# Tcell 不同疾病分组中的差异基因
# ==============================================================================
Idents(Tcell) <- "Group"
ident1 = "TCMR"
ident2 = "NR"
# Discover marker genes with the original assay and comparison settings.
CELLDEG <- FindMarkers(Tcell, ident.1 =ident1, ident.2 = ident2, logfc.threshold = 0.1, min.pct = 0.1)

noise_pattern <- "^MT-|^HB[ABDEGQZ]|^RP[0-9]+-|^AC[0-9]+\\.|^AL[0-9]+\\.|^AP[0-9]+\\.|^CT[ABCD]-|^RPS|^RPL|MALAT1|NEAT1"
CELLDEG$gene = rownames(CELLDEG)
CELLDEG <- CELLDEG[!is.na(CELLDEG$p_val),] %>%
filter(!grepl(noise_pattern, gene))

CELLDEG$gene <- rownames(CELLDEG)
write.csv(CELLDEG, paste0("./Result/Tcell/Differential_Expression/Tcell_Group_marker_",ident1,"-",ident2,".csv"))

# 火山图 & GO 富集
p_volcano <- plot_volcano(diffexp = CELLDEG,
                          top_n_labels = 10,
                          pos_group_name = ident1, 
                          neg_group_name = ident2)

ggsave(fig_file("Tcell/Differential_Expression", 1, paste0("Volcano_", ident1, "_vs_", ident2)), p_volcano, height = 5, width = 6)

# 基因表达比较
selected_gene = c("CXCR3","CXCR4","CD27", "IFNG", "GZMK", "GZMB" ,"PDCD1","CTLA4")
selected_gene = c("CXCR3", "CD27",
                  "TNF", "IFNG",
                  "HLA-A","HLA-B" ,"HLA-C","HLA-DRA","HLA-DPA","HLA-DQA",
                  "HAVCR2","TIGIT",  "PDCD1", "CTLA4")

p1 = plot_violin_comparison(seurat_obj = Tcell, group_col = "Group", 
                            comparison_groups = comparison, 
                            features_to_plot = selected_gene, plot_colors = GLOBAL_COLORS)
ggsave(fig_file("Tcell/Differential_Expression", 2, "Marker_Genes_Group_Violin"), p1, height = 9, width = 8)

# 保存
saveRDS(Tcell,file.path(getOption("hla.data_root"), "urine", "Urine_Tcell.rds"))


# ==============================================================================
# 分析CXCR3基因富集程度
# ==============================================================================
target_gene <- "CXCR3"
TCMR_subset <- subset(Tcell, subset = Group == "TCMR")

p_Vln_TCMR <- VlnPlot(TCMR_subset, features = target_gene, group.by = "Celltype", pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, linewidth = 0.6) + # 增加箱型图显式中位数
  theme_classic(base_size = 14) +
  labs(
    title = bquote(italic(.(target_gene))~"Expression in TCMR Group"), 
    x = NULL, 
    y = "Expression Level"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    legend.position = "none"
  )

# 统计 TCMR 组内，各细胞类型的 CXCR3+ 细胞比例
expr_data_tcmr <- FetchData(TCMR_subset, vars = c(target_gene, "Celltype"), layer = "data")
colnames(expr_data_tcmr) <- c("Expression", "Cluster")

cxcr3_tcmr_stats <- expr_data_tcmr %>%
  group_by(Cluster) %>%
  summarise(
    Total_Cells = n(),                                     
    Pos_Cells = sum(Expression > 0),                       
    Mean_Expression = mean(Expression),                    
    Mean_Pos_Expression = mean(Expression[Expression > 0]),
    Pos_Ratio = (Pos_Cells / Total_Cells) * 100            
  ) %>%
  arrange(desc(Pos_Ratio)) # 按表达比例从高到低排序，排在第一的就是富集最多的细胞类型

print("--- TCMR组内各细胞类型 CXCR3+ 比例统计 ---")
print(cxcr3_tcmr_stats)

p_bar_TCMR <- ggplot(cxcr3_tcmr_stats, aes(x = reorder(Cluster, -Pos_Ratio), y = Pos_Ratio, fill = Cluster)) +
  geom_col(color = "black", linewidth = 0.6, width = 0.7) + 
  geom_text(aes(label = sprintf("%.1f%%", Pos_Ratio)), vjust = -0.8, size = 3.5, fontface = "bold") +
  theme_classic(base_size = 14) +
  labs(
    title = bquote("Proportion of "~italic(.(target_gene))^"+"~" Cells in TCMR"), 
    x = "Cell Type in TCMR",
    y = "Percentage of Positive Cells (%)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

final_plot <- p_Vln_TCMR + p_bar_TCMR
print(final_plot)

ggsave(fig_file("Tcell/CXCR3", 1, "TCMR_CXCR3_Expression_by_Celltype"), final_plot, height = 9, width = 8)

# ==============================================================================
# 比较TCMR和ABMR的杀伤性基因表达
# ==============================================================================
cytox_genes <- c("GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7", "FASLG")
valid_cytox_genes <- intersect(cytox_genes, rownames(Tcell))

# 1. 提取仅包含 TCMR 和 ABMR 的细胞子集
rejection_subset <- subset(Tcell, subset = Group %in% c("TCMR", "ABMR"))
Idents(rejection_subset) <- "Group"

# 2. 强制对所有杀伤基因进行组间差异检验
cytox_deg <- FindMarkers(
  rejection_subset, 
  ident.1 = "TCMR", 
  ident.2 = "ABMR", 
  features = valid_cytox_genes,
  logfc.threshold = 0,
  min.pct = 0
) %>%
  mutate(Gene = rownames(.)) %>%
  relocate(Gene) %>%
  arrange(p_val_adj)

print("--- TCMR vs ABMR 杀伤性基因差异检验结果 (正值表示在TCMR中更高) ---")
print(cytox_deg)

# 可视化 1: 比较组别的核心杀伤基因气泡图 (DotPlot)
p_cytox_dot <- DotPlot(rejection_subset, features = valid_cytox_genes, group.by = "Group") +
  theme_classic(base_size = 14) +
  labs(
    title = "Cytotoxic Gene Expression: TCMR vs ABMR",
    x = "Cytotoxic Effectors", 
    y = "Rejection Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic", color = "black"), # 基因名斜体
    axis.text.y = element_text(color = "black", face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  ) +
  scale_color_gradientn(colours = colorRampPalette(c("navy", "white", "firebrick3"))(100))

# 可视化 2: 堆叠分裂小提琴图 (Stacked Split Violin Plot)
stacked_vln_plot = plot_violin_comparison(seurat_obj = rejection_subset, 
                                          group_col = "Group", 
                                          comparison_groups = comparison, 
                                          features_to_plot = c("GZMB", "GZMK", "PRF1","NKG7"), 
                                          plot_colors = GLOBAL_COLORS,
                                          ncol = 4)

# 终极合并与导出 (DotPlot + Stacked Violin)
combine_plot <- cowplot::plot_grid(
  p_cytox_dot, 
  stacked_vln_plot, 
  ncol = 1, 
  labels = c("A", "B"),
  rel_heights = c(2, 2)  
)

print(combine_plot)
ggsave(fig_file("Tcell/Differential_Expression", 3, "TCMR_vs_ABMR_Cytotoxic_Genes"), combine_plot, height = 6, width = 8)

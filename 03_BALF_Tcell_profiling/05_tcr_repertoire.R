# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: 04_tcr_repertoire.R
# Figure/table scope: Fig3; FigS6; Supplementary Table 13
# Data cohort: balf. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Expanded T cells: clone_ratio > 1; distinguish cell counts from clonotype counts.
# Keep within-sample TCR denominators and clonotype identifiers as defined below.

# ==============================================================================
# File: 04_tcr_repertoire.R
# Description: 完整 TCR 组库分析、克隆扩增、多样性、V/J 基因偏好、GLIPH2 抗原网络
# ==============================================================================
rm(list = ls())

source(file.path(getOption("hla.code_root"), "03_BALF_Tcell_profiling/01_setup_and_helpers.R"))

# Read the saved object or source table used by this analysis.
Tcell <- readRDS("./Result/lavage_Tcell_Scored.rds")
Sample_Info <- read.csv(file.path(getOption("hla.data_root"), "balf", "SampleInfo.csv"))

# 1. 整合TCR到seurat中 (scRepertoire & scCustomize) ----------------------------
contig_all <- data.frame()
for (i in c(1:4)){
  contig <- read.csv(paste0("./Data/TCR/Group",i,"/filtered_contig_annotations.csv"), stringsAsFactors = FALSE)
  contig$barcode <- paste0(str_sub(contig$barcode, 1, 16), "-", i)
  contig_all <- rbind(contig_all, contig)
}

contig_all <- contig_all[contig_all$barcode %in% intersect(colnames(Tcell), contig_all$barcode), ]
contig_all$Sample <- Tcell$HLA_Splitter[contig_all$barcode]
# Export the analysis source data for figure reproduction.
write.csv(contig_all, "./Result//Intergrate_TCR_contig.csv")

contig_list <- contig_all %>% group_by(Sample) %>% group_map(~ as.list(.))
combined.TCR <- combineTCR(contig_list, removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)

Tcell <- combineExpression(combined.TCR, Tcell, cloneCall = "strict", proportion = TRUE)
Tcell$clonalFrequency[is.na(Tcell$clonalFrequency)] <- 0

# 2. TCR chain质检 -------------------------------------------------------------
Tcell$TCR_chain <- "None"
Tcell$TCR_chain[grepl("TRAC", Tcell$CTgene)] <- "TRA Chain"
Tcell$TCR_chain[grepl("TRBC", Tcell$CTgene)] <- "TRB Chain"
Tcell$TCR_chain[grepl("TRAC", Tcell$CTgene) & grepl("TRBC", Tcell$CTgene)] <- "TRA+TRB Chain"

chain_proportions <- as.data.frame(table(Tcell$TCR_chain))
colnames(chain_proportions) <- c("TCR_Chain_Type", "Count")
chain_proportions$Proportion <- chain_proportions$Count / sum(chain_proportions$Count)
chain_proportions$TCR_Chain_Type <- factor(chain_proportions$TCR_Chain_Type, levels = c("None", "TRA Chain", "TRB Chain", "TRA+TRB Chain"))

p1 <- ggplot(chain_proportions, aes(x = TCR_Chain_Type, y = Proportion, fill = TCR_Chain_Type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)), vjust = -0.5, size = 4) +
  labs(title = "Overall Proportion of Cells with Different TCR Chain Types", x = "TCR Chain Type Detected", y = "Proportion of Cells") +
  theme_minimal() + scale_fill_brewer(palette = "Set2") + scale_y_continuous(expand = expansion(mult = c(0,0.2))) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, color = "black"), 
        axis.line = element_line(color = "black"))
ggsave("./Figures/Tcell_chain_percent.pdf", plot = p1, height = 5, width = 6)

# 3. 识别高扩增TCR克隆 ---------------------------------------------------------
library(scales)
df <- data.frame(clonalFrequency = Tcell$clonalFrequency, SampleName = Tcell$HLA_Splitter,
                 Group = Tcell$Group, clone_ratio = 0, 
                 CTaa = Tcell$CTaa, Celltype=Tcell$Celltype)
group <- mixedsort(unique(df$SampleName))

for (i in seq_along(group)){
  df_tmp <- df[df$SampleName == group[i], ]
  base_line <- mean(df_tmp$clonalFrequency[df_tmp$clonalFrequency != 0])
  df_tmp$clone_ratio <- df_tmp$clonalFrequency / base_line
  df_tmp$clone_ratio[is.na(df_tmp$clone_ratio)] <- 0
  df_tmp$clone_size <- ifelse(df_tmp$clone_ratio <= 1, "Un-expand",
                              ifelse(df_tmp$clone_ratio <= 2, "Small (expand:1-2)",
                                     ifelse(df_tmp$clone_ratio <= 5, "Medium (expand:2-5)",
                                            ifelse(df_tmp$clone_ratio <= 10, "Large (expand:5-10)", "Hyper (expand:>10)"))))
  df[rownames(df_tmp), "clone_size"] <- df_tmp$clone_size
  df[rownames(df_tmp), "clone_ratio"] <- df_tmp$clone_ratio
}

df$clone_size <- factor(df$clone_size, levels = c("Hyper (expand:>10)","Large (expand:5-10)","Medium (expand:2-5)",
                                                  "Small (expand:1-2)", "Un-expand"))
Tcell$clone_ratio <- df$clone_ratio
Tcell$clone_size <- df$clone_size

clone_colors <- c(
  "Un-expand" = "grey", 
  "Small (expand:1-2)" = "#f9f871", 
  "Medium (expand:2-5)" = "#ffc75f", 
  "Large (expand:5-10)" = "#ff6f91", 
  "Hyper (expand:>10)" = "#845ec2"
)

# 比较分组中不同扩增程度的TCR比例
comparison <- list(c("CAP","LC"), c("CAP","ILD"), c("ILD","LC"))
p <- ggplot(df, aes(x = Group, fill = clone_size)) +
  geom_bar(position = "fill",color = "black", width = 0.7, alpha = 0.7)+
  theme_classic() + scale_y_continuous(labels = percent) + 
  scale_fill_discrete(type =clone_colors)+
  labs(title = "TCR clone frequency by expansion levels", x = "Group", y = "Proportion of TCR clones") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./Figures/Expanded_TCR_clone_ratio.pdf", plot = p, height = 5, width = 8)

# 扩增的T细胞中不同细胞亚型比例####
df_summary <- df %>%
  dplyr::count(Group, Celltype, clone_size, name = "clone_count") %>%
  group_by(Group) %>%
  mutate(
    total_in_group = sum(clone_count),
    proportion = clone_count / total_in_group,       # 0~1 的比例
    percentage = proportion * 100                    # 0~100 的百分比
  ) %>% ungroup()

p1 <- ggplot(df_summary, aes(x = Celltype, y = proportion, fill = clone_size)) +
  geom_col(color = "black", width = 0.7, alpha = 0.7) +
  facet_grid(Group~.)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(y = "Proportion of T cells",x = NULL) +
  scale_fill_discrete(type = clone_colors)+
  theme_classic() +theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 12),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11)
  )
print(p1)
ggsave("./Figures/TCR_expanded_frequency.pdf", plot = p1, height = 4.5, width = 4.5)

#dotplot
Tcell$filter_clone_ratio = Tcell$clone_ratio[Tcell$clone_ratio>1]
Tcell$filter_clone_ratio[is.na(Tcell$filter_clone_ratio)] = 0
p_dot = plot_dotplot(seurat_object = Tcell, features = "filter_clone_ratio",
                 group_by = "Group",cluster_by = "Celltype")
p_dot
ggsave(paste0("./Figures/TCR_CloneOccupy_celltype_dotplot_group.pdf"),plot=p_dot,height = 3,width = 6)


# 4. 高扩增克隆 Heatmap 展示 (Top 1000 -> Heatmap) -----------------------------
all_clones <- Tcell@meta.data[Tcell$clonalFrequency > 0, ] %>% mutate(Clone_ID = paste(HLA_Splitter, CTaa, sep = "_"))
top_clones <- all_clones %>% arrange(desc(clone_ratio)) %>% distinct(Clone_ID) %>% head(1000) %>% pull(Clone_ID)

heatmap_data <- all_clones %>% filter(Clone_ID %in% top_clones) %>% distinct(Clone_ID, .keep_all = TRUE)
heatmap_matrix <- heatmap_data %>% select(Clone_ID, HLA_Splitter, clone_ratio) %>%
  pivot_wider(names_from = HLA_Splitter, values_from = clone_ratio, values_fill = 0) %>%
  tibble::column_to_rownames("Clone_ID") %>% as.matrix()

log_heatmap_matrix <- log10(heatmap_matrix + 1)
log_heatmap_matrix <- log_heatmap_matrix[, mixedorder(colnames(log_heatmap_matrix))]

annotation_col <- Sample_Info %>% select(Patient.ID, Group) %>% tibble::remove_rownames()  %>% tibble::column_to_rownames("Patient.ID")
annotation_row <- data.frame(Celltype = heatmap_data$Celltype, row.names = rownames(log_heatmap_matrix))
ann_colors <- list(Group = GROUP_COLORS, Celltype = setNames(GLOBAL_COLORS[1:length(levels(heatmap_data$Celltype))], levels(heatmap_data$Celltype)))

p_heatmap <- pheatmap(log_heatmap_matrix, main = "Top 1000 High-Expansion TCR Clones Across Disease Groups",
                      fontsize_row = 8, fontsize_col = 10, show_rownames = FALSE,
                      annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = ann_colors,
                      cluster_rows = TRUE, cluster_cols = FALSE, angle_col = 0, border_color = NA,
                      color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
ggsave("./Figures/TCR_clone_high_expansion.pdf", plot = p_heatmap, height = 8, width = 9)

# 5. TCR CloneSize 与 UMAP -----------------------------------------------------
p1 <- DimPlot(Tcell, group.by = "Celltype", label = FALSE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p2 <- DimPlot(Tcell, group.by = "cloneSize", alpha = 0.5) + scale_color_manual(values = clone_colors)
p3 <- DimPlot(Tcell, group.by = "Group", alpha = 0.5) + scale_color_manual(values = GROUP_COLORS)
p4 <- FeaturePlot(Tcell, features = "clone_ratio", pt.size = 0.5, cols = c("grey","#E42320"), order = TRUE)
ggsave("./Figures/Tcell_UMAP_TCR_cloneSize.pdf", p1+p2+p3+p4, height = 8, width = 10)

p5 <- FeaturePlot(Tcell, features = "clone_ratio", split.by = "Group", pt.size = 1, cols = c("grey","#E42320"), order = TRUE)
ggsave("./Figures/Tcell_UMAP_clone_ratio_Group.pdf", p5, height = 4.5, width = 13)

# 6. 不同分组/亚群的克隆占比与 Upset 图 ----------------------------------------
Tcell_recall <- subset(Tcell, subset = clonalFrequency > 0)

# Upset Plot
library(UpSetR)
df_upset <- data.frame(CDR3 = Tcell_recall$CTaa, Group = Tcell_recall$Group, Sample = Tcell_recall$HLA_Splitter)
df_upset$CDR3_TRB <- sapply(df_upset$CDR3, function(x) ifelse(grepl("_", x), unlist(str_split(x, "_"))[2], NA))
df_upset <- df_upset[!is.na(df_upset$CDR3_TRB), ]

tcr_sets <- df_upset %>% group_by(Group) %>% summarise(TCRs = list(unique(CDR3_TRB)))
tcr_sets_for_upset <- setNames(tcr_sets$TCRs, tcr_sets$Group)[mixedsort(names(setNames(tcr_sets$TCRs, tcr_sets$Group)))]

pdf("./Figures/TCR_Clone_Overlap_Upset_Plot.pdf", width = 12, height = 8)
upset(fromList(tcr_sets_for_upset), sets = names(tcr_sets_for_upset), order.by = "freq", decreasing = TRUE,
      mb.ratio = c(0.6, 0.4), point.size = 2.5, line.size = 1, text.scale = c(1.3, 1.3, 1, 1, 1.5, 1))
dev.off()

# CDR3 Length
library(purrr)
cdr3Length <- function(combined.TCR, chain ="both"){
  meta_all <- data.frame()
  chain_idx <- match(chain, c("TRA","TRB","both"))
  for (i in seq_along(combined.TCR)){
    meta <- switch(chain_idx, combined.TCR[[i]]$cdr3_aa1, combined.TCR[[i]]$cdr3_aa2, c(combined.TCR[[i]]$cdr3_aa1, combined.TCR[[i]]$cdr3_aa2)) 
    meta <- meta[meta != "NA"]
    meta <- lapply(meta, function(x) unlist(strsplit(x, split = ";"))[1]) %>% unlist()
    meta <- table(nchar(meta)) %>% prop.table() %>% data.frame()
    meta$Freq <- meta$Freq * 100
    meta$group <- names(combined.TCR)[i]
    meta_all <- rbind(meta_all, meta)
  }
  meta_all$Var1 <- factor(as.numeric(as.character(meta_all$Var1)))
  ggplot(meta_all, aes(x=Var1, y=Freq, fill=group)) + geom_col(position="dodge", width=0.9) + theme_classic() + xlab("Length") + ylab("Percent of CDR3(AA)")
}
expression2List <- function(scRNA, split.by=NULL) { split(scRNA@meta.data %>% mutate(cdr3_aa1 = map_chr(CTaa, ~str_split(.x, "_")[[1]][1]), cdr3_aa2 = map_chr(CTaa, ~str_split(.x, "_")[[1]][2])), scRNA@meta.data[[split.by]]) }
combined.TCR_grp <- expression2List(Tcell_recall, split.by = "Group")

p_len <- cdr3Length(combined.TCR_grp, chain ="both") + scale_fill_manual(values = GLOBAL_COLORS)
ggsave("./Figures/TCR_lengthContig.pdf", plot = p_len, height = 5, width = 8)

# 7. TCR 多样性 (Diversity) ----------------------------------------------------
df_list <- expression2List(Tcell_recall, split.by = "HLA_Splitter")
df_diversity <- clonalDiversity(df_list, cloneCall = "strict", exportTable = TRUE)

clonalDiversity(Tcell_recall, group.by = "HLA_Splitter", x.axis = "Group")
df_final <- df_diversity %>% melt(id.vars = "Group") %>% rename(Patient.ID = Group) %>%
  left_join(Sample_Info %>% dplyr::select(Patient.ID, Group), by = "Patient.ID") %>% rename(Diversity_Metric = variable, Value = value)
df_final$Group <- factor(df_final$Group, levels = c("CAP", "LC", "ILD"))

df_final <- left_join(df_melted, Sample_Info %>% dplyr::select(Patient.ID, Group), by = "Patient.ID") %>%
  rename(Diversity_Metric = variable, Value = value) # Clearer column names
df_final$Group <- factor(df_final$Group, levels = c("CAP", "LC","ILD"))

comparison_list <- list(c("CAP", "LC"),c("ILD", "LC"),c("CAP", "ILD"))

p_div <- ggplot(df_final, aes(x = Group, y = Value, fill = Group, color = Group)) +
  facet_wrap(~ Diversity_Metric, scales = "free_y", ncol = 3, labeller = as_labeller(function(x) paste(x, "Index"))) +
  geom_boxplot(outlier.shape = NA, fill = "white", linewidth = 0.6) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 2) +
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = GLOBAL_COLORS) + scale_color_manual(values = GLOBAL_COLORS) + theme_classic()
ggsave("./Figures/TCR_clonalDiversity.pdf", p_div, height = 8, width = 8)

# 8. V/J 基因偏好与 Circos/Venn 图 ---------------------------------------------
TCR_stats <- Tcell_recall@meta.data %>%
  dplyr::select(gene = CTgene, Group, cloneSize, celltype = Celltype, sample = HLA_Splitter, clone_ratio, CDR3 = CTaa, Freq = clonalFrequency) %>%
  tidyr::separate(gene, into = c("TRA_gene", "TRB_gene"), sep = "_", remove = FALSE, fill = "right") %>%
  dplyr::mutate(TRAV = str_extract(TRA_gene, "TRAV[^.]*"), TRBV = str_extract(TRB_gene, "TRBV[^.]*"), TRBJ = str_extract(TRB_gene, "TRBJ[0-9-]+"),
                CDR3 = purrr::map_chr(CDR3, ~ purrr::pluck(strsplit(., "_"), 1, 2), .default = NA_character_),
                CDR3 = purrr::map_chr(CDR3, ~ purrr::pluck(strsplit(., ";"), 1, 1), .default = NA_character_)) %>%
  dplyr::filter(!is.na(TRBV) & TRBV != "NA")

# Venn & Heatmap
Group_uniq <- unique(TCR_stats$Group)
for (i in seq_along(Group_uniq)){
  df_tmp <- TCR_stats[TCR_stats$Group == Group_uniq[i], ]
  Set1 <- rownames(df_tmp[df_tmp$TRAV != "NA", ])
  Set2 <- rownames(df_tmp[df_tmp$TRBV != "NA", ])
  venn.diagram(x=list(Set1, Set2), filename=paste0("./Figures/TCR_TRBV_", Group_uniq[i], "_Venn.pdf"),
               category.names = c("TRAV", "TRBV"), fill=c('#FFFFCC','#CCFFFF'), output=TRUE)
  
  df_tmp_heat <- df_tmp[df_tmp$TRBV != "NA" & !is.na(df_tmp$TRBJ), ]
  df_tmp_heat$TRBV <- factor(df_tmp_heat$TRBV, levels = mixedsort(unique(TCR_stats$TRBV)))
  mat_heat <- table(df_tmp_heat$TRBV, df_tmp_heat$TRBJ) %>% prop.table() %>% data.frame()
  p_heat <- ggplot(mat_heat, aes(Var2, Var1, color = Freq, size = Freq)) + geom_point() + scale_color_gradient(low = "white", high = "red") + theme_bw() + theme(axis.text.x = element_text(angle=270))
  ggsave(paste0("./Figures/TCR_TRBV_", Group_uniq[i], "_heatmap.pdf"), p_heat, height = 8, width = 4.5)
}

# Circos
df_circos <- Tcell@meta.data %>% filter(clone_ratio > 1) %>% left_join(TCR_stats, by = c("CTgene" = "gene"))
for (grp in unique(df_circos$Group.x)) {
  df_tmp_circ <- filter(df_circos, Group.x == grp)
  mat_circ <- table(df_tmp_circ$TRBJ, df_tmp_circ$TRBV)
  mat_circ <- mat_circ[rownames(mat_circ) != "NA", colnames(mat_circ) != "NA"]
  
  pdf(paste0("./Figures/TCell_TRBV_usage_CIRCOS_", grp, ".pdf"), width = 8, height = 8)
  circos.par(gap.after = c(rep(2, nrow(mat_circ)-1), 10, rep(2, ncol(mat_circ)-1), 10))
  chordDiagram(mat_circ, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat_circ))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + cm_h(2), CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
    circos.axis(h = "top", labels.cex = 1, major.tick.length = 0.2, sector.index = CELL_META$sector.index, track.index = 2)
  }, bg.border = NA)
  title(paste0("Pairing of TRBV with TRBJ (", grp, ")"))
  circos.clear()
  dev.off()
}

# SeqLogo
for (grp in Group_uniq){
  df_seq <- TCR_stats[TCR_stats$Group == grp & TCR_stats$CDR3 != "NA", ]
  df_seq$CDR3 <- str_sub(paste0(df_seq$CDR3, "--------"), 0, 16)
  bc_tmp <- unique(df_seq$barcode[df_seq$clone_ratio > 5 & df_seq$Freq > 10])
  if(length(bc_tmp) > 0) {
    fig_list_all <- lapply(bc_tmp, function(bc) {
      d <- df_seq[df_seq$barcode == bc, ]
      wrap_plots(lapply(1:nrow(d), function(k) {
        ggplot() + geom_logo(d$CDR3[k], method="probability") + theme_logo() + ggtitle(paste(d$sample[k], d$celltype[k], d$TRBV[k], d$TRBJ[k], sep=":"))
      }))
    })
    ggsave(paste0("./Figures/CDR3_seqlogo_", grp, ".pdf"), wrap_plots(fig_list_all, ncol = 1), height = 15, width = 15)
  }
}

# 9. GLIPH2 抗原特异性聚类 -----------------------------------------------------
gliph_data <- data.frame(chain = "TRB", CDR3b = TCR_stats$CDR3, TRBV = TCR_stats$TRBV, 
                         TRBJ = TCR_stats$TRBJ, patient = TCR_stats$sample, 
                         Group = TCR_stats$Group, counts = TCR_stats$cloneSize, celltype = TCR_stats$celltype)

# 耗时较长，请根据需要取消注释运行：
# res <- gliph2(cdr3_sequences = gliph_data, sim_depth = 1000, n_cores = 1)
# save(res, file = "./Rawdata/Gliph_result.rda")
load(file = "./Rawdata/Gliph_result.rda")

cluster_prop <- res$cluster_properties[res$cluster_properties$cluster_size >= 4, ]
nodes <- lapply(cluster_prop$members, function(x) unlist(str_split(x, " "))) %>% unlist() %>% as.data.frame() %>% distinct()
colnames(nodes) <- c("label")

meta_gliph <- do.call(rbind, res$cluster_list)
nodes <- cbind(nodes, meta_gliph[match(nodes$label, meta_gliph$CDR3b), ])
edges <- data.frame(from=res$connections$V1, to=res$connections$V2, weight=1, type=res$connections$type) %>%
  filter(from %in% nodes$label & to %in% nodes$label & from != to)

net_pc <- graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
V(net_pc)$color <- GLOBAL_COLORS[factor(nodes$Group)]

pdf("./Figures/GLIPH2_Network.pdf", width = 10, height = 10)
plot(net_pc, vertex.size=3, layout=layout_with_fr(net_pc), vertex.label.cex=0.5, vertex.label.dist=0.5, edge.color="gray50", edge.curved=.1)
legend("topleft", levels(factor(nodes$Group)), pch=21, pt.bg=GLOBAL_COLORS)
dev.off()

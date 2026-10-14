# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis.
# Original source: 06_TCR_repertoire_analysis.R
# Figure/table scope: Fig5; FigS9; Supplementary Table 18; Methods CXCR3-axis analysis
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Expanded T cells: clone_ratio > 1; distinguish cell counts from clonotype counts.
# Keep within-sample TCR denominators and clonotype identifiers as defined below.

# ==============================================================================
# File: 03_TCR_repertoire_analysis_scRepertoire2.R
# Description: TCR链分布、克隆扩增比例、V/J基因特征、组库重叠与多样性分析
# ==============================================================================
source(file.path(getOption("hla.code_root"), "04_Urine_immune_programs/01_setup_and_helpers.R"))
suppressPackageStartupMessages({
  library(Seurat); library(scRepertoire); library(dplyr); library(tidyr)
  library(tibble); library(ggplot2); library(stringr); library(purrr)
  library(scales); library(gtools); library(patchwork)
})

dir.create("./Result", showWarnings = FALSE, recursive = TRUE)
dir.create("./Figures", showWarnings = FALSE, recursive = TRUE)

# Read the saved object or source table used by this analysis.
Tcell <- readRDS("./Result/Objects/Urine_Tcell.rds")
Sample_Info <- read.csv(data_path("metadata", "Sample_Info.csv"), stringsAsFactors = FALSE)

# 1. 整合TCR到Seurat中 ---------------------------------------------------------
samplenames <- ANALYSIS_BATCHES
contig_all <- data.frame()

for (smp in samplenames) {
  contig_files <- find_tcr_contig_files(smp, full_length = FALSE)
  if (nrow(contig_files) == 0) {
    warning("No TCR contig files found for batch: ", smp)
    next
  }

  for (row_id in seq_len(nrow(contig_files))) {
    f <- contig_files$file[row_id]
    run_index <- contig_files$run_index[row_id]
    contig <- read.csv(f, stringsAsFactors = FALSE)
    contig$barcode <- paste0(stringr::str_sub(contig$barcode, 1, 16), "-", smp, "-", run_index)
    contig$sample <- smp
    contig$batch <- basename(dirname(dirname(f)))
    contig_all <- rbind(contig_all, contig)
  }
}

if (nrow(contig_all) == 0) {
  stop("No usable TCR contig files were found under Data/Group*/GEM.")
}
common_bc <- intersect(colnames(Tcell), contig_all$barcode)
contig_all <- contig_all[contig_all$barcode %in% common_bc, , drop = FALSE]
contig_all$Sample <- Tcell@meta.data[contig_all$barcode, "HLA_Splitter"]
contig_all <- contig_all[!is.na(contig_all$Sample) & contig_all$Sample != "", , drop = FALSE]
# Export the analysis source data for figure reproduction.
write.csv(contig_all, "./Result/TCR/Contigs_Metadata/Intergrate_TCR_contig.csv", row.names = FALSE)

contig_list <- split(contig_all, contig_all$Sample)
contig_list <- contig_list[gtools::mixedsort(names(contig_list))]

combined.TCR <- combineTCR(contig_list, removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)

Tcell <- combineExpression(combined.TCR, Tcell, cloneCall = "strict", proportion = TRUE)

Tcell$clonalFrequency[is.na(Tcell$clonalFrequency)] <- 0

# 2. TCR chain质检 -------------------------------------------------------------
Tcell$TCR_chain <- dplyr::case_when(
  grepl("TRAC", Tcell$CTgene) & grepl("TRBC", Tcell$CTgene) ~ "TRA+TRB Chain",
  grepl("TRAC", Tcell$CTgene) ~ "TRA Chain",
  grepl("TRBC", Tcell$CTgene) ~ "TRB Chain",
  TRUE ~ "None"
)

chain_proportions <- as.data.frame(table(Tcell$TCR_chain)) %>%
  setNames(c("TCR_Chain_Type", "Count")) %>%
  dplyr::mutate(Proportion = Count / sum(Count))

chain_proportions$TCR_Chain_Type <- forcats::fct_relevel(
  chain_proportions$TCR_Chain_Type, "None", "TRA Chain", "TRB Chain", "TRA+TRB Chain"
)

p1 <- ggplot(chain_proportions, aes(TCR_Chain_Type, Proportion, fill = TCR_Chain_Type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)), vjust = -0.5, size = 4) +
  scale_y_continuous(expand = expansion(c(0, 0.2))) +
  labs(title = "Overall Proportion of Cells with Different TCR Chain Types", x = "TCR Chain Type Detected", y = "Proportion of Cells") +
  my_theme() + scale_fill_brewer(palette = "Set2")

ggsave(fig_file("TCR/Overview", 1, "TCR_Chain_Type_Proportion"), plot = p1, height = 5, width = 5.5)

# 3. 识别高扩增TCR克隆 ---------------------------------------------------------
df <- data.frame(
  clonalFrequency = Tcell$clonalFrequency,
  SampleName = Tcell$HLA_Splitter,
  Group = Tcell$Group,
  clone_ratio = 0,
  CTaa = Tcell$CTaa,
  Celltype = Tcell$Celltype,
  row.names = colnames(Tcell)
)

for (smp in gtools::mixedsort(unique(df$SampleName))) {
  idx <- which(df$SampleName == smp)
  base_line <- mean(df$clonalFrequency[idx][df$clonalFrequency[idx] > 0], na.rm = TRUE)
  if (is.na(base_line) || base_line == 0) base_line <- 1
  df$clone_ratio[idx] <- df$clonalFrequency[idx] / base_line
}

df$clone_size <- dplyr::case_when(
  df$clone_ratio <= 1 ~ "Un-expand",
  df$clone_ratio <= 2 ~ "Small (expand:1-2)",
  df$clone_ratio <= 5 ~ "Medium (expand:2-5)",
  df$clone_ratio <= 10 ~ "Large (expand:5-10)",
  TRUE ~ "Hyper (expand:>10)"
)

df$clone_size <- factor(
  df$clone_size,
  levels = c("Hyper (expand:>10)", "Large (expand:5-10)", "Medium (expand:2-5)", "Small (expand:1-2)", "Un-expand")
)

Tcell$clone_ratio <- df$clone_ratio
Tcell$clone_size <- df$clone_size

clone_colors <- c(
  "Un-expand" = "grey",
  "Small (expand:1-2)" = "#f9f871",
  "Medium (expand:2-5)" = "#ffc75f",
  "Large (expand:5-10)" = "#ff6f91",
  "Hyper (expand:>10)" = "#845ec2"
)

comparison <- list(c("NR", "TCMR"), c("NR", "ABMR"), c("TCMR", "ABMR"))

p <- ggplot(df, aes(Group, fill = clone_size)) +
  geom_bar(position = "fill", color = "black", width = 0.7, alpha = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = clone_colors, na.value = "grey80") +
  labs(title = "TCR clone frequency by expansion levels", x = "Group", y = "Proportion of TCR clones") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(fig_file("TCR/Clonality_Expansion", 1, "Clone_Size_Proportion_by_Group"), plot = p, height = 5, width = 8)

df_summary <- df %>%
  dplyr::count(Group, Celltype, clone_size, name = "clone_count") %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(total_in_group = sum(clone_count), proportion = clone_count / total_in_group, percentage = proportion * 100) %>%
  dplyr::ungroup()

p1 <- ggplot(df_summary, aes(Celltype, proportion, fill = clone_size)) +
  geom_col(color = "black", width = 0.7, alpha = 0.7) +
  facet_grid(Group ~ .) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = clone_colors, na.value = "grey80") +
  labs(y = "Proportion of T cells", x = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 11))

ggsave(fig_file("TCR/Clonality_Expansion", 2, "Clone_Size_by_Celltype_and_Group"), plot = p1, height = 4.5, width = 4.5)

Tcell$clone_ratio <- ifelse(Tcell$clone_ratio > 1, Tcell$clone_ratio, 0)
p_dot <- plot_dotplot(Tcell, features = "clone_ratio", group_by = "Group", cluster_by = "Celltype")
ggsave(fig_file("TCR/Clonality_Expansion", 3, "Clone_Ratio_Dotplot_Celltype_Group"), plot = p_dot, height = 3, width = 6)
saveRDS(Tcell, "./Result/Objects/Urine_Tcell.rds")

# 4. 高扩增克隆 Heatmap --------------------------------------------------------
all_clones <- Tcell@meta.data[Tcell$clonalFrequency > 0, , drop = FALSE] %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::mutate(Clone_ID = paste(HLA_Splitter, CTaa, sep = "_"))

top_clones <- all_clones %>%
  dplyr::arrange(desc(clone_ratio)) %>%
  dplyr::distinct(Clone_ID) %>%
  head(1000) %>%
  dplyr::pull(Clone_ID)

heatmap_data <- all_clones %>%
  dplyr::filter(Clone_ID %in% top_clones) %>%
  dplyr::arrange(match(Clone_ID, top_clones)) %>%
  dplyr::distinct(Clone_ID, .keep_all = TRUE)

heatmap_matrix <- heatmap_data %>%
  dplyr::select(Clone_ID, HLA_Splitter, clone_ratio) %>%
  tidyr::pivot_wider(names_from = HLA_Splitter, values_from = clone_ratio, values_fill = 0) %>%
  tibble::column_to_rownames("Clone_ID") %>%
  as.matrix()

log_heatmap_matrix <- log10(heatmap_matrix + 1)
log_heatmap_matrix <- log_heatmap_matrix[, gtools::mixedorder(colnames(log_heatmap_matrix)), drop = FALSE]

annotation_col <- Sample_Info %>%
  dplyr::select(Patient.ID, Group) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("Patient.ID")
annotation_col <- annotation_col[colnames(log_heatmap_matrix), , drop = FALSE]

annotation_row <- heatmap_data %>%
  dplyr::select(Clone_ID, Celltype) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("Clone_ID")
annotation_row <- annotation_row[rownames(log_heatmap_matrix), , drop = FALSE]

cell_lv <- gtools::mixedsort(unique(as.character(annotation_row$Celltype)))
ann_colors <- list(Group = GROUP_COLORS, Celltype = setNames(GLOBAL_COLORS[seq_along(cell_lv)], cell_lv))

pdf(fig_file("TCR/Clonality_Expansion", 4, "Top1000_High_Expansion_Clone_Heatmap"), height = 8, width = 9)
pheatmap::pheatmap(
  log_heatmap_matrix,
  main = "Top 1000 High-Expansion TCR Clones Across Disease Groups",
  fontsize_row = 8, fontsize_col = 10, show_rownames = FALSE,
  annotation_col = annotation_col, annotation_row = annotation_row,
  annotation_colors = ann_colors, cluster_rows = TRUE, cluster_cols = FALSE,
  angle_col = 0, border_color = NA,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)
dev.off()

# 5. TCR CloneSize 与 UMAP -----------------------------------------------------
p1 <- DimPlot(Tcell, group.by = "Celltype", label = FALSE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p2 <- FeaturePlot(Tcell, features = "clonalProportion", pt.size = 0.5, cols = c("grey", "#E42320"), order = TRUE)
p3 <- FeaturePlot(Tcell, features = "clone_ratio", pt.size = 0.5, cols = c("grey", "#E42320"), order = TRUE)
p4 <- DimPlot(Tcell, group.by = "clone_size", alpha = 0.5) + scale_color_manual(values = clone_colors)

ggsave(fig_file("TCR/Clonality_Expansion", 5, "Tcell_UMAP_TCR_CloneSize"), p1 + p2 + p3 + p4, height = 10, width = 13)

p5 <- FeaturePlot(Tcell, features = "clone_ratio", split.by = "Group", pt.size = 1, cols = c("grey", "#E42320"), order = TRUE)
ggsave(fig_file("TCR/Clonality_Expansion", 6, "UMAP_Clone_Ratio_by_Group"), p5, height = 3, width = 12)

# 6. 克隆重叠、扩增比例与CDR3长度 ---------------------------------------------
Tcell_recall <- subset(Tcell, subset = clonalFrequency > 0)

df_upset <- data.frame(CDR3 = Tcell_recall$CTaa, Group = Tcell_recall$Group, Sample = Tcell_recall$HLA_Splitter)
df_upset$CDR3_TRB <- vapply(df_upset$CDR3, function(x) ifelse(grepl("_", x), stringr::str_split(x, "_", simplify = TRUE)[1, 2], NA), character(1))
df_upset <- df_upset[!is.na(df_upset$CDR3_TRB) & df_upset$CDR3_TRB != "NA", , drop = FALSE]

tcr_sets <- df_upset %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(TCRs = list(unique(CDR3_TRB)), .groups = "drop")

tcr_sets_for_upset <- setNames(tcr_sets$TCRs, tcr_sets$Group)
tcr_sets_for_upset <- tcr_sets_for_upset[gtools::mixedsort(names(tcr_sets_for_upset))]

pdf(fig_file("TCR/Clonality_Expansion", 7, "Clone_Overlap_Upset_Plot"), width = 12, height = 8)
UpSetR::upset(
  UpSetR::fromList(tcr_sets_for_upset),
  sets = names(tcr_sets_for_upset),
  order.by = "freq", decreasing = TRUE,
  mb.ratio = c(0.6, 0.4), point.size = 2.5,
  line.size = 1, text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)
)
dev.off()

df_expand <- df %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(total_tcr_in_group = n()) %>%
  dplyr::group_by(Group, Celltype, total_tcr_in_group) %>%
  dplyr::summarise(expanded_tcr = sum(clone_ratio > 2), .groups = "drop") %>%
  dplyr::mutate(expansion_rate = expanded_tcr / total_tcr_in_group)

p1 <- ggplot(df_expand, aes(Group, expansion_rate, fill = Celltype)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, alpha = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(y = "Proportion of T cells with High−expanded TCR", x = NULL) +
  my_theme()
p1
ggsave(fig_file("TCR/Clonality_Expansion", 8, "High_Expanded_TCR_Frequency_Overall"), plot = p1, height = 5, width = 5)

get_cdr3_len <- function(seu, group.by = "Group", chain = c("TRA", "TRB", "both")) {
  chain <- match.arg(chain)
  md <- seu@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::filter(!is.na(CTaa), CTaa != "NA", !is.na(.data[[group.by]]))
  x <- md %>%
    dplyr::transmute(group = .data[[group.by]], TRA = sub("_.*$", "", CTaa), TRB = ifelse(grepl("_", CTaa), sub("^[^_]*_", "", CTaa), NA_character_))
  x <- switch(chain,
              TRA = dplyr::select(x, group, cdr3 = TRA),
              TRB = dplyr::select(x, group, cdr3 = TRB),
              both = tidyr::pivot_longer(x, c(TRA, TRB), values_to = "cdr3") %>% dplyr::select(group, cdr3))
  x <- x %>%
    tidyr::separate_rows(cdr3, sep = ";") %>%
    dplyr::filter(!is.na(cdr3), cdr3 != "", cdr3 != "NA") %>%
    dplyr::count(group, Var1 = nchar(cdr3), name = "n") %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(Freq = n / sum(n) * 100) %>%
    dplyr::ungroup()
  x$Var1 <- factor(x$Var1, levels = gtools::mixedsort(unique(x$Var1)))
  x
}

df_cdr3 <- get_cdr3_len(Tcell_recall, group.by = "Group", chain = "both")

p_len <- ggplot(df_cdr3, aes(Var1, Freq, fill = group)) +
  geom_col(position = "dodge", width = 0.9) +
  theme_classic() + xlab("Length") + ylab("Percent of CDR3(AA)") +
  scale_fill_manual(values = GLOBAL_COLORS)
p_len
ggsave(fig_file("TCR/CDR3_Sequence", 1, "CDR3_Length_Distribution"), plot = p_len, height = 5, width = 8)

# 7. TCR 多样性 ----------------------------------------------------------------
library(vegan)
calc_diversity <- function(seu, group.by = "HLA_Splitter", clone_col = "CTstrict") {
  md <- seu@meta.data
  if (!clone_col %in% colnames(md)) clone_col <- intersect(c("CTstrict", "CTnt", "CTaa", "CTgene"), colnames(md))[1]
  if (is.na(clone_col)) stop("找不到克隆列：CTstrict / CTnt / CTaa / CTgene")
  md <- md[!is.na(md[[group.by]]) & !is.na(md[[clone_col]]) & md[[clone_col]] != "NA", , drop = FALSE]
  purrr::map_dfr(split(md, md[[group.by]]), function(x) {
    tab <- table(x[[clone_col]])
    tibble::tibble(
      Patient.ID = unique(x[[group.by]])[1],
      shannon = vegan::diversity(tab, index = "shannon"),
      inv.simpson = vegan::diversity(tab, index = "invsimpson"),
      gini.simpson = vegan::diversity(tab, index = "simpson"),
      norm.entropy = ifelse(length(tab) > 1, vegan::diversity(tab, "shannon") / log(length(tab)), 0),
      chao1 = tryCatch(vegan::estimateR(tab)["S.chao1"], error = function(e) NA_real_)
    )
  })
}

df_final <- calc_diversity(Tcell_recall, group.by = "HLA_Splitter", clone_col = "CTstrict") %>%
  tidyr::pivot_longer(-Patient.ID, names_to = "Diversity_Metric", values_to = "Value") %>%
  dplyr::left_join(Sample_Info %>% dplyr::select(Patient.ID, Group), by = "Patient.ID") %>%
  dplyr::mutate(Group = factor(Group, levels = c("NR","TCMR","ABMR","Mix")))

comparison_list <- list(c("NR","TCMR"), c("NR","ABMR"), c("NR","Mix"))

p_div <- ggplot(df_final, aes(Group, Value, fill = Group, color = Group)) +
  facet_wrap(~ Diversity_Metric, scales = "free_y", ncol = 3,
             labeller = as_labeller(function(x) paste(x, "Index"))) +
  geom_boxplot(outlier.shape = NA, fill = "white", linewidth = 0.6) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 2) +
  ggpubr::stat_compare_means(comparisons = comparison_list, method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = GROUP_COLORS) +
  scale_color_manual(values = GROUP_COLORS) +
  theme_classic()
p_div
ggsave(fig_file("TCR/Clonality_Expansion", 9, "Clonal_Diversity_by_Group"), p_div, height = 8, width = 8)

# 8. V/J 基因偏好与 Circos/Venn 图 ---------------------------------------------
TCR_stats <- Tcell_recall@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::select(barcode, gene = CTgene, Group, cloneSize, celltype = Celltype, 
                sample = HLA_Splitter, clone_ratio, CDR3 = CTaa, Freq = clonalFrequency) %>%
  tidyr::separate(gene, into = c("TRA_gene", "TRB_gene"), sep = "_", remove = FALSE, fill = "right") %>%
  dplyr::mutate(
    TRAV = stringr::str_extract(TRA_gene, "TRAV[^.]*"),
    TRBV = stringr::str_extract(TRB_gene, "TRBV[^.]*"),
    TRBJ = stringr::str_extract(TRB_gene, "TRBJ[0-9-]+"),
    CDR3 = purrr::map_chr(CDR3, ~ purrr::pluck(strsplit(.x, "_")[[1]], 2, .default = NA_character_)),
    CDR3 = purrr::map_chr(CDR3, ~ purrr::pluck(strsplit(.x, ";")[[1]], 1, .default = NA_character_))
  ) %>%
  dplyr::filter(!is.na(TRBV), TRBV != "NA")

# 选择扩增后的TCR(clone_ratio>1)
TCR_stats_expand = TCR_stats[TCR_stats$clone_ratio>1, ]
Group_uniq <- unique(TCR_stats_expand$Group)

# heatmap Figure
for (grp in Group_uniq) {
  df_tmp <- TCR_stats_expand[TCR_stats$Group == grp, , drop = FALSE]
  Set1 <- df_tmp$barcode[!is.na(df_tmp$TRAV) & df_tmp$TRAV != "NA"]
  Set2 <- df_tmp$barcode[!is.na(df_tmp$TRBV) & df_tmp$TRBV != "NA"]
  df_tmp_heat <- df_tmp[!is.na(df_tmp$TRBV) & df_tmp$TRBV != "NA" & !is.na(df_tmp$TRBJ), , drop = FALSE]
  if (nrow(df_tmp_heat) > 0) {
    df_tmp_heat$TRBV <- factor(df_tmp_heat$TRBV, levels = gtools::mixedsort(unique(TCR_stats$TRBV)))
    mat_heat <- table(df_tmp_heat$TRBV, df_tmp_heat$TRBJ) %>% prop.table() %>% as.data.frame()
    p_heat <- ggplot(mat_heat, aes(Var2, Var1, color = Freq, size = Freq)) +
      geom_point() +
      scale_color_gradient(low = "white", high = "red") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 270, hjust = 0))
    ggsave(fig_file("TCR/VJ_Gene_Usage", 1, paste0("TRBV_TRBJ_Heatmap_", grp)), p_heat, height = 8, width = 4.5)
  }
}

# Circlize Figure
for (grp in unique(TCR_stats$Group)) {
  df_tmp_circ <- TCR_stats_expand %>%
    dplyr::filter(Group == grp, !is.na(TRBJ), !is.na(TRBV), TRBJ != "NA", TRBV != "NA")
  if (nrow(df_tmp_circ) == 0) next
  mat_circ <- table(df_tmp_circ$TRBJ, df_tmp_circ$TRBV)
  if (nrow(mat_circ) > 0 && ncol(mat_circ) > 0 && sum(mat_circ) > 0) {
    pdf(fig_file("TCR/VJ_Gene_Usage", 2, paste0("TRBV_TRBJ_Circos_", grp)), width = 8, height = 8)
    circlize::circos.par(gap.after = c(rep(2, nrow(mat_circ) - 1), 10, rep(2, ncol(mat_circ) - 1), 10))
    circlize::chordDiagram(mat_circ, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat_circ))))))
    circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
      circlize::circos.text(CELL_META$xcenter, CELL_META$ylim[1] + circlize::cm_h(2), CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
      circlize::circos.axis(h = "top", labels.cex = 1, major.tick.length = 0.2, sector.index = CELL_META$sector.index, track.index = 2)
    }, bg.border = NA)
    title(paste0("Pairing of TRBV with TRBJ (", grp, ")"))
    circlize::circos.clear()
    dev.off()
  }
}

# seqlogo Figure
for (grp in Group_uniq) {
  df_seq <- TCR_stats_expand[TCR_stats_expand$Group == grp & !is.na(TCR_stats_expand$CDR3) & 
                               TCR_stats_expand$CDR3 != "NA", , drop = FALSE]
  df_seq$CDR3 <- stringr::str_sub(paste0(df_seq$CDR3, "--------"), 1, 16)
  bc_tmp <- unique(df_seq$barcode[df_seq$clone_ratio > 5 & df_seq$Freq > 10])
  if (length(bc_tmp) > 0) {
    fig_list_all <- lapply(bc_tmp, function(bc) {
      d <- df_seq[df_seq$barcode == bc, , drop = FALSE]
      patchwork::wrap_plots(lapply(seq_len(nrow(d)), function(k) {
        ggplot() +
          ggseqlogo::geom_logo(d$CDR3[k], method = "probability") +
          ggseqlogo::theme_logo() +
          ggtitle(paste(d$sample[k], d$celltype[k], d$TRBV[k], d$TRBJ[k], sep = ":"))
      }))
    })
    ggsave(fig_file("TCR/CDR3_Sequence", 2, paste0("CDR3_Seqlogo_", grp)), patchwork::wrap_plots(fig_list_all, ncol = 1), height = 15, width = 15)
  }
}

# 9. GLIPH2 抗原特异性聚类 -----------------------------------------------------
library(immGLIPH)
library(turboGliph)

gliph_data <- data.frame(
  chain = "TRB",
  CDR3b = TCR_stats$CDR3,
  TRBV = TCR_stats$TRBV,
  TRBJ = TCR_stats$TRBJ,
  patient = TCR_stats$sample,
  Group = TCR_stats$Group,
  counts = TCR_stats$Freq,
  celltype = TCR_stats$celltype
) %>%
  dplyr::filter(!is.na(CDR3b), CDR3b != "NA", !is.na(TRBV), !is.na(TRBJ))

# 耗时较长，请根据需要取消注释运行：
res <- gliph2(cdr3_sequences = gliph_data, sim_depth = 1000, n_cores = 1)
save(res, file = "./Result/TCR/GLIPH_HLA/Gliph_result.rda")

#load(file = "./Result/TCR/GLIPH_HLA/Gliph_result.rda")

table(res$cluster_properties$cluster_size)
cluster_prop <- res$cluster_properties[res$cluster_properties$cluster_size >= 2, , drop = FALSE]
nodes <- lapply(cluster_prop$members, function(x) unlist(stringr::str_split(x, " "))) %>%
  unlist() %>%
  as.data.frame() %>%
  dplyr::distinct()

colnames(nodes) <- "label"
meta_gliph <- do.call(rbind, res$cluster_list)
nodes <- cbind(nodes, meta_gliph[match(nodes$label, meta_gliph$CDR3b), ])

edges <- data.frame(from = res$connections$V1, to = res$connections$V2, weight = 1, type = res$connections$type) %>%
  dplyr::filter(from %in% nodes$label, to %in% nodes$label, from != to)

net_pc <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
grp_lv <- levels(factor(nodes$Group))
grp_cols <- setNames(GLOBAL_COLORS[seq_along(grp_lv)], grp_lv)
igraph::V(net_pc)$color <- grp_cols[as.character(nodes$Group)]

pdf(fig_file("TCR/GLIPH_HLA", 1, "GLIPH2_Network"), width = 10, height = 10)
plot(net_pc, vertex.size = 3, layout = igraph::layout_with_fr(net_pc), vertex.label.cex = 0.5, vertex.label.dist = 0.5, edge.color = "gray50", edge.curved = 0.1)
legend("topleft", grp_lv, pch = 21, pt.bg = grp_cols)
dev.off()


# ======================= 多组：TCR扩增、多样性、TRBV差异 =======================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggpubr)
  library(stringr); library(purrr); library(vegan); library(scales)
  library(gtools); library(ggrepel); library(tibble)
})

dir.create("./Result/TCR/Statistics", showWarnings = FALSE, recursive = TRUE)
dir.create("./Figures", showWarnings = FALSE, recursive = TRUE)

grp_use <- c("NR", "TCMR", "ABMR", "MIX")
ref_group <- "NR"
sample_col <- "HLA_Splitter"
group_col <- "Group"
high_cutoff <- 1
pair_comp <- combn(grp_use, 2, simplify = FALSE)

COL_GRP <- if (exists("GROUP_COLORS")) GROUP_COLORS[grp_use] else
  setNames(scales::hue_pal()(length(grp_use)), grp_use)
if (any(is.na(COL_GRP))) COL_GRP[is.na(COL_GRP)] <- scales::hue_pal()(sum(is.na(COL_GRP)))

theme_tcr <- function(base_size = 12) {
  theme_classic(base_size) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black", face = "bold"),
          strip.background = element_rect(fill = "grey95", color = "black", linewidth = .4),
          strip.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          plot.title = element_text(hjust = .5, face = "bold"))
}

sig_lab <- function(p) case_when(
  p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"
)

kw_test <- function(dat, val, by = NULL) {
  dat <- dat %>% filter(Group %in% grp_use, !is.na(.data[[val]])) %>%
    mutate(Group = factor(Group, levels = grp_use))
  f <- function(d) tibble(p = tryCatch(kruskal.test(d[[val]] ~ d$Group)$p.value, error = \(e) NA_real_))
  if (is.null(by)) f(dat) %>% mutate(p.adj = p, sig = sig_lab(p.adj)) else
    dat %>% group_by(across(all_of(by))) %>% group_modify(~f(.x)) %>% ungroup() %>%
    mutate(p.adj = p.adjust(p, "BH"), sig = sig_lab(p.adj))
}

pair_wilcox <- function(dat, val, by = NULL) {
  dat <- dat %>% filter(Group %in% grp_use, !is.na(.data[[val]])) %>%
    mutate(Group = factor(Group, levels = grp_use))
  f <- function(d) map_dfr(pair_comp, \(cc) {
    dd <- d %>% filter(Group %in% cc) %>% mutate(Group = droplevels(Group))
    tibble(group1 = cc[1], group2 = cc[2],
           p = tryCatch(wilcox.test(dd[[val]] ~ dd$Group)$p.value, error = \(e) NA_real_))
  }) %>% mutate(p.adj = p.adjust(p, "BH"), sig = sig_lab(p.adj))
  if (is.null(by)) f(dat) else dat %>% group_by(across(all_of(by))) %>% group_modify(~f(.x)) %>% ungroup()
}

# 1. 高扩增TCR比例 --------------------------------------------------------------
df_exp <- Tcell@meta.data %>%
  rownames_to_column("barcode") %>%
  filter(.data[[group_col]] %in% grp_use, !is.na(clonalFrequency), clonalFrequency > 0) %>%
  group_by(Sample = .data[[sample_col]], Group = .data[[group_col]]) %>%
  summarise(TCR_cells = n(),
            HighExpanded_cells = sum(clone_ratio > high_cutoff, na.rm = TRUE),
            HighExpanded_ratio = HighExpanded_cells / TCR_cells,
            .groups = "drop") %>%
  mutate(Group = factor(Group, grp_use))

stat_exp_kw <- kw_test(df_exp, "HighExpanded_ratio")
stat_exp_pair <- pair_wilcox(df_exp, "HighExpanded_ratio")

write.csv(df_exp, "./Result/TCR/Statistics/high_expanded_TCR_ratio_by_sample.csv", row.names = FALSE)
write.csv(stat_exp_kw, "./Result/TCR/Statistics/high_expanded_TCR_ratio_kruskal.csv", row.names = FALSE)
write.csv(stat_exp_pair, "./Result/TCR/Statistics/high_expanded_TCR_ratio_pairwise_wilcox_BH.csv", row.names = FALSE)

p_exp <- ggplot(df_exp, aes(Group, HighExpanded_ratio, fill = Group, color = Group)) +
  geom_boxplot(width = .55, outlier.shape = NA, fill = "white", linewidth = .7) +
  geom_point(position = position_jitter(.12), size = 2.6, alpha = .85) +
  stat_compare_means(comparisons = comparison, label = "p.format", size = 4) +
  scale_y_continuous(labels = percent_format(1), expand = expansion(mult = c(.03, .18))) +
  scale_fill_manual(values = COL_GRP) + scale_color_manual(values = COL_GRP) +
  labs(title = paste0("High-expanded TCR clones: clone_ratio > ", high_cutoff),
       x = NULL, y = "Proportion of TCR-positive cells") +
  theme_tcr() + theme(legend.position = "none")

ggsave(fig_file("TCR/Clonality_Expansion", 10, "High_Expanded_TCR_Ratio_by_Sample"), p_exp, width = 5, height = 5)

# 2. TCR多样性 -----------------------------------------------------------------
calc_tcr_div <- function(seu, sample_col, group_col, clone_col = "CTstrict") {
  md <- seu@meta.data
  if (!clone_col %in% colnames(md)) clone_col <- intersect(c("CTstrict", "CTaa", "CTnt", "CTgene"), colnames(md))[1]
  if (is.na(clone_col)) stop("找不到克隆列：CTstrict / CTaa / CTnt / CTgene")
  md <- md[md[[group_col]] %in% grp_use & !is.na(md[[sample_col]]) &
             !is.na(md[[clone_col]]) & md[[clone_col]] != "NA", , drop = FALSE]
  map_dfr(split(md, md[[sample_col]]), \(x) {
    tab <- table(x[[clone_col]])
    tibble(Sample = unique(x[[sample_col]])[1],
           Group = unique(x[[group_col]])[1],
           Richness = length(tab),
           Shannon = diversity(tab, "shannon"),
           InvSimpson = diversity(tab, "invsimpson"),
           GiniSimpson = diversity(tab, "simpson"),
           NormalizedEntropy = ifelse(length(tab) > 1, diversity(tab, "shannon") / log(length(tab)), 0),
           Chao1 = tryCatch(as.numeric(estimateR(tab)["S.chao1"]), error = \(e) NA_real_))
  })
}

df_div <- calc_tcr_div(Tcell_recall, sample_col, group_col) %>%
  mutate(Group = factor(Group, grp_use)) %>%
  pivot_longer(c(Richness, Shannon, InvSimpson, GiniSimpson, NormalizedEntropy, Chao1),
               names_to = "Metric", values_to = "Value")

stat_div_kw <- kw_test(df_div, "Value", "Metric")
stat_div_pair <- pair_wilcox(df_div, "Value", "Metric")

write.csv(df_div, "./Result/TCR/Statistics/TCR_diversity_by_sample.csv", row.names = FALSE)
write.csv(stat_div_kw, "./Result/TCR/Statistics/TCR_diversity_kruskal_BH.csv", row.names = FALSE)
write.csv(stat_div_pair, "./Result/TCR/Statistics/TCR_diversity_pairwise_wilcox_BH.csv", row.names = FALSE)

p_div <- ggplot(df_div, aes(Group, Value, fill = Group, color = Group)) +
  geom_boxplot(width = .55, outlier.shape = NA, fill = "white", linewidth = .65) +
  geom_point(position = position_jitter(.12), size = 2.1, alpha = .8) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3.5) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = COL_GRP) + scale_color_manual(values = COL_GRP) +
  labs(title = "TCR repertoire diversity", x = NULL, y = "Diversity index") +
  theme_tcr() + theme(legend.position = "none")

ggsave(fig_file("TCR/Clonality_Expansion", 11, "TCR_Repertoire_Diversity"), p_div, width = 9, height = 6)

# 3. TRBV使用频率多组差异 -------------------------------------------------------
if (!exists("TCR_stats")) {
  TCR_stats <- Tcell_recall@meta.data %>%
    rownames_to_column("barcode") %>%
    transmute(barcode, gene = CTgene, Group = .data[[group_col]],
              sample = .data[[sample_col]], clone_ratio,
              CDR3 = CTaa, Freq = clonalFrequency) %>%
    separate(gene, into = c("TRA_gene", "TRB_gene"), sep = "_",
             remove = FALSE, fill = "right") %>%
    mutate(TRBV = str_extract(TRB_gene, "TRBV[^.]*")) %>%
    filter(!is.na(TRBV), TRBV != "NA")
}

sample_info_tcr <- TCR_stats %>%
  filter(Group %in% grp_use, !is.na(sample), !is.na(TRBV), TRBV != "NA") %>%
  distinct(sample, Group)

trbv_usage <- TCR_stats %>%
  filter(Group %in% grp_use, !is.na(sample), !is.na(TRBV), TRBV != "NA") %>%
  dplyr::count(sample, Group, TRBV, name = "n") %>%
  group_by(sample) %>%
  mutate(total = sum(n), freq = n / total) %>%
  ungroup()

all_trbv <- mixedsort(unique(trbv_usage$TRBV))

trbv_usage_full <- expand_grid(sample = unique(sample_info_tcr$sample), TRBV = all_trbv) %>%
  left_join(sample_info_tcr, by = "sample") %>%
  left_join(select(trbv_usage, sample, TRBV, n, total, freq), by = c("sample", "TRBV")) %>%
  mutate(n = replace_na(n, 0), freq = replace_na(freq, 0),
         Group = factor(Group, grp_use))

trbv_mean <- trbv_usage_full %>%
  group_by(TRBV, Group) %>%
  summarise(mean_freq = mean(freq, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_freq, names_prefix = "mean_")

trbv_stat_kw <- kw_test(trbv_usage_full, "freq", "TRBV") %>%
  left_join(trbv_mean, by = "TRBV") %>%
  arrange(p.adj)

trbv_stat_pair <- pair_wilcox(trbv_usage_full, "freq", "TRBV")

write.csv(trbv_usage_full, "./Result/TCR/Statistics/TRBV_usage_by_sample.csv", row.names = FALSE)
write.csv(trbv_stat_kw, "./Result/TCR/Statistics/TRBV_usage_kruskal_BH.csv", row.names = FALSE)
write.csv(trbv_stat_pair, "./Result/TCR/Statistics/TRBV_usage_pairwise_wilcox_BH.csv", row.names = FALSE)

top_trbv <- trbv_stat_kw %>%
  filter(!is.na(p)) %>%
  arrange(p) %>%
  slice_head(n = 12) %>%
  pull(TRBV)

trbv_plot_df <- trbv_usage_full %>%
  filter(TRBV %in% top_trbv) %>%
  mutate(TRBV = factor(TRBV, mixedsort(top_trbv)))

p_trbv_box <- ggplot(trbv_plot_df, aes(Group, freq, fill = Group, color = Group)) +
  geom_boxplot(width = .55, outlier.shape = NA, fill = "white", linewidth = .55) +
  geom_point(position = position_jitter(.12), size = 1.8, alpha = .8) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3) +
  facet_wrap(~TRBV, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = percent_format(.1), expand = expansion(mult = c(.03, .25))) +
  scale_fill_manual(values = COL_GRP) + scale_color_manual(values = COL_GRP) +
  labs(title = "Differential TRBV usage across groups", x = NULL, y = "TRBV usage frequency") +
  theme_tcr(11) + theme(legend.position = "none")

ggsave(fig_file("TCR/VJ_Gene_Usage", 3, "Top12_TRBV_Usage_Boxplot"), p_trbv_box, width = 10, height = 8)

# 4. TRBV平均使用频率热图 -------------------------------------------------------
trbv_heat_df <- trbv_usage_full %>%
  filter(TRBV %in% top_trbv) %>%
  group_by(Group, TRBV) %>%
  summarise(mean_freq = mean(freq, na.rm = TRUE), .groups = "drop") %>%
  mutate(TRBV = factor(TRBV, rev(mixedsort(top_trbv))),
         Group = factor(Group, grp_use))

p_trbv_heat <- ggplot(trbv_heat_df, aes(Group, TRBV, fill = mean_freq)) +
  geom_tile(color = "white", linewidth = .6) +
  geom_text(aes(label = percent(mean_freq, accuracy = .1)), size = 3) +
  scale_fill_gradient(low = "white", high = "#B2182B", labels = percent_format(1)) +
  labs(title = "Mean TRBV usage frequency", x = NULL, y = NULL, fill = "Mean frequency") +
  theme_tcr() +
  theme(axis.text.x = element_text(face = "bold"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = .5))

ggsave(fig_file("TCR/VJ_Gene_Usage", 4, "Top12_TRBV_Usage_Heatmap"), p_trbv_heat, width = 5.8, height = 6.2)


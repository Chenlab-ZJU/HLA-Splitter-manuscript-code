#!/usr/bin/env Rscript
# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis.
# Original source: 13_revision_ABMR_vs_TCMR_Tcell_DE_function.R
# Figure/table scope: Fig5; FigS9; Supplementary Table 18; Methods CXCR3-axis analysis
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(clusterProfiler)
  library(Matrix)
  library(scales)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
urine_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(urine_dir, "..", ".."), winslash = "/", mustWork = TRUE)
revision_dir <- file.path(project_dir, "Revision")
out_dir <- file.path(revision_dir, "revised figure")
table_dir <- file.path(out_dir, "tables")
panel_dir <- file.path(out_dir, "panels", "ABMR_TCMR_Tcell_function")
result_dir <- file.path(urine_dir, "Result", "Reviewer_Tcell_ABMR_TCMR")
msigdb_dir <- file.path(project_dir, "Data", "MSigDB")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

GROUP_COLORS <- c("TCMR" = "#699ECA", "ABMR" = "#4DAF4A")
DIRECTION_COLORS <- c("ABMR higher" = "#4DAF4A", "TCMR higher" = "#699ECA", "Not significant" = "grey78")
SIGNATURE_LABELS <- c(
  "Cytotoxicity" = "Cytotoxicity",
  "IFNG_response" = "IFNG response",
  "Antigen_presentation_MHCII" = "Antigen presentation\n(MHC II)",
  "Exhaustion_checkpoint" = "Exhaustion/checkpoint",
  "Treg_program" = "Treg program",
  "Proliferation" = "Proliferation",
  "Naive_memory" = "Naive/memory"
)

theme_nc <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      axis.text = element_text(colour = "black", size = base_size),
      axis.title = element_text(colour = "black", size = base_size + 1),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.30, colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1.5),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.28, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(3, 3, 3, 3)
    )
}

save_panel <- function(plot, stem, width, height) {
  pdf_file <- file.path(panel_dir, paste0(stem, ".pdf"))
  png_file <- file.path(panel_dir, paste0(stem, ".png"))
  ggsave(pdf_file, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(png_file, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  invisible(c(pdf_file, png_file))
}

get_data_matrix <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "data")
  )
}

detect_fc_col <- function(tab) {
  cols <- c("avg_log2FC", "avg_logFC", "log2FC")
  hit <- intersect(cols, colnames(tab))
  if (length(hit) == 0) stop("Cannot find a log-fold-change column in marker table.")
  hit[1]
}

marker_table <- function(obj, ident.1, ident.2, group.by = NULL, min.pct = 0.05, logfc.threshold = 0) {
  if (is.null(group.by)) {
    Idents(obj) <- ident.1
    out <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct,
                       logfc.threshold = logfc.threshold, test.use = "wilcox")
  } else {
    Idents(obj) <- group.by
    out <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct,
                       logfc.threshold = logfc.threshold, test.use = "wilcox")
  }
  out %>%
    tibble::rownames_to_column("gene") %>%
    as_tibble()
}

add_marker_fields <- function(tab, fc_col = detect_fc_col(tab), fc_cutoff = 0.25, padj_cutoff = 0.05) {
  tab %>%
    mutate(
      avg_log2FC = .data[[fc_col]],
      direction = case_when(
        p_val_adj < padj_cutoff & avg_log2FC >= fc_cutoff ~ "ABMR higher",
        p_val_adj < padj_cutoff & avg_log2FC <= -fc_cutoff ~ "TCMR higher",
        TRUE ~ "Not significant"
      ),
      neg_log10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
      rank_score = avg_log2FC * neg_log10_padj
    ) %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)))
}

score_gene_set <- function(expr, gene_set) {
  genes <- intersect(gene_set, rownames(expr))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(expr)))
  }
  Matrix::colMeans(expr[genes, , drop = FALSE])
}

safe_wilcox_p <- function(values, groups, min_groups = 2) {
  keep <- !is.na(values) & !is.na(groups)
  values <- values[keep]
  groups <- droplevels(factor(groups[keep]))
  if (length(unique(groups)) != 2) return(NA_real_)
  if (any(table(groups) < min_groups)) return(NA_real_)
  wilcox.test(values ~ groups, exact = FALSE)$p.value
}

read_go_gmt <- function(msigdb_dir) {
  gmt_file <- file.path(msigdb_dir, "c5.go.bp.v2023.2.Hs.symbols.gmt")
  if (!file.exists(gmt_file)) return(NULL)
  gmt <- clusterProfiler::read.gmt(gmt_file)
  colnames(gmt) <- c("term", "gene")
  gmt
}

run_enrichment <- function(genes, universe, term2gene) {
  if (is.null(term2gene) || length(genes) < 5) return(tibble())
  enrich <- tryCatch(
    clusterProfiler::enricher(
      gene = unique(genes),
      universe = universe,
      TERM2GENE = term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = 5,
      maxGSSize = 500
    ),
    error = function(e) NULL
  )
  if (is.null(enrich)) return(tibble())
  as_tibble(enrich@result) %>%
    arrange(p.adjust) %>%
    mutate(
      Direction = if_else(grepl("^ABMR", Description), "ABMR higher", "TCMR higher")
    )
}

message("Loading T cell Seurat object...")
tcell_rds <- file.path(result_dir, "Tcell_only_ABMR_TCMR_reclustered.rds")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(tcell_rds)
DefaultAssay(tcell) <- "RNA"

needed_meta <- c("Group", "Sample", "Reviewer_Celltype")
missing_meta <- setdiff(needed_meta, colnames(tcell@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata columns: ", paste(missing_meta, collapse = ", "))
}

tcell <- subset(tcell, subset = Group %in% c("ABMR", "TCMR"))
tcell$Group <- factor(as.character(tcell$Group), levels = c("TCMR", "ABMR"))
tcell$Reviewer_Celltype <- as.character(tcell$Reviewer_Celltype)
tcell$Sample <- as.character(tcell$Sample)

qc_counts <- tcell@meta.data %>%
  count(Group, Sample, Reviewer_Celltype, name = "n_cells") %>%
  arrange(Group, Sample, desc(n_cells))
# Export the analysis source data for figure reproduction.
write_csv(qc_counts, file.path(result_dir, "ABMR_vs_TCMR_Tcell_QC_cell_counts_by_sample_celltype.csv"))
write_csv(qc_counts, file.path(table_dir, "ABMR_vs_TCMR_Tcell_QC_cell_counts_by_sample_celltype.csv"))

message("Running overall ABMR vs TCMR differential expression...")
overall_de <- marker_table(tcell, ident.1 = "ABMR", ident.2 = "TCMR", group.by = "Group") %>%
  add_marker_fields()
write_csv(overall_de, file.path(result_dir, "ABMR_vs_TCMR_Tcell_DE_full.csv"))

top_de <- bind_rows(
  overall_de %>% filter(direction == "ABMR higher") %>% arrange(p_val_adj, desc(avg_log2FC)) %>% slice_head(n = 50),
  overall_de %>% filter(direction == "TCMR higher") %>% arrange(p_val_adj, avg_log2FC) %>% slice_head(n = 50)
) %>%
  select(gene, direction, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)
write_csv(top_de, file.path(table_dir, "ABMR_vs_TCMR_Tcell_DE_top_genes.csv"))

cytotoxic_genes <- c(
  "NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH", "GZMK", "CST7",
  "CTSW", "IFNG", "FASLG", "KLRD1", "KLRK1", "KLRB1", "CCL5"
)

cytotoxic_de <- overall_de %>%
  filter(gene %in% cytotoxic_genes) %>%
  arrange(match(gene, cytotoxic_genes)) %>%
  select(gene, direction, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)
write_csv(cytotoxic_de, file.path(result_dir, "ABMR_vs_TCMR_Tcell_cytotoxic_gene_DE.csv"))
write_csv(cytotoxic_de, file.path(table_dir, "ABMR_vs_TCMR_Tcell_cytotoxic_gene_DE.csv"))

message("Calculating functional module scores...")
expr <- get_data_matrix(tcell, assay = "RNA")
gene_sets <- list(
  Cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH", "GZMK", "CST7", "CTSW", "IFNG", "FASLG", "CCL5"),
  IFNG_response = c("IFNG", "CXCL9", "CXCL10", "CXCL11", "STAT1", "IRF1", "GBP1", "GBP2", "GBP5", "TAP1", "TAP2", "PSMB8", "PSMB9"),
  Antigen_presentation_MHCII = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "CIITA"),
  Exhaustion_checkpoint = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX", "ENTPD1", "CXCL13", "LAYN"),
  Treg_program = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TNFRSF4", "TIGIT", "CCR8"),
  Proliferation = c("MKI67", "TOP2A", "TYMS", "STMN1", "HMGB2", "UBE2C", "CENPF"),
  Naive_memory = c("CCR7", "IL7R", "LEF1", "TCF7", "SELL", "LTB", "MAL")
)

cell_scores <- tcell@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, Group, Sample, Reviewer_Celltype)
for (score_name in names(gene_sets)) {
  cell_scores[[score_name]] <- score_gene_set(expr, gene_sets[[score_name]])
}

score_long <- cell_scores %>%
  pivot_longer(cols = all_of(names(gene_sets)), names_to = "Signature", values_to = "Score")
write_csv(score_long, file.path(result_dir, "ABMR_vs_TCMR_Tcell_function_scores_by_cell.csv"))

score_sample_celltype_summary <- score_long %>%
  group_by(Group, Sample, Reviewer_Celltype, Signature) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(Score, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(score_sample_celltype_summary, file.path(result_dir, "ABMR_vs_TCMR_Tcell_function_signature_sample_celltype_summary.csv"))
write_csv(score_sample_celltype_summary, file.path(table_dir, "ABMR_vs_TCMR_Tcell_function_signature_sample_celltype_summary.csv"))

score_sample_summary <- score_long %>%
  group_by(Group, Sample, Signature) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(Score, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(score_sample_summary, file.path(result_dir, "ABMR_vs_TCMR_Tcell_function_signature_sample_summary.csv"))

score_group_summary <- score_sample_summary %>%
  group_by(Group, Signature) %>%
  summarise(
    n_samples = n_distinct(Sample),
    mean_of_sample_means = mean(mean_score, na.rm = TRUE),
    median_of_sample_means = median(mean_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Group, values_from = c(n_samples, mean_of_sample_means, median_of_sample_means)) %>%
  mutate(
    ABMR_minus_TCMR_mean_score = mean_of_sample_means_ABMR - mean_of_sample_means_TCMR
  ) %>%
  arrange(desc(ABMR_minus_TCMR_mean_score))
write_csv(score_group_summary, file.path(table_dir, "ABMR_vs_TCMR_Tcell_function_signature_group_summary.csv"))

score_tests <- score_sample_summary %>%
  group_by(Signature) %>%
  summarise(
    n_TCMR = n_distinct(Sample[Group == "TCMR"]),
    n_ABMR = n_distinct(Sample[Group == "ABMR"]),
    median_TCMR = median(mean_score[Group == "TCMR"], na.rm = TRUE),
    median_ABMR = median(mean_score[Group == "ABMR"], na.rm = TRUE),
    ABMR_minus_TCMR_median = median_ABMR - median_TCMR,
    p_value = safe_wilcox_p(mean_score, Group),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj, desc(abs(ABMR_minus_TCMR_median)))
write_csv(score_tests, file.path(table_dir, "ABMR_vs_TCMR_Tcell_function_signature_wilcox_tests.csv"))

score_celltype_tests <- score_sample_celltype_summary %>%
  group_by(Reviewer_Celltype, Signature) %>%
  summarise(
    n_TCMR = n_distinct(Sample[Group == "TCMR"]),
    n_ABMR = n_distinct(Sample[Group == "ABMR"]),
    median_TCMR = median(mean_score[Group == "TCMR"], na.rm = TRUE),
    median_ABMR = median(mean_score[Group == "ABMR"], na.rm = TRUE),
    ABMR_minus_TCMR_median = median_ABMR - median_TCMR,
    p_value = safe_wilcox_p(mean_score, Group),
    .groups = "drop"
  ) %>%
  group_by(Signature) %>%
  mutate(p_adj_within_signature = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(Signature, p_adj_within_signature, desc(abs(ABMR_minus_TCMR_median)))
write_csv(score_celltype_tests, file.path(table_dir, "ABMR_vs_TCMR_Tcell_function_signature_celltype_wilcox_tests.csv"))

message("Running subset-level differential expression...")
subset_de <- list()
celltype_counts <- tcell@meta.data %>%
  count(Reviewer_Celltype, Group, name = "n_cells") %>%
  pivot_wider(names_from = Group, values_from = n_cells, values_fill = 0)

eligible_celltypes <- celltype_counts %>%
  filter(ABMR >= 20, TCMR >= 20) %>%
  pull(Reviewer_Celltype)

for (ct in eligible_celltypes) {
  message("  Cell type: ", ct)
  sub_obj <- subset(tcell, subset = Reviewer_Celltype == ct)
  tab <- tryCatch(
    marker_table(sub_obj, ident.1 = "ABMR", ident.2 = "TCMR", group.by = "Group", min.pct = 0.05, logfc.threshold = 0.1) %>%
      add_marker_fields() %>%
      mutate(Reviewer_Celltype = ct, .before = gene),
    error = function(e) {
      message("    skipped: ", e$message)
      tibble()
    }
  )
  subset_de[[ct]] <- tab
}

subset_de_all <- bind_rows(subset_de)
write_csv(subset_de_all, file.path(result_dir, "ABMR_vs_TCMR_Tcell_subset_DE_full.csv"))

subset_de_top <- bind_rows(
  subset_de_all %>% filter(direction == "ABMR higher") %>% group_by(Reviewer_Celltype) %>% arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup(),
  subset_de_all %>% filter(direction == "TCMR higher") %>% group_by(Reviewer_Celltype) %>% arrange(p_val_adj, avg_log2FC, .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup()
) %>%
  arrange(Reviewer_Celltype, direction, p_val_adj) %>%
  select(Reviewer_Celltype, gene, direction, avg_log2FC, pct.1, pct.2, p_val_adj)
write_csv(subset_de_top, file.path(table_dir, "ABMR_vs_TCMR_Tcell_subset_DE_top_genes.csv"))

message("Running GO biological process enrichment...")
go_term2gene <- read_go_gmt(msigdb_dir)
universe_genes <- rownames(expr)
abmr_up_genes <- overall_de %>% filter(direction == "ABMR higher") %>% pull(gene)
tcmr_up_genes <- overall_de %>% filter(direction == "TCMR higher") %>% pull(gene)

go_abmr <- run_enrichment(abmr_up_genes, universe_genes, go_term2gene) %>%
  mutate(Direction = "ABMR higher", .before = 1) %>%
  select(Direction, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count) %>%
  slice_head(n = 30)
go_tcmr <- run_enrichment(tcmr_up_genes, universe_genes, go_term2gene) %>%
  mutate(Direction = "TCMR higher", .before = 1) %>%
  select(Direction, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count) %>%
  slice_head(n = 30)

write_csv(go_abmr, file.path(table_dir, "ABMR_vs_TCMR_Tcell_GO_BP_enrichment_ABMR_higher.csv"))
write_csv(go_tcmr, file.path(table_dir, "ABMR_vs_TCMR_Tcell_GO_BP_enrichment_TCMR_higher.csv"))

message("Drawing reviewer-quality summary panels...")
volcano_labels <- overall_de %>%
  filter(direction != "Not significant") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 18)

p_volcano <- overall_de %>%
  mutate(direction = factor(direction, levels = c("TCMR higher", "Not significant", "ABMR higher"))) %>%
  ggplot(aes(avg_log2FC, neg_log10_padj, colour = direction)) +
  geom_point(size = 0.55, alpha = 0.55, stroke = 0) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", linewidth = 0.25, colour = "grey45") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25, colour = "grey45") +
  ggrepel::geom_text_repel(
    data = volcano_labels,
    aes(label = gene),
    size = 2.1,
    colour = "black",
    max.overlaps = Inf,
    box.padding = 0.22,
    min.segment.length = 0,
    segment.size = 0.2
  ) +
  scale_colour_manual(values = DIRECTION_COLORS, drop = FALSE) +
  labs(x = "Average log2 fold change (ABMR / TCMR)", y = expression(-log[10]("adjusted P")), colour = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.0, alpha = 1))) +
  theme_nc(8) +
  theme(legend.position = "top")
save_panel(p_volcano, "ABMR_vs_TCMR_Tcell_overall_DE_volcano", 3.7, 3.2)

cyto_plot_tab <- cytotoxic_de %>%
  mutate(
    gene = factor(gene, levels = rev(cytotoxic_genes[cytotoxic_genes %in% gene])),
    direction = factor(direction, levels = c("TCMR higher", "Not significant", "ABMR higher")),
    signif_label = case_when(
      p_val_adj < 0.001 ~ "***",
      p_val_adj < 0.01 ~ "**",
      p_val_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p_cyto <- ggplot(cyto_plot_tab, aes(avg_log2FC, gene, colour = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  geom_segment(aes(x = 0, xend = avg_log2FC, yend = gene), linewidth = 0.55, alpha = 0.8) +
  geom_point(size = 2.0) +
  geom_text(aes(label = signif_label), colour = "black", hjust = -0.1, size = 2.4) +
  scale_colour_manual(values = DIRECTION_COLORS, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(x = "Average log2 fold change (ABMR / TCMR)", y = NULL, colour = NULL) +
  theme_nc(8) +
  theme(legend.position = "top")
save_panel(p_cyto, "ABMR_vs_TCMR_Tcell_cytotoxic_gene_lollipop", 3.2, 2.8)

sig_order <- score_tests %>% arrange(desc(ABMR_minus_TCMR_median)) %>% pull(Signature)
p_scores <- score_sample_summary %>%
  mutate(Signature = factor(Signature, levels = sig_order)) %>%
  ggplot(aes(Group, mean_score, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35, alpha = 0.72) +
  geom_jitter(width = 0.12, size = 0.9, alpha = 0.75, colour = "black") +
  facet_wrap(~Signature, scales = "free_y", ncol = 4, labeller = as_labeller(SIGNATURE_LABELS)) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(x = NULL, y = "Sample-level mean module score") +
  theme_nc(7) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(p_scores, "ABMR_vs_TCMR_Tcell_function_signature_scores", 6.8, 3.8)

go_plot_tab <- bind_rows(go_abmr, go_tcmr) %>%
  group_by(Direction) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(
    Description = gsub("^GOBP_", "", Description),
    Description = gsub("_", " ", Description),
    Description = stringr::str_to_sentence(tolower(Description)),
    Description = factor(Description, levels = rev(unique(Description))),
    Direction = factor(Direction, levels = c("TCMR higher", "ABMR higher"))
  )

if (nrow(go_plot_tab) > 0) {
  p_go <- ggplot(go_plot_tab, aes(-log10(p.adjust), Description, fill = Direction)) +
    geom_col(width = 0.72) +
    facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = DIRECTION_COLORS) +
    labs(x = expression(-log[10]("adjusted P")), y = NULL) +
    theme_nc(7) +
    theme(legend.position = "none", strip.text.y = element_text(angle = 0))
  save_panel(p_go, "ABMR_vs_TCMR_Tcell_GO_BP_enrichment", 5.8, 4.4)
}

message("Done. Tables were written to: ", table_dir)
message("Panels were written to: ", panel_dir)

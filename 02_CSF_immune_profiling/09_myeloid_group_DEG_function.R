# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Required intermediate source-table / annotation generation.
# Original source: 10_revision_myeloid_group_DEG_function.R
# Figure/table scope: Fig2; FigS4; FigS5
# Data cohort: csf. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(enrichR)
})

set.seed(20260813)

env_project_dir <- Sys.getenv("CSF_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) {
    sub("^--file=", "", file_arg[1])
  } else {
    file.path(getwd(), "Code", "10_revision_myeloid_group_DEG_function.R")
  }
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}
figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
group_colors <- c(
  "Metastases" = "#E42320",
  "Hemorrhage" = "#699ECA",
  "Encephalitis" = "#4DAF4A"
)

genes_of_interest <- c(
  "TREM2", "SPP1", "APOE", "APOC1", "C1QA", "C1QB", "C1QC",
  "CXCL9", "CXCL10", "CXCL11", "IL1B", "CCL20", "HMOX1",
  "CD163", "MRC1", "MSR1", "LGALS3", "LST1", "S100A8", "S100A9"
)

theme_pub <- function(base_size = 9.5) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

save_plot_safe <- function(plot, filename, width, height) {
  target_pdf <- file.path(figure_dir, paste0(filename, ".pdf"))
  target_png <- file.path(figure_dir, paste0(filename, ".png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
  file.copy(tmp_pdf, target_pdf, overwrite = TRUE)
  file.copy(tmp_png, target_png, overwrite = TRUE)
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", signif(p, 2))
  )
}

# Read the saved object or source table used by this analysis.
Myeloid <- readRDS(file.path(getOption("hla.data_root"), "csf", "CSF_Myeloid_revision.rds"))
DefaultAssay(Myeloid) <- "RNA"
Myeloid$Group <- factor(as.character(Myeloid$Group), levels = group_levels)

genes_present <- intersect(genes_of_interest, rownames(Myeloid))
expr_mat <- FetchData(Myeloid, vars = genes_present) %>%
  as_tibble(rownames = "cell")
meta_expr <- Myeloid@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Sample = as.character(.data[["HLA_Splitter"]]),
    Group = factor(as.character(.data[["Group"]]), levels = group_levels),
    Myeloid_state = as.character(.data[["Myeloid_celltype_revision"]])
  )
expr <- left_join(meta_expr, expr_mat, by = "cell")

cell_expr_long <- expr %>%
  pivot_longer(cols = all_of(genes_present), names_to = "gene", values_to = "expression")
sample_expr <- cell_expr_long %>%
  group_by(Sample, Group, gene) %>%
  summarise(
    n_cells = n(),
    mean_expression = mean(expression, na.rm = TRUE),
    pct_expressing = mean(expression > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )
# Export the analysis source data for figure reproduction.
write_csv(sample_expr, file.path(getOption("hla.data_root"), "csf", "Myeloid_revision_key_gene_sample_expression.csv"))

gene_kw <- sample_expr %>%
  group_by(gene) %>%
  summarise(
    n_groups = n_distinct(Group),
    p.value = ifelse(n_groups >= 2, kruskal.test(mean_expression ~ Group)$p.value, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH"), label = format_p(p.value))
write_csv(gene_kw, file.path(result_dir, "Myeloid_revision_key_gene_KW_tests.csv"))

key_plot_genes <- intersect(c("TREM2", "SPP1", "APOE", "C1QA", "CXCL9", "CXCL10", "IL1B", "HMOX1"), genes_present)
p_key_expr <- sample_expr %>%
  filter(gene %in% key_plot_genes) %>%
  mutate(gene = factor(gene, levels = key_plot_genes), Group = factor(as.character(Group), levels = group_levels)) %>%
  ggplot(aes(Group, mean_expression, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.82, linewidth = 0.38) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 20260813),
    shape = 21,
    colour = "black",
    stroke = 0.20,
    size = 1.45,
    alpha = 0.90
  ) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(x = NULL, y = "Mean expression per sample") +
  theme_pub(9.6) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.0),
    strip.text = element_text(size = 9.4),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )
save_plot_safe(p_key_expr, "Myeloid_key_gene_expression_by_group", width = 6.8, height = 4.7)

Idents(Myeloid) <- "Group"
comparisons <- list(
  "Metastases_vs_Hemorrhage" = c("Metastases", "Hemorrhage"),
  "Metastases_vs_Encephalitis" = c("Metastases", "Encephalitis"),
  "Hemorrhage_vs_Encephalitis" = c("Hemorrhage", "Encephalitis")
)

deg_list <- lapply(names(comparisons), function(nm) {
  pair <- comparisons[[nm]]
  deg <- FindMarkers(
    Myeloid,
    ident.1 = pair[1],
    ident.2 = pair[2],
    logfc.threshold = 0.10,
    min.pct = 0.10,
    test.use = "wilcox",
    verbose = FALSE
  )
  deg %>%
    as_tibble(rownames = "gene") %>%
    mutate(comparison = nm, group_up = ifelse(avg_log2FC > 0, pair[1], pair[2])) %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)))
})
deg_all <- bind_rows(deg_list)
write_csv(deg_all, file.path(result_dir, "Myeloid_revision_group_DEG_all.csv"))

for (nm in names(comparisons)) {
  write_csv(filter(deg_all, comparison == nm), file.path(result_dir, paste0("Myeloid_revision_group_DEG_", nm, ".csv")))
}

myeloid_function_modules <- list(
  "Microglia/phagocytosis" = c("TREM2", "SPP1", "APOE", "APOC1", "LGALS3", "LPL", "GPNMB"),
  "Complement" = c("C1QA", "C1QB", "C1QC", "C3", "C5AR1", "FCGR3A"),
  "Inflammation/chemokine" = c("IL1B", "CXCL1", "CXCL2", "CXCL3", "CXCL9", "CXCL10", "CXCL11", "CCL20", "CCL8"),
  "Antigen presentation" = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CD74", "CIITA"),
  "Monocyte activation" = c("S100A8", "S100A9", "S100A12", "FCN1", "LST1", "VCAN", "AIF1"),
  "Repair/M2-like" = c("HMOX1", "CD163", "MRC1", "MSR1", "PPARG", "CLEC10A", "SELENOP")
)
module_df <- tibble(
  module = rep(names(myeloid_function_modules), lengths(myeloid_function_modules)),
  gene = unlist(myeloid_function_modules, use.names = FALSE)
) %>%
  filter(gene %in% rownames(Myeloid)) %>%
  distinct(module, gene) %>%
  group_by(module) %>%
  mutate(module_gene = paste(module, gene, sep = "__")) %>%
  ungroup()
top_heatmap_genes <- unique(module_df$gene)

avg_exp <- AverageExpression(Myeloid, assays = "RNA", features = top_heatmap_genes, group.by = "Group", slot = "data")$RNA %>%
  as.matrix()
avg_z <- t(scale(t(avg_exp)))
avg_z[is.na(avg_z)] <- 0
heat_df <- as.data.frame(avg_z) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = all_of(group_levels), names_to = "Group", values_to = "z_score") %>%
  left_join(module_df, by = "gene", relationship = "many-to-many") %>%
  mutate(
    module = factor(module, levels = names(myeloid_function_modules)),
    module_gene = factor(module_gene, levels = rev(module_df$module_gene)),
    Group = factor(Group, levels = group_levels)
  )
write_csv(heat_df, file.path(getOption("hla.data_root"), "csf", "Myeloid_revision_group_DEG_heatmap_source.csv"))
write_csv(module_df, file.path(result_dir, "Myeloid_revision_functional_gene_modules.csv"))

p_heat <- ggplot(heat_df, aes(Group, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#6F99C3", mid = "#F7F7F7", high = "#B7192B", midpoint = 0, name = "Z-score") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  labs(title = "Myeloid cells", x = NULL, y = NULL) +
  theme_pub(9.2) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 7.6),
    strip.text.y = element_text(size = 7.6, angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey93", colour = NA),
    panel.grid = element_blank(),
    legend.position = "right"
  )
save_plot_safe(p_heat, "Myeloid_group_DEG_key_gene_heatmap", width = 4.7, height = 6.9)

module_score_cell <- lapply(names(myeloid_function_modules), function(module_name) {
  genes <- intersect(myeloid_function_modules[[module_name]], rownames(Myeloid))
  if (length(genes) < 2) return(NULL)
  mat <- FetchData(Myeloid, vars = genes)
  tibble(
    cell = rownames(mat),
    module = module_name,
    module_score = rowMeans(as.matrix(mat), na.rm = TRUE)
  )
}) %>%
  bind_rows() %>%
  left_join(meta_expr %>% dplyr::select(cell, Sample, Group), by = "cell") %>%
  mutate(module = factor(module, levels = names(myeloid_function_modules)))

module_score_sample <- module_score_cell %>%
  group_by(Sample, Group, module) %>%
  summarise(
    n_cells = n(),
    module_score = mean(module_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Group = factor(as.character(Group), levels = group_levels))

write_csv(module_score_cell, file.path(result_dir, "Myeloid_revision_functional_module_cell_scores.csv"))
write_csv(module_score_sample, file.path(getOption("hla.data_root"), "csf", "Myeloid_revision_functional_module_sample_scores.csv"))

module_kw <- module_score_sample %>%
  group_by(module) %>%
  summarise(
    p.value = kruskal.test(module_score ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH"), label = format_p(p.value))
write_csv(module_kw, file.path(result_dir, "Myeloid_revision_functional_module_KW_tests.csv"))

module_wilcox <- module_score_sample %>%
  group_by(module) %>%
  group_modify(function(.x, .y) {
    comparisons <- t(utils::combn(group_levels, 2))
    y_max <- max(.x$module_score, na.rm = TRUE)
    y_min <- min(.x$module_score, na.rm = TRUE)
    y_range <- max(y_max - y_min, 0.08)
    tibble(
      group1 = comparisons[, 1],
      group2 = comparisons[, 2],
      x = match(comparisons[, 1], group_levels),
      xend = match(comparisons[, 2], group_levels),
      y.position = y_max + seq_len(nrow(comparisons)) * y_range * 0.18,
      tick.height = y_range * 0.035,
      label.y = y_max + seq_len(nrow(comparisons)) * y_range * 0.18 + y_range * 0.06,
      p.value = apply(comparisons, 1, function(pair) {
        values_1 <- .x$module_score[.x$Group == pair[1]]
        values_2 <- .x$module_score[.x$Group == pair[2]]
        if (length(values_1) < 2 || length(values_2) < 2) return(NA_real_)
        stats::wilcox.test(values_1, values_2, exact = FALSE)$p.value
      })
    ) %>%
      mutate(label = format_p(p.value))
  }) %>%
  ungroup()
write_csv(module_wilcox, file.path(result_dir, "Myeloid_revision_functional_module_wilcox_tests.csv"))

module_wilcox_plot <- module_wilcox %>%
  filter(!is.na(p.value), p.value < 0.1)

p_module_score <- ggplot(module_score_sample, aes(Group, module_score, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.38, alpha = 0.82) +
  geom_point(
    position = position_jitter(width = 0.11, height = 0, seed = 20260813),
    shape = 21,
    colour = "black",
    stroke = 0.20,
    size = 1.35,
    alpha = 0.90
  ) +
  geom_segment(
    data = module_wilcox_plot,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.30
  ) +
  geom_segment(
    data = module_wilcox_plot,
    aes(x = x, xend = x, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.30
  ) +
  geom_segment(
    data = module_wilcox_plot,
    aes(x = xend, xend = xend, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.30
  ) +
  geom_text(
    data = module_wilcox_plot,
    aes(x = (x + xend) / 2, y = label.y, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 2.45
  ) +
  facet_wrap(~module, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "Module score per sample") +
  theme_pub(9.2) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7.8),
    axis.text.y = element_text(size = 8.0),
    strip.text = element_text(size = 8.2),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )
save_plot_safe(p_module_score, "Myeloid_functional_module_score_compare", width = 5.4, height = 6.6)

run_go_enrichment <- function(gene_vec, comp, direction) {
  gene_vec <- unique(gene_vec)
  if (length(gene_vec) < 5) return(NULL)
  entrez <- bitr(gene_vec, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  if (nrow(entrez) < 5) return(NULL)
  ego <- tryCatch(
    enrichGO(
      gene = entrez$ENTREZID,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.2,
      qvalueCutoff = 0.5,
      readable = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  as.data.frame(ego) %>%
    as_tibble() %>%
    mutate(comparison = comp, direction = direction, database = "GO_BP")
}

enrich_res <- list()
for (nm in names(comparisons)) {
  pair <- comparisons[[nm]]
  deg <- filter(deg_all, comparison == nm, p_val_adj < 0.05, abs(avg_log2FC) >= 0.25)
  for (grp in pair) {
    genes <- deg %>% filter(group_up == grp) %>% pull(gene)
    enrich_res[[paste(nm, grp, "GO", sep = "_")]] <- run_go_enrichment(genes, nm, paste0("Up in ", grp))
  }
}
enrich_all <- bind_rows(enrich_res)
write_csv(enrich_all, file.path(result_dir, "Myeloid_revision_group_function_enrichment.csv"))

plot_enrich <- enrich_all %>%
  filter(database == "GO_BP", p.adjust < 0.25) %>%
  group_by(comparison, direction) %>%
  slice_min(n = 5, order_by = p.adjust, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    term_label = stringr::str_to_sentence(tolower(Description)),
    panel = paste(comparison, direction, sep = "\n"),
    term_label = factor(term_label, levels = rev(unique(term_label[order(p.adjust)])))
  )

if (nrow(plot_enrich) > 0) {
  p_enrich <- ggplot(plot_enrich, aes(x = -log10(p.adjust), y = term_label, size = Count)) +
    geom_point(alpha = 0.88) +
    facet_wrap(~panel, scales = "free_y", ncol = 2) +
    labs(x = "-log10 adjusted p value", y = NULL, size = "Genes") +
    theme_pub(9.8) +
    theme(
      axis.text.y = element_text(size = 8.2),
      axis.text.x = element_text(size = 8.8),
      axis.title.x = element_text(size = 10.0, face = "bold"),
      strip.text = element_text(size = 9.0),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.22),
      panel.grid.major.x = element_blank(),
      legend.title = element_text(size = 9.2, face = "bold"),
      legend.text = element_text(size = 8.6),
      legend.position = "right"
    )
  save_plot_safe(p_enrich, "Myeloid_group_function_enrichment_dotplot", width = 7.8, height = 7.0)
}

metastases_genes <- deg_all %>%
  filter(group_up == "Metastases", p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()
write_csv(tibble(gene = metastases_genes), file.path(result_dir, "Myeloid_revision_Metastases_up_genes_for_enrichR.csv"))

enrichr_dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human", "MSigDB_Hallmark_2020")
enrichr_res <- list()
if (length(metastases_genes) >= 5) {
  enrichr_raw <- enrichr(metastases_genes, enrichr_dbs)
  enrichr_res <- lapply(names(enrichr_raw), function(db) {
    enrichr_raw[[db]] %>%
      as_tibble() %>%
      mutate(Database = db)
  })
}
enrichr_all <- bind_rows(enrichr_res) %>%
  mutate(
    Adjusted.P.value = as.numeric(Adjusted.P.value),
    Combined.Score = as.numeric(Combined.Score),
    Overlap_count = as.numeric(sub("/.*$", "", Overlap)),
    Term_clean = gsub("\\s*\\([^\\)]*\\)$", "", Term),
    Term_clean = gsub("_", " ", Term_clean)
  ) %>%
  arrange(Adjusted.P.value, desc(Combined.Score))
write_csv(enrichr_all, file.path(result_dir, "Myeloid_revision_Metastases_enrichR_function_enrichment.csv"))

plot_enrichr <- enrichr_all %>%
  filter(!is.na(Adjusted.P.value), Adjusted.P.value < 0.25) %>%
  group_by(Database) %>%
  slice_min(n = 8, order_by = Adjusted.P.value, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Database = factor(Database, levels = enrichr_dbs),
    Term_clean = stringr::str_trunc(Term_clean, width = 52),
    Term_clean = factor(Term_clean, levels = rev(unique(Term_clean[order(Database, Adjusted.P.value)])))
  )

if (nrow(plot_enrichr) > 0) {
  p_enrichr <- ggplot(plot_enrichr, aes(x = -log10(Adjusted.P.value), y = Term_clean)) +
    geom_segment(aes(x = 0, xend = -log10(Adjusted.P.value), yend = Term_clean), linewidth = 0.35, colour = "grey75") +
    geom_point(aes(size = Overlap_count, fill = Database), shape = 21, colour = "black", stroke = 0.20, alpha = 0.92) +
    facet_grid(Database ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c(
      "GO_Biological_Process_2023" = "#6F99C3",
      "KEGG_2021_Human" = "#E89B1E",
      "MSigDB_Hallmark_2020" = "#B7192B"
    ), guide = "none") +
    scale_size_continuous(range = c(1.8, 4.8), name = "Genes") +
    labs(x = "-log10 adjusted p value", y = NULL) +
    theme_pub(9.8) +
    theme(
      axis.text.y = element_text(size = 8.4),
      axis.text.x = element_text(size = 8.8),
      axis.title.x = element_text(size = 10.0, face = "bold"),
      strip.text.y = element_text(size = 9.0, angle = 0, hjust = 0),
      strip.background = element_rect(fill = "grey93", colour = NA),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.22),
      panel.grid.major.x = element_blank(),
      legend.title = element_text(size = 9.2, face = "bold"),
      legend.text = element_text(size = 8.6),
      legend.position = "right"
    )
  save_plot_safe(p_enrichr, "Myeloid_group_function_enrichment_dotplot", width = 6.4, height = 6.4)
  save_plot_safe(p_enrichr, "Myeloid_Metastases_enrichR_function_enrichment_dotplot", width = 6.4, height = 6.4)
}

message("Updated myeloid group DEG and function analysis.")

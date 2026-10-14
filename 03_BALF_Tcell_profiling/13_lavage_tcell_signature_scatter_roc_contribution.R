# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Final BALF figures from saved objects and signatures.
# Original source: 21_revision_lavage_tcell_signature_scatter_roc_contribution.R
# Figure/table scope: Fig3b; Fig3c; Fig3d; Fig3e; Fig3f; Fig3g; Fig3h; Fig3i; FigS6a; FigS6b; FigS6c; FigS6d; FigS6e; FigS6f; FigS6g; FigS6h
# Data cohort: 肺泡灌洗液. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pROC)
  library(patchwork)
})

project_dir <- Sys.getenv("LAVAGE_PROJECT_DIR", unset = NA_character_)
if (is.na(project_dir) || !nzchar(project_dir)) {
  project_dir <- file.path(getOption("hla.data_root"), "\u80ba\u6ce1\u704c\u6d17\u6db2")
}
project_dir <- gsub("\\\\", "/", project_dir)
result_dir <- file.path(project_dir, "Result")
figure_dir <- file.path(project_dir, "Figures")
signature_dir <- file.path(result_dir, "Signature_score_analysis")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(signature_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("CAP", "LC", "ILD")
group_colors <- c("CAP" = "#699ECA", "LC" = "#E42320", "ILD" = "#4DAF4A")
signature_colors <- c("CAP" = "#699ECA", "LC" = "#E42320")

save_plot_safe <- function(plot, path_no_ext, width, height) {
  tmp_pdf <- file.path(tempdir(), paste0(basename(path_no_ext), "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(basename(path_no_ext), "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
  file.copy(tmp_pdf, paste0(path_no_ext, ".pdf"), overwrite = TRUE)
  file.copy(tmp_png, paste0(path_no_ext, ".png"), overwrite = TRUE)
}

format_p_label <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 2)))
  )
}

theme_pub <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      plot.margin = margin(5, 5, 5, 5)
    )
}

sample_score_path <- file.path(signature_dir, "Sample_Score.csv")
geneset_path <- file.path(signature_dir, "Geneset.csv")
if (!file.exists(sample_score_path)) {
  stop("Missing Sample_Score.csv in Signature_score_analysis.")
}
if (!file.exists(geneset_path)) {
  stop("Missing Geneset.csv in Signature_score_analysis.")
}

# Read the saved object or source table used by this analysis.
sample_score <- read.csv(sample_score_path, stringsAsFactors = FALSE, check.names = FALSE) %>%
  as_tibble(.name_repair = "unique") %>%
  rename_with(~sub("^$", "row_id", .x)) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    HLA_Splitter = as.character(HLA_Splitter),
    CAP = as.numeric(CAP),
    LC = as.numeric(LC)
  ) %>%
  filter(!is.na(Group), !is.na(CAP), !is.na(LC))

# Export the analysis source data for figure reproduction.
write.csv(sample_score, file.path(signature_dir, "Sample_Score_revision_source.csv"), row.names = FALSE)

signature_long <- sample_score %>%
  select(Group, HLA_Splitter, CAP, LC) %>%
  pivot_longer(cols = c(CAP, LC), names_to = "Signature", values_to = "Score") %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    Signature = factor(Signature, levels = c("CAP", "LC"))
  ) %>%
  filter(!is.na(Group), !is.na(Score))

comparisons <- combn(group_levels, 2, simplify = FALSE)
signature_wilcox <- signature_long %>%
  group_by(Signature) %>%
  group_modify(function(.x, .y) {
    bind_rows(lapply(comparisons, function(comp) {
      d <- .x %>% filter(as.character(Group) %in% comp) %>%
        mutate(Group_test = factor(as.character(Group), levels = comp))
      p <- if (n_distinct(as.character(d$Group_test)) == 2 && all(table(d$Group_test) > 0)) {
        suppressWarnings(wilcox.test(Score ~ Group_test, data = d, exact = FALSE)$p.value)
      } else {
        NA_real_
      }
      tibble(group1 = comp[1], group2 = comp[2], p_value = p)
    }))
  }) %>%
  ungroup() %>%
  mutate(
    p_label = format_p_label(p_value),
    group1_x = match(group1, group_levels),
    group2_x = match(group2, group_levels)
  )
write.csv(signature_wilcox, file.path(signature_dir, "TCell_signature_score_by_group_wilcox_revision.csv"), row.names = FALSE)

signature_annot <- signature_wilcox %>%
  filter(!is.na(p_value), p_value < 0.1) %>%
  left_join(
    signature_long %>%
      group_by(Signature) %>%
      summarise(y_max = max(Score, na.rm = TRUE), y_min = min(Score, na.rm = TRUE), .groups = "drop"),
    by = "Signature"
  ) %>%
  group_by(Signature) %>%
  arrange(p_value, .by_group = TRUE) %>%
  mutate(
    y_range = pmax(y_max - y_min, 0.1),
    y = y_max + y_range * (0.13 + 0.15 * (row_number() - 1)),
    y_tick = y - y_range * 0.035
  ) %>%
  ungroup()

p_signature_box <- ggplot(signature_long, aes(Group, Score, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.34, linewidth = 0.42, colour = "black") +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 11),
    shape = 21,
    size = 2.1,
    stroke = 0.52,
    alpha = 0.92
  ) +
  {
    if (nrow(signature_annot) > 0) {
      list(
        geom_segment(
          data = signature_annot,
          aes(x = group1_x, xend = group2_x, y = y, yend = y),
          inherit.aes = FALSE,
          linewidth = 0.35,
          colour = "black"
        ),
        geom_segment(
          data = signature_annot,
          aes(x = group1_x, xend = group1_x, y = y_tick, yend = y),
          inherit.aes = FALSE,
          linewidth = 0.35,
          colour = "black"
        ),
        geom_segment(
          data = signature_annot,
          aes(x = group2_x, xend = group2_x, y = y_tick, yend = y),
          inherit.aes = FALSE,
          linewidth = 0.35,
          colour = "black"
        ),
        geom_text(
          data = signature_annot,
          aes(x = (group1_x + group2_x) / 2, y = y, label = p_label),
          inherit.aes = FALSE,
          family = "Arial",
          size = 3.0,
          vjust = -0.30,
          colour = "black"
        )
      )
    }
  } +
  facet_wrap(~Signature, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "Signature score") +
  theme_pub(10.2) +
  theme(
    strip.background = element_rect(fill = "grey94", colour = NA),
    strip.text = element_text(size = 10.8, face = "bold"),
    axis.title.y = element_text(size = 11.8, face = "bold"),
    axis.text.x = element_text(size = 10.2, colour = "black"),
    axis.text.y = element_text(size = 9.8, colour = "black"),
    panel.spacing.x = unit(0.75, "lines"),
    legend.position = "none",
    plot.margin = margin(3, 4, 3, 4)
  )

save_plot_safe(p_signature_box, file.path(signature_dir, "TCell_signature_score_by_group_boxplot_wilcox_revision"), width = 4.4, height = 3.1)
save_plot_safe(p_signature_box, file.path(figure_dir, "TCell_signature_score_by_group_boxplot_wilcox_revision"), width = 4.4, height = 3.1)

scatter_df <- sample_score %>%
  filter(Group %in% c("CAP", "LC")) %>%
  mutate(Group = factor(as.character(Group), levels = c("CAP", "LC")))
write.csv(scatter_df, file.path(signature_dir, "Sample_Score_CAP_LC_scatter_revision_source.csv"), row.names = FALSE)

p_scatter <- ggplot(scatter_df, aes(CAP, LC, colour = Group, label = HLA_Splitter)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.45) +
  geom_point(shape = 21, fill = "white", size = 3.6, alpha = 0.95, stroke = 0.85) +
  ggrepel::geom_text_repel(
    size = 3.4,
    family = "Arial",
    colour = "black",
    box.padding = unit(0.30, "lines"),
    point.padding = unit(0.30, "lines"),
    segment.color = "grey65",
    segment.size = 0.25,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = group_colors[c("CAP", "LC")], drop = FALSE) +
  coord_fixed() +
  labs(
    title = "BALF T-cell signature scores",
    x = "CAP signature score",
    y = "LC signature score",
    colour = "Group"
  ) +
  theme_pub(12.0) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.45),
    axis.line = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12.0, face = "bold"),
    legend.text = element_text(size = 11.2),
    axis.title = element_text(size = 13.5, face = "bold"),
    axis.text = element_text(size = 11.5),
    plot.title = element_text(size = 14.0)
  ) +
  guides(colour = guide_legend(override.aes = list(shape = 21, fill = "white", size = 4.0, stroke = 0.85)))

save_plot_safe(p_scatter, file.path(signature_dir, "TCell_signature_scatter"), width = 5.0, height = 4.3)
save_plot_safe(p_scatter, file.path(figure_dir, "TCell_signature_scatter"), width = 5.0, height = 4.3)

make_roc_plot <- function(df, target_group) {
  score_col <- target_group
  roc_df <- df %>%
    mutate(
      target = ifelse(as.character(Group) == target_group, 1, 0),
      score = .data[[score_col]]
    )
  roc_obj <- pROC::roc(roc_df$target, roc_df$score, percent = TRUE, quiet = TRUE, direction = "auto")
  curve_df <- tibble(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  ) %>%
    arrange(desc(specificity), sensitivity)
  auc_label <- paste0("AUC = ", formatC(as.numeric(roc_obj$auc), format = "f", digits = 1), "%")
  p <- ggplot(curve_df, aes(specificity, sensitivity)) +
    geom_abline(slope = 1, intercept = 100, linetype = "dashed", colour = "grey65", linewidth = 0.45) +
    geom_line(colour = signature_colors[[target_group]], linewidth = 0.95) +
    scale_x_reverse(limits = c(100, 0), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.01, 0.03))) +
    annotate("text", x = 67, y = 18, label = auc_label, family = "Arial", fontface = "bold", size = 3.4) +
    labs(
      title = paste0(target_group, " signature ROC"),
      x = "Specificity (%)",
      y = "Sensitivity (%)"
    ) +
    theme_pub(9.5) +
    theme(plot.title = element_text(size = 10.5))
  list(plot = p, data = curve_df %>% mutate(Group = target_group, AUC = as.numeric(roc_obj$auc)))
}

roc_score <- sample_score %>%
  filter(Group %in% c("CAP", "LC")) %>%
  mutate(Group = factor(as.character(Group), levels = c("CAP", "LC")))
write.csv(roc_score, file.path(signature_dir, "Sample_Score_CAP_LC_ROC_revision_source.csv"), row.names = FALSE)

roc_cap <- make_roc_plot(roc_score, "CAP")
roc_lc <- make_roc_plot(roc_score, "LC")
roc_combined_data <- bind_rows(roc_cap$data, roc_lc$data) %>%
  mutate(
    Group = factor(Group, levels = c("CAP", "LC")),
    label = paste0(Group, " AUC = ", formatC(AUC, format = "f", digits = 1), "%")
  )
write.csv(roc_combined_data, file.path(signature_dir, "ROC_curve_CAP_LC_revision.csv"), row.names = FALSE)

auc_labels <- roc_combined_data %>%
  distinct(Group, AUC, label) %>%
  arrange(Group) %>%
  mutate(x = 62, y = c(24, 15))
p_roc_combined <- ggplot(roc_combined_data, aes(specificity, sensitivity, colour = Group)) +
  geom_abline(slope = 1, intercept = 100, linetype = "dashed", colour = "grey68", linewidth = 0.38) +
  geom_line(linewidth = 0.95, lineend = "round") +
  geom_text(
    data = auc_labels,
    aes(x = x, y = y, label = label, colour = Group),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 3.8,
    hjust = 0
  ) +
  scale_colour_manual(values = signature_colors, drop = FALSE, name = NULL) +
  scale_x_reverse(limits = c(100, 0), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Specificity (%)", y = "Sensitivity (%)") +
  theme_pub(10.8) +
  theme(
    axis.title = element_text(size = 12.4, face = "bold"),
    axis.text = element_text(size = 10.8),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.45),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = "none"
  )

save_plot_safe(roc_cap$plot, file.path(signature_dir, "ROC_Plot_CAP"), width = 3.5, height = 3.4)
save_plot_safe(roc_lc$plot, file.path(signature_dir, "ROC_Plot_LC"), width = 3.5, height = 3.4)
save_plot_safe(roc_cap$plot, file.path(figure_dir, "ROC_Plot_CAP"), width = 3.5, height = 3.4)
save_plot_safe(roc_lc$plot, file.path(figure_dir, "ROC_Plot_LC"), width = 3.5, height = 3.4)
save_plot_safe(p_roc_combined, file.path(signature_dir, "ROC_Plot_CAP_LC_combined"), width = 3.6, height = 3.2)
save_plot_safe(p_roc_combined, file.path(figure_dir, "ROC_Plot_CAP_LC_combined"), width = 3.6, height = 3.2)

geneset_all <- read.csv(geneset_path, stringsAsFactors = FALSE, check.names = FALSE) %>%
  as_tibble(.name_repair = "unique") %>%
  rename_with(~ifelse(.x == "" | grepl("^\\.\\.\\.", .x), "row_id", .x)) %>%
  mutate(
    signature = factor(as.character(cluster), levels = c("CAP", "LC")),
    gene = ifelse(!is.na(gene) & nzchar(gene), gene, X),
    avg_log2FC = as.numeric(avg_log2FC),
    pct.1 = as.numeric(pct.1),
    pct.2 = as.numeric(pct.2),
    p_val_adj = as.numeric(p_val_adj),
    contribution = avg_log2FC * abs(pct.1 - pct.2)
  ) %>%
  filter(signature %in% c("CAP", "LC"), !is.na(gene), !is.na(contribution))

signature_gene_counts <- geneset_all %>%
  distinct(signature, gene) %>%
  count(signature, name = "n_signature_genes") %>%
  mutate(strip_label = paste0(signature, " (n=", n_signature_genes, ")"))
signature_strip_labels <- setNames(signature_gene_counts$strip_label, signature_gene_counts$signature)

geneset <- geneset_all %>%
  group_by(signature) %>%
  arrange(desc(contribution), p_val_adj, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(
    gene = factor(gene, levels = rev(unique(gene))),
    signature = factor(as.character(signature), levels = c("CAP", "LC"))
  )

write.csv(geneset, file.path(signature_dir, "Top10_signature_contribution_genes_revision.csv"), row.names = FALSE)

p_contrib <- ggplot(geneset, aes(contribution, gene, colour = signature)) +
  geom_segment(aes(x = 0, xend = contribution, y = gene, yend = gene), linewidth = 0.70, alpha = 0.68) +
  geom_point(size = 2.8, alpha = 0.95) +
  facet_wrap(~signature, scales = "free_y", ncol = 2, labeller = labeller(signature = signature_strip_labels)) +
  scale_colour_manual(values = signature_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Top contribution genes for T-cell signatures",
    x = "Contribution score",
    y = NULL,
    colour = "Signature"
  ) +
  theme_pub(9.5) +
  theme(
    strip.text = element_text(face = "bold", size = 9.2),
    strip.background = element_rect(fill = "grey94", colour = NA),
    axis.title.x = element_text(size = 10.2, face = "bold"),
    axis.text.x = element_text(size = 8.8),
    axis.text.y = element_text(size = 8.0),
    legend.position = "none",
    plot.title = element_text(size = 10.2),
    panel.spacing.x = unit(0.55, "lines"),
    plot.margin = margin(3, 3, 3, 3)
  )

save_plot_safe(p_contrib, file.path(signature_dir, "Top10_signature_contribution_genes_barplot"), width = 4.9, height = 3.5)
save_plot_safe(p_contrib, file.path(figure_dir, "Top10_signature_contribution_genes_barplot"), width = 4.9, height = 3.5)

message("Saved TCell_signature_scatter, ROC_Plot_CAP, ROC_Plot_LC, and top10 contribution gene barplot.")

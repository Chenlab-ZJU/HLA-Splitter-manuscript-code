# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final robust signature construction, scoring and cross-cohort evaluation.
# Original source: 37_quick_redraw_cross_cohort_robust_boxplots_pairwise_target_vs_each_group.R
# Figure/table scope: Fig4i; Fig4j; FigS8a; FigS8b; FigS8c; FigS8d; FigS8e
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
project_dir <- if (nzchar(env_project_dir)) gsub("\\\\", "/", env_project_dir) else file.path(getOption("hla.data_root"), "\u5c3f\u6db2\u6837\u672c")

figure_dir <- file.path(project_dir, "Figures", "Signature_score_analysis", "Cross_cohort_robust_selection")
figure_dir_root <- file.path(project_dir, "Figures", "Signature_score_analysis")
result_dir <- file.path(project_dir, "Result", "Signature_score_analysis", "Cross_cohort_robust_selection")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir_root, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
signature_levels <- c("TCMR", "ABMR", "BK")
group_colors <- c(
  "NR" = "#E42320",
  "TCMR" = "#699ECA",
  "ABMR" = "#4DAF4A",
  "Mix" = "#FF8C00",
  "BK" = "#F898CB"
)

write_csv_safe <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}

format_p_value <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    p < 0.01 ~ formatC(p, format = "f", digits = 3),
    TRUE ~ formatC(p, format = "f", digits = 2)
  )
}

theme_pub <- function(base_size = 9.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      strip.background = element_rect(fill = "grey93", colour = NA),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.margin = margin(4, 6, 4, 6)
    )
}

save_plot_safe <- function(plot, filename, width, height) {
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(filename, ".png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
  file.copy(tmp_pdf, pdf_path, overwrite = TRUE)
  file.copy(tmp_png, png_path, overwrite = TRUE)
  if (grepl("^Cross_cohort_robust_sample_score_boxplots_", filename)) {
    file.copy(tmp_pdf, file.path(figure_dir_root, paste0(filename, ".pdf")), overwrite = TRUE)
    file.copy(tmp_png, file.path(figure_dir_root, paste0(filename, ".png")), overwrite = TRUE)
  }
  invisible(c(pdf = pdf_path, png = png_path))
}

make_pairwise_target_stats <- function(df, signatures, groups) {
  bind_rows(lapply(signatures, function(sig) {
    tmp <- df %>%
      filter(as.character(signature) == sig, Group %in% groups) %>%
      mutate(Group = factor(as.character(Group), levels = groups))
    target_group <- sig
    other_groups <- setdiff(groups[groups %in% as.character(tmp$Group)], target_group)
    other_groups <- other_groups[vapply(other_groups, function(g) any(as.character(tmp$Group) == g), logical(1))]
    target_values <- tmp$mean_score[as.character(tmp$Group) == target_group]
    max_y <- max(tmp$mean_score, na.rm = TRUE)
    min_y <- min(tmp$mean_score, na.rm = TRUE)
    y_range <- max(max_y - min_y, 0.001)
    step_y <- max(0.135 * y_range, 0.00135)
    tick_height <- max(0.030 * y_range, 0.00035)
    bind_rows(lapply(seq_along(other_groups), function(i) {
      other_group <- other_groups[i]
      other_values <- tmp$mean_score[as.character(tmp$Group) == other_group]
      p <- if (length(target_values) > 0 && length(other_values) > 0 && length(unique(c(target_values, other_values))) > 1) {
        suppressWarnings(wilcox.test(target_values, other_values, exact = FALSE)$p.value)
      } else {
        NA_real_
      }
      x_target <- as.numeric(factor(target_group, levels = groups))
      x_other <- as.numeric(factor(other_group, levels = groups))
      tibble(
        dataset = unique(as.character(tmp$dataset)),
        signature = factor(sig, levels = signatures),
        target_group = target_group,
        compared_group = other_group,
        group1 = target_group,
        group2 = other_group,
        x = min(x_target, x_other),
        xend = max(x_target, x_other),
        y = max_y + step_y * i,
        y_text = max_y + step_y * i + tick_height * 0.75,
        tick_height = tick_height,
        p_value = p,
        p_label = format_p_value(p)
      )
    }))
  }))
}

plot_pairwise_boxplot <- function(df, signatures, groups, filename, width, height, y_label) {
  plot_df <- df %>%
    filter(as.character(signature) %in% signatures, Group %in% groups, is.finite(mean_score)) %>%
    mutate(
      signature = factor(as.character(signature), levels = signatures),
      Group = factor(as.character(Group), levels = groups)
    )
  stat_df <- make_pairwise_target_stats(plot_df, signatures, groups)
# Export the analysis source data for figure reproduction.
  write_csv_safe(stat_df, file.path(result_dir, paste0(filename, "_wilcox.csv")))
  stat_df_plot <- stat_df %>%
    filter(!is.na(p_value), p_value < 0.1)

  p <- ggplot(plot_df, aes(Group, mean_score, fill = Group, colour = Group)) +
    geom_boxplot(width = 0.56, outlier.shape = NA, alpha = 0.38, linewidth = 0.34) +
    geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.55, alpha = 0.88) +
    geom_segment(data = stat_df_plot, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.30, colour = "black") +
    geom_segment(data = stat_df_plot, aes(x = x, xend = x, y = y, yend = y - tick_height), inherit.aes = FALSE, linewidth = 0.30, colour = "black") +
    geom_segment(data = stat_df_plot, aes(x = xend, xend = xend, y = y, yend = y - tick_height), inherit.aes = FALSE, linewidth = 0.30, colour = "black") +
    geom_text(data = stat_df_plot, aes(x = (x + xend) / 2, y = y_text, label = p_label), inherit.aes = FALSE, size = 3.45, family = "Arial", colour = "black") +
    facet_wrap(~ signature, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = group_colors, guide = "none", drop = FALSE) +
    scale_colour_manual(values = group_colors, guide = "none", drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.035, 0.10))) +
    labs(x = NULL, y = y_label) +
    theme_pub(9.1) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 8.8),
      axis.title.y = element_text(size = 10.8),
      strip.text = element_text(size = 10.2),
      panel.spacing.x = unit(0.72, "lines")
    )
  save_plot_safe(p, filename, width, height)
  invisible(p)
}

score_file <- file.path(result_dir, "Cross_cohort_robust_sample_scores.csv")
if (!file.exists(score_file)) stop("Missing sample score table: ", score_file)

# Read the saved object or source table used by this analysis.
all_scores <- read.csv(score_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(Group = as.character(Group), signature = as.character(signature))

plot_pairwise_boxplot(
  all_scores %>% filter(dataset == "Urine"),
  c("TCMR", "ABMR", "BK"),
  c("NR", "TCMR", "ABMR", "Mix", "BK"),
  "Cross_cohort_robust_sample_score_boxplots_urine",
  7.8,
  3.7,
  "Mean AUCell score per urine sample"
)

plot_pairwise_boxplot(
  all_scores %>% filter(dataset == "GSE189536"),
  c("TCMR", "ABMR"),
  c("NR", "TCMR", "ABMR", "Mix"),
  "Cross_cohort_robust_sample_score_boxplots_GSE189536",
  5.7,
  3.7,
  "Mean AUCell score per biopsy"
)

message("Done. Redrawn target-vs-each-group robust score boxplots.")

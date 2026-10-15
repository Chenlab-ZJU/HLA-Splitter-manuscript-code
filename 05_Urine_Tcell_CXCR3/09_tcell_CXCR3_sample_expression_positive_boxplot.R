# Manuscript section: Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR
# Purpose: Final urine T-cell, CXCR3 and expanded-clone figures.
# Original source: 32_revision_tcell_CXCR3_sample_expression_positive_boxplot.R
# Figure/table scope: Fig5a; Fig5b; Fig5c; Fig5d; Fig5e; Fig5f; Fig5g; Fig5h; Fig5i; FigS9a; FigS9b; FigS9c; FigS9d; FigS9e; FigS9f; FigS9g
# Data cohort: urine. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "32_revision_tcell_CXCR3_sample_expression_positive_boxplot.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

object_file <- file.path(getOption("hla.data_root"), "urine", "Urine_Tcell.rds")
figure_dir <- file.path(project_dir, "Figures", "Tcell", "Marker_Frequency")
result_dir <- file.path(project_dir, "Result", "Tcell", "Marker_Frequency")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("NR", "TCMR", "ABMR", "Mix", "BK")
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

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "p=NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
  )
}

theme_pub <- function(base_size = 8.8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.36),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      strip.background = element_rect(fill = "grey93", colour = NA),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 6, 4, 6)
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
  invisible(c(pdf = target_pdf, png = target_png))
}

safe_wilcox <- function(v1, v2) {
  v1 <- v1[is.finite(v1)]
  v2 <- v2[is.finite(v2)]
  if (length(v1) < 2 || length(v2) < 2 || dplyr::n_distinct(c(v1, v2)) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(v1, v2, exact = FALSE)$p.value)
}

if (!file.exists(object_file)) {
  stop("Cannot find urine T-cell object: ", object_file)
}

message("Loading urine T-cell object and calculating sample-level CXCR3 metrics...")
# Read the saved object or source table used by this analysis.
tcell <- readRDS(object_file)
DefaultAssay(tcell) <- "RNA"
if (!"CXCR3" %in% rownames(tcell)) stop("CXCR3 was not found in Urine_Tcell.rds.")

cxcr3_cell <- FetchData(tcell, vars = c("CXCR3")) %>%
  rownames_to_column("cell") %>%
  as_tibble() %>%
  left_join(
    tcell@meta.data %>%
      rownames_to_column("cell") %>%
      transmute(
        cell,
        Sample = dplyr::coalesce(as.character(Sample), as.character(HLA_Splitter)),
        Group = as.character(Group),
        T_celltype = as.character(revised_celltype)
      ),
    by = "cell"
  ) %>%
  filter(!is.na(Sample), Sample != "", !is.na(Group), Group %in% group_levels) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    CXCR3_positive = CXCR3 > 0
  )

sample_summary <- cxcr3_cell %>%
  group_by(Group, Sample) %>%
  summarise(
    n_t_cells = n(),
    mean_CXCR3_expression = mean(CXCR3, na.rm = TRUE),
    median_CXCR3_expression = median(CXCR3, na.rm = TRUE),
    pct_CXCR3_positive = mean(CXCR3_positive, na.rm = TRUE) * 100,
    n_CXCR3_positive = sum(CXCR3_positive, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Group = factor(as.character(Group), levels = group_levels))

group_summary <- sample_summary %>%
  group_by(Group) %>%
  summarise(
    n_samples = n_distinct(Sample),
    total_t_cells = sum(n_t_cells),
    mean_CXCR3_expression = mean(mean_CXCR3_expression, na.rm = TRUE),
    median_CXCR3_expression = median(mean_CXCR3_expression, na.rm = TRUE),
    mean_pct_CXCR3_positive = mean(pct_CXCR3_positive, na.rm = TRUE),
    median_pct_CXCR3_positive = median(pct_CXCR3_positive, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- sample_summary %>%
  select(Group, Sample, n_t_cells, mean_CXCR3_expression, pct_CXCR3_positive) %>%
  pivot_longer(
    cols = c(mean_CXCR3_expression, pct_CXCR3_positive),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      mean_CXCR3_expression = "Mean CXCR3 expression",
      pct_CXCR3_positive = "CXCR3+ T cells (%)"
    ),
    Metric = factor(Metric, levels = c("Mean CXCR3 expression", "CXCR3+ T cells (%)"))
  )

pair_grid <- tibble::tibble(
  group1 = c("TCMR", "TCMR", "TCMR", "TCMR"),
  group2 = c("NR", "ABMR", "Mix", "BK")
)

pairwise_stats <- plot_df %>%
  group_by(Metric) %>%
  group_modify(function(dat, key) {
    pair_grid %>%
      mutate(
        n_group1 = vapply(group1, function(g) sum(dat$Group == g), numeric(1)),
        n_group2 = vapply(group2, function(g) sum(dat$Group == g), numeric(1)),
        p_value = mapply(
          function(g1, g2) safe_wilcox(dat$Value[dat$Group == g1], dat$Value[dat$Group == g2]),
          group1,
          group2
        )
      )
  }) %>%
  ungroup() %>%
  mutate(
    p_adjust = p.adjust(p_value, method = "BH"),
    label = format_p(p_value)
  )

all_pairwise_stats <- plot_df %>%
  group_by(Metric) %>%
  group_modify(function(dat, key) {
    all_pairs <- as_tibble(t(utils::combn(group_levels, 2)), .name_repair = "minimal") %>%
      setNames(c("group1", "group2"))
    all_pairs %>%
      mutate(
        p_value = mapply(
          function(g1, g2) safe_wilcox(dat$Value[dat$Group == g1], dat$Value[dat$Group == g2]),
          group1,
          group2
        )
      )
  }) %>%
  ungroup() %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH"), label = format_p(p_value))

# Export the analysis source data for figure reproduction.
write_csv_safe(cxcr3_cell, file.path(result_dir, "Tcell_CXCR3_cell_expression_revision_source.csv"))
write_csv_safe(sample_summary, file.path(result_dir, "Tcell_CXCR3_sample_expression_positive_summary_revision.csv"))
write_csv_safe(group_summary, file.path(result_dir, "Tcell_CXCR3_group_summary_revision.csv"))
write_csv_safe(pairwise_stats, file.path(result_dir, "Tcell_CXCR3_sample_level_TCMR_vs_other_wilcox_revision.csv"))
write_csv_safe(all_pairwise_stats, file.path(result_dir, "Tcell_CXCR3_sample_level_all_pairwise_wilcox_revision.csv"))

sig_stats <- pairwise_stats %>%
  filter(!is.na(p_value)) %>%
  mutate(
    group1 = factor(group1, levels = group_levels),
    group2 = factor(group2, levels = group_levels),
    x = as.numeric(group1),
    xend = as.numeric(group2)
  )

y_pos <- plot_df %>%
  group_by(Metric) %>%
  summarise(
    y_max = max(Value, na.rm = TRUE),
    y_min = min(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_step = pmax((y_max - y_min) * 0.22, ifelse(Metric == "CXCR3+ T cells (%)", 10, 0.11)))

sig_stats <- sig_stats %>%
  left_join(y_pos, by = "Metric") %>%
  group_by(Metric) %>%
  arrange(p_value, .by_group = TRUE) %>%
  mutate(
    bracket_id = row_number(),
    y = y_max + y_step * bracket_id,
    y_text = y + y_step * 0.22
  ) %>%
  ungroup()

p <- ggplot(plot_df, aes(Group, Value, fill = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.36, alpha = 0.82) +
  geom_point(
    position = position_jitter(width = 0.11, height = 0, seed = 20260826),
    shape = 21,
    colour = "black",
    stroke = 0.18,
    size = 1.45,
    alpha = 0.88
  ) +
  geom_segment(
    data = sig_stats,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = sig_stats,
    aes(x = x, xend = x, y = y - y_step * 0.08, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_segment(
    data = sig_stats,
    aes(x = xend, xend = xend, y = y - y_step * 0.08, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_text(
    data = sig_stats,
    aes(x = (x + xend) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 2.65
  ) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_colors, drop = FALSE, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.30))) +
  labs(x = NULL, y = NULL) +
  theme_pub(8.8) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7.4),
    strip.text = element_text(size = 8.6),
    legend.position = "right",
    legend.text = element_text(size = 7.6),
    legend.title = element_text(size = 8.0, face = "bold"),
    panel.spacing.x = unit(0.85, "lines")
  )

save_plot_safe(p, "Tcell_CXCR3_sample_mean_expression_positive_fraction_boxplot_revision", width = 7.0, height = 3.6)

message("Saved sample-level urine T-cell CXCR3 expression and positivity boxplots with Wilcoxon tests.")

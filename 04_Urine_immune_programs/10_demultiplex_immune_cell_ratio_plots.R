# Manuscript section: Urinary immune-cell states distinguish renal allograft rejection phenotypes
# Purpose: Final urine immune-cell composition, module and demultiplexing figures.
# Original source: 27_revision_demultiplex_immune_cell_ratio_plots.R
# Figure/table scope: Fig4b; Fig4c; Fig4d; Fig4e; Fig4f; Fig4g; Fig4h; FigS7a; FigS7b; FigS7c; FigS7d; FigS7e
# Data cohort: 尿液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

env_project_dir <- Sys.getenv("URINE_PROJECT_DIR", unset = "")
if (nzchar(env_project_dir)) {
  stopifnot(dir.exists(env_project_dir))
  project_dir <- gsub("\\\\", "/", env_project_dir)
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "Code", "27_revision_demultiplex_immune_cell_ratio_plots.R")
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
project_dir <- getOption("hla.cohort_dir")
}

figure_dir <- file.path(project_dir, "Figures", "Demultiplex")
result_dir <- file.path(project_dir, "Result", "Demultiplex")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

rejection_colors <- c("NR" = "#E42320", "RE" = "#699ECA", "BK" = "#F898CB")
rejection_levels <- c("NR", "RE", "BK")

read_csv_safe <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  x <- x[, nzchar(names(x)), drop = FALSE]
  x
}

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

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.margin = margin(5, 7, 5, 7)
    )
}

save_plot <- function(plot, filename, width, height) {
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
  invisible(c(pdf = pdf_path, png = png_path))
}

input_file <- file.path(result_dir, "Comparision_immune_cell_ratio.csv")
if (!file.exists(input_file)) {
  stop("Missing input table: ", input_file)
}

ratio_df <- read_csv_safe(input_file) %>%
  transmute(
    Sample = as.character(Sample),
    Group = as.character(Group),
    Rejection = as.character(Rejection),
    TotalCell = suppressWarnings(as.numeric(TotalCell)),
    PTD = suppressWarnings(as.numeric(PTD)),
    ImmuneCell = suppressWarnings(as.numeric(ImmuneCell)),
    Immune_cell_ratio_percent = suppressWarnings(as.numeric(CellRatio))
  ) %>%
  filter(!is.na(Sample), Sample != "/", Rejection %in% rejection_levels) %>%
  mutate(
    Rejection = factor(Rejection, levels = rejection_levels),
    detected = ImmuneCell > 0,
    CD45_status = ifelse(ImmuneCell > 0, "CD45+ cells > 0", "CD45+ cells = 0")
  )

summary_by_rejection <- ratio_df %>%
  group_by(Rejection) %>%
  summarise(
    n_samples = dplyr::n(),
    median_ratio_percent = median(Immune_cell_ratio_percent, na.rm = TRUE),
    mean_ratio_percent = mean(Immune_cell_ratio_percent, na.rm = TRUE),
    detection_rate_percent = mean(detected, na.rm = TRUE) * 100,
    median_immune_cells = median(ImmuneCell, na.rm = TRUE),
    .groups = "drop"
  )

pairwise_tests <- bind_rows(lapply(list(c("NR", "RE"), c("NR", "BK"), c("RE", "BK")), function(comp) {
  sub <- ratio_df %>% filter(Rejection %in% comp)
  tibble(
    group1 = comp[1],
    group2 = comp[2],
    p_value = if (n_distinct(sub$Rejection) == 2 && dplyr::n_distinct(sub$Immune_cell_ratio_percent) > 1) {
      suppressWarnings(wilcox.test(Immune_cell_ratio_percent ~ Rejection, data = sub, exact = FALSE)$p.value)
    } else {
      NA_real_
    }
  )
})) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = format_p(p_value)
  )

ptd_tests <- ratio_df %>%
  filter(!is.na(PTD)) %>%
  group_by(Rejection) %>%
  summarise(
    n_samples = dplyr::n(),
    spearman_rho = if (dplyr::n() >= 3 && dplyr::n_distinct(PTD) > 1 && dplyr::n_distinct(Immune_cell_ratio_percent) > 1) {
      suppressWarnings(cor(PTD, Immune_cell_ratio_percent, method = "spearman", use = "complete.obs"))
    } else {
      NA_real_
    },
    spearman_p = if (dplyr::n() >= 3 && dplyr::n_distinct(PTD) > 1 && dplyr::n_distinct(Immune_cell_ratio_percent) > 1) {
      suppressWarnings(cor.test(PTD, Immune_cell_ratio_percent, method = "spearman", exact = FALSE)$p.value)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(
      n_samples >= 3 & !is.na(spearman_rho),
      paste0(as.character(Rejection), ": rho=", sprintf("%.2f", spearman_rho), ", ", format_p(spearman_p)),
      paste0(as.character(Rejection), ": n=", n_samples, ", not tested")
    )
  )

ptd_inclusion_summary <- ratio_df %>%
  mutate(PTD_available = !is.na(PTD)) %>%
  group_by(Rejection, CD45_status, PTD_available) %>%
  summarise(n_samples = dplyr::n(), .groups = "drop")

# Export the analysis source data for figure reproduction.
write_csv_safe(summary_by_rejection, file.path(result_dir, "Comparision_immune_cell_ratio_Rejection_summary_revision.csv"))
write_csv_safe(pairwise_tests, file.path(result_dir, "Comparision_immune_cell_ratio_Rejection_wilcox_revision.csv"))
write_csv_safe(ptd_tests, file.path(result_dir, "Comparision_immune_cell_ratio_PTD_spearman_revision.csv"))
write_csv_safe(ptd_inclusion_summary, file.path(result_dir, "Comparision_immune_cell_ratio_PTD_included_zero_CD45_samples_revision.csv"))

y_max <- max(ratio_df$Immune_cell_ratio_percent, na.rm = TRUE)
pairwise_plot <- pairwise_tests %>%
  mutate(
    y.position = y_max + c(0.08, 0.16, 0.24),
    group1 = factor(group1, levels = rejection_levels),
    group2 = factor(group2, levels = rejection_levels)
  )

p_rejection <- ggplot(ratio_df, aes(Rejection, Immune_cell_ratio_percent, fill = Rejection)) +
  geom_violin(
    width = 0.74,
    alpha = 0.34,
    colour = "black",
    linewidth = 0.34,
    trim = FALSE,
    scale = "width"
  ) +
  geom_boxplot(width = 0.17, outlier.shape = NA, alpha = 0.78, linewidth = 0.32, colour = "black") +
  geom_point(
    aes(colour = Rejection),
    position = position_jitter(width = 0.10, height = 0, seed = 20260826),
    size = 1.75,
    alpha = 0.90,
    stroke = 0
  ) +
  stat_pvalue_manual(
    pairwise_plot,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.32,
    size = 2.45,
    family = "Arial"
  ) +
  scale_fill_manual(values = rejection_colors, drop = FALSE, name = "Clinical status") +
  scale_colour_manual(values = rejection_colors, drop = FALSE, name = "Clinical status") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(x = "Clinical status", y = "Immune-cell ratio (%)") +
  guides(fill = "none", colour = guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
  theme_pub(8.8) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 9.0, face = "bold"),
    axis.text = element_text(size = 8.0)
  )

ptd_df <- ratio_df %>%
  filter(!is.na(PTD)) %>%
  mutate(CD45_status = factor(CD45_status, levels = c("CD45+ cells > 0", "CD45+ cells = 0")))
label_x <- quantile(ptd_df$PTD, 0.04, na.rm = TRUE)
label_y_start <- max(ptd_df$Immune_cell_ratio_percent, na.rm = TRUE) * 0.96
label_y_step <- max(ptd_df$Immune_cell_ratio_percent, na.rm = TRUE) * 0.08
ptd_label_df <- ptd_tests %>%
  mutate(
    label_x = label_x,
    label_y = label_y_start - (as.integer(factor(Rejection, levels = rejection_levels)) - 1) * label_y_step,
    label_colour = rejection_colors[as.character(Rejection)]
  )

p_ptd <- ggplot(ptd_df, aes(PTD, Immune_cell_ratio_percent, colour = Rejection)) +
  geom_point(aes(shape = CD45_status), size = 2.0, alpha = 0.90, stroke = 0.55) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.48, alpha = 0.95) +
  geom_text(
    data = ptd_label_df,
    aes(x = label_x, y = label_y, label = label, colour = Rejection),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    family = "Arial",
    size = 2.55,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = rejection_colors, drop = FALSE, name = "Clinical status") +
  scale_shape_manual(values = c("CD45+ cells > 0" = 16, "CD45+ cells = 0" = 1), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = "PTD (post-transplant day)", y = "Immune-cell ratio (%)") +
  guides(
    colour = guide_legend(override.aes = list(size = 2.8, alpha = 1), order = 1),
    shape = guide_legend(override.aes = list(size = 2.8, alpha = 1, colour = "black"), order = 2)
  ) +
  theme_pub(8.8) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7.8),
    legend.title = element_text(size = 8.4, face = "bold"),
    axis.title = element_text(size = 9.0, face = "bold"),
    axis.text = element_text(size = 8.0)
  )

save_plot(p_rejection, "Comparision_immune_cell_ratio_Rejection", width = 3.7, height = 3.4)
save_plot(p_ptd, "Comparision_immune_cell_ratio_PTD", width = 4.8, height = 3.4)

message("Saved revised demultiplex immune-cell ratio plots.")

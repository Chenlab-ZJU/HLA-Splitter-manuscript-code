#!/usr/bin/env Rscript
# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Final benchmark / HLA-signal figures.
# Original source: 12_revision_ext_data_fig2_method_limits.R
# Figure/table scope: Fig1e; Fig1f; Fig1g; Fig1h; FigS2a; FigS2b; FigS2c
# Data cohort: Syn_data. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(readr)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
syn_dir <- getOption("hla.cohort_dir")
project_dir <- getOption("hla.project_root")
revision_dir <- if (dir.exists(file.path(project_dir, "Revison"))) {
  file.path(project_dir, "Revison")
} else {
  file.path(project_dir, "Revision")
}
out_dir <- file.path(revision_dir, "revised figure")
panel_dir <- file.path(out_dir, "panels", "Ext_Data_Fig2")
table_dir <- file.path(out_dir, "tables")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(file) readr::read_csv(file.path(syn_dir, "Result", file), show_col_types = FALSE)

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
      legend.spacing.y = unit(0.04, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

panel_tag_theme <- theme(plot.tag = element_text(face = "bold", size = 12, family = "Arial"))
my_color_palette <- c("#CE3D32FF", "#F0E685FF", "#749B58FF", "#5050FFFF", "#466983FF", "#845EC2", "#FF9E4A")
pal_tools <- c(
  "Demuxlet" = my_color_palette[1],
  "Vireo_Genotype" = my_color_palette[2],
  "Vireo-HLA" = my_color_palette[3],
  "Vireo_GeneFree" = my_color_palette[4],
  "Souporcell" = my_color_palette[5],
  "HLA-Splitter" = my_color_palette[6],
  "Demuxlet-HLA" = my_color_palette[1]
)
pal_status <- c("Correct" = "#CEDE94", "Incorrect" = "firebrick", "Undefined" = "grey80")
pal_depth <- setNames(my_color_palette[seq_len(6)], c("10", "25", "50", "100", "500", "1000"))
pal_restricted <- c(
  "Demuxlet-HLA" = "#0072B2",
  "Vireo-HLA" = "#E69F00",
  "HLA-Splitter" = "#CC79A7"
)
pale_restricted <- c(
  "Demuxlet-HLA" = "#D8EEF9",
  "Vireo-HLA" = "#F8E5A7",
  "HLA-Splitter" = "#F1D4E5"
)

set_rect_fill <- function(grob, fill_vec, idx_env = new.env(parent = emptyenv())) {
  if (is.null(idx_env$i)) idx_env$i <- 0L
  if (inherits(grob, "rect")) {
    idx_env$i <- idx_env$i + 1L
    grob$gp$fill <- fill_vec[((idx_env$i - 1L) %% length(fill_vec)) + 1L]
    grob$gp$col <- "grey35"
    grob$gp$lwd <- 0.45
  }
  if (!is.null(grob$grobs)) {
    grob$grobs <- lapply(grob$grobs, set_rect_fill, fill_vec = fill_vec, idx_env = idx_env)
  }
  if (!is.null(grob$children)) {
    grob$children <- lapply(grob$children, set_rect_fill, fill_vec = fill_vec, idx_env = idx_env)
  }
  grob
}

color_left_method_strips <- function(plot, fills) {
  gt <- ggplotGrob(plot)
  strip_idx <- grep("strip-l", gt$layout$name)
  strip_idx <- strip_idx[order(gt$layout$t[strip_idx])]
  for (i in seq_along(strip_idx)) {
    gt$grobs[[strip_idx[i]]] <- set_rect_fill(gt$grobs[[strip_idx[i]]], fills[i])
  }
  gt
}

# a. HLA-restricted Vireo/Demuxlet benchmark, paired with HLA-Splitter from the same
# downsampling benchmark for a direct operating-limit comparison.
restricted <- read_result("HLA_Restricted_benchmark.csv") %>%
  mutate(
    Cell_Number = as.numeric(sub("cells-.*", "", Run_ID)),
    Tool = recode(Tool_Name, Vireo = "Vireo-HLA", Demuxlet = "Demuxlet-HLA"),
    Accuracy_pct = Accuracy * 100
  )

splitter_accuracy <- read_result("Demultiplex_performance.csv") %>%
  group_by(Cell_Number) %>%
  summarise(
    Accuracy_pct = mean(Correct_Ratio) * 100,
    Accuracy_sd = sd(Correct_Ratio) * 100,
    .groups = "drop"
  ) %>%
  filter(Cell_Number %in% unique(restricted$Cell_Number)) %>%
  mutate(Tool = "HLA-Splitter")

restricted_plot_df <- bind_rows(
  restricted %>% select(Cell_Number, Tool, Accuracy_pct),
  splitter_accuracy %>% select(Cell_Number, Tool, Accuracy_pct)
) %>%
  mutate(
    Tool = factor(Tool, levels = c("Demuxlet-HLA", "Vireo-HLA", "HLA-Splitter")),
    Tool_Label = recode(as.character(Tool), "Demuxlet-HLA" = "Demuxlet-\nHLA", "Vireo-HLA" = "Vireo-HLA", "HLA-Splitter" = "HLA-Splitter"),
    Tool_Label = factor(Tool_Label, levels = c("Demuxlet-\nHLA", "Vireo-HLA", "HLA-Splitter")),
    Cell_Number = factor(as.character(Cell_Number), levels = rev(c("10", "25", "50", "100", "500", "1000")))
  )

p_a_base <- ggplot(restricted_plot_df, aes(x = Accuracy_pct, y = Cell_Number, fill = Tool)) +
  geom_col(width = 0.68, colour = "grey30", linewidth = 0.20, alpha = 0.88) +
  facet_grid(Tool_Label ~ ., switch = "y") +
  scale_fill_manual(values = pal_restricted) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(limits = c(0, 105), breaks = c(0, 50, 100), labels = function(x) as.character(x), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "HLA-restricted benchmark", x = "Correct assignment (%)", y = "Cells per donor") +
  theme_nc(8.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
    axis.title = element_text(face = "bold", size = 9.8),
    axis.text = element_text(size = 9.2),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 9.0),
    strip.background = element_rect(fill = "grey96", colour = "grey35", linewidth = 0.35),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.25),
    panel.spacing.y = unit(0.55, "lines"),
    axis.line.y = element_line(linewidth = 0.35, colour = "black"),
    plot.margin = margin(3, 3, 3, 3)
  )

p_a <- color_left_method_strips(p_a_base, pale_restricted[levels(restricted_plot_df$Tool)])

# b. Cell recovery versus full-genotyping tools.
full_benchmark <- read_result("Cellnum_benchmark.csv") %>%
  mutate(
    Average_cell = ifelse(is.na(Average_cell), Total_Cells / 40, Average_cell),
    Recovery_pct = Recall * 100,
    Tool_Name = factor(Tool_Name, levels = c("HLA-Splitter", "Demuxlet", "Vireo_Genotype", "Vireo_GeneFree", "Souporcell"))
  ) %>%
  filter(Tool_Name %in% levels(Tool_Name))

p_b <- ggplot(full_benchmark, aes(Average_cell, Recovery_pct, colour = Tool_Name, group = Tool_Name)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.45) +
  scale_x_log10(breaks = c(10, 25, 50, 100, 250, 500, 1000, 2500)) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 105)) +
  scale_colour_manual(values = pal_tools, na.value = "grey50") +
  labs(title = "Assigned-cell recovery", x = "Cells per donor", y = "Recovered singlets", colour = NULL) +
  theme_nc(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

# c. Cell-number downsampling, shown as mean +/- SD across three replicates.
downsampling <- read_result("Demultiplex_performance.csv") %>%
  pivot_longer(c(Correct_Ratio, Incorrect_Ratio, Undefined_Ratio), names_to = "Metric", values_to = "Ratio") %>%
  mutate(
    Metric = recode(Metric, Correct_Ratio = "Correct", Incorrect_Ratio = "Incorrect", Undefined_Ratio = "Undefined")
  ) %>%
  group_by(Cell_Number, Metric) %>%
  summarise(Mean = mean(Ratio) * 100, SD = sd(Ratio) * 100, .groups = "drop") %>%
  mutate(Metric = factor(Metric, levels = c("Correct", "Incorrect", "Undefined")))

p_c <- ggplot(downsampling, aes(Cell_Number, Mean, colour = Metric, group = Metric)) +
  geom_ribbon(aes(ymin = pmax(0, Mean - SD), ymax = pmin(100, Mean + SD), fill = Metric), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.58) +
  geom_point(size = 1.45) +
  scale_x_log10(breaks = c(10, 25, 50, 100, 250, 500, 1000, 2500)) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 100)) +
  scale_colour_manual(values = pal_status) +
  scale_fill_manual(values = pal_status) +
  labs(title = "Cell-number downsampling", x = "Cells per donor", y = "Fraction of cells", colour = NULL, fill = NULL) +
  theme_nc(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# d. Reads-per-cell rarefaction.
depth <- read_result("SeqDepth_Benchmark.csv") %>%
  mutate(
    Cells_per_donor = as.numeric(sub("cells/Sample", "", Sample)),
    Accuracy_pct = Accuracy * 100
  )

p_d <- ggplot(depth, aes(Depth_RPC, Accuracy_pct, colour = factor(Cells_per_donor), group = Cells_per_donor)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.35) +
  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000, 20000, 40000), labels = comma) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(35, 100)) +
  scale_colour_manual(values = pal_depth) +
  labs(title = "Reads-per-cell rarefaction", x = "Reads per cell", y = "Correct assignment", colour = "Cells per donor") +
  theme_nc(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

depth_summary <- depth %>%
  filter(Depth_RPC %in% c(10000, 20000)) %>%
  group_by(Depth_RPC) %>%
  summarise(Mean_accuracy_pct = mean(Accuracy_pct), SD_accuracy_pct = sd(Accuracy_pct), .groups = "drop")

# e. Practical QC decision rule derived from the downsampling and rarefaction operating limits.
qc_rule <- tibble::tibble(
  Step = factor(c("Input cells", "Pilot depth", "Proceed"), levels = c("Input cells", "Pilot depth", "Proceed")),
  Rule = c(">=100 cells per donor", ">=10,000 reads per cell", "Expected accuracy >=90%"),
  Evidence = c(
    sprintf("100-cell mean correct %.1f%%", downsampling %>% filter(Cell_Number == 100, Metric == "Correct") %>% pull(Mean)),
    sprintf("10k reads mean correct %.1f%%", depth_summary %>% filter(Depth_RPC == 10000) %>% pull(Mean_accuracy_pct)),
    "Flag samples below either threshold"
  ),
  x = c(1, 2, 3),
  y = c(1, 1, 1)
)

p_e <- ggplot(qc_rule, aes(x, y)) +
  annotate("segment", x = 1.25, xend = 1.75, y = 1, yend = 1, linewidth = 0.32, arrow = arrow(length = unit(0.10, "inches"))) +
  annotate("segment", x = 2.25, xend = 2.75, y = 1, yend = 1, linewidth = 0.32, arrow = arrow(length = unit(0.10, "inches"))) +
  geom_label(aes(label = paste0(Rule, "\n", Evidence)), size = 2.35, linewidth = 0.22, fill = "white", label.padding = unit(0.16, "lines"), family = "Arial") +
  scale_x_continuous(limits = c(0.55, 3.45)) +
  scale_y_continuous(limits = c(0.65, 1.35)) +
  labs(title = "Predictive QC rule") +
  theme_void(base_size = 8) +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0, size = 8.5),
    plot.margin = margin(4, 4, 4, 4)
  )

save_panel <- function(plot, name, width = 3.2, height = 2.45) {
  ggsave(file.path(panel_dir, paste0(name, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(panel_dir, paste0(name, ".png")), plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

safe_write_csv <- function(x, path) {
  tryCatch(
# Export the analysis source data for figure reproduction.
    write_csv(x, path),
    error = function(e) {
      alt <- file.path(dirname(path), paste0(tools::file_path_sans_ext(basename(path)), "_new.csv"))
      warning("Could not overwrite ", path, "; writing ", alt, " instead. Original error: ", conditionMessage(e))
      write_csv(x, alt)
    }
  )
}

save_panel(p_a, "Ext_Data_Fig2a_HLA_restricted_benchmark", width = 3.55, height = 4.6)
save_panel(p_b, "Ext_Data_Fig2b_cell_recovery")
save_panel(p_c, "Ext_Data_Fig2c_cell_number_downsampling")
save_panel(p_d, "Ext_Data_Fig2d_reads_per_cell_rarefaction")
save_panel(p_e, "Ext_Data_Fig2e_predictive_QC_rule", width = 6.8, height = 1.35)

combined <- (wrap_elements(p_a) | p_b) / (p_c | p_d) / p_e +
  plot_layout(heights = c(1, 1, 0.45)) +
  plot_annotation(tag_levels = "a") & panel_tag_theme

ggsave(file.path(out_dir, "Ext_Data_Fig2_method_operating_limits.pdf"), combined, width = 7.2, height = 7.1, units = "in", device = cairo_pdf)
ggsave(file.path(out_dir, "Ext_Data_Fig2_method_operating_limits.png"), combined, width = 7.2, height = 7.1, units = "in", dpi = 600, bg = "white")

safe_write_csv(restricted_plot_df, file.path(table_dir, "Ext_Data_Fig2a_HLA_restricted_benchmark_source.csv"))
safe_write_csv(full_benchmark, file.path(table_dir, "Ext_Data_Fig2b_full_tool_recovery_source.csv"))
safe_write_csv(downsampling, file.path(table_dir, "Ext_Data_Fig2c_downsampling_summary.csv"))
safe_write_csv(depth_summary, file.path(table_dir, "Ext_Data_Fig2d_depth_summary_10k_20k.csv"))
safe_write_csv(qc_rule %>% select(Step, Rule, Evidence), file.path(table_dir, "Ext_Data_Fig2e_predictive_QC_rule.csv"))

message("Saved Ext. Data Fig.2 outputs to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

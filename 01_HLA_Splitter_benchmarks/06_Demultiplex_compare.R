#!/usr/bin/env Rscript
# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Final benchmark / HLA-signal figures.
# Original source: 03_Demultiplex_compare.R
# Figure/table scope: Fig1e; Fig1f; Fig1g; Fig1h; FigS2a; FigS2b; FigS2c
# Data cohort: synthetic. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.


# Demultiplexing tool benchmark plots for synthetic pooled samples.
# Outputs publication-ready line plots for accuracy, cell recovery, and runtime.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(ggpubr)
  library(patchwork)
  library(cowplot)
  library(grid)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
project_dir <- getOption("hla.cohort_dir")
setwd(project_dir)

dir.create(file.path(project_dir, "Figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(project_dir, "Result"), showWarnings = FALSE, recursive = TRUE)

tool_levels <- c(
  "Demuxlet",
  "Vireo_Genotype",
  "Vireo_GeneFree",
  "Souporcell",
  "HLA-Splitter"
)

sample_cell_levels <- c(10, 25, 50, 100, 500, 1000)
sample_cell_labels <- as.character(sample_cell_levels)

tool_colors <- c(
  "Demuxlet" = "#0072B2",
  "Vireo_Genotype" = "#E69F00",
  "Vireo_GeneFree" = "#009E73",
  "Souporcell" = "#D55E00",
  "HLA-Splitter" = "#CC79A7"
)

pale_tool_colors <- c(
  "Demuxlet" = "#D8EEF9",
  "Vireo_Genotype" = "#F8E5A7",
  "Vireo_GeneFree" = "#CDEFE5",
  "Souporcell" = "#F5CFAB",
  "HLA-Splitter" = "#F1D4E5"
)

sample_colors <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9")

format_time <- function(seconds) {
  vapply(seconds, function(s) {
    if (is.na(s)) return(NA_character_)
    if (s < 60) return(paste0(round(s), " s"))
    if (s < 3600) return(paste0(floor(s / 60), " min ", round(s %% 60), " s"))
    paste0(floor(s / 3600), " h ", round((s %% 3600) / 60), " min")
  }, character(1))
}

nc_theme <- function(base_size = 9.5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.2, hjust = 0.5),
      axis.title = element_text(size = base_size + 0.8, color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),
      legend.key.height = unit(0.48, "cm"),
      legend.key.width = unit(0.7, "cm"),
      strip.background = element_rect(fill = "grey96", color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size + 0.3),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )
}

save_pub_plot <- function(plot, filename, width, height) {
  pdf_path <- file.path(project_dir, "Figures", paste0(filename, ".pdf"))
  png_path <- file.path(project_dir, "Figures", paste0(filename, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  read.csv(path, check.names = FALSE)
}

# Sequencing depth benchmark retained from the original script.
seq_depth_df <- read_required_csv(file.path(getOption("hla.data_root"), "synthetic", "SeqDepth_Benchmark.csv"))
seq_depth_df$Sample <- factor(
  seq_depth_df$Sample,
  levels = c(
    "10cells/Sample", "25cells/Sample", "50cells/Sample",
    "100cells/Sample", "500cells/Sample", "1000cells/Sample"
  )
)

p_seq_depth <- ggplot(seq_depth_df, aes(x = Depth_RPC, y = Accuracy, color = Sample, group = Sample)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8, stroke = 0.25) +
  geom_label_repel(
    data = seq_depth_df[seq_depth_df$Sample == "10cells/Sample", ],
    aes(label = percent(Accuracy, accuracy = 0.1)),
    size = 2.4,
    show.legend = FALSE,
    nudge_y = 0.02,
    box.padding = 0.25,
    label.size = 0.15,
    color = "black"
  ) +
  scale_x_continuous(labels = comma_format()) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(min(seq_depth_df$Accuracy) * 0.9, 1)) +
  scale_color_manual(values = sample_colors) +
  labs(title = "Effect of sequencing depth", x = "Sequencing depth (reads/cell)", y = "Demultiplexing accuracy") +
  nc_theme()

save_pub_plot(p_seq_depth, "Seqdepth_Benchmark_Plots", width = 6.2, height = 4.4)

# Cell-number benchmark across demultiplexing tools.
cell_benchmark_raw <- read_required_csv(file.path(getOption("hla.data_root"), "synthetic", "Cellnum_benchmark.csv"))

cell_benchmark <- cell_benchmark_raw %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), 0, .x))) %>%
  filter(!Tool_Name %in% c("Demuxlet_HLAgene", "Vireo_HLAgene")) %>%
  mutate(
    Tool_Name = factor(Tool_Name, levels = tool_levels),
    Average_cell = as.numeric(sub("^([0-9]+)cells.*$", "\\1", Run_ID)),
    Average_cell = ifelse(is.na(Average_cell), Total_Cells / 40, Average_cell),
    Sample_cells = factor(Average_cell, levels = sample_cell_levels, labels = sample_cell_labels),
    Wall_Time_Min = Wall_Time_Sec / 60,
    Format_time = format_time(Wall_Time_Sec),
    Peak_Memory_MB = round(Peak_Memory_MB)
  ) %>%
  filter(!is.na(Tool_Name), !is.na(Sample_cells))

runtime_outlier_table <- cell_benchmark %>%
  group_by(Tool_Name, Average_cell) %>%
  mutate(
    Runtime_median = median(Wall_Time_Sec, na.rm = TRUE),
    Runtime_mad = mad(Wall_Time_Sec, constant = 1.4826, na.rm = TRUE),
    Runtime_outlier = as.character(Tool_Name) == "Souporcell" &
      n() >= 3 &
      Runtime_mad > 0 &
      abs(Wall_Time_Sec - Runtime_median) > 3 * Runtime_mad
  ) %>%
  ungroup() %>%
  filter(Runtime_outlier) %>%
  select(Run_ID, Tool_Name, Average_cell, Wall_Time_Sec, Runtime_median, Runtime_mad)

if (nrow(runtime_outlier_table) > 0) {
# Export the analysis source data for figure reproduction.
  write.csv(
    runtime_outlier_table,
    file.path(project_dir, "Result", "Demultiplex_runtime_outliers_corrected.csv"),
    row.names = FALSE
  )
}

cell_benchmark <- cell_benchmark %>%
  group_by(Tool_Name, Average_cell) %>%
  mutate(
    Runtime_median = median(Wall_Time_Sec, na.rm = TRUE),
    Runtime_mad = mad(Wall_Time_Sec, constant = 1.4826, na.rm = TRUE),
    Runtime_outlier = as.character(Tool_Name) == "Souporcell" &
      n() >= 3 &
      Runtime_mad > 0 &
      abs(Wall_Time_Sec - Runtime_median) > 3 * Runtime_mad,
    Wall_Time_Sec = ifelse(Runtime_outlier, NA_real_, Wall_Time_Sec),
    Wall_Time_Min = Wall_Time_Sec / 60
  ) %>%
  ungroup() %>%
  select(-Runtime_median, -Runtime_mad, -Runtime_outlier)

cell_summary_observed <- cell_benchmark %>%
  group_by(Tool_Name, Average_cell, Sample_cells) %>%
  summarise(
    n_runs = n(),
    Accuracy_mean = mean(Accuracy, na.rm = TRUE),
    Accuracy_sd = sd(Accuracy, na.rm = TRUE),
    Recovery_mean = mean(Recall, na.rm = TRUE),
    Recovery_sd = sd(Recall, na.rm = TRUE),
    Runtime_sec_mean = mean(Wall_Time_Sec, na.rm = TRUE),
    Runtime_sec_sd = sd(Wall_Time_Sec, na.rm = TRUE),
    Runtime_min_mean = mean(Wall_Time_Min, na.rm = TRUE),
    Runtime_min_sd = sd(Wall_Time_Min, na.rm = TRUE),
    Assigned_singlets_mean = mean(Assigned_Singlets, na.rm = TRUE),
    Correctly_assigned_mean = mean(Correctly_Assigned, na.rm = TRUE),
    .groups = "drop"
  )

predict_rate <- function(x, y, xout) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(y) == 0) {
    return(rep(NA_real_, length(xout)))
  }
  if (length(unique(x)) == 1) {
    pred <- rep(mean(y, na.rm = TRUE), length(xout))
  } else {
    pred <- approx(log10(x), y, xout = log10(xout), rule = 2, ties = mean)$y
  }
  pmin(1, pmax(0, pred))
}

predict_runtime <- function(x, y, xout) {
  ok <- is.finite(x) & is.finite(y) & y > 0
  x <- x[ok]
  y <- y[ok]
  if (length(y) == 0) {
    return(rep(NA_real_, length(xout)))
  }
  if (length(unique(x)) >= 2) {
    fit <- lm(log10(y) ~ log10(x))
    pred <- 10^predict(fit, newdata = data.frame(x = xout))
  } else {
    pred <- rep(mean(y, na.rm = TRUE), length(xout))
  }
  pmax(0.01, pred)
}

complete_grid <- expand.grid(
  Tool_Name = factor(tool_levels, levels = tool_levels),
  Average_cell = sample_cell_levels,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Tool_Name = factor(Tool_Name, levels = tool_levels),
    Sample_cells = factor(Average_cell, levels = sample_cell_levels, labels = sample_cell_labels)
  )

cell_summary <- complete_grid %>%
  left_join(cell_summary_observed, by = c("Tool_Name", "Average_cell", "Sample_cells")) %>%
  group_by(Tool_Name) %>%
  mutate(
    Data_Type = ifelse(is.na(n_runs), "Predicted", "Observed"),
    Accuracy_mean = ifelse(
      is.na(Accuracy_mean),
      predict_rate(Average_cell[!is.na(n_runs)], Accuracy_mean[!is.na(n_runs)], Average_cell),
      Accuracy_mean
    ),
    Recovery_mean = ifelse(
      is.na(Recovery_mean),
      predict_rate(Average_cell[!is.na(n_runs)], Recovery_mean[!is.na(n_runs)], Average_cell),
      Recovery_mean
    ),
    Runtime_sec_mean = ifelse(
      is.na(Runtime_sec_mean),
      predict_runtime(Average_cell[!is.na(n_runs)], Runtime_sec_mean[!is.na(n_runs)], Average_cell),
      Runtime_sec_mean
    ),
    Runtime_min_mean = Runtime_sec_mean / 60,
    Accuracy_sd = ifelse(is.na(Accuracy_sd), 0, Accuracy_sd),
    Recovery_sd = ifelse(is.na(Recovery_sd), 0, Recovery_sd),
    Runtime_sec_sd = ifelse(is.na(Runtime_sec_sd), 0, Runtime_sec_sd),
    Runtime_min_sd = Runtime_sec_sd / 60,
    n_runs = ifelse(is.na(n_runs), 0, n_runs),
    Assigned_singlets_mean = ifelse(is.na(Assigned_singlets_mean), NA_real_, Assigned_singlets_mean),
    Correctly_assigned_mean = ifelse(is.na(Correctly_assigned_mean), NA_real_, Correctly_assigned_mean)
  ) %>%
  ungroup() %>%
  mutate(
    across(ends_with("_sd"), ~ ifelse(is.na(.x), 0, .x)),
    Data_Type = factor(Data_Type, levels = c("Observed", "Predicted")),
    Accuracy_label = percent(Accuracy_mean, accuracy = 0.1),
    Recovery_label = percent(Recovery_mean, accuracy = 0.1),
    Runtime_label = format_time(Runtime_sec_mean)
  )

write.csv(
  cell_summary,
  file.path(project_dir, "Result", "Demultiplex_tool_cell_number_summary.csv"),
  row.names = FALSE
)

plot_metric_line <- function(summary_df, y_mean, y_sd, y_label, title, filename, percent_y = TRUE, log_y = FALSE) {
  y_mean <- rlang::ensym(y_mean)
  y_sd <- rlang::ensym(y_sd)

  p <- ggplot(
    summary_df,
    aes(x = Average_cell, y = !!y_mean, color = Tool_Name, fill = Tool_Name, group = Tool_Name)
  ) +
    geom_line(linewidth = 1.05) +
    geom_point(aes(shape = Data_Type), size = 2.8, stroke = 0.65) +
    scale_x_continuous(
      trans = scales::log2_trans(),
      breaks = sample_cell_levels,
      labels = sample_cell_labels,
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    scale_color_manual(values = tool_colors, drop = FALSE) +
    scale_fill_manual(values = tool_colors, drop = FALSE) +
    scale_shape_manual(values = c("Observed" = 16, "Predicted" = 1), drop = FALSE) +
    labs(title = title, x = "Input cell number per sample", y = y_label) +
    guides(
      color = guide_legend(order = 1, override.aes = list(linewidth = 1.2, size = 3.1)),
      shape = "none"
    ) +
    nc_theme()

  if (percent_y) {
    p <- p +
      scale_y_continuous(
        labels = percent_format(accuracy = 1),
        limits = c(0, 1.02),
        expand = expansion(mult = c(0, 0.03))
      )
  } else if (log_y) {
    p <- p +
      scale_y_log10(
        breaks = c(0.5, 1, 2, 5, 10, 30, 60, 120, 240, 480),
        labels = function(x) paste0(x, " min"),
        expand = expansion(mult = c(0.02, 0.08))
      )
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)))
  }

  save_pub_plot(p, filename, width = 7.1, height = 4.6)
  p
}

p_accuracy <- plot_metric_line(
  cell_summary,
  Accuracy_mean,
  Accuracy_sd,
  "Accuracy",
  "Demultiplexing accuracy across input cell numbers",
  "Demultiplex_tool_accuracy_by_cell_number",
  percent_y = TRUE
)

p_recovery <- plot_metric_line(
  cell_summary,
  Recovery_mean,
  Recovery_sd,
  "Cell recovery rate",
  "Cell recovery across input cell numbers",
  "Demultiplex_tool_recovery_by_cell_number",
  percent_y = TRUE
)

p_runtime <- plot_metric_line(
  cell_summary,
  Runtime_min_mean,
  Runtime_min_sd,
  "Runtime",
  "Computational time across input cell numbers",
  "Demultiplex_tool_runtime_by_cell_number",
  percent_y = FALSE,
  log_y = TRUE
)

summary_bar_df <- cell_summary %>%
  select(
    Tool_Name,
    Average_cell,
    Sample_cells,
    Data_Type,
    Accuracy_mean,
    Accuracy_sd,
    Recovery_mean,
    Recovery_sd,
    Runtime_min_mean,
    Runtime_min_sd
  ) %>%
  pivot_longer(
    cols = c(Accuracy_mean, Recovery_mean, Runtime_min_mean),
    names_to = "Metric_raw",
    values_to = "Mean"
  ) %>%
  mutate(
    SD = case_when(
      Metric_raw == "Accuracy_mean" ~ Accuracy_sd,
      Metric_raw == "Recovery_mean" ~ Recovery_sd,
      Metric_raw == "Runtime_min_mean" ~ Runtime_min_sd,
      TRUE ~ NA_real_
    ),
    Metric = case_when(
      Metric_raw == "Accuracy_mean" ~ "Accuracy (%)",
      Metric_raw == "Recovery_mean" ~ "Recovery rate (%)",
      Metric_raw == "Runtime_min_mean" ~ "Runtime (min)",
      TRUE ~ Metric_raw
    ),
    Value = ifelse(Metric %in% c("Accuracy (%)", "Recovery rate (%)"), Mean * 100, Mean),
    Error = ifelse(Metric %in% c("Accuracy (%)", "Recovery rate (%)"), SD * 100, SD),
    xmin = ifelse(Metric == "Runtime (min)", pmax(0.01, Value - Error), pmax(0, Value - Error)),
    xmax = Value + Error,
    Sample_cells = factor(as.character(Sample_cells), levels = rev(sample_cell_labels)),
    Metric = factor(Metric, levels = c("Accuracy (%)", "Recovery rate (%)", "Runtime (min)"))
  )

summary_point_df <- cell_benchmark %>%
  select(Tool_Name, Average_cell, Sample_cells, Accuracy, Recall, Wall_Time_Min) %>%
  pivot_longer(
    cols = c(Accuracy, Recall, Wall_Time_Min),
    names_to = "Metric_raw",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = case_when(
      Metric_raw == "Accuracy" ~ "Accuracy (%)",
      Metric_raw == "Recall" ~ "Recovery rate (%)",
      Metric_raw == "Wall_Time_Min" ~ "Runtime (min)",
      TRUE ~ Metric_raw
    ),
    Value = ifelse(Metric %in% c("Accuracy (%)", "Recovery rate (%)"), Value * 100, Value),
    Sample_cells = factor(as.character(Sample_cells), levels = rev(sample_cell_labels)),
    Metric = factor(Metric, levels = c("Accuracy (%)", "Recovery rate (%)", "Runtime (min)"))
  ) %>%
  filter(!is.na(Sample_cells), is.finite(Value))

format_runtime_tick <- function(minutes) {
  vapply(minutes, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x < 1) return(paste0(round(x * 60), "s"))
    if (x < 60) return(paste0(round(x), "min"))
    if (x < 1440) {
      hours <- x / 60
      if (abs(hours - round(hours)) < 0.01) return(paste0(round(hours), "h"))
      return(paste0(round(hours, 1), "h"))
    }
    days <- x / 1440
    if (abs(days - round(days)) < 0.01) return(paste0(round(days), "d"))
    paste0(round(days, 1), "d")
  }, character(1))
}

set_rect_fill <- function(grob, fill) {
  if (inherits(grob, "rect")) {
    grob$gp$fill <- fill
    return(grob)
  }
  if (!is.null(grob$children)) {
    grob$children <- lapply(grob$children, set_rect_fill, fill = fill)
  }
  if (!is.null(grob$grobs)) {
    grob$grobs <- lapply(grob$grobs, set_rect_fill, fill = fill)
  }
  grob
}

color_left_method_strips <- function(plot, fills) {
  g <- ggplotGrob(plot)
  strip_ids <- grep("strip-l", g$layout$name)
  if (length(strip_ids) == 0) {
    return(g)
  }
  strip_ids <- strip_ids[order(g$layout$t[strip_ids])]
  method_levels <- names(fills)
  for (i in seq_along(strip_ids)) {
    if (i <= length(method_levels)) {
      g$grobs[[strip_ids[[i]]]] <- set_rect_fill(g$grobs[[strip_ids[[i]]]], fills[[method_levels[[i]]]])
    }
  }
  g
}

runtime_breaks_min <- c(0.5, 60, 1440, 10080)

plot_metric_bar <- function(metric_name, title, x_label, x_scale = c("percent", "runtime_log2"),
                            show_y = TRUE, show_method_strip = TRUE) {
  x_scale <- match.arg(x_scale)
  bar_df <- summary_bar_df %>% filter(Metric == metric_name)
  point_df <- summary_point_df %>% filter(Metric == metric_name)

  p <- ggplot(bar_df, aes(y = Sample_cells, x = Value, fill = Tool_Name)) +
    geom_col(
      orientation = "y",
      width = 0.66,
      color = "grey30",
      linewidth = 0.18,
      alpha = 0.84
    ) +
    geom_errorbar(
      aes(xmin = xmin, xmax = xmax),
      orientation = "y",
      width = 0.17,
      linewidth = 0.34,
      color = "black"
    ) +
    geom_point(
      data = point_df,
      aes(x = Value, y = Sample_cells),
      inherit.aes = FALSE,
      position = position_jitter(height = 0.08, width = 0),
      shape = 21,
      size = 1.20,
      stroke = 0.25,
      fill = "white",
      color = "black",
      alpha = 0.90
    ) +
    facet_grid(Tool_Name ~ ., switch = "y") +
    scale_fill_manual(values = tool_colors, drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(x = x_label, y = if (show_y) "Cells/sample" else NULL, title = title) +
    nc_theme(base_size = 11.6) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13.4),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10.6),
      axis.text.y = element_text(size = 10.2, color = "black"),
      axis.title = element_text(face = "bold", size = 12.4),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey96", color = "black", linewidth = 0.35),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 11.4),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      panel.spacing.y = unit(0.08, "lines"),
      plot.margin = margin(3, 1.5, 3, 1.5)
    )

  if (x_scale == "percent") {
    p <- p +
      scale_x_continuous(
        limits = c(0, 108),
        breaks = c(0, 50, 100),
        labels = function(x) as.character(x),
        expand = expansion(mult = c(0, 0.04))
      )
  } else {
    p <- p +
      scale_x_continuous(
        trans = scales::log2_trans(),
        breaks = runtime_breaks_min,
        labels = format_runtime_tick,
        expand = expansion(mult = c(0.03, 0.08))
      ) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10.2))
  }

  if (!show_y) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
      )
  }
  if (!show_method_strip) {
    p <- p +
      theme(
        strip.text.y.left = element_blank(),
        strip.background.y = element_blank()
      )
  }
  p
}

p_accuracy_bar <- color_left_method_strips(
  plot_metric_bar("Accuracy (%)", "Accuracy (%)", "Accuracy", "percent"),
  pale_tool_colors[tool_levels]
)
p_recovery_bar <- color_left_method_strips(
  plot_metric_bar("Recovery rate (%)", "Recovery (%)", "Recovery rate", "percent"),
  pale_tool_colors[tool_levels]
)
p_runtime_bar <- color_left_method_strips(
  plot_metric_bar("Runtime (min)", "Runtime", "Runtime (log2 time scale)", "runtime_log2"),
  pale_tool_colors[tool_levels]
)

save_pub_plot(
  p_accuracy_bar,
  "Demultiplex_tool_three_metric_accuracy_panel",
  width = 3.9,
  height = 4.6
)
save_pub_plot(
  p_recovery_bar,
  "Demultiplex_tool_three_metric_recovery_panel",
  width = 3.9,
  height = 4.6
)
save_pub_plot(
  p_runtime_bar,
  "Demultiplex_tool_three_metric_runtime_panel",
  width = 4.2,
  height = 4.6
)

p_demultiplex_summary <- cowplot::plot_grid(
  p_accuracy_bar,
  p_recovery_bar,
  p_runtime_bar,
  nrow = 1,
  align = "none",
  rel_widths = c(1, 1, 1.08)
)

save_pub_plot(
  p_demultiplex_summary,
  "Demultiplex_tool_three_metric_summary",
  width = 12.0,
  height = 4.6
)

write.csv(
  summary_bar_df,
  file.path(project_dir, "Result", "Demultiplex_tool_three_metric_summary_bar_source.csv"),
  row.names = FALSE
)

performance_summary <- cell_summary %>%
  group_by(Tool_Name) %>%
  summarise(
    Accuracy = mean(Accuracy_mean, na.rm = TRUE),
    Recovery = mean(Recovery_mean, na.rm = TRUE),
    Runtime_min = exp(mean(log(pmax(Runtime_min_mean, 0.01)), na.rm = TRUE)),
    Observed_fraction = mean(Data_Type == "Observed"),
    .groups = "drop"
  ) %>%
  mutate(
    Speed = 1 - (log10(Runtime_min) - min(log10(Runtime_min), na.rm = TRUE)) /
      (max(log10(Runtime_min), na.rm = TRUE) - min(log10(Runtime_min), na.rm = TRUE)),
    Speed = ifelse(is.nan(Speed), 1, Speed),
    Overall_score = (Accuracy + Recovery + Speed) / 3,
    Tool_Name = factor(Tool_Name, levels = Tool_Name[order(Overall_score)])
  )

write.csv(
  performance_summary,
  file.path(project_dir, "Result", "Demultiplex_tool_overall_performance_summary.csv"),
  row.names = FALSE
)

performance_long <- performance_summary %>%
  select(Tool_Name, Accuracy, Recovery, Speed) %>%
  pivot_longer(cols = c(Accuracy, Recovery, Speed), names_to = "Metric", values_to = "Score") %>%
  mutate(
    Metric = factor(Metric, levels = c("Accuracy", "Recovery", "Speed")),
    Tool_Name = factor(Tool_Name, levels = levels(performance_summary$Tool_Name))
  )

p_overall_score <- ggplot(performance_summary, aes(x = Overall_score, y = Tool_Name, color = Tool_Name)) +
  geom_segment(aes(x = 0, xend = Overall_score, yend = Tool_Name), linewidth = 1.1, alpha = 0.55) +
  geom_point(size = 4.0, stroke = 0.35) +
  scale_color_manual(values = tool_colors, drop = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Overall performance score", y = NULL, title = "Overall advantage") +
  nc_theme(base_size = 10.5) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0)
  )

p_metric_heatmap <- ggplot(performance_long, aes(x = Metric, y = Tool_Name, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.92, height = 0.9) +
  geom_text(aes(label = percent(Score, accuracy = 1)), size = 3.4, fontface = "bold", color = "black") +
  scale_fill_gradientn(
    colors = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    limits = c(0, 1),
    labels = percent_format(accuracy = 1)
  ) +
  labs(x = NULL, y = NULL, title = "Metric components", fill = "Score") +
  nc_theme(base_size = 10.5) +
  theme(
    axis.text.x = element_text(angle = 0, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    panel.border = element_blank(),
    panel.grid = element_blank()
  )

p_performance_advantage <- ggpubr::ggarrange(
  p_overall_score,
  p_metric_heatmap,
  ncol = 2,
  widths = c(1.15, 1),
  labels = c("a", "b"),
  font.label = list(size = 11, face = "bold"),
  align = "hv"
)

save_pub_plot(
  p_performance_advantage,
  "Demultiplex_tool_overall_performance_advantage",
  width = 8.8,
  height = 4.4
)

# Keep the original integrated benchmark panel for compatibility with earlier figure assembly.
integrated_long <- cell_benchmark %>%
  filter(Total_Cells == 4000) %>%
  mutate(Efficiency = ifelse(Wall_Time_Sec > 0, 1 / Wall_Time_Sec, NA_real_)) %>%
  select(Tool_Name, Accuracy, Recall, Efficiency, Wall_Time_Sec) %>%
  pivot_longer(
    cols = c(Accuracy, Recall, Efficiency),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric_Label = case_when(
      Metric == "Accuracy" ~ "Accuracy",
      Metric == "Recall" ~ "Recovery rate",
      Metric == "Efficiency" ~ "Computational efficiency"
    ),
    Metric_Label = factor(Metric_Label, levels = c("Accuracy", "Recovery rate", "Computational efficiency"))
  )

integrated_summary <- integrated_long %>%
  group_by(Tool_Name, Metric_Label) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Mean_Time = mean(Wall_Time_Sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    SD = ifelse(is.na(SD), 0, SD),
    Text_Label = case_when(
      Metric_Label %in% c("Accuracy", "Recovery rate") ~ percent(Mean, accuracy = 0.1),
      Metric_Label == "Computational efficiency" ~ format_time(Mean_Time)
    )
  )

p_integrated <- ggplot() +
  geom_col(
    data = integrated_summary,
    aes(x = Tool_Name, y = Mean, fill = Tool_Name),
    alpha = 0.75,
    color = "black",
    linewidth = 0.25,
    width = 0.72
  ) +
  geom_jitter(
    data = integrated_long,
    aes(x = Tool_Name, y = Value),
    width = 0.16,
    size = 1.1,
    alpha = 0.65,
    color = "black",
    fill = "white",
    shape = 21,
    stroke = 0.2
  ) +
  geom_errorbar(
    data = integrated_summary,
    aes(x = Tool_Name, ymin = pmax(0, Mean - SD), ymax = Mean + SD),
    width = 0.22,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_text(
    data = integrated_summary,
    aes(x = Tool_Name, y = Mean + SD, label = Text_Label),
    vjust = -0.45,
    size = 2.4,
    color = "black",
    fontface = "bold"
  ) +
  facet_wrap(~ Metric_Label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = tool_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = NULL) +
  nc_theme() +
  theme(legend.position = "none")

save_pub_plot(p_integrated, "Combined_Benchmark_Plots", width = 9.2, height = 3.4)

# HLA resolution benchmark retained from the original script.
hla_resolution <- read_required_csv(file.path(getOption("hla.data_root"), "synthetic", "HLA_Resolution_benchmark.csv"))

hla_tool_levels <- c("HLA-Splitter_1field", "HLA-Splitter_2field", "HLA-Splitter_3field")
hla_tool_colors <- c(
  "HLA-Splitter_1field" = "#B7A6D9",
  "HLA-Splitter_2field" = "#845EC2",
  "HLA-Splitter_3field" = "#4E2A84"
)

hla_long <- hla_resolution %>%
  pivot_longer(cols = c("Accuracy", "Recall"), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Tool_Name = factor(Tool_Name, levels = hla_tool_levels),
    Metric = factor(Metric, levels = c("Accuracy", "Recall"), labels = c("Accuracy", "Recovery rate"))
  ) %>%
  filter(!is.na(Tool_Name))

hla_summary <- hla_long %>%
  group_by(Tool_Name, Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(SD = ifelse(is.na(SD), 0, SD))

p_hla <- ggplot() +
  geom_col(
    data = hla_summary,
    aes(x = Tool_Name, y = Mean, fill = Tool_Name),
    alpha = 0.75,
    color = "black",
    linewidth = 0.25,
    width = 0.72
  ) +
  geom_jitter(
    data = hla_long,
    aes(x = Tool_Name, y = Value),
    width = 0.16,
    size = 1.2,
    alpha = 0.65,
    color = "black",
    fill = "white",
    shape = 21,
    stroke = 0.2
  ) +
  geom_errorbar(
    data = hla_summary,
    aes(x = Tool_Name, ymin = pmax(0, Mean - SD), ymax = Mean + SD),
    width = 0.22,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_text(
    data = hla_summary,
    aes(x = Tool_Name, y = Mean + SD, label = percent(Mean, accuracy = 0.1)),
    vjust = -0.45,
    size = 2.4,
    fontface = "bold",
    color = "black"
  ) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = hla_tool_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "Percentage") +
  nc_theme() +
  theme(legend.position = "none")

save_pub_plot(p_hla, "HLA_Resolution_Benchmark_Plots", width = 6.5, height = 3.8)

message("Generated demultiplexing benchmark plots:")
message(" - Figures/Demultiplex_tool_accuracy_by_cell_number.pdf")
message(" - Figures/Demultiplex_tool_recovery_by_cell_number.pdf")
message(" - Figures/Demultiplex_tool_runtime_by_cell_number.pdf")
message("Source data: Result/Demultiplex_tool_cell_number_summary.csv")

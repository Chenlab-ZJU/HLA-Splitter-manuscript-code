# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: plot_cell_ratio_comparison.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

#' Plot Cell Type Ratios with Statistical Comparison and Optional CSV Export
#'
#' This function calculates the ratio of each cell type per sample and
#' generates a faceted barplot with error bars and jittered points to compare 
#' these ratios across different groups. It performs statistical tests (default: t-test for small n), 
#' displays significance levels on the plot, and can optionally export the calculated
#' percentages and statistical test results to CSV files.
#'
#' @param cell_subset A Seurat object or a data.frame containing 'Sample' and 'clusters' columns.
#' @param sample_col The name of the column in `cell_subset` for sample IDs (default: "Sample").
#' @param cluster_col The name of the column in `cell_subset` for cell type clusters (default: "Celltype").
#' @param group_col The name of the column in `cell_subset` or its metadata that contains the experimental groups (default: "Group").
#' @param comparison_groups A list of character vectors specifying the pairwise comparisons to perform.
#' @param stat_method Statistical test to use: "t.test" (recommended for small n) or "wilcox.test". Default is "t.test".
#' @param facet A logical value (TRUE/FALSE) to control whether to use `facet_wrap` for plotting. Default is FALSE.
#' @param plot_colors A named character vector of colors for the groups. If NULL, default colors will be used.
#' @param write_to_csv A character string specifying a file path prefix for CSV output.
#' @return A list containing:
#'   - `plot`: A ggplot object.
#'   - `cell_percentages`: A data.frame of calculated cell type percentages per sample.
#'   - `statistical_results`: A data.frame of statistical test results.
#'
#' @importFrom ggplot2 ggplot aes stat_summary geom_jitter facet_wrap ylab scale_y_continuous theme_classic theme element_text expansion scale_color_manual scale_fill_manual facet_grid position_jitter
#' @importFrom ggpubr stat_compare_means
#' @importFrom stats aggregate t.test wilcox.test
#' @importFrom dplyr left_join
#' @importFrom scales hue_pal
#' @importFrom utils write.csv
#'
plot_cell_ratio_comparison <- function(cell_subset,
                                       sample_col = "Sample",
                                       cluster_col = "Celltype",
                                       group_col = "Group",
                                       comparison_groups,
                                       stat_method = "t.test", # 默认改用 t.test 适应小样本
                                       facet = FALSE,
                                       plot_colors = NULL,
                                       write_to_csv = NULL) {
  
  # 1. Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Package 'ggpubr' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  
  # 2. Prepare the data frame
  if (inherits(cell_subset, "Seurat")) {
    meta_data <- cell_subset@meta.data
    if (!(sample_col %in% colnames(meta_data) && cluster_col %in% colnames(meta_data) && group_col %in% colnames(meta_data))) {
      stop("Specified columns not found in Seurat object metadata.")
    }
    df <- data.frame(
      Sample = meta_data[[sample_col]],
      clusters = meta_data[[cluster_col]],
      group = meta_data[[group_col]]
    )
  } else if (is.data.frame(cell_subset)) {
    if (!(sample_col %in% colnames(cell_subset) && cluster_col %in% colnames(cell_subset) && group_col %in% colnames(cell_subset))) {
      stop("Specified columns not found in the data.frame.")
    }
    df <- data.frame(
      Sample = cell_subset[[sample_col]],
      clusters = cell_subset[[cluster_col]],
      group = cell_subset[[group_col]]
    )
  } else {
    stop("Input 'cell_subset' must be a Seurat object or a data.frame.")
  }
  
  df$group <- factor(df$group)
  
  # 3. Calculate cell type ratios per sample
  df_counts <- as.data.frame(table(df$Sample, df$clusters))
  colnames(df_counts) <- c("Sample", "clusters", "Count")
  
  df_total_counts <- aggregate(Count ~ Sample, data = df_counts, sum)
  colnames(df_total_counts) <- c("Sample", "Total_Count")
  
  df_ratio <- dplyr::left_join(df_counts, df_total_counts, by = "Sample")
  df_sample_group <- unique(df[, c("Sample", "group")])
  df_ratio <- dplyr::left_join(df_ratio, df_sample_group, by = "Sample")
  
  df_ratio$Freq <- (df_ratio$Count / df_ratio$Total_Count) * 100
  df_ratio <- df_ratio[df_ratio$Total_Count > 0 & !is.na(df_ratio$group), ]
  
  # 4. Handle plot colors
  all_groups <- unique(df_ratio$group)
  if (is.null(plot_colors)) {
    warning("plot_colors not provided. Using default colors generated by scales::hue_pal().")
    default_colors <- scales::hue_pal()(length(all_groups))
    plot_colors <- stats::setNames(default_colors, all_groups)
  }
  df_ratio$group <- factor(df_ratio$group, levels = levels(all_groups))
  
  # 5. Optional: Write cell percentages to CSV
  if (!is.null(write_to_csv)) {
    percentage_output_path <- paste0(write_to_csv, "_cell_percentages.csv")
    utils::write.csv(df_ratio, percentage_output_path, row.names = FALSE)
    message(paste("Cell type percentages written to:", percentage_output_path))
  }
  
  # 6. Perform statistical tests and store results
  statistical_results <- data.frame()
  unique_clusters <- unique(df_ratio$clusters)
  
  for (cluster_name in unique_clusters) {
    cluster_data <- df_ratio[df_ratio$clusters == cluster_name, ]
    
    for (comp in comparison_groups) {
      group1_data <- cluster_data[cluster_data$group == comp[1], "Freq"]
      group2_data <- cluster_data[cluster_data$group == comp[2], "Freq"]
      
      if (length(group1_data) > 1 && length(group2_data) > 1) { 
        tryCatch({
          # 根据用户选择的方法进行检验
          if (stat_method == "t.test") {
            test_result <- stats::t.test(group1_data, group2_data)
          } else {
            test_result <- stats::wilcox.test(group1_data, group2_data, exact = FALSE)
          }
          
          # Convert p-value to significance label
          p_label <- ggpubr::symnum.args$cutpoints
          p_symbols <- ggpubr::symnum.args$symbols
          signif_label <- ""
          for (i in seq_along(p_label)) {
            if (test_result$p.value <= p_label[i]) {
              signif_label <- p_symbols[i]
              break
            }
          }
          if (signif_label == "" && test_result$p.value > tail(p_label, 1)) {
            signif_label <- ggpubr::symnum.args$nsymbol
          }
          
          statistical_results <- rbind(statistical_results, data.frame(
            Cluster = cluster_name, Group1 = comp[1], Group2 = comp[2],
            P_value = test_result$p.value, Significance = signif_label, Method = stat_method
          ))
        }, error = function(e) {
          warning(paste("Could not perform", stat_method, "for cluster", cluster_name, ":", e$message))
        })
      }
    }
  }
  
  # Optional: Write statistical test results to CSV
  if (!is.null(write_to_csv) && nrow(statistical_results) > 0) {
    stats_output_path <- paste0(write_to_csv, "_statistical_tests.csv")
    utils::write.csv(statistical_results, stats_output_path, row.names = FALSE)
  }
  
  # 7. Create the plot (Barplot + Error bars + Jitter points)
  p <- ggplot(df_ratio, aes(x = group, y = Freq, fill = group)) +
    # 绘制柱状图 (均值)
    stat_summary(fun = mean, geom = "bar", color = "black", linewidth = 0.6, width = 0.7, alpha = 0.6) +
    # 绘制误差线 (均值的标准误 SE)
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.6, color = "black") +
    # 绘制抖动点 (Jitter)
    geom_jitter(aes(color = group), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.9, shape = 21, color = "black") +
    ylab("Percentage (%)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    ggpubr::stat_compare_means(
      comparisons = comparison_groups,
      method = stat_method, # 使用统一的检验方法
      label = "p.signif",
      step.increase = 0.1,
      hide.ns = FALSE,
      size = 4,
      tip.length = 0.02
    ) +
    scale_fill_manual(values = plot_colors) +
    scale_color_manual(values = plot_colors) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      legend.position = "none", # 柱状图X轴已有说明，通常不再需要图例
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  # 8. Apply faceting based on the 'facet' parameter
  if (isTRUE(facet)) {
    p <- p + facet_wrap(~clusters, scales = "free_y")
  } else {
    p <- p + facet_grid(~clusters, scales = "fixed")
  }
  
  return(p)
}

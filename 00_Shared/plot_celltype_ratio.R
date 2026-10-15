# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: plot_celltype_ratio.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

#' Perform Cell Type Ratio Analysis and Plotting
#'
#' This function analyzes cell type proportions from a Seurat object,
#' generates stacked bar plots, heatmaps, and comparison plots with statistical significance.
#'
#' @param seurat_obj A Seurat object containing cell type and sample information in its metadata.
#' @param sample_col The metadata column name for sample identifiers (e.g., "Sample").
#' @param cluster_col The metadata column name for cell type clusters (e.g., "revised_celltype").
#' @param group_col The metadata column name for disease/experimental groups (e.g., "Group").
#' @param comparison_groups A list of character vectors for pairwise comparisons (e.g., list(c("NR","TCMR"))).
#' @param group_colors A named character vector for group colors (e.g., c("NR" = "red", "TCMR" = "blue")).
#' @param cluster_colors A named character vector for cluster colors.
#' @param output_dir The directory to save the plots.
#'
#' @return A list containing the generated ggplot objects.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom dplyr as_tibble group_by mutate ungroup select distinct arrange summarise pull
#' @importFrom ggplot2 ggplot aes geom_col geom_tile geom_text scale_fill_manual scale_fill_gradient2 labs theme_classic theme_void element_text element_blank expansion unit
#' @importFrom patchwork plot_layout
#' @importFrom ggpubr ggboxplot stat_compare_means
#'
plot_celltype_ratio <- function(seurat_obj,
                                sample_col,
                                group_col,
                                cluster_col,
                                comparison_groups,
                                group_colors,
                                cluster_colors) {
  
  # Ensure required packages are installed
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package 'patchwork' is required.")
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Package 'ggpubr' is required.")
  
  # --- Data Preparation ---
  df <- seurat_obj@meta.data %>%
    dplyr::select(!!sym(sample_col), !!sym(cluster_col), !!sym(group_col))
  
  df_ratio <- df %>%
    dplyr::group_by(!!sym(sample_col), !!sym(cluster_col)) %>%
    dplyr::summarise(Freq = n(), .groups = "drop_last") %>%
    dplyr::mutate(Freq = Freq / sum(Freq)) %>%
    dplyr::ungroup()
  
  colnames(df_ratio) <- c("Sample", "clusters", "Freq")
  
  # --- 1. Combined Stacked Bar Plot with Group Bar ---
  # Create the group color bar
  df_group_bar <- df %>%
    dplyr::select(Sample = !!sym(sample_col), Group = !!sym(group_col)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(factor(Sample, levels = unique(df_ratio$Sample)))
  
  p_group_bar <- ggplot2::ggplot(df_group_bar, aes(x = Sample, y = 1, fill = Group)) +
    ggplot2::geom_tile(color = "black", linewidth = 0.2, alpha = 0.5) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = unit(c(0.1, 0.2, 0, 0.2), "cm")
    )
  
  # Create the main stacked bar plot
  p_main_bars <- ggplot2::ggplot(df_ratio, aes(x = Sample, y = Freq, fill = clusters)) +
    ggplot2::geom_col(position = "stack", alpha = 0.6, color = "black", linewidth = 0.2) +
    ggplot2::ylab("Frequency in CD45+ cells") +
    ggplot2::xlab(NULL) +
    ggplot2::scale_fill_manual(values = cluster_colors) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      text = element_text(size = 12),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  # Combine plots using patchwork
  p_stacked_bar <- p_group_bar / p_main_bars +
    patchwork::plot_layout(heights = c(0.1, 1))
  
  # --- 2. Cell Type Heatmap with Significance by Group ---
  df_comp <- df_ratio %>%
    dplyr::left_join(df_group_bar, by = "Sample")
  
  # Calculate mean frequencies and p-values for heatmap
  df_heatmap <- df_comp %>%
    dplyr::group_by(clusters, Group) %>%
    dplyr::summarise(MeanFreq = mean(Freq, na.rm = TRUE), .groups = "drop") %>%
    dplyr::ungroup() %>%
    tidyr::complete(clusters, Group, fill = list(MeanFreq = 0))
  
  df_p_values <- data.frame()
  for (comp in comparison_groups) {
    if (length(comp) == 2) {
      group1 <- comp[1]
      group2 <- comp[2]
      
      p_vals <- df_comp %>%
        dplyr::filter(Group %in% c(group1, group2)) %>%
        dplyr::group_by(clusters) %>%
        dplyr::summarise(
          p_value = tryCatch({
            t.test(Freq ~ Group, data = .)$p.value
          }, error = function(e) NA)
        ) %>%
        dplyr::ungroup()
      
      df_p_values <- rbind(df_p_values,
                           p_vals %>% mutate(Group = group1, Comparison_Group = group2))
      df_p_values <- rbind(df_p_values,
                           p_vals %>% mutate(Group = group2, Comparison_Group = group1))
    }
  }
  
  # Map p-values to stars
  df_p_values <- df_p_values %>%
    dplyr::mutate(
      stars = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(clusters, Group, stars)
  
  # Merge mean frequencies with stars
  df_heatmap_plot <- dplyr::left_join(df_heatmap, df_p_values, by = c("clusters", "Group"))
  
  # --- MODIFICATION START ---
  # Get the factor levels from the original data
  group_levels <- levels(df_comp$Group)
  
  # Set the order of the Group factor in the plotting data frame
  # This is the key step to ensure the heatmap x-axis is ordered correctly
  df_heatmap_plot$Group <- factor(df_heatmap_plot$Group, levels = group_levels)
  # --- MODIFICATION END ---
  
  p_heatmap <- ggplot2::ggplot(df_heatmap_plot, aes(x = Group, y = clusters, fill = MeanFreq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(aes(label = stars), color = "black", size = 5, fontface = "bold") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = mean(df_heatmap_plot$MeanFreq),
                                  name = "Mean Frequency") +
    ggplot2::labs(
      title = "Mean Frequencies of Cell Types",
      x = "Disease Group",
      y = "Cell Cluster"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  return(list(stacked_bar = p_stacked_bar, heatmap = p_heatmap))
}

# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: plot_dotplot.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

#' Create a dot plot for features expression in scRNA-seq data.
#'
#' This function featuresrates a dot plot to visualize the average expression and
#' percentage of cells expressing a features of interest across different groups
#' and cell subtypes from a Seurat object.
#'
#' @param seurat_object A Seurat object containing scRNA-seq data.
#' @param features A character string specifying the features to plot.
#' @param group_by A character string specifying the metadata column for
#'   grouping the cells (e.g., "Group"). Defaults to "Group".
#' @param cluster_by A character string specifying the metadata column for
#'   cell subtypes (e.g., "revised_celltype"). Defaults to "revised_celltype".
#' @param color_low A character string for the low end of the color gradient.
#'   Defaults to "lightblue".
#' @param color_high A character string for the high end of the color gradient.
#'   Defaults to "red".
#'
#' @return A ggplot2 object representing the dot plot.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom dplyr %>% mutate group_by across all_of summarise select
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient
#' @importFrom ggplot2 scale_size_continuous labs theme_classic theme
#' @importFrom ggplot2 element_text
#' @importFrom tidyr drop_na
#' @export
#'
plot_dotplot <- function(seurat_object,
                         features,
                         group_by = NULL,
                         cluster_by = NULL
                         ) {
  
  # 1. 数据准备和计算
  # Check if the feature(s) exist as genes in the Seurat object.
  is_gene_feature <- base::all(features %in% base::rownames(seurat_object))
  # Check if the feature(s) exist in the object's metadata.
  is_meta_feature <- base::all(features %in% base::colnames(seurat_object@meta.data))
  
  if (is_gene_feature) {
    # If the feature is a gene, extract the normalized expression data.
    expression_data <- Seurat::GetAssayData(seurat_object, assay = "RNA", layer = "data")[features, , drop = FALSE]
    base::cat("Extracting gene expression data.\n")
    
  } else if (is_meta_feature) {
    if (!base::is.numeric(seurat_object@meta.data[[features]])) {
      base::stop("The metadata feature is not a numeric value. It cannot be processed for quantitative analysis.")
    }
    expression_data <- base::as.data.frame(seurat_object@meta.data[, features, drop = FALSE])
    base::cat("Extracting numerical metadata feature(s).\n")
    
  } else {
    # Stop the function if the feature is not found in either genes or metadata.
    base::stop("Feature(s) not found in either the Seurat object's gene expression data or its metadata.")
  }
  
  # 提取元数据
  meta_data <- seurat_object@meta.data
  expression_vector <- base::as.numeric(unlist(expression_data))
  
  # 将表达数据和元数据合并到一个数据框中
  df <- meta_data %>%
    dplyr::mutate(featuresExpression = expression_vector)
  
  # 计算每个细胞亚群和分组中，该基因的平均表达量和细胞占比
  dotplot_data <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_by, cluster_by)))) %>%
    dplyr::summarise(
      AvgExpression = base::mean(featuresExpression),
      PctExpressed = base::sum(featuresExpression > 0) / dplyr::n() * 100,
      .groups = "drop"
    ) %>%
    # 移除 AvgExpression 为 NA 的行
    tidyr::drop_na(AvgExpression)
  
  # 2. 绘制 DotPlot
  p = ggplot2::ggplot(dotplot_data, ggplot2::aes(x = .data[[cluster_by]], y = .data[[group_by]])) +
    ggplot2::geom_point(ggplot2::aes(size = PctExpressed, color = AvgExpression)) +
    
    # 美化颜色
    ggplot2::scale_color_gradientn(colours = c("navy","white","firebrick3")) +
    
    # 调整点的大小范围
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    
    # 添加标题和标签
    ggplot2::labs(
      title = base::paste0(features, " Expression Dot Plot"),
      x = "Cell Subtype",
      y = "Group",
      size = "Cell percent (%)",
      color = "Average Expression"
    ) +
    
    # 调整图表主题和细节
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
  
  return(p)
}

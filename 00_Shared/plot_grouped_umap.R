# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: plot_grouped_umap.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

#' Plot faceted UMAP density with background for different groups
#'
#' This function generates UMAP density plots faceted by a specified group,
#' with all cells displayed as a gray background to show global context.
#'
#' @param scRNA Seurat object.
#' @param split_by The column name in scRNA@meta.data to be used for faceting (e.g., "Group").
#' @param bg_color The color for the background cells (default is "grey85").
#' @param pt_size Size of the background points in the UMAP plot (default is 0.1).
#' @param density_alpha Transparency of the density layer (default is 0.6).
#' @return A list of ggplot objects, one for each group.
#' @importFrom ggplot2 ggplot geom_point stat_density_2d scale_fill_viridis_c theme_void theme ggtitle coord_fixed unit element_text
#'
plot_grouped_umap <- function(scRNA, split_by = "Group", 
                                      bg_color = "grey85", pt_size = 0.1, 
                                      density_alpha = 0.6) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  
  # 检查分组列是否存在
  if (!split_by %in% colnames(scRNA@meta.data)) {
    stop(paste("Column", split_by, "not found in Seurat metadata."))
  }
  
  groups <- levels(factor(scRNA@meta.data[[split_by]]))
  
  # 提取 UMAP 坐标并与 meta.data 合并
  if (!"umap" %in% names(scRNA@reductions)) {
    stop("UMAP reduction not found in the Seurat object. Please run RunUMAP first.")
  }
  
  umap_coords <- as.data.frame(scRNA@reductions$umap@cell.embeddings)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2") # 确保列名标准
  plot_data <- cbind(umap_coords, scRNA@meta.data)
  
  plot_list <- list()
  
  for (current_group in groups) {
    # 提取当前 group 的数据用于绘制密度
    group_data <- plot_data[plot_data[[split_by]] == current_group, ]
    
    # 构建 ggplot
    p <- ggplot() +
      # 1. 背景层：将所有细胞画成灰色散点
      geom_point(data = plot_data, aes(x = UMAP_1, y = UMAP_2), 
                 color = bg_color, size = pt_size) +
      
      # 2. 密度层：针对当前 group 绘制密度多边形
      stat_density_2d(data = group_data, 
                      aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level)), 
                      geom = "polygon", 
                      alpha = density_alpha, 
                      bins = 100, 
                      n = 100, 
                      h = c(0.5, 0.5)) +# 根据你的 UMAP 坐标范围微调这个数值
      
      # 3. 颜色映射：使用连续型 (Continuous) 调色板
      scale_fill_viridis_c(option = "plasma", name = "Density") +
      
      # 4. 主题与细节调整
      ggtitle(current_group) +
      theme_void() +
      theme(
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
      ) +
      coord_fixed()
    
    plot_list[[current_group]] <- p
  }
  
  return(plot_list)
}

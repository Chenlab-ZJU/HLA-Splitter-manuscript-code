# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: plot_pseudobulk_heatmap.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

#' @title Plot Enhanced Pseudo-bulk Gene Expression Heatmap (Simplified)
#' @description This function generates a publication-ready heatmap of pseudo-bulk gene
#' expression from a Seurat object, with advanced customization for colors,
#' annotations, and specific gene labeling.
#'
#' @param seurat_obj A Seurat object containing expression data and metadata.
#' @param selected_genes Character vector. A list of all genes to be included
#'   in the heatmap rows. Only genes present in this list and in the Seurat object
#'   will be plotted.
#' @param markered_genes A named list (optional). If provided, names are the desired
#'   cluster/group names (e.g., "Brain Tumor", "Encephalitis") and values are
#'   character vectors of gene symbols to be explicitly highlighted and marked
#'   on the heatmap (with lines and labels). These genes should ideally be a subset
#'   of `selected_genes`. If `NULL` (default), all `selected_genes` will be marked,
#'   and their corresponding group for marking will be determined by which group
#    they have the highest average expression in the pseudo-bulk matrix.
#' @param group_by_col Character string. The name of the metadata column in
#'   `seurat_obj` used to group cells for pseudo-bulk aggregation (e.g., "Group").
#' @param output_dir Character string. Directory to save the heatmap PDF.
#'   Defaults to "./Figures/".
#' @param pdf_filename Character string. The name of the output PDF file.
#'   Defaults to "scRNA_all_Heatmap_typical_markers_Enhanced.pdf".
#' @param group_colors A named character vector for the fill colors of
#'   the `Group` column annotation. Names should correspond to the group names
#'   in `group_by_col`. If `NULL` (default), colors will be generated using `scales::hue_pal()`.
#' @param heatmap_main_title Character string. The main title of the heatmap.
#'   Defaults to "Pseudo-bulk Gene Expression Heatmap by Group".
#' @param heatmap_legend_title Character string. Title for the expression color legend.
#'   Defaults to "Scaled\nExpression".
#' @param cell_width_val Numeric. Desired width of each cell in the heatmap (in cm).
#'   This will be used to calculate the overall heatmap width. Defaults to 0.6 cm.
#' @param cell_height_val Numeric. Desired height of each cell in the heatmap (in cm).
#'   This will be used to calculate the overall heatmap height. Defaults to 0.3 cm.
#'
#' @return A `ComplexHeatmap` object (invisibly). The heatmap is saved as a PDF.
#'
#' @importFrom Seurat AverageExpression GetAssayData
#' @importFrom dplyr group_by slice_max ungroup mutate arrange rename filter distinct pull bind_rows setdiff
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation anno_mark draw max_text_width width height
#' @importFrom circlize colorRamp2
#' @importFrom grid unit gpar grid.grabExpr grid.draw convertWidth convertHeight
#' @importFrom grDevices pdf dev.off
#' @importFrom methods is
#' @importFrom scales hue_pal
#'
plot_pseudobulk_heatmap <- function(seurat_obj,
                                    selected_genes,
                                    markered_genes = NULL,
                                    group_by_col = "Group",
                                    group_colors = NULL,
                                    pdf_filename = NULL,
                                    cell_width_val = 2,
                                    cell_height_val = 0.3) {
  
  # --- 0. 检查所需包是否已安装 ---
  required_packages <- c("Seurat", "dplyr", "tidyr", "tibble", "ComplexHeatmap", "circlize", "grid", "grDevices", "methods", "scales")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required but not installed. Please install it."))
    }
  }
  
  # --- 1. 数据准备与验证 ---
  message("Step 1: Preparing data for heatmap...")
  
  if (!is(seurat_obj, "Seurat")) {
    stop("Input 'seurat_obj' must be a Seurat object.")
  }
  if (!group_by_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("Metadata column '", group_by_col, "' not found in Seurat object."))
  }
  if (is.null(pdf_filename)) {
    stop(paste0("Input 'pdf_filename' not found."))
  }
  
  # 计算伪批量表达矩阵并进行标准化
  avg_expr <- AverageExpression(seurat_obj, group.by = group_by_col, assays = "RNA", slot = "data")
  gene_matrix <- as.matrix(avg_expr$RNA)
  gene_matrix_scaled <- t(scale(t(gene_matrix)))
  gene_matrix_scaled[is.nan(gene_matrix_scaled)] <- 0 # 将标准化产生的NaN替换为0
  
  # 筛选出在矩阵中存在的基因
  plot_genes_final <- intersect(unique(selected_genes), rownames(gene_matrix_scaled))
  if (length(plot_genes_final) == 0) {
    stop("None of the 'selected_genes' were found in the Seurat object. Cannot plot heatmap.")
  }
  gene_matrix_filtered <- gene_matrix_scaled[plot_genes_final, , drop = FALSE]
  
  # 确定列的顺序
  valid_col_order <- if (is.factor(seurat_obj@meta.data[[group_by_col]])) {
    levels(seurat_obj@meta.data[[group_by_col]])
  } else {
    sort(unique(seurat_obj@meta.data[[group_by_col]]))
  }
  valid_col_order <- intersect(valid_col_order, colnames(gene_matrix_filtered))
  gene_matrix_ordered <- gene_matrix_filtered[, valid_col_order, drop = FALSE]
  
  # 找到每个基因表达量最高的组
  gene_order_df <- as.data.frame(gene_matrix_ordered) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(cols = -gene, names_to = "group", values_to = "expression") %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(order_by = expression, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group = factor(group, levels = valid_col_order)) %>% # 按预设顺序排列分组
    dplyr::arrange(group)
  
  ordered_genes_for_heatmap <- gene_order_df$gene
  row_split_factor_labels <- gene_order_df$group
  
  if (length(ordered_genes_for_heatmap) == 0) {
    stop("No genes found for heatmap plotting after processing.")
  }
  
  # 最终的热图矩阵和行分割因子
  gene_matrix_final <- gene_matrix_ordered[ordered_genes_for_heatmap, , drop = FALSE]
  row_split_factor <- factor(row_split_factor_labels, levels = unique(row_split_factor_labels))
  
  # 准备行注释的标签
  labels_to_mark <- if (is.null(markered_genes)) {
    ordered_genes_for_heatmap
  } else {
    intersect(unlist(markered_genes, use.names = FALSE), ordered_genes_for_heatmap)
  }
  
  # --- 2. 定义热图颜色和注释 ---
  message("Step 2: Defining heatmap colors and annotations...")
  
  abs_max <- max(abs(range(gene_matrix_final, na.rm = TRUE)))
  col_fun <- circlize::colorRamp2(c(-abs_max, 0, abs_max), c("#546de5", "#d2dae2", "#ff4757"))
  
  col_anno_df <- data.frame(Group = colnames(gene_matrix_final))
  rownames(col_anno_df) <- colnames(gene_matrix_final)
  
  if (is.null(group_colors)) {
    group_colors_final <- scales::hue_pal()(length(valid_col_order))
    names(group_colors_final) <- valid_col_order
  } else {
    group_colors_final <- group_colors[1:length(valid_col_order)]
    names(group_colors_final) = valid_col_order
  }
  
  col_ha <- ComplexHeatmap::HeatmapAnnotation(
    Group = col_anno_df$Group,
    col = list(Group = group_colors_final),
    gp = grid::gpar(col = "black"),
    annotation_name_side = "left"
  )
  
  row_ha_right <- ComplexHeatmap::rowAnnotation(
    gene_label = ComplexHeatmap::anno_mark(
      at = match(labels_to_mark, rownames(gene_matrix_final)),
      labels = labels_to_mark,
      labels_gp = grid::gpar(fontsize = 12, fontface = "italic"),
      link_width = grid::unit(4, "mm")
    )
  )
  
  # --- 3. 绘制热图 ---
  message("Step 3: Drawing heatmap...")
  ht <- ComplexHeatmap::Heatmap(
    gene_matrix_final,
    name = "Scaled\nExpression",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 12),
    column_names_gp = grid::gpar(fontsize = 12),
    width = grid::unit(ncol(gene_matrix_final) * cell_width_val, "cm"),
    height = grid::unit(nrow(gene_matrix_final) * cell_height_val, "cm"),
    row_split = row_split_factor,
    row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    row_title_rot = 90,
    column_title = "Pseudo-bulk Gene Expression Heatmap by Group",
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    column_names_rot = 30,
    heatmap_legend_param = list(title = "Scaled\nExpression"),
    top_annotation = col_ha,
    right_annotation = row_ha_right
  )
  
  # --- 4. 保存图表 ---
  message(paste0("Step 4: Saving heatmap to ", file.path(pdf_filename), "..."))
  
  ht_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right"))
  pdf_width_in <- ncol(gene_matrix_final)*cell_width_val+2
  pdf_height_in <- nrow(gene_matrix_final)*cell_height_val/2+3
  
  grDevices::pdf(file.path(pdf_filename), width = pdf_width_in, height = pdf_height_in )
  grid::grid.draw(ht_grob)
  grDevices::dev.off()
  
  message("Heatmap generation complete.")
  invisible(ht)
}

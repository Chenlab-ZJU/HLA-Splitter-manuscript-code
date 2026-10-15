# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: Enrichment_analysis.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

#' Perform Gene Enrichment Analysis for a Two-Group Comparison
#'
#' This function takes a data.frame of differentially expressed genes from a two-group
#' comparison (e.g., Seurat FindMarkers), performs enrichment analysis for the up-regulated
#' and down-regulated genes separately, and generates a faceted bar plot.
#'
#' @param diffexp A data.frame of differential expression results containing 'p_val_adj', 'avg_log2FC'.
#' @param pos_group_name Name of the group with positive log2FC (Up-regulated).
#' @param neg_group_name Name of the group with negative log2FC (Down-regulated).
#' @param dbs A character vector of database names for enrichR.
#' @param top_n_pathways The number of top enriched pathways to display per group (default: 6).
#' @param plot_colors Colors for the positive and negative groups.
#' @param p_val_adj_threshold The adjusted p-value threshold (default: 0.05).
#' @param log2fc_threshold The log2 fold change threshold (default: 1).
#' @param output_csv_path File path to save the full enrichment results.
#' @return A ggplot object of the enrichment results.
#' @importFrom enrichR enrichr
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual coord_flip labs theme_minimal theme element_blank element_text element_line margin facet_wrap
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_replace_all

Enrichment_analysis<- function(diffexp,
                               pos_group_name = "Group A",
                               neg_group_name = "Group B",
                               dbs = c("GO_Biological_Process_2023"),
                               top_n_pathways = 6, 
                               plot_colors = c("#e74c3c", "#3498db"),
                               p_val_adj_threshold = 0.05,
                               log2fc_threshold = 1,
                               output_csv_path = "./TwoGroup_Enrichment_Results.csv") {
  
  if (!requireNamespace("enrichR", quietly = TRUE)) stop("Package 'enrichR' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  
  # 1. 确保 gene 列存在 (修复 FindMarkers 默认把基因名放在行名的问题)
  if (!"gene" %in% colnames(diffexp)) {
    diffexp$gene <- rownames(diffexp)
  }
  
  # 2. 提取上调和下调的基因列表
  genes_up <- diffexp$gene[diffexp$p_val_adj < p_val_adj_threshold & diffexp$avg_log2FC >= log2fc_threshold]
  genes_down <- diffexp$gene[diffexp$p_val_adj < p_val_adj_threshold & diffexp$avg_log2FC <= -log2fc_threshold]
  
  # 3. 定义一个内部辅助函数来运行 enrichR 并整理格式
  run_enrichment <- function(gene_list, group_label) {
    if (length(gene_list) == 0) return(NULL)
    
    res <- enrichr(gene_list, dbs)
    
    # 合并多数据库结果
    df_list <- lapply(names(res), function(db) {
      if (nrow(res[[db]]) > 0) {
        temp <- res[[db]]
        temp$Database <- db
        return(temp)
      }
      return(NULL)
    })
    
    enriched_df <- do.call(rbind, df_list[!sapply(df_list, is.null)])
    if (is.null(enriched_df) || nrow(enriched_df) == 0) return(NULL)
    
    # 提取需要的列并计算 -log10(P_value)
    df_processed <- data.frame(
      Term = enriched_df$Term,
      P_value = enriched_df$P.value,
      Genes = enriched_df$Genes,
      Group = group_label,
      Database = enriched_df$Database
    )
    
    df_processed <- df_processed[order(df_processed$P_value, decreasing = FALSE), ]
    df_processed$minus_log_pvalue <- -log10(df_processed$P_value)
    return(df_processed)
  }
  
  # 4. 执行富集分析
  message("Analyzing enriched pathways for: ", pos_group_name)
  res_up <- run_enrichment(genes_up, pos_group_name)
  
  message("Analyzing enriched pathways for: ", neg_group_name)
  res_down <- run_enrichment(genes_down, neg_group_name)
  
  # 合并保存完整结果
  combined_all_results <- bind_rows(res_up, res_down)
  if (nrow(combined_all_results) == 0) {
    warning("No enriched pathways found for either group.")
    return(NULL)
  }
# Export the analysis source data for figure reproduction.
  write.csv(combined_all_results, output_csv_path, row.names = FALSE)
  message("Full enrichment results saved to: ", output_csv_path)
  
  # 5. 提取每组 Top N 绘制图表
  top_results_list <- list()
  if (!is.null(res_up)) top_results_list[[pos_group_name]] <- head(res_up, top_n_pathways)
  if (!is.null(res_down)) top_results_list[[neg_group_name]] <- head(res_down, top_n_pathways)
  
  combined_top_results <- bind_rows(top_results_list)
  
  # 设置因子水平以控制绘图顺序
  combined_top_results$Group <- factor(combined_top_results$Group, levels = c(pos_group_name, neg_group_name))
  combined_top_results$Term <- factor(combined_top_results$Term, levels = rev(unique(combined_top_results$Term)))
  
  # 分配颜色
  names(plot_colors) <- c(pos_group_name, neg_group_name)
  
  # 6. 绘制分面柱状图
  p <- ggplot(combined_top_results, aes(x = Term, y = minus_log_pvalue, fill = Group)) +
    geom_col(linewidth = 0.5, alpha = 0.7, color = "black") +
    scale_fill_manual(values = plot_colors) +
    coord_flip() +
    labs(
      x = "Enriched Pathway",
      y = expression(-Log[10](P~value)),
      fill = "Enriched in"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_blank(), # 隐藏原生Y轴文字，因为后面会加在柱子里
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_blank(), 
      strip.text = element_text(size = 14, face = "bold"), # 分面标题
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    # 将通路名称内嵌到柱状图左侧 (移除 GO 编号使名称更整洁)
    geom_text(aes(label = sub("\\s*\\(.*", "", Term), y = 0.2),
              hjust = 0, 
              color = "black", 
              size = 4.5, 
              fontface = "bold") + 
    # 使用 free_y 分面，使两组不互相挤占空间
    facet_wrap(~Group, scales = "free_y", ncol = 1, strip.position = "top")
  
  return(p)
}

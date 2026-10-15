# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: plot_volcano.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

#' Plot a highly optimized volcano plot for two-group comparison
#'
#' @param diffexp A data.frame of differential expression results containing: 'p_val_adj', 'avg_log2FC', and 'gene'.
#' @param pos_group_name Character, name of the group with positive log2FC (e.g., "Treatment").
#' @param neg_group_name Character, name of the group with negative log2FC (e.g., "Control").
#' @param p_val_adj_threshold The adjusted p-value threshold (default: 0.05).
#' @param log2fc_threshold The log2 fold change threshold (default: 1).
#' @param colors A named character vector for coloring points.
#' @param top_n_labels Number of top genes to label per direction (default: 10, means 10 Up and 10 Down).
#' @param title Title of the plot.
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_size_manual scale_alpha_manual geom_vline geom_hline labs theme element_text element_blank element_rect ggtitle
#' @importFrom ggrepel geom_text_repel

plot_volcano <- function(diffexp,
                        pos_group_name = "Group A",
                        neg_group_name = "Group B",
                        p_val_adj_threshold = 0.05,
                        log2fc_threshold = 1,
                        colors = c('Down' = "#3498db", 'Un-sigif' = "grey80", 'Up' = "#e74c3c"),
                        top_n_labels = 10,
                        title = "Differential Expression") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' is required.")
  
  # 1. 处理 P 值边界问题 (将 p_val_adj = 0 的值替换为 R 中可表示的最小正数)
  diffexp$p_val_adj_capped <- ifelse(diffexp$p_val_adj == 0, .Machine$double.xmin, diffexp$p_val_adj)
  diffexp$logP <- -log10(diffexp$p_val_adj_capped)
  
  # 2. 标记上下调状态
  diffexp$change <- "Un-sigif"
  diffexp$change[diffexp$p_val_adj < p_val_adj_threshold & diffexp$avg_log2FC >= log2fc_threshold] <- "Up"
  diffexp$change[diffexp$p_val_adj < p_val_adj_threshold & diffexp$avg_log2FC <= -log2fc_threshold] <- "Down"
  
  # 设定因子水平，确保画图时 Un-sigif 在底层，显著点在上层
  diffexp$change <- factor(diffexp$change, levels = c("Un-sigif", "Down", "Up"))
  
  # 3. 筛选差异性最强的 Top N 基因用于标注
  # 排序逻辑：综合基础表达量、差异倍数 (|log2FC|) 和 显著性 (-log10P)
  if(all(c("pct.1", "pct.2") %in% colnames(diffexp))) {
    diffexp$base_exp <- (diffexp$pct.1 + diffexp$pct.2) / 2
  } else {
    diffexp$base_exp <- 1 
  }
  
  # 计算标准化综合打分 (Composite Score)
  diffexp$composite_score <- 0
  sig_idx <- diffexp$change %in% c("Up", "Down")
  
  if(sum(sig_idx) > 0) {
    diffexp$composite_score[sig_idx] <- 
      scale(-log10(diffexp$p_val_adj_capped[sig_idx])) +  # 显著性得分
      scale(abs(diffexp$avg_log2FC[sig_idx])) +           # 差异幅度得分
      scale(diffexp$base_exp[sig_idx])                    # 基础表达量得分
  }
  
  # 分别提取上调和下调，并按综合打分降序排列
  top_up_genes <- diffexp[diffexp$change == "Up", ]
  top_up_genes <- top_up_genes[order(top_up_genes$composite_score, decreasing = TRUE), ]
  
  top_down_genes <- diffexp[diffexp$change == "Down", ]
  top_down_genes <- top_down_genes[order(top_down_genes$composite_score, decreasing = TRUE), ]
  
  # 提取 Top N 基因进行标注
  label_genes_df <- rbind(head(top_up_genes, top_n_labels), head(top_down_genes, top_n_labels))
  diffexp$label <- ifelse(diffexp$gene %in% label_genes_df$gene, as.character(diffexp$gene), "")
  
  # 4. 构建图例标签
  legend_labels <- c("Un-sigif" = "Not Significant",
                     "Down" = paste0("Enriched in ", neg_group_name),
                     "Up" = paste0("Enriched in ", pos_group_name))
  
  # 5. 绘图
  p <- ggplot(diffexp, aes(x = avg_log2FC, y = logP, color = change)) +
    # 根据状态动态调整点的大小和透明度，突出显著基因
    geom_point(aes(size = change, alpha = change)) +
    scale_size_manual(values = c("Un-sigif" = 0.5, "Down" = 1.5, "Up" = 1.5)) +
    scale_alpha_manual(values = c("Un-sigif" = 0.5, "Down" = 0.9, "Up" = 0.9)) +
    scale_color_manual(values = colors, labels = legend_labels) +
    
    # 添加阈值虚线
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black", linewidth = 0.6, alpha = 0.7) +
    geom_hline(yintercept = -log10(p_val_adj_threshold), linetype = "dashed", color = "black", linewidth = 0.6, alpha = 0.7) +
    
    # 轴标签
    labs(x = expression(Log[2](Fold~Change)),
         y = expression(-Log[10](Adjusted~P~value)),
         title = title,
         subtitle = paste0("Thresholds: |log2FC| >= ", log2fc_threshold, " & adj.P < ", p_val_adj_threshold)) +
    
    # 优化基因标注 (斜体、防止重叠、带小圆点指示)
    geom_text_repel(aes(label = label),
                    size = 4,
                    fontface = "italic",
                    box.padding = unit(0.6, "lines"),
                    point.padding = unit(0.4, "lines"),
                    max.overlaps = 50,
                    segment.color = 'grey50',
                    segment.size = 0.4,
                    min.segment.length = 0,
                    show.legend = FALSE) +
    
    # 顶部和底部留白
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    
    # 主题美化
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.line = element_blank() # 被 panel.border 替代
    )
  
  return(p)
}

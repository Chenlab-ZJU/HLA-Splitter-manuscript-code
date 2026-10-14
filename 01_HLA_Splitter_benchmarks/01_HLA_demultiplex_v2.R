# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Shared algorithm / evaluation function.
# Original source: HLA_demultiplex_v2.R
# Figure/table scope: Fig1c; Fig1d; FigS3
# Data cohort: Syn_data. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Normalize cell-by-allele counts by the total HLA library size and log-transform.
# Scale positive allele expression without centering; retain zero observations.
# Combine locus-specific scores and refine assignments by neighborhood voting.
# v2 is retained because the cohort workflows explicitly source this implementation.

HLA_demultiplex = function(HLA_list, HLA_matrix, graph_corr = TRUE, outdir) {
  namelist = names(HLA_list)
  library(ggsci) # 确保加载调色板包
  color = pal_igv("default")(length(namelist) + 1)
  
  # Step 1: normalize each cell by its total HLA counts and log-transform.
  library_sizes <- rowSums(HLA_matrix)
  HLA_matrix_norm <- HLA_matrix / library_sizes * 1e2
  HLA_matrix_norm <- log(HLA_matrix_norm + 1)
  HLA_matrix_norm[is.na(HLA_matrix_norm)] = 0
  
  # Step 2: scale positive allele counts without centering.
  for (i in 1:ncol(HLA_matrix_norm)) {
    exp = HLA_matrix_norm[, i]
    names(exp) = rownames(HLA_matrix_norm)
    exp_positive = exp[exp > 0]
    if(length(exp_positive) > 0) {
      exp_positive = scale(exp_positive, center = F, scale = T)[, 1]
      HLA_matrix_norm[names(exp_positive), i] = exp_positive
    }
  }
  
  # Step 3: calculate donor-specific scores across the six HLA loci.
  HLA_type = c("A", "B", "C", "DRB1", "DQB1", "DPB1")
  score = data.frame(matrix(nrow = nrow(HLA_matrix_norm), ncol = 0))
  rownames(score) = rownames(HLA_matrix_norm)
  
  for (i in 1:length(namelist)) {
    score_tmp = rep(0, nrow(HLA_matrix_norm))
    for (j in 1:length(HLA_type)) {
      # Apply the original homozygous-locus correction factor.
      HLA_tmp = grep(paste0("^", HLA_type[j]), HLA_list[[i]], value = T, fixed = F)
      if (length(HLA_tmp) > 0) {
        homo_factor = ifelse(length(HLA_tmp) == 1, 1.8, 1)
        for (k in 1:length(HLA_tmp)) {
          if (HLA_tmp[k] %in% colnames(HLA_matrix_norm)) {
            tmp = unlist(HLA_matrix_norm[, HLA_tmp[k]])
            if (sum(tmp) != 0) { 
              score_tmp = score_tmp + tmp * homo_factor 
            }
          }
        }
      }
    }
    score_tmp = as.data.frame(score_tmp)
    colnames(score_tmp) = paste0(namelist[i], ".score")
    score = cbind(score, score_tmp)
  }
  
  # --- 4. 注释单细胞样本来源 ---
  sample_predict = apply(score, 1, function(x) {
    if (sum(x) == 0) {
      return("Undefined")
    } else {
      return(namelist[which.max(x)])
    }
  })
  names(sample_predict) = rownames(HLA_matrix_norm)
  table(sample_predict)
  
  # --- 5. HLA TSNE聚类 ---
  library(Rtsne)
  score_mat = data.frame(predict = sample_predict)
  score_mat$predict = factor(score_mat$predict, levels = c(namelist, "Undefined"))
  score_mat$HLA_count = log(rowSums(HLA_matrix) + 1, 10)
  set.seed(321)
  
  # 动态计算 perplexity 对齐 Python
  perplexity_val = min(max(5, floor((nrow(HLA_matrix_norm) - 1) / 3)), 50)
  
  tsne_out = Rtsne(HLA_matrix_norm[rownames(score_mat), ], dims = 2, pca = T,
                   max_iter = 1000, theta = 0.4, perplexity = perplexity_val, 
                   verbose = F, check_duplicates = F)
  tsne_result = as.data.frame(tsne_out$Y)
  colnames(tsne_result) = c("tSNE1", "tSNE2")
  tsne_result = cbind(tsne_result, score_mat)
  
  # --- 6. 动态 KNN 置信度计算与图关联校正 (对齐 analysis.py 逻辑) ---
  library(RANN)
  message("Adjusting predictions based on dynamic neighbors...")
  
  # 预先统计每个预测样本簇的总细胞数
  cluster_sizes = table(tsne_result$predict)
  # 找到全局需要的最大邻居数
  max_k = floor(max(cluster_sizes) / 2)
  
  confidence = rep(0, nrow(tsne_result))
  tsne_result$adjust_predict = as.character(tsne_result$predict)
  
  if (max_k >= 1) {
    point_set = as.matrix(tsne_result[, c("tSNE1", "tSNE2")])
    message(paste0("Querying up to ", max_k, " neighbors..."))
    
    # 批量查询最大邻居数 (包含自身，所以 k = max_k + 1)
    d = nn2(point_set, point_set, k = max_k + 1)
    all_neighbour_indices = d$nn.idx
    
    message("Adjusting labels based on sample-specific neighbors...")
    for (i in 1:nrow(tsne_result)) {
      current_label = as.character(tsne_result$predict[i])
      
      # 动态获取当前细胞所需的邻居数 k_i
      if (current_label %in% names(cluster_sizes)) {
        k_i = floor(cluster_sizes[current_label] / 2)
      } else {
        k_i = 0
      }
      
      if (k_i < 1) next
      
      # 截取 k_i + 1 个邻居
      cell_neighbors = all_neighbour_indices[i, 1:(k_i + 1)]
      # 排除自身索引
      valid_neighbor_indices = cell_neighbors[cell_neighbors != i]
      # 确保切片后长度严格等于 k_i
      if (length(valid_neighbor_indices) > k_i) {
        valid_neighbor_indices = valid_neighbor_indices[1:k_i]
      }
      
      if (length(valid_neighbor_indices) == 0) next
      
      neighbor_labels = as.character(tsne_result$predict[valid_neighbor_indices])
      neighbor_counts = sort(table(neighbor_labels), decreasing = TRUE)
      
      if (length(neighbor_counts) == 0) next
      
      dominate_type = names(neighbor_counts)[1]
      dominate_count = as.numeric(neighbor_counts[1])
      cert = dominate_count / length(valid_neighbor_indices) # 置信度
      confidence[i] = cert
      
      # 如果开启了图关联校正且周围优势细胞不同于自身
      if (graph_corr) {
        if (dominate_type != current_label) {
          if (cert > 0.5) {
            tsne_result$adjust_predict[i] = dominate_type
          } else {
            tsne_result$adjust_predict[i] = "Undefined"
          }
        }
      }
    }
    
    if (!graph_corr) {
      message("Graph-based correlation is disabled. Skipping label adjustment.")
      tsne_result$adjust_predict = as.character(tsne_result$predict)
    }
    
    sample_predict[rownames(tsne_result)] = as.character(tsne_result$adjust_predict)
    sample_predict = factor(sample_predict, levels = c(namelist, "Undefined"))
    
  } else {
    message("Max calculated n_neighbor < 1. Skipping adjustment.")
    tsne_result$adjust_predict = as.character(tsne_result$predict)
  }
  
  # --- 7. TSNE Plot ---
  tsne_result_tmp = tsne_result[!tsne_result$adjust_predict %in% c("Undefined"), ]
  tsne_result_tmp$adjust_predict = factor(tsne_result_tmp$adjust_predict, levels = c(namelist, "Undefined"))
  
  df = aggregate(tsne_result_tmp$tSNE1, by = list(tsne_result_tmp$adjust_predict), mean)
  df$y = aggregate(tsne_result_tmp$tSNE2, by = list(tsne_result_tmp$adjust_predict), mean)$x
  colnames(df) = c("predict", "tSNE1", "tSNE2")
  
  library(ggplot2)
  p1 = ggplot(tsne_result_tmp, aes(tSNE1, tSNE2, color = adjust_predict)) +
    geom_point(size = 0.3) + 
    theme_classic() +
    theme(text = element_text(size = 15)) +
    geom_text(data = df, aes(label = predict), size = 5, color = "black") +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = color)
  
  ggsave(paste0(outdir, "/HLA_TSNE_Adjusted_Prediction.pdf"), p1, width = 10, height = 8)
  
  # --- 8. Count plot ---
  df_count = data.frame(predict_label = sample_predict)
  df_tmp = data.frame(table(df_count$predict_label))
  colnames(df_tmp) = c("Sample", "Cellnum")
  
  p2 = ggplot(df_tmp, aes(x = Sample, y = Cellnum, fill = Sample)) +
    geom_col(alpha = 0.7) + 
    theme_classic() + 
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = color)
  
  ggsave(paste0(outdir, "/Countplot_Adjusted_Prediction.pdf"), p2, width = 10, height = 5)
  
  # --- 9. write demultiplex result ---
  out_csv = data.frame("Cell barcode" = rownames(HLA_matrix_norm),
                       "Sample_prediction" = sample_predict,
                       "Certainty" = confidence)
# Export the analysis source data for figure reproduction.
  write.csv(out_csv, paste0(outdir, "/HLA_demultiplex.csv"), row.names = F)
  
  return(sample_predict)
}

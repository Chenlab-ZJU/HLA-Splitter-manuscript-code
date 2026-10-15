# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: plot_HLA_tSNE.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(stringr)
library(ggsci)
library(Rtsne)

#' Perform t-SNE on HLA matrix and generate a plot
#'
#' @param seurat_obj A Seurat object containing the single-cell data.
#' @param hla_dir Character. The base directory containing the HLA count folders.
#' @param filenames Character vector. The names of the group folders (e.g., c("Group1", "Group2")).
#' @param out_dir Character. Directory to save the output plot.
#' @param seed Numeric. Random seed for t-SNE reproducibility.
#' 
#' @return A list containing the generated ggplot object and the tSNE result dataframe.
plot_HLA_tSNE <- function(seurat_obj, 
                          hla_dir = "./HLA count", 
                          samplename = NA,
                          filenames = NA,
                          out_dir = "./Demultiplex",
                          base_colors = pal_igv("default")(50),
                          thresh = 50,
                          seed = 321) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # 1. Read and combine HLA counts using a list for better performance
  message("Reading HLA counts...")
  hla_list <- list()
  for (i in seq_along(filenames)) {
    grp <- filenames[i]
    label_file <- file.path(hla_dir, grp, "sorted_genotypes.txt")
    barcode_file <- file.path(hla_dir, grp, "output.barcodes.txt")
    mtx_file <- file.path(hla_dir, grp, "HLA_type.mtx")
    
    HLA_label <- read.delim(label_file, header = FALSE)
    barcodes <- read.delim(barcode_file, header = FALSE)
    HLA_matrix <- as.data.frame(as.matrix(readMM(mtx_file)))
    
    colnames(HLA_matrix) <- HLA_label$V1
    
    if (is.na(samplename)){
      rownames(HLA_matrix) <- paste0(barcodes$V1, "-", i)
    }else{
      rownames(HLA_matrix) <- paste0(barcodes$V1, "-", samplename, "-", i)
    }
    
    hla_list[[i]] <- HLA_matrix
  }
  HLA_matrix_All <- do.call(rbind, hla_list)
  
  # 2. Filter Seurat object and match barcodes
  message("Filtering Seurat object and matching barcodes...")
  # Subset to remove Doublets (assuming 'Group' is in meta.data)
  scRNA_filter <- subset(seurat_obj, subset = Group != "Doublets")
  
  common_barcodes <- intersect(rownames(HLA_matrix_All), colnames(scRNA_filter))
  HLA_matrix <- HLA_matrix_All[common_barcodes, ]
  scRNA_filter <- scRNA_filter[, common_barcodes]
  
  message(sprintf("Matched %d barcodes.", length(common_barcodes)))
  
  # 3. Normalize HLA expression
  message("Normalizing and scaling HLA counts...")
  library_sizes <- rowSums(HLA_matrix)
  HLA_matrix_norm <- HLA_matrix / library_sizes * 1e2
  HLA_matrix_norm <- log(HLA_matrix_norm + 1)
  HLA_matrix_norm[is.na(HLA_matrix_norm)] <- 0
  
  # 4. Scale HLA expression (positive values only per column)
  for (i in seq_len(ncol(HLA_matrix_norm))) {
    exp_vals <- HLA_matrix_norm[, i]
    positive_idx <- exp_vals > 0
    if (sum(positive_idx) > 1) {
      exp_positive <- scale(exp_vals[positive_idx], center = FALSE, scale = TRUE)[, 1]
      HLA_matrix_norm[names(exp_positive), i] <- exp_positive
    }
  }
  
  # 5. Prepare TSNE inputs
  message("Running Rtsne...")
  score_mat <- data.frame(predict = scRNA_filter$HLA_Splitter)
  namelist = mixedsort(unique(scRNA_filter$HLA_Splitter))
  score_mat$predict <- factor(score_mat$predict, levels = c(namelist))
  score_mat$HLA_count <- log10(rowSums(HLA_matrix) + 1)
  
  # Run TSNE
  set.seed(seed)
  perplexity_val = min(max(10, floor((nrow(HLA_matrix_norm) - 1) / 3)), 50)
  
  tsne_out = Rtsne(as.matrix(HLA_matrix_norm), dims = 2, pca = TRUE,
                   max_iter = 1000, theta = 0.4, perplexity = perplexity_val, 
                   verbose = F, check_duplicates = F)
  
  tsne_result <- as.data.frame(tsne_out$Y)
  colnames(tsne_result) <- c("tSNE1", "tSNE2")
  tsne_result <- cbind(tsne_result, score_mat)
  
  # 6. TSNE plot
  message("Generating plot...")
  
  # Safely generate colors for the namelist, plus a color for "Undefined" and "Doublets"
  plot_colors <- c(base_colors[1:(length(namelist)-1)], "grey") 
  
  # Calculate centroids for labels
  df_centroids <- aggregate(cbind(tSNE1, tSNE2) ~ predict, data = tsne_result, FUN = mean)
  
  p2 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = predict)) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_text(data = df_centroids, aes(label = predict), size = 4, color = "black", position = position_nudge(y = 2)) +
    theme_void() +
    theme(text = element_text(size = 15)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = plot_colors)
  
  # 7. #统计拆分样本的细胞数量
  df_plot <- data.frame(predict_label = scRNA_filter$HLA_Splitter, Group = scRNA_filter$Group) %>%
    dplyr::count(Group, predict_label, name = "Cellnum")
  group_list = levels(scRNA_filter$Group)
  group_colors <- c(base_colors[1:(length(group_list)-1)], "grey")  
  
  p3 <- ggplot(data = df_plot, aes(x = predict_label, y = Cellnum, fill = Group)) +
    geom_col(alpha = 0.8, width = 0.5) + 
    geom_text(aes(label = Cellnum), vjust = -0.5, size = 4, color = "black") + 
    geom_hline(yintercept = thresh, linetype = "dashed", color = "darkblue", linewidth = 0.8) + # Added dashed line at y=100
    theme_classic() + scale_y_log10() +
    theme(
      text = element_text(size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), # Rotate x-axis labels for readability
      axis.text.y = element_text(size = 12)) +
    scale_fill_manual(values = group_colors) + 
    ylab("Number of Cells") + # Clearer y-axis label
    xlab("Predicted Sample Label")
  
  # Save plot
  out_path <- file.path(out_dir, "HLA_TSNE_Adjusted_Prediction.pdf")
  ggsave(out_path, plot = p2, width = 8, height = 6)
  message(paste("Plot successfully saved to:", out_path))
  
  out_path <- file.path(out_dir, "Predict_sample_number.pdf")
  ggsave(out_path, plot = p3, width = length(namelist)/2+3, height = 4)
  message(paste("Plot successfully saved to:", out_path))
  
}

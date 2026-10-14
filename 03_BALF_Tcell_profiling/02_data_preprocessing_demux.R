# Manuscript section: BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: 01_data_preprocessing_demux.R
# Figure/table scope: Fig3; FigS6; Supplementary Table 13
# Data cohort: 肺泡灌洗液. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

# ==============================================================================
# File: 01_data_preprocessing_demux.R
# Description: 数据合并、质控、HLA-Splitter拆分可视化及 Scrublet 去双细胞
# ==============================================================================

rm(list = ls())
source(file.path(getOption("hla.code_root"), "03_BALF_Tcell_profiling/01_setup_and_helpers.R"))
dir.create("Figures/QC", recursive = TRUE, showWarnings = FALSE)

# 1. 读取单细胞表达矩阵与质控 --------------------------------------------------
filenames <- c("Group1","Group2","Group3","Group4")
scRNA_list <- lapply(seq_along(filenames), function(i) {
  Counts <- Read10X(data.dir = paste0("./Data/GEM/", filenames[i], "/filtered_feature_bc_matrix/"))
  obj <- CreateSeuratObject(Counts, project = filenames[i], assay = "RNA")
  obj <- RenameCells(obj, new.names = paste0(str_sub(colnames(obj), 1, 16), "-", i))
  return(obj)
})
scRNA_Intergrate <- merge(scRNA_list[[1]], y = scRNA_list[-1], project = "lavage_Total")


# 2. 单细胞 QC 与过滤 ----------------------------------------------------------
scRNA_Intergrate[["percent.mt"]] <- PercentageFeatureSet(scRNA_Intergrate, pattern = "^MT-")
scRNA_Intergrate <- subset(scRNA_Intergrate, subset = nCount_RNA > 200 & percent.mt < 20)

# 批量绘制QC小提琴图
qc_metrics <- list(nFeature = c("nFeature_RNA", 5000), 
                   nCount = c("nCount_RNA", 10000), 
                   MT = c("percent.mt", 20))

dir.create("./Figures/QC", recursive = TRUE, showWarnings = FALSE)
for(name in names(qc_metrics)) {
  p <- ggplot(scRNA_Intergrate@meta.data, aes(x = orig.ident, y = .data[[qc_metrics[[name]][1]]])) +
    geom_violin(aes(fill = orig.ident), width = 1, alpha = 0.5) +
    geom_jitter(aes(color = orig.ident), width = 0.2, size = 0.01) +
    geom_boxplot(fill = "white", width = 0.3, outlier.size = 0) +
    ylim(c(0, as.numeric(qc_metrics[[name]][2]))) + theme_bw()
  ggsave(paste0("./Figures/QC/Vln_", name, ".pdf"), plot = p, height = 5, width = 8)
}

# Export the analysis source data for figure reproduction.
saveRDS(scRNA_Intergrate, "./Result/lavage_scRNA_Intergrate.rds")

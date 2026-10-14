# Manuscript section: CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs
# Purpose: Upstream cohort processing and reusable helper definitions.
# Original source: 02_demultiplexing_and_doublets.R
# Figure/table scope: Fig2; FigS4; FigS5
# Data cohort: 脑脊液样本. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.

# ==============================================================================
# File: 02_demultiplexing_and_doublets.R
# Description: 整合 HLA-Splitter 拆分结果，去除 Doublets
# ==============================================================================
rm(list = ls())
source(file.path(getOption("hla.code_root"), "02_CSF_immune_profiling/01_setup_and_helpers.R"))

# Read the saved object or source table used by this analysis.
scRNA_Intergrate <- readRDS("./Result/CSF_scRNA_Intergrate.rds")
filenames <- c("Batch1","Batch2","Batch3")

# 1. 整合 HLA-Splitter 拆分结果 ------------------------------------------------
source(file.path(getOption("hla.code_root"), "01_HLA_Splitter_benchmarks/01_HLA_demultiplex_v2.R"))
sampleInfo <- read.csv("./Data/HLA_list/HLA-genotypes.csv", header = TRUE, row.names = 1)
HLA_list <- lapply(sampleInfo, unique)

for (i in seq_along(filenames)) {
  base_path <- paste0("./Data/HLA_counts/", filenames[i])
  HLA_label <- read.delim(paste0(base_path, "/sorted_genotypes.txt"), header = FALSE)
  barcodes <- read.delim(paste0(base_path, "/output.barcodes.txt"), header = FALSE)
  
  HLA_matrix <- readMM(paste0(base_path, "/HLA_type.mtx")) %>% as.matrix() %>% as.data.frame()
  colnames(HLA_matrix) <- HLA_label$V1
  rownames(HLA_matrix) <- barcodes$V1
  
  outdir <- paste0("./Result/Demultiplex/", filenames[i])
  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  sample_predict <- HLA_demultiplex(HLA_list = HLA_list, HLA_matrix = HLA_matrix, 
                                    graph_corr = TRUE, outdir = outdir)
}

# 提取拆分结果回注 Metadata
demultiplex_result <- do.call(rbind, lapply(seq_along(filenames), function(i) {
  tmp <- read.csv(paste0("./Result/Demultiplex/", filenames[i], "/HLA_demultiplex.csv"))
  tmp$Cell.barcode <- paste0(tmp$Cell.barcode, "-", i)
  return(tmp)
}))
rownames(demultiplex_result) <- demultiplex_result$Cell.barcode

scRNA_Intergrate$HLA_Splitter <- demultiplex_result[colnames(scRNA_Intergrate), "Sample_prediction"]
scRNA_Intergrate$HLA_Splitter[is.na(scRNA_Intergrate$HLA_Splitter)] <- "Undefined"

# 2. 去除 Doublets (scrubletR) -------------------------------------------------
options(future.globals.maxSize = 10000 * 1024^2)

scRNA_Intergrate <- NormalizeData(scRNA_Intergrate) %>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData() %>% RunPCA() %>% 
  RunHarmony(group.by.vars = "orig.ident") %>% 
  JoinLayers()

scRNA_Intergrate <- scrublet(scRNA_Intergrate)
scRNA_Intergrate <- subset(scRNA_Intergrate, subset = doublet_scores < 0.2)

# 3. 添加分组信息与保存 --------------------------------------------------------
SampleInfo <- read.csv("./Data/SampleInfo.csv")
scRNA_Intergrate$Group <- SampleInfo$Group[match(scRNA_Intergrate$HLA_Splitter, SampleInfo$Patient.ID)]
scRNA_Intergrate$Group[is.na(scRNA_Intergrate$Group)] <- 'Undefined'
Group_levels <- c("Metastases", "Hemorrhage", "Encephalitis", "Undefined")
scRNA_Intergrate$Group <- factor(scRNA_Intergrate$Group, levels = Group_levels)
scRNA_Intergrate$HLA_Splitter = factor(scRNA_Intergrate$HLA_Splitter, 
                                       levels = mixedsort(unique(scRNA_Intergrate$HLA_Splitter)))

# Export the analysis source data for figure reproduction.
saveRDS(scRNA_Intergrate, "Result/CSF_scRNA_Intergrate.rds")


# 绘制总体的HLA t-SNE 图 ------------------------------------------------
#scRNA_Intergrate = readRDS("Result/CSF_scRNA_Intergrate.rds")
scRNA_Intergrate = scRNA_Intergrate[, scRNA_Intergrate$HLA_Splitter != "Undefined"]

results <- plot_HLA_tSNE(seurat_obj = scRNA_Intergrate, 
                         hla_dir = "./Data/HLA_counts/",
                         filenames = c("Batch1","Batch2","Batch3"),
                         out_dir = "./Result/Demultiplex",
                         base_colors = GLOBAL_COLORS)

df <- data.frame(HLA_demutliplex = scRNA_Intergrate$HLA_Splitter)
df <- table(df) %>% prop.table(1)*100
write.csv(df, paste0("./Result/Demultiplex/demultiplex_res.csv"))

## 比较样本性别基因与拆分结果一致性 ------------------------------------------------
#scRNA_Intergrate = readRDS("Result/CSF_scRNA_Intergrate.rds")
scRNA_Intergrate = scRNA_Intergrate[, scRNA_Intergrate$HLA_Splitter != "Undefined"]

gene_matrix <- AverageExpression(scRNA_Intergrate, group.by = c("HLA_Splitter"), assays = "RNA", return.seurat = TRUE)
gene_matrix <- gene_matrix@assays$RNA$scale.data %>% as.matrix()

Sample_Info = read.csv("./Data/SampleInfo.csv")
annotation_col = data.frame(Gender = Sample_Info$Gender[match(colnames(gene_matrix), Sample_Info$Patient.ID)])
 
cell_counts <- table(scRNA_Intergrate$HLA_Splitter)
annotation_col$Cell_num = as.numeric(cell_counts[colnames(gene_matrix)])
rownames(annotation_col) = colnames(gene_matrix)

# XIST: 女性特异表达(X不激活); RPS4Y1, EIF1AY, DDX3Y, KDM5D: 男性Y染色体特异表达
gender_genes <- c("XIST", "RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "UTY")

gender_genes <- intersect(gender_genes, rownames(gene_matrix))

sample_order <- order(annotation_col$Gender)
mat <- gene_matrix[gender_genes, sample_order, drop = FALSE]
annotation_col <- annotation_col[sample_order, , drop = FALSE]

library(pheatmap)
plot_heatmap = pheatmap(mat, 
                        cluster_rows = TRUE,    # 开启行聚类，X和Y染色体的基因会自动分块
                        cluster_cols = FALSE,   # 关闭列聚类，保持我们上面手动按性别排序的顺序
                        border_color = "black", 
                        fontsize = 10, 
                        cellwidth = 20, 
                        cellheight = 20,
                        annotation_col = annotation_col)

ggsave(paste0("./Result/demultiplex/Gender_related_GENEs.pdf"), plot_heatmap, dpi = 300, height = 6, width = 10)

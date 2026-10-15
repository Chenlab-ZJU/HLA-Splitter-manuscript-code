# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: NMF analysis.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

GetData = function(
    srt,
    genes = NULL,
    slot = 'scale.data',
    assay = NULL
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (is.null(genes)){
    if ('RNA' %in% names(srt@assays)){
      genes = rownames(GetAssayData(srt, assay = 'RNA', slot = 'counts'))
    } else if  ('Spatial' %in% names(srt@assays)){
      genes = rownames(GetAssayData(srt, assay = 'Spatial', slot = 'counts'))
    } else {
      genes = rownames(GetAssayData(srt, assay = assay, slot = 'counts'))
    }
  }
  data = GetAssayData(srt, assay = assay, slot = slot)
  missing = setdiff(genes, rownames(data))
  add = matrix(0, nrow = length(missing), ncol = ncol(data))
  rownames(add) = missing
  data = rbind(data, add)
  data = data[genes, ]
  return(data)
}
NMFToModules = function(
    res,
    gmin = 5
){
  
  scores = basis(res)
  coefs = coefficients(res)
  
  # Remove if fewer than gmin genes
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  l = sapply(modules, length)
  keep = (l >= gmin)
  scores = scores[, keep]
  coefs = coefs[keep, ]
  
  # Find modules
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  
  names(modules) = sapply(modules, '[', 1)
  names(modules) = paste('m', names(modules), sep = '_')
  names(modules) = gsub('-','_',names(modules))
  
  return(modules)
}
BuildEnrichmentMatrix = function(
    genes,
    type = 'GO',
    db = NULL
){
  if (is.null(db)){
    db = FindMSigDB(type)
  }
  terms = names(db)
  enrichment.matrix = sapply(db, function(term){
    genes %in% term
  })
  rownames(enrichment.matrix) = genes
  #enrichment.matrix = enrichment.matrix[, colSums(enrichment.matrix) > 0]
  enrichment.matrix[, colSums(enrichment.matrix) == 0] = NA
  return(enrichment.matrix)
}
FindMSigDB = function(
    type
){
  library(msigdbr)
  if (type == 'GO'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C5')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
  } else if (type == 'HALLMARK'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'H')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
  } else if (type == 'MOTIF'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C3')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[-grep('UNKNOWN|MIR',names(db))]
  } else if (type == 'PATHWAYS'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C2')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[grep('BIOCARTA|REACTOME|KEGG',names(db))]
  } else if (type == 'BIOCARTA'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C2')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[grep('BIOCARTA',names(db))]
  } else if (type == 'KEGG'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C2')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[grep('KEGG',names(db))]
  } else if (type == 'REACTOME'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C2')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[grep('REACTOME',names(db))]
  } else {
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C5')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[grep(type, names(db), ignore.case = TRUE, value = TRUE)]
  }
  return(db)
}
loadRData = function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
Heatmap = function(
    matrix,
    colors = NULL,
    breaks = NULL,
    col = NULL,
    clustering_method_rows = 'ward.D2', 
    clustering_method_columns = 'ward.D2',
    clustering_distance_rows = 'pearson',
    clustering_distance_columns = 'pearson',
    show_row_names = FALSE,
    show_column_names = FALSE,
    is.symmetric = FALSE,
    ...
){
  if (is.null(breaks)){
    mi = min(matrix, na.rm = TRUE)
    ma = max(matrix, na.rm = TRUE)
    if (mi < 0 & ma > 0){
      breaks = c(mi, 0, ma)
    } else {
      breaks = c(mi, (mi + ma)/2, ma)
    }
  }
  if (is.null(colors)){
    colors = c('blue','white','red')
  }
  if (is.null(col)){
    col = colorRamp2(breaks = breaks, colors = colors)
  }
  if (is.symmetric == TRUE){
    h = ComplexHeatmap::Heatmap(matrix, 
                                col = col,
                                clustering_method_rows = clustering_method_rows, 
                                clustering_distance_rows = clustering_distance_rows,
                                cluster_columns = FALSE,
                                show_row_names = show_row_names,
                                show_column_names = show_column_names,
                                ...)
    o = unlist(row_order(h))
    return(ComplexHeatmap::Heatmap(matrix, 
                                   col = col,
                                   clustering_method_rows = clustering_method_rows, 
                                   clustering_distance_rows = clustering_distance_rows,
                                   cluster_columns = FALSE,
                                   column_order = o,
                                   show_row_names = show_row_names,
                                   show_column_names = show_column_names, 
                                   ...))
  } else {
    return(ComplexHeatmap::Heatmap(matrix, 
                                   col = col,
                                   clustering_method_rows = clustering_method_rows, 
                                   clustering_method_columns = clustering_method_columns,
                                   clustering_distance_rows = clustering_distance_rows,
                                   clustering_distance_columns = clustering_distance_columns,
                                   show_row_names = show_row_names,
                                   show_column_names = show_column_names,
                                   ...))
  }
  
}
GeneToEnrichment <- function (srt, type = "GO", db = NULL, method = "rand", genes = NULL, 
                              assay = NULL, do.rescale = FALSE, min.cells = 0, min.genes = 0, 
                              min.var = 0, min.var.rescaled = 0, auc_percentile = 0.05, 
                              db_rand = NULL, nrand = 4, nbin = 25, ...) 
{
  if (is.null(assay)) {
    assay = DefaultAssay(srt)
  }
  if (is.null(db)) {
    db = FindMSigDB(type)
  }
  counts = as.matrix(GetData(srt, assay = assay, slot = "counts"))
  genes = rownames(counts)
  genes.expr = rownames(counts)[rowSums(counts) > min.cells]
  if (method == "metagene") {
    data = as.matrix(GetAssayData(srt, assay = assay, slot = "scale.data"))
    db = lapply(db, intersect, genes.expr)
    enrichment.profile = t(sapply(names(db), function(m) {
      colMeans(data[db[[m]], ], na.rm = TRUE)
    }))
    enrichment.profile = enrichment.profile[sapply(names(db), 
                                                   function(x) {
                                                     v = var(enrichment.profile[x, ])
                                                     l = length(db[[x]])
                                                     return(l > min.genes && v > min.var && v * l^2 > 
                                                              min.var.rescaled)
                                                   }), ]
    if (do.rescale) {
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn)/sqrt(v)
    }
    srt = AddMetaData(srt, t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  if (method == "auc") {
    data = as.matrix(GetData(srt, assay = assay, slot = "data"))
    cells_rankings = AUCell_buildRankings(data)
    cells_AUC = AUCell_calcAUC(db, cells_rankings, aucMaxRank = nrow(cells_rankings) * 
                                 auc_percentile)
    enrichment.profile = getAUC(cells_AUC)
    if (do.rescale) {
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn)/sqrt(v)
    }
    srt = AddMetaData(srt, t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  if (method == "score") {
    temp = AddModuleScore(srt, features = db, assay = assay, 
                          name = names(db), nbin = nbin, ...)
    enrichment.profile = t(temp@meta.data[, names(db)])
    if (do.rescale) {
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn)/sqrt(v)
    }
    srt = AddMetaData(srt, t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  if (method == "rand") {
    data = as.matrix(GetData(srt, assay = assay, slot = "scale.data"))
    db = lapply(db, intersect, genes)
    if (is.null(db_rand)) {
      db_rand = MakeRand(srt, db, nrand = nrand, nbin = nbin)
    }
    else {
      nrand = log10(length(db_rand[[1]]))
    }
    enrichment.profile = t(sapply(names(db), function(m) {
      ra = sapply(db_rand[[m]], function(i) {
        colMeans(data[i, ], na.rm = TRUE)
      })
      re = colMeans(data[db[[m]], ], na.rm = TRUE)
      p = rowMeans(ra >= re)
      p = -log10(p)
      return(p)
    }))
    enrichment.profile[is.infinite(enrichment.profile)] = nrand
    enrichment.profile = enrichment.profile/nrand
    srt = AddMetaData(srt, t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  return(srt)
}
MakeRand <-function (srt, db, assay = NULL, nrand = 3, nbin = 25) 
{
  if (is.null(assay)) {
    assay = DefaultAssay(srt)
  }
  data = GetData(srt, slot = "data")
  db = lapply(db, intersect, rownames(data))
  data.avg = sort(rowMeans(x = data))
  data.cut = cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                        n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) = names(x = data.avg)
  binned = split(names(data.cut), data.cut)
  db_rand = lapply(names(db), function(m) {
    lapply(1:10^nrand, function(i) {
      used = vector()
      unused = binned
      for (g in db[[m]]) {
        pool = data.cut[g]
        new = sample(unused[[pool]], 1)
        used = c(used, new)
        unused[[pool]] = setdiff(unused[[pool]], new)
      }
      return(used)
    })
  })
  names(db_rand) = names(db)
  return(db_rand)
}
#读取数据####

#设置运行参数####
range = 4:12 # 迭代次数
gmin = 5 #最低关联数量
ncores = 1 #运行核数量

#拆分单样本数据集####
cell_subset = Myeloid
DefaultAssay(cell_subset) <- 'RNA'
patient_list <- SplitObject(object = cell_subset, split.by = "HLA_predict")
# Export the analysis source data for figure reproduction.
saveRDS(patient_list ,file = "./NMF/patient_list.rds")
# Read the saved object or source table used by this analysis.
patient_list <- readRDS(file = "./NMF/patient_list.rds")
patient_list_names <- levels(cell_subset$HLA_predict)
for(i in patient_list_names){
  saveRDS(patient_list[[i]],file = paste0('./NMF/patient_RDS/',i,'_Seurat.RDS'))
}
patient_list_names = patient_list_names[table(cell_subset$HLA_predict)>50  ]

#install.packages("NMF")
library(NMF)

for(i in patient_list_names) {
  obj <- readRDS(file = paste0('./NMF/patient_RDS/',i,'_Seurat.RDS'))
  obj <- SCTransform(obj, return.only.var.genes = TRUE)
  saveRDS(obj,file = paste0('./NMF/patient_RDS/',i,'_Seurat.RDS'))
  data <- GetAssayData(obj, slot = 'scale.data') %>% as.matrix()
  
  data <- data[VariableFeatures(obj),]
  rm(obj)
  gc()
  data[data < 0] <- 0
  data <- data[apply(data, 1, var) > 0, ]
  res.list <- parallel::mclapply(range, function(r){
    nmf(data, rank = r, nrun = 1, seed = 'ica', method = 'nsNMF')
  }, mc.cores = ncores)
  names(res.list) = range
  saveRDS(res.list,file = paste0("./NMF/result/",i,"res.list.rds"))
  print(paste0(i, " is already over!"))
}

#拆分后重新读取NMF结果####
res.list_L<- list()
for(i in patient_list_names){
  a <- readRDS(paste0("./NMF/result/",i ,"res.list.rds"))
  res.list_L[[i]] <- a
}
names(res.list_L) <- patient_list_names

#### 获取最优的NMF聚类rank，提取 ####
modules_L <- list()
res_optima_L <- list()
for (i in 1:length(patient_list_names)){ 
    modules.list = lapply(res.list_L[[i]], NMFToModules, gmin = gmin)
    print(sapply(modules.list,length))
    comp = as.numeric(names(modules.list)) - sapply(modules.list, length)
    mi = min(comp)
    r = names(which(comp == mi))
    r = r[length(r)]
    res = res.list_L[[i]][[r]]
    modules = NMFToModules(res, gmin = gmin)   
    modules_L[[i]] <- modules
    res_optima_L[[i]] <- res  
}
saveRDS(object = modules_L,file = "./NMF/result/modules_L.rds")
saveRDS(object = res_optima_L,file = "./NMF/result/res_optima_L.rds")

names(res_optima_L) <- patient_list_names
modules.list <- lapply(res_optima_L, NMFToModules)

all = unlist(modules.list, recursive = FALSE, use.names = FALSE)
names(all) = unlist(sapply(modules.list, names))
ta = table(unlist(all))
genes.use = names(ta)[ta > 1]
for (i in 1:5){
  all = unlist(modules.list, recursive = FALSE, use.names = TRUE)
  all = lapply(all, intersect, genes.use)
  sim = sapply(all, function(x){
    sapply(all, function(y){
      length(intersect(x,y))/length(union(x,y))
    })
  })
  keep = rownames(sim)[apply(sim, 1, function(x){
    sum(x > 0.05) >= 3 #大于等于模块数相关的情况，我这里2或3就行
  })]
  print(keep)
  all = all[keep]
  modules.list = lapply(names(modules.list), function(x){
    li = modules.list[[x]]
    li[names(li)[paste(x,names(li),sep='.') %in% keep]]
  })
  names(modules.list) =  names(res_optima_L)
  ta = table(unlist(all))
  genes.use = names(ta)[ta > 1] 
  print(length(all))
}
options(repr.plot.height = 7.5, repr.plot.width = 8)
pheatmap(sim, show_rownames = F, show_colnames = F)

sub = matrix(0, nrow = length(genes.use), ncol = length(genes.use))
rownames(sub) = genes.use
colnames(sub) = genes.use
for (s in names(modules.list)){
  for (mod in modules.list[[s]]){
    mod = intersect(mod, genes.use)
    for (x in mod){
      for (y in mod){
        sub[x,y] = sub[x,y] + 1
      }
    }
  }
}
diag(sub) = 0
adj_keep = sub
# Remove low connections
library(igraph)
adj = adj_keep
adj[] = (adj >= 3)
#adj[adj <= 1] = 0
for (i in 1:5){
  keep = names(which(rowSums(adj) >= 5))
  adj = adj[keep,keep]
  print(dim(adj))
}
g = graph_from_adjacency_matrix(adj, diag = FALSE, mode = 'undirected', weighted = TRUE)
modules = communities(cluster_infomap(g, nb.trials = 100))
names(modules) = paste0('m_', sapply(modules, '[', 1))
save(modules,file = "./NMF/result/modules.Rdata")


#差异基因GO功能富集分析####
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
#基因ID转换为ENTREZID
for(i in 1:length(modules)){
  modules[[i]] <- bitr(modules[[i]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  modules[[i]] <- modules[[i]][,2]
}
#GO富集分析
ego_GO <- list()
for(i in 1:length(modules)){
  ego_GO[[i]] <- enrichGO(gene = modules[[i]],
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "ALL", #CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程
                          pAdjustMethod = "BH",
                          minGSSize = 1,
                          pvalueCutoff = 0.01, #通常为0.01或0.05
                          qvalueCutoff = 0.05,
                          readable = TRUE)
}
library(Hmisc)
#根据GO富集的信息最终注释为以下模块
ego_GO[[9]]@result[["Description"]]
module_GO <- c("lipid transport","response to chemokine","NA.1","response to metal ion",
               "Antigen processing and presentation","response to unfolded protein",
               "response to virus","NA.2")
load("./NMF/result/modules.Rdata")
for(i in 1:length(modules)){
  names(modules)[[i]] <- module_GO[i]
}
sim = sapply(all, function(x){
  sapply(modules, function(y){
    pval = phyper(length(intersect(x, y)), length(x), 2*10^4 - length(x), length(y), lower.tail = FALSE)
    return(-log10(pval))
  })
})
sim[is.infinite(sim)] = max(sim[is.finite(sim)])
sim = sim[,apply(sim, 2, function(x){any(x > 3)})]
df = data.frame('sample' = sapply(colnames(sim), function(x){
  y = sapply(strsplit(x, '.', fixed = TRUE), '[', 1)
}))
df$top = apply(sim, 2, which.max)
df$top = factor(names(modules)[df$top], levels = names(modules))


cell_subset.meta <- cell_subset@meta.data
cell_subset.meta$Molecular.typing[cell_subset.meta$Molecular.typing == "Triple negative tumor"] <- "TNBC"
Epi.meta$Molecular.typing[Epi.meta$Molecular.typing == "Her-2"] <- "HER2+"
Epi.meta$Molecular.typing[Epi.meta$Molecular.typing == "ER+ tumour_primary"] <- "ER+"
Epi.meta$Molecular.typing[Epi.meta$Molecular.typing == "ER+ tumour_lymph"] <- "ER+"
Epi.meta$Molecular.typing[Epi.meta$Molecular.typing == "Triple negative BRCA1 tumour"] <- "TNBC"
Epi.meta$Molecular.typing[Epi.meta$Molecular.typing == "HER2+ tumour"] <- "HER2+"
Epi.meta$Molecular.typing[Epi.meta$Molecular.typing == "PR+ tumour"] <- "PR+"
Epi.meta$Molecular.typing[Epi.meta$Molecular.typing == "ER+ tumour"] <- "ER+"
names(Epi.meta)[2] <- "sample"
Epi.meta <- distinct(Epi.meta)
df$subtype <- "non"
for(i in 1:dim(df)[1]){
  df$subtype[i] <- Epi.meta$Molecular.typing[which(Epi.meta$sample == df$sample[i])]
}
side_ann = HeatmapAnnotation(df = df[, c('subtype','top')], 
                             which = 'row')
library(ComplexHeatmap)
options(repr.plot.height = 7.5, repr.plot.width = 8)
h = Heatmap(name = 'Module overlap p-value', sim, 
            cluster_columns = FALSE, cluster_rows = FALSE,
            column_order = order(df$top),
            show_row_names = TRUE,
            border_gp = gpar(lty = 2, col = "red"), show_heatmap_legend = T,
            heatmap_legend_param = list(
              title = "Overlap p-value",
              legend_height = unit(4, "cm"))
)
print(h)

ovlp = sapply(all, function(x){
  sapply(all, function(y){
    pval = phyper(length(intersect(x, y)), length(x), 2*10^4 - length(x), length(y), lower.tail = FALSE)
    return(-log10(pval))
  })
})
ovlp[is.infinite(ovlp)] = max(ovlp[is.finite(ovlp)])
ovlp = ovlp[colnames(sim),colnames(sim)]
save(ovlp,file = "./ovlp.Rdata")
bk <-  seq(0, 6,length = 7)
mycolors<-rev(brewer.pal(7, "RdBu"))
set.seed(123)
h = Heatmap(name = 'Module overlp p-value', ovlp,column_order = order(df$top),  row_order = order(df$top),
            breaks = bk,colors =mycolors
) + side_ann
print(h)
fi = data.frame(sapply(modules,'[', 1:max(sapply(modules, length))), stringsAsFactors = FALSE)
fi[is.na(fi)] = ''
write.table(fi, file = './modules.csv', sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)

#根据细胞可能富集的模块计算模块活性####
library(AUCell)
calAUCell <- function(mode,x){
  cells_rankings <- AUCell_buildRankings(x) #x为表达矩阵
  cells_AUC <- AUCell_calcAUC(mode, cells_rankings, aucMaxRank = nrow(cells_rankings)*0.1)
  return(cells_AUC)
}
filterAUCell <- function(x){
  cells_assignment <-AUCell_exploreThresholds(x, plotHist = TRUE, nCores = 1, assign = TRUE) #x为cells_AUC
  return(cells_assignment)
}
cells_rankings <- AUCell_buildRankings(cell_subset@assays$RNA@data)
modules_list2 <- list()
for(i in names(modules)){
  modules_list2[[i]] <- intersect(modules[[i]],all_gene)
}
cells_AUC <- AUCell_calcAUC(modules_list2, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
warningMsg <- sapply(cells_assignment, function(x){x$aucThr$comment})
warningMsg[which(warningMsg!="")]
AUC <- getAUC(cells_AUC)
sample.meta <- unique(cell_subset@meta.data$orig.ident)
for(i in 1:length(sample.meta)){
  sample.data <- AUC[,grep(sample.meta[i],colnames(AUC))]
  medians <- apply(sample.data,1,median)
  if(i == 1){
    AUC.sample <- data.frame(medians)
    colnames(AUC.sample)[i] <- sample.meta[i]
  } else {
    AUC.sample <- cbind(AUC.sample, data.frame(medians))
    colnames(AUC.sample)[i] <- sample.meta[i]
  }
}
sample.meta <- cell_subset@meta.data[,c("orig.ident","TumorType")]
sample.meta <- distinct(sample.meta)
p.ind <- c()
l.ind <- c()
for(i in 1:dim(AUC.sample)[2]){
  type <- sample.meta[which(colnames(AUC.sample)[i] == sample.meta[,1]),2]
  if(type == "Primary"){
    p.ind <- c(p.ind,i)
  } else if(type == "Lymph"){
    l.ind <- c(l.ind,i)
  }
}
AUC.L <- AUC.sample[,l.ind]
AUC.P <- AUC.sample[,p.ind]
sample.data <- data.frame(value = 0, id = "GSM", module = "mode", type = "tumor")
for(i in 1:dim(AUC.L)[2]){
  data <- data.frame(value = AUC.L[,i], id = colnames(AUC.L)[i],
                     module = rownames(AUC.L), type = "Lymph")
  sample.data <- rbind(sample.data, data)
}
for(i in 1:dim(AUC.P)[2]){
  data <- data.frame(value = AUC.P[,i], id = colnames(AUC.P)[i],
                     module = rownames(AUC.P), type = "Primary")
  sample.data <- rbind(sample.data, data)
}
sample.data <- sample.data[-1,]
library(ggpubr)
library(ggsignif)
options(repr.plot.height = 8, repr.plot.width = 10)
ggboxplot(sample.data, x = "module", y = "AUC value",
          color = "type", add = "jitter")+ #palette可以按照期刊选择相应的配色，如"npg"等
  stat_compare_means(aes(group = type, label = paste0("p = ",signif(as.numeric(..p.format..), 2))))+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size = 12, family = "serif"),
        axis.text.x = element_text(angle = 90))

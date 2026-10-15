# Manuscript section: Shared analysis and plotting functions
# Purpose: Transitive shared function dependency.
# Original source: HLA_TSNE.R
# Figure/table scope: Referenced by analysis scripts
# Data cohort: shared. Inputs use the compact English data layout.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

#HLA TSNE聚类 ####
library(Rtsne)
score_mat = data.frame(predict=sample_predict,sample_label=sample_label)
score_mat$sample_label = factor(score_mat$sample_label,levels=namelist)
score_mat$predict = factor(score_mat$predict,levels=c(namelist,"undefined"))
set.seed(321)
tsne_out = Rtsne(HLA_matrix_norm,dims = 2,pca = T,max_iter = 1000,theta = 0.4,perplexity = 20,verbose = F,check_duplicates=F)
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) = c("tSNE1","tSNE2")
tsne_result = cbind(tsne_result,score_mat)

color = pal_igv("default")(51)
#show_col(color)

p1 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=sample_label)) +
  geom_point(size=0.3)+theme_classic()+
  theme(text = element_text (size = 15))+
  guides(color=guide_legend(override.aes =list(size = 5)))+
  scale_color_discrete(type = c(color))

p2 = ggplot(tsne_result,aes(tSNE1,tSNE2,color=predict)) +
  geom_point(size=0.3)+theme_classic()+
  theme(text = element_text (size = 15))+
  guides(color=guide_legend(override.aes =list(size = 5)))+
  scale_color_discrete(type = c(color[1:length(namelist)],"grey"))
p1+p2
ggsave(paste0("./Figures/HLA_TSNE.pdf"),p1+p2,dpi=300,height=3,width = 8)

# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Shared algorithm / evaluation function.
# Original source: Demultiplex_evaluation.R
# Figure/table scope: Fig1c; Fig1d; FigS3
# Data cohort: Syn_data. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

demultiplex_evaluation = function(sample_label,sample_predict,outdir){
  namelist = unique(sample_label)
  color = pal_igv("default")(length(namelist)+1)
  #label proportion####
  df = data.frame(label = sample_label,predict = sample_predict)
  df$class = "Correct"
  df$class[as.character(df$label)!=as.character(df$predict)] = "Incorrect"
  df$class[df$predict=="Undefined"] = "Undefined"
  
  df$class = factor(df$class, levels = c("Correct","Undefined","Incorrect"))
  df_tmp = data.frame(prop.table(table(df$class))*100)
  df_tmp$label = paste0(df_tmp$Var1,": ",round(df_tmp$Freq,2),"%")
  
  p = ggplot(df_tmp,aes(x= "ratio",y=Freq,fill=label))+
    geom_col(position="stack",alpha=0.5)+coord_polar(theta = 'y')+
    theme_void()+theme(text = element_text (size =15))+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    scale_fill_discrete(type = c("#CEDE94","firebrick","grey"))
  p
  ggsave(paste0(outdir,"/label_proportion.pdf"),p,height=5,width = 8)
  
  #计算回收样本比例####
  df_tmp = data.frame(prop.table(table(df[,c(1,3)]),1)*100)
  df_tmp$label = factor(df_tmp$label,levels = namelist)
  df_tmp = df_tmp[!is.na(df_tmp$label),]
  p = ggplot(df_tmp ,aes(x = class,y=Freq))+
    geom_boxplot(fill="white",color="black",outliers = F)+
    geom_point(aes(color=label), position = position_jitter(width = 0.3),size=2)+
    theme_bw()+ylab("Percentage in each sample(%)")+
    scale_color_discrete(type = c(color))+
    theme_classic()+theme(text = element_text (size =15),
                          plot.title = element_text(hjust = 0.5),
                          axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75))
  p
  ggsave(paste0(outdir,"/predict_fraction.pdf"),p,dpi=300,height=6,width=6)
# Export the analysis source data for figure reproduction.
  write.csv(df_tmp,paste0(outdir,"/Sample_recover_ratio.csv"))
  
  #计算预测准确率
  predict_result = data.frame(sample_label=sample_label,predict=sample_predict)
  predict_result$sample_label = as.character(predict_result$sample_label)
  predict_result$predict = as.character(predict_result$predict)
  Accuracy =c()
  Precision=c()
  Recall =c()
  for (i in 1:length(namelist)){
    TP = nrow(predict_result[predict_result$sample_label==namelist[i] & predict_result$predict==namelist[i],])/nrow(predict_result)
    TN = nrow(predict_result[predict_result$sample_label!=namelist[i] & predict_result$predict!=namelist[i],])/nrow(predict_result)
    FP = nrow(predict_result[predict_result$sample_label!=namelist[i] & predict_result$predict==namelist[i],])/nrow(predict_result)
    FN = nrow(predict_result[predict_result$sample_label==namelist[i] & predict_result$predict!=namelist[i],])/nrow(predict_result)
    Accuracy = append(Accuracy,(TP+TN))
    Precision = append(Precision,TP/(TP+FP))
    Recall = append(Recall,TP/(TP+FN))
  }
  fraction = table(predict_result) %>% data.frame()
  
  #计算Accuracy,Precision,Recall,F1等参数
  F1 = 2*Precision*Recall/(Precision+Recall)
  df = data.frame(sample=namelist,Accuracy=Accuracy,Precision=Precision,Recall=Recall,F1=F1)
  write.csv(df,paste0(outdir,"/performance.csv"))
}

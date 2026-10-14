# Urinary immune-cell states distinguish renal allograft rejection phenotypes

Manuscript scope: Fig. 4; Supplementary Figs. S7-S8.

Run individual scripts using `Rscript run_analysis.R <relative-script-path>`. Numbers organize the manuscript section; they do not imply that all scripts should be rerun consecutively. Use saved objects and source tables for final figure regeneration. Upstream steps can overwrite objects and require annotation review.

- `01_setup_and_helpers.R`: Upstream urine processing and sample-level signature analysis. Original name: `00_setup_and_helpers.R`.
- `02_data_integration_and_QC.R`: Upstream urine processing and sample-level signature analysis. Original name: `01_data_integration_and_QC.R`.
- `03_demultiplexing.R`: Upstream urine processing and sample-level signature analysis. Original name: `02_demultiplexing.R`.
- `04_global_clustering.R`: Upstream urine processing and sample-level signature analysis. Original name: `03_global_clustering.R`.
- `05_AUCell_scoring.R`: Upstream urine processing and sample-level signature analysis. Original name: `08_AUCell_scoring.R`.
- `06_Pseudobulk_analysis.R`: Upstream urine processing and sample-level signature analysis. Original name: `09_Pseudobulk_analysis.R`.
- `07_scRNA_UMAP_celltype_and_TCMR_ABMR_NR_proportion.R`: Final urine immune-cell composition, module and demultiplexing figures. Original name: `26_revision_scRNA_UMAP_celltype_and_TCMR_ABMR_NR_proportion.R`.
- `08_scRNA_UMAP_Rejection_Major_celltype.R`: Final urine immune-cell composition, module and demultiplexing figures. Original name: `42_revision_scRNA_UMAP_Rejection_Major_celltype.R`.
- `09_scRNA_celltype_propration_group_compare.R`: Final urine immune-cell composition, module and demultiplexing figures. Original name: `41_revision_scRNA_celltype_propration_group_compare.R`.
- `10_demultiplex_immune_cell_ratio_plots.R`: Final urine immune-cell composition, module and demultiplexing figures. Original name: `27_revision_demultiplex_immune_cell_ratio_plots.R`.
- `11_all_cells_rejection_module_heatmap_score.R`: Final urine immune-cell composition, module and demultiplexing figures. Original name: `40_revision_all_cells_rejection_module_heatmap_score.R`.
- `12_demultiplex_HLA_tsne_predict_sample_number.R`: Final urine immune-cell composition, module and demultiplexing figures. Original name: `30_revision_demultiplex_HLA_tsne_predict_sample_number.R`.
- `13_scRNA_celltype_composition_CSF_style.R`: Final urine immune-cell composition, module and demultiplexing figures. Original name: `25_revision_scRNA_celltype_composition_CSF_style.R`.
- `14_major_celltype_HLA_cytotoxic_exhaustion_violin.R`: Final urine immune-cell composition, module and demultiplexing figures. Original name: `29_revision_major_celltype_HLA_cytotoxic_exhaustion_violin.R`.
- `15_signature_score_analysis_gene_contribution.R`: Final robust signature construction, scoring and cross-cohort evaluation. Original name: `30_revision_signature_score_analysis_gene_contribution.R`.
- `16_signature_score_advanced_optimization_GSE189536_validation.R`: Final robust signature construction, scoring and cross-cohort evaluation. Original name: `31_revision_signature_score_advanced_optimization_GSE189536_validation.R`.
- `17_signature_cross_cohort_robust_selection.R`: Final robust signature construction, scoring and cross-cohort evaluation. Original name: `32_revision_signature_cross_cohort_robust_selection.R`.
- `18_cross_cohort_robust_signature_distribution_plots.R`: Final robust signature construction, scoring and cross-cohort evaluation. Original name: `34_quick_cross_cohort_robust_signature_distribution_plots.R`.
- `19_cross_cohort_signature_gene_expression_correlation.R`: Final robust signature construction, scoring and cross-cohort evaluation. Original name: `35_cross_cohort_signature_gene_expression_correlation.R`.
- `20_cross_cohort_robust_boxplots_pairwise_target_vs_each_group.R`: Final robust signature construction, scoring and cross-cohort evaluation. Original name: `37_quick_redraw_cross_cohort_robust_boxplots_pairwise_target_vs_each_group.R`.

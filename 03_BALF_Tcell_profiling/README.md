# BALF T cell profiling reveals T cell heterogeneity and TCR clonal alterations in lung diseases

Manuscript scope: Fig. 3; Supplementary Fig. S6.

Run individual scripts using `Rscript run_analysis.R <relative-script-path>`. Numbers organize the manuscript section; they do not imply that all scripts should be rerun consecutively. Use saved objects and source tables for final figure regeneration. Upstream steps can overwrite objects and require annotation review.

- `01_setup_and_helpers.R`: Upstream BALF processing, repertoire and disease-signature computation. Original name: `00_setup_and_helpers.R`.
- `02_data_preprocessing_demux.R`: Upstream BALF processing, repertoire and disease-signature computation. Original name: `01_data_preprocessing_demux.R`.
- `03_demultiplexing_and_doublets.R`: Upstream BALF processing, repertoire and disease-signature computation. Original name: `02_demultiplexing_and_doublets.R`.
- `04_global_clustering.R`: Upstream BALF processing, repertoire and disease-signature computation. Original name: `03_global_clustering.R`.
- `05_tcr_repertoire.R`: Upstream BALF processing, repertoire and disease-signature computation. Original name: `04_tcr_repertoire.R`.
- `06_AUCell_scoring.R`: Upstream BALF processing, repertoire and disease-signature computation. Original name: `05_AUCell_scoring.R`.
- `07_run_sc_pseudobulk_functional_analysis.R`: Upstream BALF processing, repertoire and disease-signature computation. Original name: `run_sc_pseudobulk_functional_analysis.R`.
- `08_lavage_find_signature.R`: Upstream BALF processing, repertoire and disease-signature computation. Original name: `lavage_find_signature.R`.
- `09_lavage_tcell_group_facet_umap.R`: Final BALF figures from saved objects and signatures. Original name: `17_revision_lavage_tcell_group_facet_umap.R`.
- `10_lavage_tcell_subtype_proportion_boxplot_heatmap.R`: Final BALF figures from saved objects and signatures. Original name: `18_revision_lavage_tcell_subtype_proportion_boxplot_heatmap.R`.
- `11_lavage_tcell_subtype_deg_module_heatmap.R`: Final BALF figures from saved objects and signatures. Original name: `19_revision_lavage_tcell_subtype_deg_module_heatmap.R`.
- `12_lavage_tcell_functional_module_score_compare.R`: Final BALF figures from saved objects and signatures. Original name: `20_revision_lavage_tcell_functional_module_score_compare.R`.
- `13_lavage_tcell_signature_scatter_roc_contribution.R`: Final BALF figures from saved objects and signatures. Original name: `21_revision_lavage_tcell_signature_scatter_roc_contribution.R`.
- `14_lavage_tcell_selected_module_score_umap.R`: Final BALF figures from saved objects and signatures. Original name: `23_revision_lavage_tcell_selected_module_score_umap.R`.
- `15_lavage_demultiplex_hla_tsne_predict_sample_number.R`: Final BALF figures from saved objects and signatures. Original name: `23_revision_lavage_demultiplex_hla_tsne_predict_sample_number.R`.
- `16_lavage_tcell_CAP_LC_DEG_volcano.R`: Final BALF figures from saved objects and signatures. Original name: `22_revision_lavage_tcell_CAP_LC_DEG_volcano.R`.
- `17_lavage_tcell_tcr_clone_ratio_umap_expanded_summary.R`: Final BALF figures from saved objects and signatures. Original name: `22_revision_lavage_tcell_tcr_clone_ratio_umap_expanded_summary.R`.

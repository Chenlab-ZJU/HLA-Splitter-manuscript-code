# CSF immune profiling identifies leptomeningeal metastasis-associated myeloid programs

Manuscript scope: Fig. 2; Supplementary Figs. S4-S5.

Run individual scripts using `Rscript run_analysis.R <relative-script-path>`. Numbers organize the manuscript section; they do not imply that all scripts should be rerun consecutively. Use saved objects and source tables for final figure regeneration. Upstream steps can overwrite objects and require annotation review.

- `01_setup_and_helpers.R`: Upstream cohort processing and reusable helper definitions. Original name: `00_setup_and_helpers.R`.
- `02_data_preprocessing.R`: Upstream cohort processing and reusable helper definitions. Original name: `01_data_preprocessing.R`.
- `03_demultiplexing_and_doublets.R`: Upstream cohort processing and reusable helper definitions. Original name: `02_demultiplexing_and_doublets.R`.
- `04_global_clustering.R`: Upstream cohort processing and reusable helper definitions. Original name: `03_global_clustering.R`.
- `05_myeloid_analysis.R`: Upstream cohort processing and reusable helper definitions. Original name: `04_myeloid_analysis.R`.
- `06_tcell_analysis.R`: Upstream cohort processing and reusable helper definitions. Original name: `05_tcell_analysis.R`.
- `07_scRNA_celltype_composition_plots.R`: Required intermediate source-table / annotation generation. Original name: `08_revision_scRNA_celltype_composition_plots.R`.
- `08_myeloid_reanalysis.R`: Required intermediate source-table / annotation generation. Original name: `09_revision_myeloid_reanalysis.R`.
- `09_myeloid_group_DEG_function.R`: Required intermediate source-table / annotation generation. Original name: `10_revision_myeloid_group_DEG_function.R`.
- `10_tcell_subtype_function_tcr.R`: Required intermediate source-table / annotation generation. Original name: `11_revision_tcell_subtype_function_tcr.R`.
- `11_scRNA_all_celltype_facet.R`: Final CSF figures from saved objects or source tables. Original name: `07_revision_scRNA_all_celltype_facet.R`.
- `12_selected_celltype_compare_width.R`: Final CSF figures from saved objects or source tables. Original name: `08b_revision_selected_celltype_compare_width.R`.
- `13_myeloid_umap_cluster_merge_bc_legend.R`: Final CSF figures from saved objects or source tables. Original name: `09b_revision_myeloid_umap_cluster_merge_bc_legend.R`.
- `14_myeloid_heatmap_label.R`: Final CSF figures from saved objects or source tables. Original name: `10c_revision_myeloid_heatmap_label.R`.
- `15_myeloid_module_score_compare_layout.R`: Final CSF figures from saved objects or source tables. Original name: `10b_revision_myeloid_module_score_compare_layout.R`.
- `16_scRNA_all_major_celltype_heatmap.R`: Final CSF figures from saved objects or source tables. Original name: `16_revision_scRNA_all_major_celltype_heatmap.R`.
- `17_myeloid_key_gene_expression_wilcox.R`: Final CSF figures from saved objects or source tables. Original name: `10d_revision_myeloid_key_gene_expression_wilcox.R`.
- `18_myeloid_key_gene_feature_umap.R`: Final CSF figures from saved objects or source tables. Original name: `10e_revision_myeloid_key_gene_feature_umap.R`.
- `19_redraw_tcell_umap_from_csf_tcell.R`: Final CSF figures from saved objects or source tables. Original name: `11b_revision_redraw_tcell_umap_from_csf_tcell.R`.
- `20_tcell_celltype_proportion_heatmap.R`: Final CSF figures from saved objects or source tables. Original name: `13b_revision_tcell_celltype_proportion_heatmap.R`.
- `21_tcell_celltype_heatmap_align_umap.R`: Final CSF figures from saved objects or source tables. Original name: `17_revision_tcell_celltype_heatmap_align_umap.R`.
- `22_tcell_subtype_proportion_compare.R`: Final CSF figures from saved objects or source tables. Original name: `13_revision_tcell_subtype_proportion_compare.R`.
- `23_tcell_key_gene_violin_compare.R`: Final CSF figures from saved objects or source tables. Original name: `15_revision_tcell_key_gene_violin_compare.R`.
- `24_tcell_functional_module_score_compare.R`: Final CSF figures from saved objects or source tables. Original name: `11c_revision_tcell_functional_module_score_compare.R`.

# Urinary T-cell profiling identifies CXCR3 as a candidate biomarker for TCMR

Manuscript scope: Fig. 5; Supplementary Figs. S9-S10.

Run individual scripts using `Rscript run_analysis.R <relative-script-path>`. Numbers organize the manuscript section; they do not imply that all scripts should be rerun consecutively. Use saved objects and source tables for final figure regeneration. Upstream steps can overwrite objects and require annotation review.

- `01_Tcell_subclustering.R`: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis. Original name: `05_Tcell_subclustering.R`.
- `02_TCR_repertoire_analysis.R`: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis. Original name: `06_TCR_repertoire_analysis.R`.
- `03_Cell_communication.R`: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis. Original name: `07_Cell_communication.R`.
- `04_Tcell_recluster_subtype_frequency.R`: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis. Original name: `19_revision_Tcell_recluster_subtype_frequency.R`.
- `05_ABMR_vs_TCMR_Tcell_DE_function.R`: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis. Original name: `13_revision_ABMR_vs_TCMR_Tcell_DE_function.R`.
- `06_balanced_CXCR3_CXCL9_CXCL10_LR_communication.R`: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis. Original name: `17_revision_balanced_CXCR3_CXCL9_CXCL10_LR_communication.R`.
- `07_tcell_celltype_umap_TCMR_ABMR_facet.R`: Final urine T-cell, CXCR3 and expanded-clone figures. Original name: `31_revision_tcell_celltype_umap_TCMR_ABMR_facet.R`.
- `08_tcell_group_DEG_module_heatmap_score.R`: Final urine T-cell, CXCR3 and expanded-clone figures. Original name: `28_revision_tcell_group_DEG_module_heatmap_score.R`.
- `09_tcell_CXCR3_sample_expression_positive_boxplot.R`: Final urine T-cell, CXCR3 and expanded-clone figures. Original name: `32_revision_tcell_CXCR3_sample_expression_positive_boxplot.R`.
- `10_tcell_CXCR3_positive_ratio_by_subtype_barplot.R`: Final urine T-cell, CXCR3 and expanded-clone figures. Original name: `39_revision_tcell_CXCR3_positive_ratio_by_subtype_barplot.R`.
- `11_tcell_expanded_clone_fraction_by_group_celltype.R`: Final urine T-cell, CXCR3 and expanded-clone figures. Original name: `33_revision_tcell_expanded_clone_fraction_by_group_celltype.R`.
- `12_CXCR3_positive_negative_Tcell_expanded_clone_fraction.R`: Final urine T-cell, CXCR3 and expanded-clone figures. Original name: `34_revision_CXCR3_positive_negative_Tcell_expanded_clone_fraction.R`.
- `13_TCR_clone_ratio_UMAP_histogram_group_fraction.R`: Final urine T-cell, CXCR3 and expanded-clone figures. Original name: `35_revision_TCR_clone_ratio_UMAP_histogram_group_fraction.R`.
- `14_expanded_TCR_clone_diversity_by_group.R`: Final urine T-cell, CXCR3 and expanded-clone figures. Original name: `36_revision_expanded_TCR_clone_diversity_by_group.R`.
- `15_annotate_GSE189536_immune_like_urine.R`: External Harmony annotation and final CXCR3/module validation. Original name: `annotate_GSE189536_immune_like_urine.R`.
- `16_run_GSE189536_revision_rejection_axis_validation.R`: External Harmony annotation and final CXCR3/module validation. Original name: `run_GSE189536_revision_rejection_axis_validation.R`.
- `17_plot_GSE189536_urine_style_Tcell_CXCR3_marker_dotplot.R`: External Harmony annotation and final CXCR3/module validation. Original name: `plot_GSE189536_urine_style_Tcell_CXCR3_marker_dotplot.R`.
- `18_plot_GSE189536_urine_style_Tcell_DEG_module_heatmap_score.R`: External Harmony annotation and final CXCR3/module validation. Original name: `plot_GSE189536_urine_style_Tcell_DEG_module_heatmap_score.R`.

# HLA-Splitter enables accurate demultiplexing across low-input settings

Manuscript scope: Fig. 1; Supplementary Figs. S1-S3.

Run individual scripts using `Rscript run_analysis.R <relative-script-path>`. Numbers organize the manuscript section; they do not imply that all scripts should be rerun consecutively. Use saved objects and source tables for final figure regeneration. Upstream steps can overwrite objects and require annotation review.

- `01_HLA_demultiplex_v2.R`: Shared algorithm / evaluation function. Original name: `HLA_demultiplex_v2.R`.
- `02_Demultiplex_evaluation.R`: Shared algorithm / evaluation function. Original name: `Demultiplex_evaluation.R`.
- `03_scRNA_analysis.R`: Upstream synthetic-data preparation and benchmark computation. Original name: `00_scRNA_analysis.R`.
- `04_Syn_data.R`: Upstream synthetic-data preparation and benchmark computation. Original name: `01_Syn_data.R`.
- `05_ext_data_fig2_method_limits.R`: Final benchmark / HLA-signal figures. Original name: `12_revision_ext_data_fig2_method_limits.R`.
- `06_Demultiplex_compare.R`: Final benchmark / HLA-signal figures. Original name: `03_Demultiplex_compare.R`.
- `07_syn_figures_publication_style.R`: Final benchmark / HLA-signal figures. Original name: `13_revision_syn_figures_publication_style.R`.
- `08_figure_s1e_hla_gene_abundance.R`: Final benchmark / HLA-signal figures. Original name: `14_revision_figure_s1e_hla_gene_abundance.R`.
- `09_Nanopore.R`: Upstream HLA typing, resolution and PBMC demultiplexing. Original name: `Nanopore.R`.
- `10_30PBMC.R`: Upstream HLA typing, resolution and PBMC demultiplexing. Original name: `30PBMC.R`.
- `11_Compare_HLA_resolution.R`: Upstream HLA typing, resolution and PBMC demultiplexing. Original name: `Compare_HLA_resolution.R`.
- `12_nanopore_hla_publication_plots.R`: Final nanopore typing and reference-concordance figures. Original name: `14_revision_nanopore_hla_publication_plots.R`.
- `13_nanopore_hla_locus_coverage_violin_result.R`: Final nanopore typing and reference-concordance figures. Original name: `14b_revision_nanopore_hla_locus_coverage_violin_result.R`.
- `14_HLA_TSNE_ONT_SBT_publication.R`: Final nanopore typing and reference-concordance figures. Original name: `15_revision_HLA_TSNE_ONT_SBT_publication.R`.
- `15_4mix_PBMC.R`: Experimental sorted-cell mixture validation. Original name: `4mix_PBMC.R`.
- `16_10X-hash.R`: Experimental cell-hashing validation. Original name: `10X-hash.R`.

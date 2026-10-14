# HLA-Splitter manuscript analysis code

This code collection follows the Results sections of `HLA-Splitter_manuscript_v3.docx`. It contains the retained analysis workflows, their required helper functions, and final figure-generation scripts. Historical copies and exploratory alternatives are archived outside this directory and should not be uploaded as manuscript code.

## Manuscript order

1. `01_HLA_Splitter_benchmarks`: algorithm, synthetic benchmarking, nanopore typing, sorted-cell mixtures and cell hashing (Fig. 1; Supplementary Figs. S1-S3).
2. `02_CSF_immune_profiling`: CSF composition, myeloid programs and T-cell support analyses (Fig. 2; Supplementary Figs. S4-S5).
3. `03_BALF_Tcell_profiling`: BALF T-cell states, signatures and TCR expansion (Fig. 3; Supplementary Fig. S6).
4. `04_Urine_immune_programs`: urine immune composition, rejection programs and cross-cohort signatures (Fig. 4; Supplementary Figs. S7-S8).
5. `05_Urine_Tcell_CXCR3`: urine T-cell states, CXCR3, expanded clones and GSE189536 evaluation (Fig. 5; Supplementary Figs. S9-S10).
6. `06_Supplementary_tables`: complete signature-gene source metrics.

`00_Shared` contains only transitive shared function dependencies. `manuscript_code_manifest.csv` records original filenames, retained roles and checksums. `manuscript_section_index.csv` follows the Results and Methods headings. `manuscript_figure_code_index.csv` maps the v3 panel order to retained workflows. `figure_panel_code_index.csv` retains the previous panel mapping as provenance; caption inconsistencies are documented separately. `data_file_catalog.csv` lists static data-file references and should be interpreted together with the directory construction in each script.

## Inputs and setup

Use R and install the dependencies listed in `validation/r_package_inventory.csv`. The installed versions in that file describe the validation environment, not a verified locked environment for every historical analysis. Packages from CRAN, Bioconductor and external repositories must be installed from their respective distributions. The scripts do not install packages automatically.

Place the project `Data` directory outside the code repository. Keep its existing cohort layout: `Syn_data`, `30PBMC`, `4mix_PBMC`, the cell-hashing cohort, CSF, BALF, urine and `GSE189536`. Set its location once:

```r
Sys.setenv(HLA_SPLITTER_DATA_ROOT = "/absolute/path/to/Data")
```

For command-line runs, set the same variable in your shell. Optional variables `HLA_SPLITTER_PROJECT_ROOT` and `HLA_SPLITTER_REVISION_ROOT` configure non-cohort table and manuscript output locations. By default the runner assumes the original layout `project/Revison/Code` and `project/Data`. When uploading this directory as a standalone GitHub repository, configure all three roots explicitly.

Analyses require normalized/annotated Seurat objects, metadata, HLA count matrices, V(D)J contig annotations, intermediate result tables and, for upstream reconstruction, sequencing-derived count matrices. Large primary inputs are not duplicated in this code directory. Compact source data remain in `Revison/SourceData`; availability of those tables alone does not imply that every upstream workflow can run.

## Running and checking

From this code directory:

```sh
Rscript validate_code.R
Rscript run_analysis.R --list
Rscript run_analysis.R --self-test
Rscript run_analysis.R --check 01_HLA_Splitter_benchmarks/05_ext_data_fig2_method_limits.R
Rscript run_analysis.R 01_HLA_Splitter_benchmarks/05_ext_data_fig2_method_limits.R
```

Use the exact relative filename from `--list`. Run one script per R process through the runner so configured paths survive the original cleanup calls. Final plotting scripts reuse saved results. Scripts marked upstream or source-table generation are retained to explain how those inputs were produced; they can overwrite intermediate annotations and older plots and should not be rerun solely to change figure styling. Follow each section README and inspect the script's input filenames first. Numbering follows manuscript organization, not an automatic end-to-end execution plan.

Outputs generally retain the original cohort `Figures` and `Result` paths; scripts handling supplementary tables may write under the configured revision root. Vector PDF and raster PNG exports retain existing plotting settings and named color mappings.

## Statistical and scientific interpretation

Patients/samples/biopsies, rather than cells, are the biological replicates for clinical comparisons. Pooled-cell summaries are descriptive unless the manuscript specifies otherwise. Expanded clones use `clone_ratio > 1`; expanded-cell fractions and expanded-clonotype fractions have different denominators. Preserve the code's clonotype metadata and sample identifiers.

Read `manuscript_consistency_review.md` before public release. The v3 manuscript and some retained historical scripts differ in CSF inferential testing and certain panel mappings. These differences are explicitly recorded rather than silently changing scientific methods during packaging. A successful syntax check is not an end-to-end validation of the manuscript results.

## GitHub submission

Upload this code directory after resolving the manuscript-consistency items. Keep protected clinical data, raw reads, large R objects and the history archive outside the public repository. The included GitHub Actions workflow checks R syntax and literal internal source dependencies on Linux; it does not run cohort analyses. Add the authors' approved software license before publication; no license or redistribution rights have been assumed during preparation.

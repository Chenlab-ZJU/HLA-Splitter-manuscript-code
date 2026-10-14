# HLA-Splitter analysis workflows

Analysis code accompanying the article:

**Label-free HLA-guided demultiplexing enables single-cell analysis of scarce clinical samples to reveal disease-specific immune signatures**

R workflows for HLA-guided sample demultiplexing, immune-cell profiling, disease signatures and TCR analysis in benchmark, CSF, BALF, urine and GSE189536 datasets.

## Code structure

| Directory | Main analyses |
| --- | --- |
| `00_Shared` | Shared plotting, enrichment and pseudobulk functions. |
| `01_HLA_Splitter_benchmarks` | HLA assignment, accuracy/recovery/runtime benchmarks, HLA typing and experimental validation. |
| `02_CSF_immune_profiling` | CSF immune composition, myeloid/T-cell annotation and functional programs. |
| `03_BALF_Tcell_profiling` | BALF T-cell states, disease signatures, ROC analysis and TCR expansion. |
| `04_Urine_immune_programs` | Urine immune composition, rejection modules and cross-cohort signatures. |
| `05_Urine_Tcell_CXCR3` | T-cell CXCR3, differential expression, ligand-receptor analysis, TCR expansion/diversity and external validation. |

See each directory's README and script headers for individual workflows.

## Environment

Reference R version: **4.5.3**. Install common dependencies in R:

```r
install.packages(c(
  "Seurat", "harmony", "Matrix", "dplyr", "tidyr", "tibble",
  "readr", "readxl", "stringr", "purrr", "data.table",
  "ggplot2", "ggpubr", "patchwork", "cowplot", "ggrepel",
  "ggsci", "scales", "pheatmap", "RColorBrewer", "circlize",
  "ragg", "pROC", "RANN", "Rtsne", "enrichR", "BiocManager"
), repos = "https://cloud.r-project.org")

BiocManager::install(c(
  "AUCell", "ComplexHeatmap", "clusterProfiler", "GSVA",
  "limma", "MAST", "scRepertoire"
), ask = FALSE, update = FALSE)
```

Install additional packages listed by the selected script; CellChat installation is described [here](https://github.com/jinworks/CellChat). Use a UTF-8-capable environment and Arial/Cairo for the corresponding figure exports.

## Data and paths

Keep cohort folders under `Data/`: `Syn_data`, `30PBMC`, `4mix_PBMC`, `10X-dataset`, `脑脊液样本`, `肺泡灌洗液`, `尿液样本` and `GSE189536`. Preserve their existing input/result hierarchy.

Inputs include HLA genotypes/counts, clinical metadata, annotated Seurat objects, VDJ annotations and summary tables. Required files are specified in each script and are not included in the code directory.

Set these variables in your shell environment or `.Renviron`, replacing the example paths:

```text
HLA_SPLITTER_PROJECT_ROOT=/path/to/project
HLA_SPLITTER_DATA_ROOT=/path/to/project/Data
HLA_SPLITTER_REVISION_ROOT=/path/to/project/Revison
```

Without explicit configuration, the runner assumes `project/Data` and `project/Revison/Code`.

## Running

From the code directory:

```sh
Rscript run_analysis.R --list
Rscript run_analysis.R --self-test
Rscript run_analysis.R --check 03_BALF_Tcell_profiling/09_lavage_tcell_group_facet_umap.R
Rscript run_analysis.R 03_BALF_Tcell_profiling/09_lavage_tcell_group_facet_umap.R
```

Select an exact script path from `--list`. Run scripts individually and provide their upstream inputs; numbering is not an automatic execution order. Figure scripts can reuse saved objects/tables. Preflight checks verify function loading or syntax/directories, not complete analysis execution.

Outputs include PDF/PNG figures and result tables, usually under the cohort's `Figures` and `Result` directories. Some workflows use nested `Tables` or revision-level output directories.

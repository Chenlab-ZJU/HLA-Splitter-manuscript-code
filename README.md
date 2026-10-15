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

The compact processed-data bundle is stored beside the code:

```text
Revison/
|-- Code/
`-- Data/
    |-- synthetic/                    # Benchmark source tables and HLA metadata
    |-- pbmc30/                       # HLA typing, coverage and embedding source
    |-- pbmc4/                        # Annotated mixture object and HLA metadata
    |-- hash10x/                      # Annotated hashing object and HLA metadata
    |-- csf/                          # Three annotated objects, metadata and source tables
    |-- balf/                         # T-cell object, metadata, signatures and embedding source
    |-- urine/                        # Three annotated objects, metadata and signature tables
    |-- gse189536/                    # Immune_filtered, Immune_annotated and Tcell_harmony
    |-- references/                   # GO biological-process gene sets
    `-- data_manifest.csv
```

Each dataset has one flat directory containing its required input files. Objects retain expression data, annotations, embeddings and embedded metadata; separate sample/HLA metadata are also included. `data_manifest.csv` lists every input and its SHA256 checksum. Original-to-current paths are recorded in `validation/data_path_migration.csv`. Raw reads, large workspace backups and general intermediate files are excluded. Upstream preprocessing and signature optimization require additional original inputs.

With `Code/` and `Data/` as sibling directories, the runner automatically uses this local data bundle. To use a relocated copy of the same compact layout, set the following variable in your shell environment or `.Renviron`:

```text
HLA_SPLITTER_DATA_ROOT=/path/to/Revison/Data
```

## Running

From the code directory:

```sh
Rscript run_analysis.R --list
Rscript validate_data.R
Rscript run_analysis.R --self-test
Rscript run_analysis.R --check 03_BALF_Tcell_profiling/09_lavage_tcell_group_facet_umap.R
Rscript run_analysis.R 03_BALF_Tcell_profiling/09_lavage_tcell_group_facet_umap.R
```

Select an exact script path from `--list`. Run scripts individually and provide their upstream inputs; numbering is not an automatic execution order. Figure scripts can reuse saved objects/tables. Preflight checks verify function loading or syntax/directories, not complete analysis execution.

Outputs include PDF/PNG figures and result tables, usually under the cohort's `Figures` and `Result` directories. Some workflows use nested `Tables` or revision-level output directories.

`validation/compact_layout_runtime_check.csv` records runtime checks using the English input layout; `runtime_check.csv` retains the earlier validation. `workflow_coverage.csv` distinguishes downstream workflows from function/setup components and upstream scripts requiring additional inputs. `validate_data.R` verifies the bundled files and R syntax; add `--objects` to inspect single-cell object structure.

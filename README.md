# HLA-Splitter analysis workflows

R workflows for HLA-guided single-cell sample demultiplexing and downstream immune profiling. The analyses cover synthetic and experimental benchmarks, cerebrospinal fluid (CSF), bronchoalveolar lavage fluid (BALF), urine, and the external kidney-biopsy dataset GSE189536.

The code includes sample-assignment functions, expression-data processing, cell-type annotation, sample-level comparisons, functional-module and disease-signature analyses, T-cell receptor (TCR) analyses, and figure generation.

## Code structure

```text
Code/
|-- README.md
|-- run_analysis.R
|-- 00_Shared/
|-- 01_HLA_Splitter_benchmarks/
|-- 02_CSF_immune_profiling/
|-- 03_BALF_Tcell_profiling/
|-- 04_Urine_immune_programs/
`-- 05_Urine_Tcell_CXCR3/
```

| Directory | Functions |
| --- | --- |
| `00_Shared` | Shared plotting, enrichment, pseudobulk and matrix-analysis functions used by the cohort workflows. |
| `01_HLA_Splitter_benchmarks` | HLA-based donor assignment, synthetic-data construction, accuracy/recovery/runtime comparisons, HLA-expression and typing-resolution analyses, nanopore HLA typing, and validation using sorted-cell mixtures and cell hashing. |
| `02_CSF_immune_profiling` | CSF quality control and demultiplexing, immune-cell composition, myeloid and T-cell analysis, marker expression, functional-module scoring, heatmaps and UMAP visualization across metastases, hemorrhage and encephalitis groups. |
| `03_BALF_Tcell_profiling` | BALF T-cell processing and annotation, subtype composition, functional programs, CAP/LC disease signatures, ROC analysis and TCR clonal expansion across CAP, LC and ILD. |
| `04_Urine_immune_programs` | Urine quality control and demultiplexing, major immune-cell composition, immune-cell ratios, rejection-associated programs, HLA/module scores, and disease-signature construction and evaluation across urine and GSE189536. |
| `05_Urine_Tcell_CXCR3` | Urine T-cell annotation, ABMR-versus-TCMR differential expression, CXCR3 expression and positivity, composition-balanced CXCR3-CXCL9/CXCL10 analyses, TCR expansion/diversity and external T-cell validation in GSE189536. |

Each analysis directory has a README listing its scripts. Script headers describe their purpose, data cohort and related analyses. Numbering organizes the code by analysis topic; it does not define a single automatic pipeline.

## Environment

### R and packages

The current local environment uses R 4.5.3. Selected installed versions are listed as an environment reference. This table is not a complete environment lock or a statement that every workflow has been executed with these versions.

| Package | Installed version |
| --- | --- |
| Seurat | 5.5.0 |
| SeuratObject | 5.4.0 |
| harmony | 2.0.2 |
| ggplot2 | 4.0.3 |
| AUCell | 1.32.0 |
| scRepertoire | 2.6.2 |
| CellChat | 2.2.0.9001 |

Install common data-processing and plotting dependencies in an R session:

```r
install.packages(c(
  "Seurat", "harmony", "Matrix", "dplyr", "tidyr", "tibble",
  "readr", "readxl", "stringr", "purrr", "data.table",
  "ggplot2", "ggpubr", "patchwork", "cowplot", "ggrepel",
  "ggsci", "scales", "pheatmap", "RColorBrewer", "circlize",
  "ragg", "pROC", "RANN", "Rtsne", "enrichR", "BiocManager"
), repos = "https://cloud.r-project.org")
```

Use Bioconductor for packages distributed there:

```r
BiocManager::install(c(
  "AUCell", "ComplexHeatmap", "clusterProfiler", "GSVA",
  "limma", "MAST", "scRepertoire"
), ask = FALSE, update = FALSE)
```

These commands provide a starting environment. Individual workflows may additionally require packages such as `scrubletR`, `SCIPAC`, `immGLIPH`, `turboGliph`, `immunarch`, `presto` or organism annotation packages. Check the package-loading statements in the selected script and its sourced setup/helpers. Install workflow-specific packages using their maintainers' instructions. See the official [Seurat installation guide](https://satijalab.org/seurat/articles/install_v5.html) and [Bioconductor installation guide](https://bioconductor.org/install/).

CellChat workflows require CellChat and its dependencies; follow the [CellChat installation instructions](https://github.com/jinworks/CellChat). Installing current packages does not recreate the development package versions listed above.

RStudio is optional. Use a UTF-8-capable environment for Chinese input directory names and source annotations. PDF export uses Cairo where specified. Install Arial or apply a consistent font substitution if Arial is unavailable.

### External tools

The R workflows consume sequencing-derived inputs. Reconstructing these inputs from raw reads may require Cell Ranger for expression/VDJ processing, kallisto-bustools for HLA allele counts, and the HLA-typing or demultiplexing tools used by the selected benchmark. These external workflows are not automatically installed or invoked by `run_analysis.R`.

## Input data

Keep code and data in separate directories. The default layout is:

```text
project/
|-- Data/
|   |-- Syn_data/
|   |-- 30PBMC/
|   |-- 4mix_PBMC/
|   |-- 10X-dataset/
|   |-- 脑脊液样本/
|   |-- 肺泡灌洗液/
|   |-- 尿液样本/
|   `-- GSE189536/
`-- Revison/
    `-- Code/
```

Within each cohort, preserve the `Data`, `Result` and `Figures` hierarchy referenced by the scripts. Cohorts do not all use the same subdirectories. Example saved-object inputs are:

| Analysis | Input relative to `Data/` |
| --- | --- |
| CSF immune composition | `脑脊液样本/Result/CSF_scRNA.rds` |
| CSF T-cell figures | `脑脊液样本/Result/CSF_Tcell.rds` |
| BALF T-cell figures | `肺泡灌洗液/Result/lavage_Tcell.rds` |
| Urine immune-cell module figures | `尿液样本/Result/Objects/Urine_ALL.rds` |
| Urine T-cell figures | `尿液样本/Result/Objects/Urine_Tcell.rds` |
| GSE189536 annotation | `GSE189536/Harmony_immune_cells_manual_annotation_results/GSE189536_manual_annotation_immune_filtered_seurat.rds` |

Additional inputs depend on the analysis: HLA genotypes, cell-by-allele count matrices, barcodes, clinical metadata, annotated Seurat objects, VDJ contig annotations and saved summary tables. Required filenames are defined in the script's input-loading statements. These large biological datasets are not included in the code directory.

Saved Seurat objects must retain the expected metadata, assays and reductions. Plotting scripts that read an existing UMAP reuse its coordinates; clustering and annotation scripts have separate input requirements.

## Path configuration

| Environment variable | Meaning | Default |
| --- | --- | --- |
| `HLA_SPLITTER_DATA_ROOT` | Directory containing the cohort data folders. | `<project>/Data` |
| `HLA_SPLITTER_PROJECT_ROOT` | Project directory used for project-level paths. | Parent of the revision directory. |
| `HLA_SPLITTER_REVISION_ROOT` | Directory for revision-level output paths used by some scripts. | Parent of the code directory. |

Set explicit absolute paths when using another layout. The runner configures the cohort-specific variables and persistent R path options used by the analysis scripts.

Windows PowerShell:

```powershell
$env:HLA_SPLITTER_PROJECT_ROOT = "D:/HLA-Splitter-project"
$env:HLA_SPLITTER_DATA_ROOT = "D:/HLA-Splitter-project/Data"
$env:HLA_SPLITTER_REVISION_ROOT = "D:/HLA-Splitter-project/Revison"
Set-Location "D:/HLA-Splitter-project/Revison/Code"
```

Linux/macOS:

```bash
export HLA_SPLITTER_PROJECT_ROOT="/path/to/project"
export HLA_SPLITTER_DATA_ROOT="/path/to/project/Data"
export HLA_SPLITTER_REVISION_ROOT="/path/to/project/Revison"
cd /path/to/project/Revison/Code
```

Replace the example paths with the locations on your computer. Cohort folder names must match the code. When setting paths in R with `Sys.setenv()`, launch the runner from the same R session so the child process inherits the settings.

## Running analyses

### List available scripts

From the code directory:

```sh
Rscript run_analysis.R --list
```

The runner lists available relative filenames and their purposes using the script metadata. Select a script using its exact listed path.

### Check the entry point and configured paths

```sh
Rscript run_analysis.R --self-test
Rscript run_analysis.R --check 03_BALF_Tcell_profiling/09_lavage_tcell_group_facet_umap.R
```

`--self-test` checks that the retained HLA-Splitter function can be loaded. `--check` confirms that the selected script parses and that the configured data root and cohort directory exist. These checks do not load biological data, verify every dynamic input, check all package dependencies or execute the analysis.

### Run one analysis

Regenerate the BALF T-cell group-faceted UMAP from the saved `lavage_Tcell.rds` object:

```sh
Rscript run_analysis.R 03_BALF_Tcell_profiling/09_lavage_tcell_group_facet_umap.R
```

Other examples:

```sh
# Synthetic demultiplexing performance plots
Rscript run_analysis.R 01_HLA_Splitter_benchmarks/06_Demultiplex_compare.R

# CSF selected cell-type frequency comparisons
Rscript run_analysis.R 02_CSF_immune_profiling/12_selected_celltype_compare_width.R

# Urine immune-cell rejection modules
Rscript run_analysis.R 04_Urine_immune_programs/11_all_cells_rejection_module_heatmap_score.R

# Urine T-cell UMAP and subtype-frequency summaries
Rscript run_analysis.R 05_Urine_Tcell_CXCR3/07_tcell_celltype_umap_TCMR_ABMR_facet.R
```

Run one analysis per R process through `run_analysis.R`. Some scripts clear their local workspace; the runner's persistent path options preserve the configured directories.

For a new dataset, prepare count matrices and sample metadata, perform quality control/demultiplexing, create annotated expression objects, generate analysis summaries, and run the relevant figure scripts. Some workflows require marker-based manual annotation before downstream analysis. Consult the directory README and input-loading statements for step dependencies.

For figure regeneration from existing results, use the relevant saved-object or source-table plotting script. Upstream processing can overwrite intermediate objects and annotations. Do not execute every numbered script consecutively as an assumed pipeline.

## Outputs and interpretation

Outputs are written primarily to the cohort's `Figures` and `Result` directories. GSE189536 workflows can use nested `Figures` and `Tables` directories under their annotation results. Some project-level outputs use the configured revision root. Each script specifies its output filenames and export settings.

Typical outputs include PDF/PNG figures, cell-type frequencies, differential-expression tables, module scores, disease-signature scores and ROC results, ligand-receptor summaries, and TCR expansion/diversity summaries. Most figure scripts export vector PDF and high-resolution PNG files.

Clinical comparisons commonly use sample-, patient- or biopsy-level summaries. Distinguish pooled-cell summaries from sample-level inference. The selected script defines the statistical test, multiple-testing correction and displayed significance threshold.

CXCR3 positivity is defined by expression greater than zero where specified. Expanded T cells use `clone_ratio > 1`; expanded-cell and expanded-clonotype fractions have different denominators. Preserve sample identifiers and clonotype metadata during TCR analysis.

## Troubleshooting

| Message or symptom | Action |
| --- | --- |
| `Rscript` is not recognized | Add the R `bin` directory to the executable search path or use the absolute path to `Rscript`. |
| Missing data root or cohort directory | Check the path variables and exact cohort folder names. |
| Missing package | Install the reported package using its CRAN, Bioconductor or maintainer instructions. Setup/helper scripts may require additional packages. |
| Missing RDS or source table | Check the input-loading statement and provide or generate its upstream object/table. |
| Missing metadata column or reduction | Use the matching annotated object and confirm its metadata, assay and UMAP/PCA/Harmony reductions. |
| Different fonts or PDF export failure | Check Arial availability and Cairo support in the R installation. |
| Memory allocation error | Use sufficient memory for the selected Seurat object and matrix operations; requirements depend on dataset size. |

Record the analysis environment with `sessionInfo()` and retain the input versions, sample metadata and script parameters with the outputs.

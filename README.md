# AorticStenOmics

Analysis code for the AorticStenOmics study of aortic stenosis, including bulk RNA-seq, single-nucleus RNA-seq, single-cell RNA-seq, Xenium spatial transcriptomics, imaging-linked statistics, manuscript figure generation, and the companion Shiny app.

## Repository Structure

``` text
00_study_design/          Study metadata and sequencing design files
01_bulk_seq/              Bulk RNA-seq analysis notebooks
02_UBC_snRNA_seq/         UBC single-nucleus RNA-seq analysis notebooks
03_SALTIRE3_scRNA_seq/    SALTIRE3 single-cell RNA-seq analysis notebooks
04_SALTIRE3_Xenium/       Xenium spatial transcriptomics analysis notebooks
05_statistics/            mIF, FAPI PET-CT, and follow-up statistical analyses
06_data_deposit/          Object preparation for data deposition
07_shiny_app/             Shiny app and app data-preparation script
99_figures/               Manuscript figure and supplemental table scripts
utils.R                   Shared functions, plotting helpers, and colour palettes
```

Each analysis module follows the same broad layout:

``` text
01_data/       Local input data
02_notebooks/  Source notebooks and scripts
03_output/     Generated analysis outputs; not tracked by Git
```

## Data Availability

Large input data are not included in this repository. Raw input files are available on The University of Edinburgh's DataShare (LINK provided upon publication). After downloading, place the data in the following directories:

- `01_bulk_seq/01_data/featureCounts_output/`
- `02_UBC_snRNA_seq/01_data/cellbender_output/`
- `03_SALTIRE3_scRNA_seq/01_data/cellbender_output/`
- `04_SALTIRE3_Xenium/01_data/20251029__151959__MBSV20251027/`

The notebooks expect these files to be available in the local directory structure shown above. To rerun the analyses, download the deposited data and regenerate the required intermediate files, then place them in the corresponding `01_data/` or `03_output/` directories.

## Environment

The R package environment is managed with `renv`. After cloning the repository, restore the R environment from `renv.lock`:

``` r
renv::restore()
```

CellBender and souporcell outputs are non-deterministic; pre-computed outputs are provided and should not be rerun to reproduce results. Analyses outside of R are not strictly environment-tracked and their outputs have been included for reproducibility.

## Analysis Workflow

The main analysis notebooks are numbered in execution order within each module. In general, run the notebooks in this order:

1.  Study design and sequencing metadata preparation.
2.  Bulk RNA-seq processing and downstream pathway analysis.
3.  UBC snRNA-seq preprocessing, annotation, VIC subclustering, and trajectory analyses.
4.  SALTIRE3 scRNA-seq preprocessing, annotation, CellPhoneDB, and Drug2Cell analyses.
5.  Xenium spatial preprocessing, annotation transfer, histology alignment, and scNiche analysis.
6.  Statistical analyses for mIF, FAPI PET/CT, and follow-up data.
7.  Data-deposition object preparation.
8.  Manuscript figure and supplemental table generation.

Many notebooks include checksum blocks to validate expected inputs and outputs.

## Shiny App

The Shiny app is in `07_shiny_app/`. It uses curated, preprocessed objects rather than loading the full analysis objects at runtime.

Prepare app data with:

``` r
old_wd <- setwd("07_shiny_app")
source("prepare_data.R")
setwd(old_wd)
```

Then run the app with:

``` r
shiny::runApp("07_shiny_app")
```

## Figures and Supplemental Tables

All manuscript figures and supplemental tables are generated directly from individual scripts in `99_figures/02_notebooks/`.

## License

The content of this project itself is licensed under the [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/), and the underlying source code used to format and display that content is licensed under the [MIT license](LICENSE).
 
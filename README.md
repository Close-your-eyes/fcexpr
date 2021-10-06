
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fcexpr

<!-- badges: start -->
<!-- badges: end -->

A set functions to organize and analyze flow cytometry experiments.

## Installation

Please install the necessary packages from bioconductor first:

``` r
install.packages("BiocManager")
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("CytoML")
```

Please install the development version of fcexpr from GitHub. This
requires devtools to be insalled first.

``` r
install.packages("devtools")
devtools::install_github("Close-your-eyes/fcexpr")
```

## Vignettes

Please see vignette(s) for workflows:  
[synchronize FCS files and
sampledescription](https://close-your-eyes.github.io/fcexpr/articles/synchronizing_FCS_files_with_an_xlsx_file.html)  
[import data from flowjo
workspaces](https://close-your-eyes.github.io/fcexpr/articles/import_data_from_fj_workspaces.html)

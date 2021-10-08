
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fcexpr

<!-- badges: start -->
<!-- badges: end -->

A set of functions to organize and analyze flow cytometry experiments.

## Installation

Please install the necessary packages from bioconductor first:

``` r
install.packages("BiocManager")
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("CytoML")
```

Please install the development version of fcexpr from GitHub. This
requires devtools to be installed first.

``` r
install.packages("devtools")
devtools::install_github("Close-your-eyes/fcexpr")
```

If you want you can start by creating a default folder for your
experiment:

``` r
fcecpr::new_exp(path = 'your system directory', name = 'your folder name')
```

Please see vignette(s) for further workflows:

## Vignettes (tutorials)

[synchronize FCS files and
sampledescription](https://close-your-eyes.github.io/fcexpr/articles/synchronizing_FCS_files_with_an_xlsx_file.html)  
[import data from flowjo
workspaces](https://close-your-eyes.github.io/fcexpr/articles/import_data_from_fj_workspaces.html)

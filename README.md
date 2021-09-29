
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fcexpr

<!-- badges: start -->
<!-- badges: end -->

A set of beginner-friendly functions to organize flow cytometry
experiments.

## Installation

Please install the necessary packages from bioconductor first:

``` r
install.packages("BiocManager")
BiocManager::install(version = "3.13")
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

A release on CRAN is not intended.

## Example

``` r
library(fcexpr)
## To create a new experiment template folder on your disk run
new_exp(path = 'full path to the parent directory', name = 'name of the folder (e.g. CD3_titration)', date_prefix = T)
## as date_prefix is set to TRUE (or short 'T' only) by default the folder will be prefixed by the current date 
```

## Vignettes

Please see vignette(s) for other workflows:

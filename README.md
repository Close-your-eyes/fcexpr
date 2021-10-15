
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fcexpr

<!-- badges: start -->
<!-- badges: end -->

A set of functions to organize and analyze flow cytometry experiments in
R and to overcome a few limitations from FlowJo.

## Installation

Please install optional packages from bioconductor. Some functions can
be used without them though.

``` r
install.packages("BiocManager")
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("CytoML")
```

Please install fcexpr from GitHub. This requires devtools to be
installed first.

``` r
install.packages("devtools")
devtools::install_github("Close-your-eyes/fcexpr")
```

## The idea

FlowJo is, to my mind, still the best tool to explore and gate flow
cytometric data. I am not aware of a purely R-based procedure which is
as convenient and elaborated and which hence could replace this
commercial software. Other commercial software may also be convenient
though. When it comes to statistical analysis of gated data scientists
usually export aggregated values from FlowJo (MFIs or cell counts) to
other software (see workflow (i) below). Not only that the subsequent
software is commercial as well, it is also  
1) often inconvenient to manually copy and arrange the data  
2) hence inefficient with respect to labor time  
3) problematic with respect to reproducibility (unambiguous links
between meta-data and the phenotypic data may be lost in the process of
clicking, dragging and dropping)  
4) not scalable to many samples and to problems where the exact read-out
parameter is unknown yet.  
So, I want to make a plea for workflow (ii) and provide tools to easily
accomplish it. Once you have a data.frame with aggregated values from
your cytometric experiment you can apply common procedures to tidy and
plot the data with tools from the
[tidyverse](https://www.tidyverse.org).

![alt
text](https://github.com/Close-your-eyes/fcexpr/blob/main/inst/extdata/workflows.png)

## How to start

You may adopt the whole concept of organizing your flow cytrometric
experiment by creating a folder as follows. Then paste your fcs files
and see the vignette below of how to synchronize and document them with
a sampledescription file.

``` r
fcecpr::new_exp(path = 'your system directory', name = 'your folder name')
```

Alternatively, you may only be interested in pulling out the aggregated
data from your flowjo workspace (cell counts and/or MFIs), then use:

``` r
fcexpr::wsx_get_popstats(ws)
```

## Vignettes (tutorials)

[synchronize FCS files and
sampledescription](https://close-your-eyes.github.io/fcexpr/articles/synchronizing_FCS_files_with_an_xlsx_file.html)  
[import data from flowjo
workspaces](https://close-your-eyes.github.io/fcexpr/articles/import_data_from_fj_workspaces.html)

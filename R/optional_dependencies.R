.fcexpr_bioconductor_packages <- c(
  "BiocGenerics",
  "CytoML",
  "flowCore",
  "flowWorkspace",
  "ggcyto"
)

.fcexpr_pak_packages <- c(
  brathering = "close-your-eyes/brathering",
  colrr = "close-your-eyes/colrr",
  EmbedSOM = "exaexa/EmbedSOM",
  presto = "immunogenomics/presto",
  scattermore = "exaexa/scattermore"
)

.fcexpr_cran_packages <- c(
  "BiocManager",
  "cowplot",
  "digest",
  "diptest",
  "ggplot2",
  "ggpointdensity",
  "ggrepel",
  "ggtext",
  "ggraph",
  "Gmisc",
  "igraph",
  "knitr",
  "lubridate",
  "magick",
  "matrixStats",
  "matrixTests",
  "mclust",
  "netstat",
  "officer",
  "openxlsx",
  "pak",
  "parallel",
  "parallelDist",
  "rmarkdown",
  "RSelenium",
  "rstudioapi",
  "Rtsne",
  "scales",
  "Seurat",
  "stringdist",
  "uwot",
  "waldo"
)

.fcexpr_optional_packages <- c(
  .fcexpr_bioconductor_packages,
  names(.fcexpr_pak_packages),
  .fcexpr_cran_packages
)

.ensure_package <- function(package) {
  if (!package %in% .fcexpr_optional_packages) {
    stop("Unknown optional package: ", package, call. = FALSE)
  }
  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible(TRUE))
  }

  message("Installing optional package '", package, "'.")
  if (package %in% .fcexpr_bioconductor_packages) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager")
    }
    BiocManager::install(package, ask = FALSE, update = FALSE)
  } else if (package %in% names(.fcexpr_pak_packages)) {
    if (!requireNamespace("pak", quietly = TRUE)) {
      utils::install.packages("pak")
    }
    pak::pak(unname(.fcexpr_pak_packages[[package]]))
  } else {
    utils::install.packages(package)
  }

  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package '", package, "' is required for this function but could not be installed.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.ensure_packages <- function(packages) {
  invisible(lapply(unique(packages), .ensure_package))
}

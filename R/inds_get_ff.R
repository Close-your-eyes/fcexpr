#' Title
#'
#' @param ind_mat
#' @param population
#' @param alias_attr_name
#' @param path_attr_name
#' @param downsample
#' @param inverse_transform
#' @param lapply_fun
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
inds_get_ff <- function(ind_mat,
                        population,
                        alias_attr_name = "short_names",
                        path_attr_name = "FilePath",
                        downsample = 1,
                        inverse_transform = c(T,F),
                        lapply_fun = lapply,
                        ...) {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }
  lapply_fun <- match.fun(lapply_fun)
  check_in(wsp = "wsp", samples = NULL, groups = NULL, FCS.file.folder = NULL, inverse_transform = inverse_transform)
  ff.list <- lapply_fun(ind_mat,
                        get_ff2,
                        downsample = downsample,
                        population = population,
                        inverse_transform = inverse_transform)

  ind_mat <- ind_mat[which(!sapply(ff.list, is.null))]
  ff.list <- ff.list[which(!sapply(ff.list, is.null))]

  names(ff.list) <- unname(sapply(ind_mat, function(x) basename(attr(x, path_attr_name))))
  ffs <- lapply(seq_along(ff.list[[1]]), function(x) {
    sapply(ff.list, "[", x, simplify = T)
  })
  names(ffs) <- stats::setNames(c("inverse", "logicle"), c(T,F))[as.character(inverse_transform)]

  return(ffs)
}

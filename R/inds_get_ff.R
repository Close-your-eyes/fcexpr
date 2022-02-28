#' Get (subsetted) flowFrames from FCS files
#'
#' @param ind_mat a list of indices matrices, preferentially generated with fcexpr::wsp_get_indices
#' @param population which population (=node, =gate) to subset flowFrames on; must be a column name of ind_mat or an alias stored in alias_attr_name
#' @param alias_attr_name the name of attribute containing aliases (shortest unique names) of node-names (gating paths)
#' @param path_attr_name the name of attribute containing the path (url) to the fcs file on which to apply the subsetting
#' @param downsample numeric, if < 0 then a fraction of events is sampled, if > 0 an absolute number of events is sampled
#' @param inverse_transform return inverse- (T) or logicle- (F) transform or both (c(T,F))
#' @param lapply_fun lapply function name, unquoted; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen as lapply_fun
#'
#' @return list of flow frames, one for each ind_mat
#' @export
#'
#' @examples
#' \dontrun{
#' ind_mat <- fcexpr::wsp_get_indices("mypath/my.wsp")
#' ff <- inds_get_ff(ind_mat = ind_mat, population = "CD8+")
#' }
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
  if (length(population) > 1) {
    stop("Only provide one population.")
  }
  lapply_fun <- match.fun(lapply_fun)
  check_in(wsp = "wsp", samples = NULL, groups = NULL, FCS.file.folder = NULL, inverse_transform = inverse_transform)
  ff.list <- lapply_fun(ind_mat,
                        get_ff2,
                        downsample = downsample,
                        population = population,
                        inverse_transform = inverse_transform,
                        alias_attr_name = alias_attr_name,
                        path_attr_name = path_attr_name,
                        ...)

  ind_mat <- ind_mat[which(!sapply(ff.list, is.null))]
  ff.list <- ff.list[which(!sapply(ff.list, is.null))]

  # maybe not necessary
  names(ff.list) <- unname(sapply(ind_mat, function(x) basename(attr(x, path_attr_name))))

  ffs <- lapply(seq_along(inverse_transform), function(x) {
    sapply(ff.list, "[", x, simplify = T)
  })
  names(ffs) <- stats::setNames(c("inverse", "logicle"), c(T,F))[as.character(inverse_transform)]

  return(ffs)
}

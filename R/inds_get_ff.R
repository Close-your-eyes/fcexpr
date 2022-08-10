#' Get (subsetted) flowFrames from FCS files
#'
#' From a matrix of indices generated with fcexpr::wsp_get_indices flowframes are generated.
#' Only those events within the selected population are included in flowframes. By default
#' no compensation will be applied on fluorescence intensites. This can be done afterwards
#' though with in appropriate compensation matrix.
#'
#' @param ind_mat a list of indices matrices, preferentially generated with fcexpr::wsp_get_indices
#' @param population which population (=node, =gate) to subset flowFrames on; must be a column name of ind_mat or an alias stored in alias_attr_name
#' @param alias_attr_name the name of attribute containing aliases (shortest unique names) of node-names (gating paths)
#' @param path_attr_name the name of attribute containing the path (url) to the fcs file on which to apply the subsetting
#' @param downsample numeric, if < 0 then a fraction of events is sampled, if > 0 an absolute number of events is sampled; or set to "min"
#' which will lead to downsampling each flowframe to the number of events in the flowframe with lowest number of events; can be a single value to treat all
#' FCS files equally or can be a vector of same length as FCS files
#' @param return_untransformed logical; do return untransformed (inverse) data
#' @param return_logicle_transformed logical; do return logicle-transformed data
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
                        return_untransformed = T,
                        return_logicle_transformed = T,
                        lapply_fun = lapply,
                        ...) {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }
  if (missing(population)) {
    stop("Plesae provide a population to get flowframes for. To get all events, set population = 'root'.")
  }
  if (length(population) > 1) {
    stop("Only provide one population.")
  }

  lapply_fun <- match.fun(lapply_fun)

  if (is.numeric(downsample)) {
    ds <- downsample
  } else if (all(downsample == "min")) {
    ds <- 1
  } else {
    stop("downsample has to be numeric of 'min'. With min all flowframes will be downsampled to that flowframe with the lowest number of events.")
  }
  # check length of downsample equal to length of ind_mat or equal to 1
  if (length(ds) != 1 && length(ds) != length(ind_mat)) {
    stop("downsample has to have length 1 or length of ind_mat (one value for each FCS file).")
  }

  if (length(ds) != 1) {
    for (x in seq_along(ind_mat)) {
      attr(ind_mat[[x]], "downsample") <- ds[x]
    }
  }

  check_in(wsp = "wsp", samples = NULL, groups = NULL, FCS.file.folder = NULL, return_untransformed = return_untransformed,
           return_logicle_transformed = return_logicle_transformed)

  ## loop over ind_mat_indices = loop over fcs files
  ff.list <- lapply_fun(ind_mat,
                        get_ff2,
                        downsample = ds,
                        population = population,
                        return_untransformed = return_untransformed,
                        return_logicle_transformed = return_logicle_transformed,
                        alias_attr_name = alias_attr_name,
                        path_attr_name = path_attr_name,
                        ...)

  if (all(downsample == "min")) {
    min <- min(unlist(lapply(sapply(sapply(ff.list, "[", 1), "[", 1), nrow)))
    for (x in seq_along(ff.list)) {
      inds <- c(rep(T, min), rep(F, nrow(ff.list[[x]][[1]][[1]])-min))
      ff.list[[x]][[1]][[1]] <- subset(ff.list[[x]][[1]][[1]], sample(inds))
      ff.list[[x]][[1]][[2]] <- subset(ff.list[[x]][[1]][[2]], sample(inds))
    }
  }

  ind_mat <- ind_mat[which(!sapply(ff.list, is.null))]
  ff.list <- ff.list[which(!sapply(ff.list, is.null))]

  # maybe not necessary
  names(ff.list) <- unname(sapply(ind_mat, function(x) basename(attr(x, path_attr_name))))

  ff.list <- stats::setNames(lapply(names(ff.list[[1]]), function(x) {
    sapply(ff.list, "[", x, simplify = T)
  }), nm = names(ff.list[[1]]))

  return(ff.list)
}

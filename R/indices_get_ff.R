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
#' @param lapply_fun lapply function name, unquoted; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen as lapply_fun
#' @param seed set a seed to reproduce downsampling
#' @param channels channels to use for leverage score calculation; use wsx_get_keywords to retrieve channel names/descriptions
#' @param leverage_score_for_sampling logical whether to use leverage scores for downsampling
#'
#' @return list of flow frames, one for each ind_mat
#' @export
#'
#' @examples
#' \dontrun{
#' ind_mat <- fcexpr::wsp_get_indices("mypath/my.wsp")
#' ff <- indices_get_ff(ind_mat = ind_mat, population = "CD8+")
#' }
indices_get_ff <- function(ind_mat,
                           population,
                           alias_attr_name = "short_names",
                           path_attr_name = "FilePath",
                           downsample = 1,
                           lapply_fun = lapply,
                           seed = 42,
                           channels = NULL,
                           leverage_score_for_sampling = F,
                           ...) {

  ## check and update


  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }
  if (missing(population)) {
    stop("Plesae provide a population to get flowframes for. To get all events, set population = 'root'.")
  }

  if (is.matrix(ind_mat)) {
    ind_mat <- list(ind_mat)
  }
  lapply_fun <- match.fun(lapply_fun)


  if (is.numeric(downsample)) {
    ds <- downsample
  } else if (all(downsample == "min")) {
    ds <- 1
  } else {
    stop("downsample has to be numeric or 'min'. With min, all flowframes will be downsampled to that flowframe with the lowest number of events.")
  }
  if (length(ds) != 1 && length(ds) != length(ind_mat)) {
    stop("downsample has to have length 1 or length of ind_mats.")
  } else if (length(ds) == 1) {
    ds <- rep(ds, length(ind_mat))
  }

  for (x in seq_along(ind_mat)) {
    attr(ind_mat[[x]], "downsample") <- ds[x]
  }

  ## loop over ind_mat_indices = loop over fcs files
  ff.list <- lapply_fun(ind_mat,
                        get_ff2,
                        population = population,
                        alias_attr_name = alias_attr_name,
                        path_attr_name = path_attr_name,
                        seed = seed,
                        channels = channels,
                        leverage_score_for_sampling = leverage_score_for_sampling,
                        ...)
  if (is.null(names(ff.list))) {
    names(ff.list) <- purrr::map_chr(ind_mat, function(x) basename(attr(x, "FilePath")))
  }

  if (all(downsample == "min") && length(ff.list) > 1) {
    pops <- unique(unlist(purrr::map(ff.list, names)))
    names(pops) <- pops
    minrows <- purrr::map_int(pops, function(pop) {
      min(purrr::map_int(purrr::map(ff.list, ~ purrr::pluck(.x, pop)), purrr::compose(nrow, flowCore::exprs)))
    })
    for (pop in pops) {
      for (z in seq_along(ff.list)) {
        if (pop %in% names(ff.list[[z]])) {
          set.seed(seed)
          inds <- sample(c(rep(T, minrows[pop]), rep(F, nrow(ff.list[[z]][[pop]])-minrows[pop])))
          ff.list[[z]][[pop]] <- flowCore::Subset(ff.list[[z]][[pop]], inds)
        }
      }
    }
  }

  return(ff.list)
}

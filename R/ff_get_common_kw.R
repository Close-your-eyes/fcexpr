#' Get intersecting keywords from flow frames
#'
#' @param ff_list list of flow frames
#'
#' @return list of keywords
#' @export
#'
#' @examples
ff_get_common_kw <- function(ff_list) {
  # get common (intersecting keywords)
  # a bit unhandy but vectorized version (below) did not work
  kw_list <- lapply(ff_list, flowCore::keyword)
  common_kw <- as.list(unlist(lapply(names(kw_list[[1]]), function(x) {
    if (length(unique(unlist(sapply(kw_list, "[", x)))) == 1) {
      return(stats::setNames(unique(unlist(sapply(kw_list, "[", x))), nm = x))
    }
  })))
  #common_kw <- Reduce(intersect, lapply(ff_list, flowCore::keyword))
  return(common_kw)
}

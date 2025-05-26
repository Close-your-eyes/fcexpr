#' Get intersecting channel description from flow frames
#'
#' @param ff_list list of flow frames
#'
#' @return character vector
#' @export
#'
#' @examples
ff_get_common_desc <- function(ff_list) {
  # flowCore::markernames(ff_list[[1]])
  desclist <- purrr::map(ff_list, ~flowCore::parameters(.x)@data$desc)
  unique_descs <- purrr::pmap_lgl(desclist, ~length(unique(.x)) == 1)
  common_desc <- character(length(unique_descs))
  for (i in seq_along(unique_descs)) {
    if (!unique_descs[i]) {
      common_desc[i] <- NA
    } else {
      common_desc[i] <- desclist[[1]][i]
    }
  }
  return(common_desc)
}

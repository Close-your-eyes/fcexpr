#' Update params and keywords when exprs was changed
#'
#' @param ff flowframe
#'
#' @returns flowframe
#' @export
#'
#' @examples
#' \dontrun{
#' ff <- ff_update(ff)
#' }
ff_update <- function(ff) {

  if (!is.list(ff)) {
    ff <- list(ff)
  }

  ff <- purrr::map(ff, function(f) {

    f@parameters@data <- f@parameters@data[which(f@parameters@data$name %in% colnames(f@exprs)),]
    keypar <- get_kw_and_pars(f@exprs, ff = f)
    f@parameters <- keypar$params
    f@description <- keypar$keywrd
    return(f)
  })
  return(ff)
}

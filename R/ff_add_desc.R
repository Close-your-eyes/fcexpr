#' Add channel descriptions to flow frame
#'
#' @param ff one flowframe or a list of such
#' @param desc named vector of channel descriptions, names must be channel names
#'
#' @returns ff list
#' @export
#'
#' @examples
#' \dontrun{
#'   ps <- fcexpr::wsx_get_popstats(w1, strip_data = F)[[1]]
#'   channels <- fcexpr::guess_marker_channels(ps) |>
#'     dplyr::filter(!ind %in% c(3,10,7))
#'   channels <- stats::setNames(channels$marker, channels$value)
#'
#'   fflist <- fcexpr::wsp_get_ff(w1,
#'                                FCS.file.folder = file.path(wd, "FCS_files"),
#'                                samples = grep("full", basename(f1[[1]]), value = T),
#'                                population = "Time")
#'   ffs <- purrr::map(fflist[[1]], ~.x[[1]][[1]])
#'   ffs <- ff_add_desc(ffs, desc = channels)
#'
#'   dr <- fcexpr::dimred_to_fcs(ffs, channels = channels)
#' }
ff_add_desc <- function(ff, desc) {

  if (!is.list(ff)) {
    ff <- list(ff)
  }
  if (missing(ff) || missing(desc)) {
    stop("ff or desc missing.")
  }
  if (is.null(names(desc))){
    stop("desc needs channel names as names.")
  }
  ff <- purrr::map(ff, function(f) {
    f@parameters@data$desc <- desc[f@parameters@data$name]
    kwpar <- get_kw_and_pars(exprs = f@exprs, ff = f)
    f@description <- kwpar[["keywrd"]]
    return(f)
  })
  return(ff)
}

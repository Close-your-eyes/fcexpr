#' Title
#'
#' @param n
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(base::pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

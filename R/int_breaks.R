#' Make axis labels that are numeric but should only integer values
#'
#' @param n integer giving the desired number of intervals. Non-integer values are rounded down. (passed to base::pretty)
#' @param ... further arguments passed to base::pretty
#'
#' @return a vector or axis labels of length n
#' @export
#'
#' @examples
#' int_breaks(n=3)(c(1,4,5,6,7,10))
#' int_breaks(n=6)(c(1.1,4.9,5.8,6,7,10.78))
#' # or pass to ggplot2: ggplot2::scale_x_continuous(breaks = int_breaks(n = 4))
int_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(base::pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}



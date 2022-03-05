#' Geometric mean
#'
#' In flow cytometry the geometric mean is a common metric of central tendency.
#' It is less susceptible to outliers than the arithmetic mean, but more susceptible
#' as compared to the median. If you have untransformed values of fluorescence
#' intensities (original output of the machines) this function may not work as you
#' may have negative values in your vector x. Hence, the median has to be used or
#' the data needs transformation (e.g. logicle transformation).
#'
#' @param x numeric vector
#' @param rm.na logical if to remove na or not
#'
#' @return a numeric representing the geometric mean of values in x
#' @export
#'
#' @examples
#' geo_mean(c(1,2,3,10,100))
geo_mean <- function(x, rm.na = F) {
    if (min(x) < 0) {
        stop("No values below zero allowed.")
    }
    if (NA %in% x && !rm.na) {
        stop("NAs found, please fix.")
    } else if (NA %in% x && !rm.na) {
        warning("NAs removed.")
        x <- x[which(!is.na(x))]
    }
    return(exp(mean(log(x))))
}

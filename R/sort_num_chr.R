#' Sort numerical character according to the numbers that the character values represent.
#'
#' Just a convenience function often used to sort patient IDs by numerical value
#' while being considered as factor though.
#'
#' @param x vector of character or factor which represent numbers
#' @param decreasing logical, sort decreasing
#'
#' @return character vector of sorted values
#' @export
#'
#' @examples
#' sort_num_chr(c("3", "4", "2"))
#' sort_num_chr(as.factor(c("3", "4", "2")))
sort_num_chr <- function(x, decreasing = F) {

  if (anyNA(suppressWarnings(as.numeric(as.character(x))))) {
    stop("Not all values in x are numeric. Please check.")
  }

  # as.character first to rescue factors from coercion to numeric factor levels
  return(as.character(sort(as.numeric(as.character(x)), decreasing = decreasing)))
}

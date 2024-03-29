#' Make appealing axis labels for log axis
#'
#' Any numbers below 100 remain as they are. So 1 remain 1 and 10 remains 10 etc. 100 becomes 10^2, 1000 become 10^3 etc.
#' If all numbers in the vector start with a 1 (like c(1, 10, 100, 1000)) then the leading 1x in 1x10^2 is omitted.
#' If any number does not start with a 1, then the leading 1x or 3x or so are added.
#'
#' @param x vector of axis numbers to convert
#'
#' @return vector of modified numbers
#' @export
#'
#' @examples
#' sci10(c(NA,1,10,100,1000,NA))
#' sci10(c(NA,1,10,30))
#' # or pass to ggplot2::scale_y_log10(label = sci10)
sci10 <- function(x) {

  if (!requireNamespace("scales", quietly = T)){
    utils::install.packages("scales")
  }

  if (all(substr(formatC(x[which(!is.na(x))], format = "e"),1,1) %in% c("0", "1"))) {
    x <- sapply(x, function(y) {
      if (is.na(y)) {
        return(NA)
      }
      if (y >= 100) {
        return(gsub("1e\\+", "10^", scales::scientific_format()(y)))
      } else if (y < 1 && y != 0) {
        return(gsub("1e", "10^", scales::scientific_format()(y)))
      } else {
        return(y)
      }
    })
  } else {
    x <- sapply(x, function(y) {
      if (is.na(y)) {
        return(NA)
      }
      if (y >= 100) {
        return(gsub("e\\+", "%*%10^", scales::scientific_format()(y)))
      } else if (y < 1 && y != 0) {
        return(gsub("e", "%*%10^", scales::scientific_format()(y)))
      } else {
        return(y)
      }
    })
  }
  parse(text = x)
}





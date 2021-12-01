#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
sci10 <- function(x) {
  if (all(substr(as.character(x[which(!is.na(x))]),1,1) == "1")) {
    x[which(x>=100)] <- gsub("1e\\+", "10^", scales::scientific_format()(x[which(x>=100)]))
  } else if (any(substr(as.character(x[which(!is.na(x))]),1,1) != "1")) {
    x[which(x>=100)] <- gsub("e\\+", "%*%10^", scales::scientific_format()(x[which(x>=100)]))
  }
  parse(text = x)
}

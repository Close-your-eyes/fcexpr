#' Calculate PCA based on data in flow frame
#'
#' @param ff flow frame
#' @param channels which channels to use, see ff_get_channels for selection
#' options
#' @param exprs expression matrix, either provide ff or exprs
#' @param args arguments to stats::prcomp
#' @param return return pca object or flow frame with PC-channels added
#' @param ... arguments to ff_get_channels
#'
#' @return flow frame or pca object
#' @export
#'
#' @examples
ff_calc_pca <- function(ff,
                        channels = NULL,
                        exprs = NULL,
                        args = list(scale. = F, center = F),
                        return = c("pca", "ff"),
                        ...) {

  return <- rlang::arg_match(return)
  if (return == "ff" && missing(ff)) {
    return <- "pca"
  }

  if (is.null(exprs)) {
    channels <- ff_get_channels(ff, channels = channels, ...)
    stopifnot("pca: less than 3 channels provided." = length(channels) > 2)
    exprs <- flowCore::exprs(ff)[,channels]
    message("Channels for pca: ", paste(channels, collapse = ", "))
  } else {
    stopifnot("exprs must be matrix." = is.matrix(exprs))
  }

  pca.result <- Gmisc::fastDoCall(stats::prcomp, args = c(list(x = exprs), args))

  if (return == "pca") {
    return(pca.result)
  }
  if (return == "ff") {
    return(ff_add_columns(ff, pca.result[["x"]]))
  }
}

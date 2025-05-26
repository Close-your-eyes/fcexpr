#' Calculate UMAP based on flow frame
#'
#' @param ff flow frame
#' @param channels which channels to use, see ff_get_channels for selection
#' options
#' @param exprs expression matrix, either provide ff or exprs
#' @param args arguments to uwot::umap
#' @param seed random seed
#' @param names names of new dimensions to add as channels
#' @param return return UMAP object or flow frame with umap channels added
#' @param ... arguments to ff_get_channels
#'
#' @return flow frame or umap object
#' @export
#'
#' @examples
ff_calc_umap <- function(ff,
                         channels = NULL,
                         exprs = NULL,
                         args = list(n_neighbors = 20, metric = "cosine", verbose = T, scale = T),
                         seed = 42,
                         names = c("UMAP_1", "UMAP_2"),
                         return = c("umap", "ff"),
                         ...) {

  return <- rlang::arg_match(return)
  if (return == "ff" && missing(ff)) {
    return <- "umap"
  }

  if (is.null(exprs)) {
    channels <- ff_get_channels(ff, channels = channels, ...)
    stopifnot("umap: less than 3 channels provided." = length(channels) > 2)
    exprs <- flowCore::exprs(ff)[,channels]
    message("Channels for umap: ", paste(channels, collapse = ", "))
  } else {
    stopifnot("exprs must be matrix." = is.matrix(exprs))
  }

  set.seed(seed)
  umap.dims <- Gmisc::fastDoCall(uwot::umap, args = c(list(X = exprs), args))
  colnames(umap.dims) <- names

  if (return == "umap") {
    return(umap.dims)
  }
  if (return == "ff") {
    return(ff_add_columns(ff, umap.dims))
  }
}

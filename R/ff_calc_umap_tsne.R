#' Calculate UMAP or tSNE dimension reduction based on flow frame
#'
#' @param ff flow frame
#' @param channels which channels to use, see ff_get_channels for selection
#' options
#' @param exprs expression matrix, either provide ff or exprs
#' @param args arguments to fun
#' @param seed random seed
#' @param names names of new dimensions to add as channels
#' @param return return UMAP object or flow frame with umap channels added
#' @param ... arguments to ff_get_channels
#' @param fun uwot::umap or Rtsne::Rtsne
#'
#' @return flow frame or umap object
#' @export
#'
#' @examples
ff_calc_umap_tsne <- function(ff,
                              channels = NULL,
                              exprs = NULL,
                              args = list(n_neighbors = 20, metric = "cosine", verbose = T, scale = T),
                              seed = 42,
                              names = c("UMAP_1", "UMAP_2"),
                              return = c("coords", "ff"),
                              fun = uwot::umap,
                              ...) {

  return <- rlang::arg_match(return)
  if (return == "ff" && missing(ff)) {
    return <- "coords"
  }

  if (is.null(exprs)) {
    channels <- ff_get_channels(ff, channels = channels, ...)
    stopifnot("umap/tSNE: less than 3 channels provided." = length(channels) > 2)
    exprs <- flowCore::exprs(ff)[,channels]
    message("Channels for umap/tSNE: ", paste(channels, collapse = ", "))
  } else {
    stopifnot("exprs must be matrix." = is.matrix(exprs))
  }

  set.seed(seed)
  coords <- Gmisc::fastDoCall(fun, args = c(list(X = exprs), args))
  colnames(coords) <- names

  if (return == "coords") {
    return(coords)
  }
  if (return == "ff") {
    return(ff_add_columns(ff, coords))
  }
}

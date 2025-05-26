#' Calculate self-organizing map based on flow frame
#'
#' @param ff flow frame
#' @param channels which channels to use, see ff_get_channels for selection
#' options
#' @param exprs expression matrix, either provide ff or exprs
#' @param som_args arguments to som_fun
#' @param embedsom_args arguments to EmbedSOM::EmbedSOM
#' @param seed random seed
#' @param names names for SOM channels
#' @param return return SOM object or flow frame with SOM channels added
#' @param ... arguments to ff_get_channels
#' @param som_fun EmbedSOM::SOM or EmbedSOM::GQTSOM
#'
#' @return flow frame or SOM
#' @export
#'
#' @examples
ff_calc_som <- function(ff,
                        channels = NULL,
                        exprs = NULL,
                        som_fun = EmbedSOM::SOM,
                        som_args = list(),
                        embedsom_args = list(),
                        seed = 42,
                        names = c("SOM_1", "SOM_2"),
                        return = c("som", "ff"),
                        ...) {

  return <- rlang::arg_match(return)
  if (return == "ff" && missing(ff)) {
    return <- "som"
  }
  som_fun <- match.fun(som_fun)

  if (is.null(exprs)) {
    channels <- ff_get_channels(ff, channels = channels, ...)
    stopifnot("som: less than 3 channels provided." = length(channels) > 2)
    exprs <- flowCore::exprs(ff)[,channels]
    message("Channels for som: ", paste(channels, collapse = ", "))
  } else {
    stopifnot("exprs must be matrix." = is.matrix(exprs))
  }

  set.seed(seed)
  som.map.dr <- Gmisc::fastDoCall(som_fun, args = c(list(data = exprs), som_args))
  set.seed(seed)
  som.dims <- Gmisc::fastDoCall(EmbedSOM::EmbedSOM, args = c(list(data = exprs, map = som.map.dr), embedsom_args))
  colnames(som.dims) <- names

  if (return == "som") {
    return(som.dims)
  }
  if (return == "ff") {
    return(ff_add_columns(ff, som.dims))
  }
}

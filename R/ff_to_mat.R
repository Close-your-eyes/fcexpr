#' Get expression matrix from a flow frame
#'
#' Channel descriptions from parameters are used as colnames.
#'
#' @param ff flowframe
#' @param grep character to select channel
#' @param channels in addition to grep, pass channels to ff_get_channels
#'
#' @returns matrix or data frame
#' @export
#'
#' @examples
#' \dontrun{
#' mat <- ff_to_mat(ff)
#' }
ff_to_mat <- function(ff,
                      grep = c("trans", "umap", "som", "tsne", "ident"),
                      channels = NULL) {

  # keywords are too complicated. e.g non-numeric
  # could become df then by default

  mat <- ff@exprs

  if (is.null(grep) && is.null(channels)) {
    grep <- colnames(mat)
  }

  if (!is.null(grep)) {
    colinds <- which(grepl(paste(grep, collapse = "|"), colnames(mat), ignore.case = T))
  }

  if (!is.null(channels)) {
    channels <- ff_get_channels(ff, channels = channels)
    channinds <- match(channels, colnames(mat))
    colinds <- unique(c(colinds, channinds))
  }

  if (!length(colinds)) {
    stop("No colnames according to grep found.")
  }
  mat <- mat[, colinds]

  conv <- stats::setNames(ff@parameters@data$desc, ff@parameters@data$name)
  conv <- conv[colnames(mat)]
  conv[which(is.na(conv))] <- names(conv[which(is.na(conv))]) # may not exist
  conv <- stats::setNames(make.unique(conv), names(conv))

  colnames(mat) <- conv[colnames(mat)]

  # put dr cols first if they are there
  um <- which(grepl("umap", colnames(mat), ignore.case = T))
  ts <- which(grepl("tsne", colnames(mat), ignore.case = T))
  so <- which(grepl("som", colnames(mat), ignore.case = T))
  drind <- unique(c(um, ts, so))
  remain <- seq_along(colnames(mat))
  remain <- remain[which(!remain %in% drind)]

  mat <- mat[,c(drind, remain)]

  return(mat)
}

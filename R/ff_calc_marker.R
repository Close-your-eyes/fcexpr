#' Calculate cluster markers in data in a flow frame
#'
#' @param ff flow frame
#' @param exprs expression matrix, if provided is preferred over ff
#' @param channels which channels to use, see ff_get_channels for selection
#' options
#' @param cluster_cols name of columns that contain cluster identities
#' @param mc.cores number of cores to use for calculation
#' @param global do calculate global cluster markers
#' @param pairwise do calculate pairwise cluster markers
#' @param ... arguments to ff_get_channels
#'
#' @return list with data frames
#' @export
#'
#' @examples
ff_calc_marker <- function(ff,
                           exprs = NULL,
                           channels = NULL,
                           cluster_cols,
                           mc.cores = 1,
                           global = T,
                           pairwise = F,
                           ...) {

  # messages if chnalle not found

  # exprs over ff
  if (!is.null(exprs)) {
    if (!is.matrix(exprs)) {
      stop("exprs has to be matrix.")
    }
    if (!is.numeric(exprs)) {
      stop("exprs is not numeric. e.g. use exprs <- matrix(as.numeric(exprs), nrow = nrow(exprs), ncol = ncol(exprs))")
    }
    if (!is.null(channels)) {
      channels <- channels[which(channels %in% colnames(exprs))]
    }

  } else {
    channels <- ff_get_channels(ff, channels = channels, ...)
    exprs <- flowCore::exprs(ff)
  }

  cluster_cols <- cluster_cols[which(cluster_cols %in% colnames(exprs))]
  if (is.null(channels)) {
    # wehn exprs provided but no channels
    channels <- colnames(exprs)[which(!colnames(exprs) %in% cluster_cols)]
  }

  if (length(cluster_cols) == 0) {
    stop("cluster_cols not found.")
  }
  if (length(channels) == 0) {
    stop("channels not found.")
  }

  marker <- lapply(cluster_cols, function(clust_col) {
    # global markers
    if (global) {
      message("Global markers.\nStart: ", Sys.time())
      marker_table <- .calc.global.cluster.marker(
        dat = exprs[,channels],
        cluster = exprs[,clust_col,drop = T],
        mc.cores = mc.cores
      )
      #marker_table[,"channel_desc"] <- channel.desc_augment[marker_table[,"channel"]]
      message("End: ", Sys.time())
    } else {
      marker_table <- NULL
    }

    ## pairwise markers
    if (pairwise) {
      message("Pairwise markers.\nStart: ", Sys.time())
      pair_marker_table <- .calc.pairwise.cluster.marker(
        dat = exprs[,channels],
        cluster = exprs[,clust_col,drop = T],
        mc.cores = mc.cores
      )
      message("End: ", Sys.time())
    } else {
      pair_marker_table <- NULL
    }
    return(list(
      global = marker_table,
      pairwise = pair_marker_table,
      cluster_sizes = table(exprs[,clust_col,drop = T])
    ))
  })
  names(marker) <- cluster_cols
  return(marker)
}

#' Turn subset of ff exprs slot into data frame
#'
#' If there columns like ident or other grouplist they are recoded back
#' from numeric to character. attributes(ff) has to contain grouplist.
#'
#' @param ff flow frame
#' @param channels channels to include, passed to ff_get_channels
#' @param ... arguments to ff_get_channels
#'
#' @return data frame
#' @export
#'
#' @examples
ff_get_exprsdf <- function(ff,
                           channels = c("umap", "som", "res\\.", "ident"),
                           ...) {


  # convert numeric groups back to chars via grouplist attr
  if ("grouplist" %in% names(attributes(ff))) {
    channels <- unique(c(channels, names(attr(ff, "grouplist"))))
    groupcols <- names(attr(ff, "grouplist"))
  }
  # replace_NA_desc: if desc is NA then name of channelis used
  channels <- ff_get_channels(ff, channels = channels, replace_NA_desc = T, ...)
  df <- as.data.frame(flowCore::exprs(ff)[,channels])
  for (i in groupcols) {
    df[,i] <- names(attr(ff, "grouplist")[[i]][ff@exprs[,i]])
  }

  # use desc as colnames
  channelsinv <- stats::setNames(names(channels), channels)
  names(df) <- channelsinv[names(df)]

  return(df)
}

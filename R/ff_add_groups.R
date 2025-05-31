#' Add sample group identifier to a flow frame
#'
#' @param ff flow frame
#' @param ident_col column name with sample identifier
#' @param grouplist named list of sample information
#' @param overwrite do overwrite exisitng columns in ff
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' ff <- ff_merge(ff_list = list(ff1, ff2))
#' ff_add_groups(ff, grouplist = list(cond2 = c("d","k"), dis2 = c("e", "r")))
#' }
ff_add_groups <- function(ff,
                          ident_col = "ident",
                          grouplist = NULL,
                          overwrite = F) {


  if (is.null(grouplist)) {
    return(ff) #see ff_merge
  }
  if (!is.list(grouplist)) {
    stop("grouplist has to be a list.")
  }
  if (is.null(names(grouplist))) {
    stop("grouplist has to have names. These names will become channel names in the FCS file.")
  }
  if (any(unlist(lapply(grouplist, function(x) is.na(x))))) {
    stop("NA found in sample infos.")
  }
  if (!ident_col %in% colnames(flowCore::exprs(ff))) {
    stop("ident_col not found in exprs.")
  }

  # convert non numeric to numeric
  if (any(non_num <- !purrr::map_lgl(grouplist, is.numeric))) {
    grouplist[non_num] <- purrr::map(grouplist[non_num], ~stats::setNames(as.numeric(as.factor(.x)), .x))
  }

  relident <- rle(flowCore::exprs(ff)[,ident_col])
  nident <- length(unique(relident[["values"]]))

  if (any(lengths(grouplist) != nident)) {
    stop("Length of each grouplist entry has to match the number of samples which is: ", nident,".")
  } else {
    consecutive_ident <- nident == length(relident[["values"]])
    if (!consecutive_ident) {
      stop("idents are not consecutive. please reorder to guarantee that.")
    }
  }

  # expand group idents accordingly
  groupcols <- sapply(grouplist, rep, times = relident[["lengths"]], simplify = T, USE.NAMES = F)
  rownames(groupcols) <- NULL

  # run the actual addition of columns
  ff <- ff_add_columns(
    ff = ff,
    mat = groupcols,
    overwrite = overwrite
  )

  # join lists if grouplist is already in attr
  if (!is.null(attr(ff, "grouplist"))) {
    grouplist <- c(attr(ff, "grouplist"), grouplist)
  }
  attr(ff, "grouplist") <- grouplist
  return(ff)

}




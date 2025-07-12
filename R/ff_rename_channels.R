#' Rename one or many channels of a flowframe
#'
#' Alter all slots: exprs, params, desc (keys).
#'
#' @param ff flowframe
#' @param rename_vec named vector with current channel name as name
#' and new names as values (see example)
#'
#' @return updated flowframe
#' @export
#'
#' @examples
#' \dontrun{
#' ff <- flowCore::read.FCS("path/to/file/Specimen_001_Tube_001_001.fcs")
#' colnames(flowCore::exprs(ff)) # inspect
#' rename_vec <- setNames(c("Pe", "APC"), c("yg-586/15-E-A", "r-670/30-C-A"))
#' ffnew <- ff_rename_channels(ff, rename_vec)
#' # optionally: overwrite fcs file on disk
#' flowCore::write.FCS(ffnew, filename = "path/to/file/Specimen_001_Tube_001_001.fcs")
#' }
ff_rename_channels <- function(ff, rename_vec) {

  rename_vec <- rename_vec[which(names(rename_vec) %in% colnames(ff@exprs))]
  if (!length(rename_vec)) {
    message("names of rename_vec not found in ff.")
    return(ff)
  }
  # change exprs
  colnames(ff@exprs)[match(names(rename_vec), colnames(ff@exprs))] <- rename_vec

  # change params
  channelnums <- names(ff@parameters@data$name[match(names(rename_vec), ff@parameters@data$name)])
  ff@parameters@data$name[match(names(rename), ff@parameters@data$name)] <- rename_vec

  # change desc
  for (i in channelnums) {
    ff@description[[i]] <- ff@parameters@data[gsub("N", "", i),"name"]
  }
  return(ff)
}

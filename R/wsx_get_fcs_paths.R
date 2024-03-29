#' Get file paths to the fcs files in a workspace
#'
#' If one has a large experiment with many fcs files and many workspaces with various
#' gatings, finding out which fcs files are gated in which workspace may be helpful.
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#' @param split logical whether to split file paths by groups; if yes a list is returned, if no a data frame
#'
#' @return a list or data frame of paths to fcs files
#' @export
#'
#' @examples
#' \dontrun{
#' wsx_get_fcs_paths(ws)
#' }
wsx_get_fcs_paths <- function(ws, split = T) {

  ws <- check_ws(ws)

  FilePath <- gsub("^file:", "", xml2::xml_attr(xml2::xml_find_all(ws, ".//DataSet"), "uri"))
  sampleID <- gsub("^file:", "", xml2::xml_attr(xml2::xml_find_all(ws, ".//DataSet"), "sampleID"))
  paths <- dplyr::full_join(wsx_get_groups(ws, filter_AllSamples = F, collapse_groups = F), data.frame(FilePath, sampleID), by = "sampleID")

  paths$FileName <- basename(paths$FilePath)
  if (split) {
    paths <- split(paths[,"FilePath",drop=T], paths[,"FlowJoGroup",drop=T])
  }

  return(paths)
}

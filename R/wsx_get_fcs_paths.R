#' Get file paths to the fcs files in a workspace
#'
#' If one has a large experiment with many fcs files and many workspaces with various
#' gatings, finding out which fcs files are gated in which workspace may be helpful.
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#' @param split logical whether to split file paths by groups; if TRUE a list is returned, if FALSE a data frame
#' @param filter_AllSamples logical whether to filter the All Samples Group in case the fcs file is also part of another group
#'
#' @return a list or data frame of paths to fcs files
#' @export
#'
#' @examples
#' \dontrun{
#' wsx_get_fcs_paths(ws)
#' }
wsx_get_fcs_paths <- function(ws,
                              split = T,
                              filter_AllSamples = F) {

  ws <- check_ws(ws)
  if (!is.logical(split)) {
    stop("split has be TRUE or FALSE.")
  }
  paths <- dplyr::full_join(wsx_get_groups(ws,
                                           filter_AllSamples = filter_AllSamples,
                                           collapse = "nest",
                                           force_collapse = T),
                            data.frame(FilePath = gsub("^file:", "", xml2::xml_attr(xml2::xml_find_all(ws, ".//DataSet"), "uri")),
                                       sampleID = gsub("^file:", "", xml2::xml_attr(xml2::xml_find_all(ws, ".//DataSet"), "sampleID"))),
                            by = "sampleID")

  paths$FileName <- basename(paths$FilePath)
  paths <- tidyr::unnest(paths, "FlowJoGroup")

  if (split) {
    paths <- split(paths[,"FilePath",drop=T], paths[,"FlowJoGroup",drop=T])
  }
  return(paths)
}

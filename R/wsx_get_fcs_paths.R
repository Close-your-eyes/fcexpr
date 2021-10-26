#' Title
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#'
#' @return a list of paths to fcs files grouped by groups in the ws
#' @export
#'
#' @examples
#' \dontrun{
#' wsx_get_fcs_paths(ws)
#' }
wsx_get_fcs_paths <- function(ws) {

  if (is.character(ws)) {
    if (!file.exists(ws)) {
      stop("ws file not found.")
    }
    if (length(ws) > 1) {
      stop("Only one ws at a time.")
    }
    ws <- xml2::read_xml(ws)
  }
  if (!any(class(ws) == "xml_document")) {
    stop("ws must be a xml-document or a character path to its location on disk")
  }

  if (xml2::xml_attr(ws, "flowJoVersion") != "10.7.1") {
    warning("This function was tested with a FlowJo wsp from version 10.7.1. Other version may lead to unexpected results.")
  }

  FilePath <- gsub("^file:", "", xml2::xml_attr(xml2::xml_find_all(ws, ".//DataSet"), "uri"))
  sampleID <- gsub("^file:", "", xml2::xml_attr(xml2::xml_find_all(ws, ".//DataSet"), "sampleID"))
  df <- dplyr::full_join(wsx_get_groups(ws, filter_AllSamples = F), data.frame(FilePath, sampleID))[-2]
  paths <- split(df[,"FilePath"], df$group)

  return(paths)
}

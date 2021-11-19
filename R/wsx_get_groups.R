#' Retrieve groups within a flowjo workspace and associated samples (sampleID)
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#' @param filter_AllSamples logical whether to filter the All Samples Group in case the fcs file is also part of another group
#' @param collapse_groups logical whether to collapse multiple group-belongings of samples into a list-column in the data frame
#'
#' @return a data frame
#' @export
#'
#' @examples
#' \dontrun{
#' wsx_get_groups(ws)
#' }
wsx_get_groups <- function(ws, filter_AllSamples = T, collapse_groups = T) {

  ws <- check_ws(ws)
  if (!any(class(ws) == "xml_document")) {
    stop("ws must be a xml-document or a character path to its location on disk")
  }

  if (xml2::xml_attr(ws, "flowJoVersion") != "10.7.1") {
    warning("This function was tested with a FlowJo wsp from version 10.7.1. Other version may lead to unexpected results.")
  }

  g <- sapply(xml2::xml_children(xml2::xml_child(ws, "Groups")), function(y) {
    xml2::xml_attrs(y)[["name"]]
  })
  gs <- lapply(xml2::xml_children(xml2::xml_child(ws, "Groups")), function(y) {
    unlist(xml2::xml_attrs(xml2::xml_children(xml2::xml_child(xml2::xml_child(y, "Group"), "SampleRefs"))))
  })
  gr <- data.frame(group = rep(g, lengths(gs)),  sampleID = unlist(gs))

  if (filter_AllSamples) {
    gr <- do.call(rbind, lapply(unique(gr$sampleID), function(y) {
      if (length(gr[which(gr$sampleID == y),"group"]) > 1) {
        gr[base::intersect(which(gr$sampleID == y), which(gr$group != "All Samples")), ]
      } else if (gr[which(gr$sampleID == y),"group"] == "All Samples") {
        gr[which(gr$sampleID == y), ]
      } else {
        stop("group error occured.")
      }
    }))
  }

  if (collapse_groups) {
    gr <- do.call(rbind, lapply(unique(gr$sampleID), function(y) {
      data.frame(group = I(list(gr[which(gr$sampleID == y),"group"])), sampleID = y)
    }))
  }

  return(gr)
}



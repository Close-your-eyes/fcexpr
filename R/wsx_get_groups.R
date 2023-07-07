#' Retrieve groups within a flowjo workspace and associated samples (sampleID)
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#' @param filter_AllSamples logical whether to filter the All Samples Group in case the fcs file is also part of another group
#' @param collapse_groups logical whether to collapse multiple group-belongings of samples into a list-column in the data frame
#' @param collapse_to string how to collapse groups; to collapse to a list-column pass 'list'; to collapse to a string provide any separator string like ';', ',' or '_-_'
#'
#' @return a data frame
#' @export
#'
#' @examples
#' \dontrun{
#' wsx_get_groups(ws)
#' }
wsx_get_groups <- function(ws,
                           filter_AllSamples = T,
                           collapse_groups = T,
                           collapse_to = "list") {

  ws <- check_ws(ws)

  g <- sapply(xml2::xml_children(xml2::xml_child(ws, "Groups")), function(y) {
    xml2::xml_attrs(y)[["name"]]
  })
  gs <- lapply(xml2::xml_children(xml2::xml_child(ws, "Groups")), function(y) {
    unlist(xml2::xml_attrs(xml2::xml_children(xml2::xml_child(xml2::xml_child(y, "Group"), "SampleRefs"))))
  })
  gr <- data.frame(FlowJoGroup = rep(g, lengths(gs)),
                   sampleID = unlist(gs))

  if (filter_AllSamples) {
    gr <- do.call(rbind, lapply(unique(gr$sampleID), function(y) {
      if (length(gr[which(gr$sampleID == y),"FlowJoGroup"]) > 1) {
        gr[base::intersect(which(gr$sampleID == y), which(gr$FlowJoGroup != "All Samples")), ]
      } else if (gr[which(gr$sampleID == y),"FlowJoGroup"] == "All Samples") {
        gr[which(gr$sampleID == y), ]
      } else {
        stop("FlowJoGroup error occured.")
      }
    }))
  }

  if (collapse_groups && any(duplicated(gr$sampleID))) {
    if (collapse_to == "list") {
      gr <- do.call(rbind, lapply(unique(gr$sampleID), function(y) {
        data.frame(FlowJoGroup = I(list(gr[which(gr$sampleID == y),"FlowJoGroup"])), sampleID = y)
      }))
    }
    if (collapse_to != "list") {
      gr <- do.call(rbind, lapply(unique(gr$sampleID), function(y) {
        data.frame(FlowJoGroup = paste(gr[which(gr$sampleID == y),"FlowJoGroup"], collapse = collapse_to), sampleID = y)
      }))
    }
  }

  return(gr)
}

#' Retrieve groups within a flowjo workspace and associated samples (sampleID)
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#' @param filter_AllSamples logical whether to filter the All Samples Group in case the fcs file is also part of another group
#' @param collapse string how to collapse groups; to collapse to a list-column pass 'list'; to collapse to a string provide any separator string like ';', ',' or '_-_'
#' @param force_collapse collapse in any case, even though there is only one group per sample
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
                           collapse = "nest",
                           force_collapse = F) {

  if (!is.null(collapse) && collapse %in% c("nest", "list") && !is.character(collapse)) {
    stop("collapse has to be NULL, 'nest', 'list' or a character.")
  }

  if (!is.logical(filter_AllSamples)) {
    stop("filter_AllSamples has be TRUE or FALSE.")
  }

  ws <- check_ws(ws)
  g <- sapply(xml2::xml_attrs(xml2::xml_children(xml2::xml_child(ws, "Groups"))), "[[", "name")
  gs <- sapply(xml2::xml_child(xml2::xml_child(xml2::xml_children(xml2::xml_child(ws, "Groups")), "Group"), "SampleRefs"), function(x) unlist(xml2::xml_attrs(xml2::xml_children(x))))
  gr <- data.frame(FlowJoGroup = rep(g, lengths(gs)),
                   sampleID = unlist(gs))

  if (filter_AllSamples) {
    gr <- do.call(rbind, lapply(split(gr, gr$sampleID), function(x) {
      if (nrow(x) > 1) {
        return(x[which(x$FlowJoGroup != "All Samples"),])
      } else {
        return(x)
      }
    }))
  }

  if (!is.null(collapse) && (anyDuplicated(gr$sampleID) != 0 || force_collapse)) {
    if (collapse == "nest") {
      gr <- tidyr::nest(gr, FlowJoGroup = FlowJoGroup)
    } else if (collapse == "list") {
      gr <- do.call(rbind, lapply(split(gr, gr$sampleID), function(x) {
        data.frame(FlowJoGroup = I(list(x$FlowJoGroup)), sampleID = x$sampleID[1])
      }))
    } else {
      gr <- do.call(rbind, lapply(split(gr, gr$sampleID), function(x) {
        data.frame(FlowJoGroup = paste(x$FlowJoGroup, collapse = collapse), sampleID = x$sampleID[1])
      }))
    }
  }
  return(gr)
}

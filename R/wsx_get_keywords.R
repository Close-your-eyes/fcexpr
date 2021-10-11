#' Title
#'
#' @param ws
#'
#' @return
#' @export
#'
#' @examples
wsx_get_keywords <- function(ws) {

  if (is.character(ws)) {
    if (!file.exists(ws)) {
      stop("wsp file not found.")
    }
    ws <- xml2::read_xml(ws)
  } else if (!any(class(ws) == "xml_document")) {
    stop("x must be a xml-document or a character path to its location on disk")
  }

  keywords <- lapply(xml2::xml_children(xml2::xml_child(f, "SampleList")), function(x) {
    k <- xml2::xml_attrs(xml2::xml_contents(xml2::xml_child(x, "Keywords")))
    do.call(rbind, lapply(k, function(y) data.frame(as.list(y))))
  })

  return(keywords)
}



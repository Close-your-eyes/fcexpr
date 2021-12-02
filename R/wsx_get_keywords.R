#' Get keywords from a flowjo workspace
#'
#' Keywords that are annotated in a flowjo workspace only exist in that very workspace. They are not hard-coded into the respective FCS files, as long as you do not re-export them.
#' This behavior is part of flowjos philosophy not to alter the FCS files. So, the keywords that you add in workspace are like a mask on top of entries in the FCS files.
#' When the workspace is lost though or the connection of samples in the workspace and FCS files gets corrupted your keyword-annotation is lost.
#' Also, when you decide to put the same FCS files into another workspace the keywords will not be transferred.
#'
#'
#' @param ws a path to a workspace or a the parsed xml document (xml2::read_xml(ws))
#'
#' @return a list of data.frames. one list entry for each sample, each row of data.frame representing one keyword
#' @export
#'
#' @examples
#' \dontrun{
#' # ws is the path to a flowjo workspace to extract keywords from
#' k <- fcexpr:::wsx_get_keywords(ws)
#' # make it a data.frame
#' kk <- do.call(rbind, k)
#' # create a FileName-column
#' kk$FileName <- rep(names(k), sapply(k,nrow))
#' # do not use kk$FileName <- rownames(kk) as rownames have a suffix to make them unique
#' # select keywords and make it a wider data.frame for subsequent joining with the sampledescription
#' kk <- kk[which(kk$name %in% c("$CYT", "$OP")),]
#' kk <- tidyr::pivot_wider(kk, names_from = name, values_from = value)
#' # make valid column names (general)
#' names(kk) <- make.names(names(kk))
#' # or in this special case
#' names(kk) <- gsub("\\$", "", names(kk))
#' # sd is your sampledescription
#' sd <- dplyr::left_join(sd, kk, by = "FileName")
#' }
wsx_get_keywords <- function(ws) {

  ws <- check_ws(ws)

  keywords <- lapply(xml2::xml_children(xml2::xml_child(ws, "SampleList")), function(x) {
    k <- xml2::xml_attrs(xml2::xml_contents(xml2::xml_child(x, "Keywords")))
    df <- data.frame(name = sapply(k, "[", 1), value = sapply(k, "[", 2))
    rownames(df) <- NULL
    return(df)
  })

  names(keywords) <- wsp_xml_get_samples(ws)[,"FileName"]
  return(keywords)
}



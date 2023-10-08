#' Get keywords from a flowjo workspace
#'
#' Keywords that are annotated in a flowjo workspace only exist in that very workspace. They are not hard-coded into the respective FCS files, as long as you do not re-export them.
#' This behavior is part of flowjos philosophy not to alter the FCS files. So, the keywords that you add in workspace are like a mask on top of entries in the FCS files.
#' When the workspace is lost though or the connection of samples in the workspace and FCS files gets corrupted your keyword-annotation is lost.
#' Also, when you decide to put the same FCS files into another workspace the keywords will not be transferred.
#'
#'
#' @param ws a path to a workspace or a the parsed xml document (xml2::read_xml(ws))
#' @param return how to return keywords: data.frame with 2 column or named vector?
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param keywords which keywords to return
#' @param ... ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
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
wsx_get_keywords <- function(ws,
                             keywords = NULL,
                             return = c("data.frame", "vector"),
                             lapply_fun = lapply,
                             ...) {

  ws <- check_ws(ws)
  lapply_fun <- match.fun(lapply_fun)
  return <- match.arg(arg = return, choices = c("data.frame", "vector"))

  k_return <- lapply_fun(xml2::xml_children(xml2::xml_child(ws, "SampleList")), function(x) {
    keys <- xml2::xml_attrs(xml2::xml_contents(xml2::xml_child(x, "Keywords")))
    keys <- stats::setNames(sapply(keys, "[", 2), sapply(keys, "[", 1))

    if (!is.null(keywords)) {
      keys <- keys[which(names(keys) %in% keywords)]
    }

    if (length(keys) == 0) {
      keys <- NULL
    }

    if (return == "data.frame" && !is.null(keys)) {
      keys <- utils::stack(keys)
      names(keys) <- c("value", "name")
      keys <- keys[,c(2,1)]
      keys$name <- as.character(keys$name)
    }
    return(keys)
  }, ...)

  names(k_return) <- basename(xml2::xml_attr(xml2::xml_child(xml2::xml_children(xml2::xml_child(ws, "SampleList")), "DataSet"), "uri"))
  return(k_return)
}



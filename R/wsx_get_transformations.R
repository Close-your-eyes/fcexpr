#' Get transformation parameters from a flowjo workspace
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param groups which flowjo groups to include
#' @param group_filenames create a list column of filenames that have equal transformation pararmeters
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#'
#' @return data frame with transformation parameters for each channel
#' @export
#'
#' @examples
#' \dontrun{
#' trafo_df <- wsx_get_transformations(ws = ws)
#' }
wsx_get_transformations <- function(ws,
                                    lapply_fun = lapply,
                                    groups = NULL,
                                    group_filenames = F,
                                    ...) {


  lapply_fun <- match.fun(lapply_fun)

  ws <- check_ws(ws)

  ids <- wsx_get_groups(ws)
  if (is.null(groups)) {
    groups <- unique(ids[,"group", drop=T])
  }
  ids <- ids[which(ids$group %in% groups),"sampleID"]
  rel_nodes <- xml2::xml_children(xml2::xml_child(ws, "SampleList"))
  rel_nodes <- rel_nodes[which(sapply(seq_along(rel_nodes), function(x) xml2::xml_attrs(xml2::xml_child(rel_nodes[[x]], "DataSet"))[["sampleID"]]) %in% ids)]

  trans_df <- dplyr::bind_rows(lapply_fun(seq_along(rel_nodes), function(n) {

    FilePath <- gsub("^file:", "", xml2::xml_attr(xml2::xml_child(rel_nodes[n], "DataSet")[[1]], "uri"))
    FileName <- basename(FilePath)

    tr <- xml2::xml_child(rel_nodes[[n]], "Transformations")

    dplyr::bind_rows(lapply(seq_along(xml2::xml_children(tr)), function(nn) {
      tr_attr <- xml2::xml_attrs(xml2::xml_child(tr, nn))
      tr_name <- xml2::xml_name(xml2::xml_child(tr, nn))
      channel_name <- xml2::xml_attrs(xml2::xml_child(xml2::xml_child(tr, nn), 1))

      tibble::tibble(FileName = FileName,
                     channel = channel_name,
                     transformation = tr_name,
                     transformation_pars = list(tr_attr))
    }))
  }, ...))

  if (group_filenames) {
    trans_df <- dplyr::summarise(dplyr::group_by(trans_df, channel, transformation, transformation_pars), FileName = list(FileName), .groups = "drop")
  }

  return(trans_df)
}

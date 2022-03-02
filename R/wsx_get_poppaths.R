#' Obtain paths to populations/nodes/gates in a flowjo workspace
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#' @param collapse logical whether to collapse FileNames to a list column
#'
#' @return data frame
#' @export
#'
#' @examples
#' \dontrun{
#' pp <- wsx_get_poppaths(ws, collapse = F)
#' # check which files have equal gating trees
#' # and get the node (or path, or population) names
#' pp <- pp %>%
#' dplyr::group_by(PopulationFullPath, Population, ws) %>%
#' dplyr::summarise(FileName = list(FileName), .groups = "drop")
#' }
wsx_get_poppaths <- function(ws, collapse = T) {

  ws <- check_ws(ws)

  gg <- xml2::xml_find_all(xml2::xml_child(ws, "SampleList"), ".//Gate|.//Dependents")

  gates <- lapply(seq_along(gg), function(n) {

    prnts <- xml2::xml_parents(gg[n])

    s_node <- prnts[which(xml2::xml_name(prnts) == "Sample")]
    sampleID <- xml2::xml_attr(xml2::xml_child(s_node, "DataSet")[[1]], "sampleID")
    FilePath <- gsub("^file:", "", xml2::xml_attr(xml2::xml_child(s_node, "DataSet")[[1]], "uri"))
    FileName <- basename(FilePath)

    p_nodes <- prnts[which(xml2::xml_name(prnts) %in% c("AndNode", "OrNode", "NotNode", "Population"))]
    PopulationFullPath <- paste(rev(xml2::xml_attr(p_nodes, "name")), collapse = "/")
    Population <- basename(PopulationFullPath)

    if (xml2::xml_name(gg[n]) == "Dependents") {
      origin <- "Dependents"
    } else {
      origin <- "Gate"
    }


    return(data.frame(FileName = FileName,
                      PopulationFullPath = PopulationFullPath,
                      Population = Population,
                      sampleID = sampleID,
                      origin = origin,
                      stringsAsFactors = F))
  })


  gates_df <- do.call(rbind, gates)
  gates_list <- split(gates_df, gates_df$sampleID)
  # remove duplicate rows from gate+dependents
  gates_list <- lapply(gates_list, function(y) {
    ex <- base::intersect(c(which(duplicated(y$PopulationFullPath)),
                            which(duplicated(y$PopulationFullPath, fromLast=T))),
                          which(y$origin == "Dependents"))
    if (length(ex) > 0) {
      y <- y[-ex,]
    }
    return(y)
  })

  full_paths <- unique(lapply(gates_list, function(y) {
    if (length(unique(y$PopulationFullPath)) != length(y$PopulationFullPath)) {
      stop("PopulationFullPaths not unique which cannot or should not be. Check.")
    }
    y$PopulationFullPath
  }))
  auto_paths <- lapply(full_paths, function(y) shortest_unique_path(y))

  for (y in seq_along(gates_list)) {
    gates_list[[y]][["Population"]] <- auto_paths[[which(sapply(full_paths, function(z) identical(z,  gates_list[[y]][["PopulationFullPath"]])))]]
  }
  gates_out <- do.call(rbind, gates_list)
  gates_out <- gates_out[order(gates_out$FileName),]
  rownames(gates_out) = seq(1,nrow(gates_out),1)
  gates_out[,"ws"] <- basename(xml2::xml_attr(ws, "nonAutoSaveFileName"))
  gates_out <- gates_out[,-which(names(gates_out) == "sampleID")]

  if (collapse) {
    gates_out <- gates_out %>%
      dplyr::group_by(PopulationFullPath, Population, ws) %>%
      dplyr::summarise(FileName = list(FileName), .groups = "drop")
  }

  return(as.data.frame(gates_out))
}



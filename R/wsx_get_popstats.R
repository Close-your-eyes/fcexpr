#' Import data from a wsp file
#'
#' Flowjo wsp files contain many information like gated event counts, statistics and keywords
#' from FCS files. These may be accessed without a dongle and can be read completely independent of
#' the respective FCS files once the gating has been conducted.
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#' @param return_stats return statistics next to cells counts
#' @param groups which flowjo groups to include
#' @param ... arguments passed to wsx_get_groups
#'
#' @return data frame with cells counts or a list with counts and statistics if return_stats = T
#' @export
#'
#' @examples
#' \dontrun{
#' # When the script is saved to R_scripts in the experiment folder,
#' # get the absolute path to the folder
#' wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
#' # find workspaces
#' ws <- list.files(path = wd, pattern = '\\.wsp$', recursive = T, full.names = T)
#' # import the population counts:
#' wsx_get_popstats(ws = ws[[1]])
#' }
wsx_get_popstats <- function(ws,
                             groups = NULL,
                             return_stats = T,
                             ...) {

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
    Parent <- if (PopulationFullPath == basename(PopulationFullPath)) {"root"} else {dirname(PopulationFullPath)}
    Population <- basename(PopulationFullPath)
    Count <- xml2::xml_attr(p_nodes[1], "count")
    ParentCount <- if (length(p_nodes) > 1) {xml2::xml_attr(p_nodes[2], "count")} else {xml2::xml_attr(xml2::xml_child(s_node, "SampleNode"), "count")}
    gate_level <- length(p_nodes)

    xDim <- tryCatch({
      xml2::xml_attr(xml2::xml_child(xml2::xml_child(xml2::xml_child(gg[n]), 1)), "name")
    }, error = function(e) {
      NA
    })

    yDim <- tryCatch({
      xml2::xml_attr(xml2::xml_child(xml2::xml_child(xml2::xml_child(gg[n]), 2)), "name")
    }, error = function(e) {
      NA
    })

    if (xml2::xml_name(gg[n]) == "Dependents") {
      origin <- "Dependents"
    } else {
      origin <- "Gate"
    }

    gate_id <- xml2::xml_attr(gg[n], "id")
    parentgate_id <- xml2::xml_attr(gg[n], "parent_id")
    eventsInside <- xml2::xml_attr(xml2::xml_child(gg[n]), "eventsInside")

    return(data.frame(FileName = FileName,
                      PopulationFullPath = PopulationFullPath,
                      Parent = Parent,
                      Population = Population,
                      Count = as.numeric(Count),
                      ParentCount = as.numeric(ParentCount),
                      FractionOfParent = as.numeric(Count)/as.numeric(ParentCount),
                      xDim = xDim,
                      yDim = yDim,
                      gate_id = gate_id,
                      parentgate_id = parentgate_id,
                      eventsInside = eventsInside,
                      sampleID = sampleID,
                      FilePath = FilePath,
                      gate_level = gate_level,
                      origin = origin,
                      n = n,
                      stringsAsFactors = F)
    )
  })

  roots <- do.call(rbind, lapply(xml2::xml_children(xml2::xml_child(ws, "SampleList")), function(y) {
    data.frame(FileName = basename(xml2::xml_attr(xml2::xml_child(y, "DataSet"), "uri")),
               PopulationFullPath = "root",
               Parent = NA,
               Population = "root",
               Count = as.numeric(xml2::xml_attr(xml2::xml_child(y, "SampleNode"), "count")),
               ParentCount = NA,
               FractionOfParent = NA,
               xDim = NA,
               yDim = NA,
               gate_id = NA,
               parentgate_id = NA,
               eventsInside = NA,
               sampleID = xml2::xml_attr(xml2::xml_child(y, "DataSet"), "sampleID"),
               FilePath = gsub("^file:", "", xml2::xml_attr(xml2::xml_child(y, "DataSet"), "uri")),
               gate_level = 0,
               origin = "root",
               n = 0,
               stringsAsFactors = F)
  }))

  gates_df <- do.call(rbind, gates)
  gates_df <- rbind(roots,gates_df)
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
  gates_out <- dplyr::left_join(gates_out, wsx_get_groups(ws, ...), by = "sampleID")
  if (!is.null(groups)) {
    gates_out <- gates_out[which(sapply(gates_out$group, function(x) length(intersect(groups, x)) > 0)),]
    gates_out$group <- sapply(gates_out$group, function(x) intersect(groups, x))
  }
  gates_out[,"ws"] <- basename(xml2::xml_attr(ws, "nonAutoSaveFileName"))
  gates_out <- gates_out[order(gates_out$FileName),]
  rownames(gates_out) = seq(1,nrow(gates_out),1)

  gates_out <- gates_out[,which(!names(gates_out) %in% c("gate_id", "parentgate_id", "sampleID", "origin", "n", "gate_level"))]

  if (return_stats) {
    stats_out <- do.call(rbind, lapply(seq_along(xml2::xml_children(xml2::xml_child(ws, "SampleList"))), function(n) {

      node <- xml2::xml_children(xml2::xml_child(ws, "SampleList"))[n]
      stats <- xml2::xml_find_all(node, ".//Statistic")

      stats_df <- do.call(rbind, lapply(stats, function(x) {
        prnts <- xml2::xml_parents(x)
        p_nodes <- prnts[which(xml2::xml_name(prnts) %in% c("AndNode", "OrNode", "NotNode", "Population"))]

        sampleID <- xml2::xml_attr(xml2::xml_child(x, "DataSet"), "sampleID")
        FilePath <- gsub("^file:", "", xml2::xml_attr(xml2::xml_child(node, "DataSet"), "uri"))
        FileName <- basename(FilePath)
        PopulationFullPath <- if (length(p_nodes) == 0) {"root"} else {paste(rev(xml2::xml_attr(p_nodes, "name")), collapse = "/")}


        data.frame(FileName = FileName,
                   PopulationFullPath = PopulationFullPath,
                   statistic = xml2::xml_attr(x, "name"),
                   channel = xml2::xml_attr(x, "id"),
                   value = as.numeric(xml2::xml_attr(x, "value")),
                   FilePath = FilePath,
                   stringsAsFactors = F)
      }))

      return(stats_df)
    }))
    return(list(counts = gates_out, stats = stats_out))
  }
  return(gates_out)
}

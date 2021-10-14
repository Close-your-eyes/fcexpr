#' Convenient function to import population counts from a flowjo workspace
#'
#' In comparison to ws_get_popstats another backend (xml2) is used to read data from a flowjo .wsp-file (under the hood an .xml-file).
#' No groups or similar have to be selected - filtering has to be conducted afterwards.
#' Also, no restriction with respect to inconsistent gating trees are made. Every sample in handled invididually.
#' Since original FCS files are not touch (all infoamtion are stored in the .wsp file) no FCS folder path has to be provided.
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#'
#' @return returns a data.frame with population counts
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
#' wsx_get_popstats(ws = ws[1])
#' }
wsx_get_popstats <- function(ws) {

  ## graph in xml refers to the current (or last) channels looked at in FJ; not where the Gate is set
  if (is.character(ws)) {
    if (!file.exists(ws)) {
      stop("wsp file not found.")
    }
    ws <- xml2::read_xml(ws)
  }
  if (!any(class(ws) == "xml_document")) {
    stop("x must be a xml-document or a character path to its location on disk")
  }

  ## check FJ version
  if (xml2::xml_attr(ws, "flowJoVersion") != "10.7.1") {
    warning("This function was tested with a FlowJo wsp from version 10.7.1. Other version may lead to unexpected results.")
  }

  gg <- xml2::xml_find_all(xml2::xml_child(ws, "SampleList"), ".//Gate|.//Dependents")
  gates <- lapply(seq_along(gg), function(n) {

    prnts <- xml2::xml_parents(gg[n])

    s_node <- prnts[which(xml2::xml_name(prnts) == "Sample")]
    sampleID <- xml2::xml_attr(xml2::xml_child(s_node, "DataSet")[[1]], "sampleID")
    FilePath <- xml2::xml_attr(xml2::xml_child(s_node, "DataSet")[[1]], "uri")
    FileName <- basename(FilePath)

    p_nodes <- prnts[which(xml2::xml_name(prnts) %in% c("AndNode", "OrNode", "NotNode", "Population"))]
    PopulationFullPath <- paste(rev(xml2::xml_attr(p_nodes, "name")), collapse = "/")
    Parent <- if (PopulationFullPath == basename(PopulationFullPath)) {"root"} else {dirname(PopulationFullPath)}
    Population <- basename(PopulationFullPath)
    Count <- xml2::xml_attr(p_nodes[1], "count")
    ParentCount <- if (length(p_nodes) > 1) {xml2::xml_attr(p_nodes[2], "count")} else {xml2::xml_attr(xml2::xml_child(s_node, "SampleNode"), "count")}
    gate_level <- length(p_nodes)

    if (xml2::xml_name(gg[n]) == "Dependents") {
      xDim <- NA
      yDim <- NA
    } else {
      xDim <- xml2::xml_attr(xml2::xml_child(xml2::xml_child(xml2::xml_child(gg[n]), 1)), "name")
      yDim <- xml2::xml_attr(xml2::xml_child(xml2::xml_child(xml2::xml_child(gg[n]), 2)), "name")
    }

    if (xml2::xml_name(gg[n]) == "Dependents") {
      origin <- "Dependents"
    } else {
      origin <- "Gate"
    }

    # handle different gatetpyes??
    #xGateLim <- as.numeric(xml2::xml_attrs(xml2::xml_child(xml2::xml_child(g), 1)))
    #yGateLim <- as.numeric(xml2::xml_attrs(xml2::xml_child(xml2::xml_child(g), 2)))

    gate_id <- xml2::xml_attr(gg[n], "id")
    parentgate_id <- xml2::xml_attr(gg[n], "parent_id")
    eventsInside <- xml2::xml_attr(xml2::xml_child(gg[n]), "eventsInside")

    return(data.frame(FileName = FileName,
                      PopulationFullPath = PopulationFullPath,
                      Parent = Parent,
                      Population = Population,
                      Count = Count,
                      ParentCount = ParentCount,
                      xDim = xDim,
                      yDim = yDim,
                      gate_id = gate_id,
                      parentgate_id = parentgate_id,
                      eventsInside = eventsInside,
                      sampleID = sampleID,
                      FilePath = gsub("^file:", "", FilePath),
                      gate_level = gate_level,
                      origin = origin,
                      n = n)
    )
  })

  roots <- do.call(rbind, lapply(xml2::xml_children(xml2::xml_child(ws, "SampleList")), function(y) {
    data.frame(FileName = basename(xml2::xml_attr(xml2::xml_child(y, "DataSet"), "uri")),
               PopulationFullPath = "root",
               Parent = NA,
               Population = "root",
               Count = xml2::xml_attr(xml2::xml_child(y, "SampleNode"), "count"),
               ParentCount = NA,
               xDim = NA,
               yDim = NA,
               gate_id = NA,
               parentgate_id = NA,
               eventsInside = NA,
               sampleID = xml2::xml_attr(xml2::xml_child(y, "DataSet"), "sampleID"),
               FilePath = gsub("^file:", "", xml2::xml_attr(xml2::xml_child(y, "DataSet"), "uri")),
               gate_level = 0,
               origin = "root",
               n = 0)
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
      browser()
      stop("PopulationFullPaths not unique which cannot or should not be. Check.")
    }
    y$PopulationFullPath
  }))
  auto_paths <- lapply(full_paths, function(y) shortest_unique_path(y))

  for (y in seq_along(gates_list)) {
    gates_list[[y]][["Population"]] <- auto_paths[[which(sapply(full_paths, function(z) identical(z,  gates_list[[y]][["PopulationFullPath"]])))]]
  }
  gates_out <- do.call(rbind, gates_list)

  gates_out <- dplyr::left_join(gates_out, wsp_xml_get_groups(ws), by = "sampleID")
  gates_out[,"ws"] <- basename(xml2::xml_attr(ws, "nonAutoSaveFileName"))
  gates_out <- gates_out[order(gates_out$FileName),]
  rownames(gates_out) = seq(1,nrow(gates_out),1)
  return(gates_out)
}


wsp_xml_get_samples <- function(x) {
  if (is.character(x)) {
    x <- xml2::read_xml(x)
  }
  if (!any(class(x) == "xml_document")) {
    stop("x must be a xml-document or a character path to its location on disk")
  }
  s <- as.data.frame(t(sapply(xml2::xml_children(xml2::xml_child(x, "SampleList")), function(y) {
    xml2::xml_attrs(xml2::xml_child(y, "DataSet"))
  })))
  names(s) <- c("FilePath", "sampleID")
  s$FilePath <- gsub("file:", "", s$FilePath)
  s$FileName <- basename(s$FilePath)
  return(s)
}

wsp_xml_get_groups <- function(x) {
  if (is.character(x)) {
    x <- xml2::read_xml(x)
  }
  if (!any(class(x) == "xml_document")) {
    stop("x must be a xml-document or a character path to its location on disk")
  }

  g <- sapply(xml2::xml_children(xml2::xml_child(x, "Groups")), function(y) {
    xml2::xml_attrs(y)[["name"]]
  })
  gs <- lapply(xml2::xml_children(xml2::xml_child(x, "Groups")), function(y) {
    unlist(xml2::xml_attrs(xml2::xml_children(xml2::xml_child(xml2::xml_child(y, "Group"), "SampleRefs"))))
  })
  gr <- data.frame(group = rep(g, lengths(gs)),  sampleID = unlist(gs))
  gr <- do.call(rbind, lapply(unique(gr$sampleID), function(y) {
    if (nrow(gr[which(gr$sampleID == y),]) > 1) {
      g <- gr[base::intersect(which(gr$sampleID == y), which(gr$group != "All Samples")), ]
      g <- data.frame(group = paste(g$group, collapse = ", "), sampleID = g$sampleID)
    }
  }))
  return(gr)
}

shortest_unique_path <- function(p) {
  p_rev <- sapply(strsplit(p, "/"), rev)
  p_rev <- lapply(seq_along(p_rev), function(x) {
    i<-1
    while (any(sapply(p_rev[-x], function(y) {
      identical(p_rev[[x]][1:i], y[1:i])
    }))) {
      i<-i+1
    }
    return(p_rev[[x]][1:i])
  })
  p <- sapply(sapply(p_rev, rev), function(x) paste(x, collapse = "/"))
  return(p)
}

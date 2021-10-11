#' Convenient function to import population counts from a flowjo workspace
#'
#' In comparison to ws_get_popstats another backend (xml2) is used to read data from a flowjo .wsp-file (under the hood an .xml-file).
#' As is should be quick to extract data no groups have to be selected. Filtering has to be conducted afterwards.
#' Also, no restriction to gating trees within a group are made and no folder of FCS has to be provided.
#'
#' @param ws character of path to flowjo workspace
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

  if (class(ws) != "character") {
    stop("Provide a path to wsp file.")
  }
  if (!file.exists(ws)) {
    stop("wsp file not found.")
  }

  samples <- wsp_xml_get_samples(ws)
  groups <- wsp_xml_get_groups(ws)
  gates <- wsp_xml_get_gates(ws)

  ##fix
  table <- do.call(rbind, lapply(1:length(s), function(x) {

    data.frame(FileName = unname(s0[[x]]["name"]),
               PopulationFullPath = c("root", gates[[x]]$PopulationFullPath),
               Parent = stats::lag(gates[[x]]$PopulationFullPath, k=1),
               Population = unname(shortest_unique_path(gates[[x]]$PopulationFullPath)),
               count = as.numeric(unname(c(s0[[x]]["count"], gates[[x]]$count))),
               ParentCount = as.numeric(stats::lag(unname(c(s0[[x]]["count"], gates[[x]]$count)), k=1)),
               FractionOfParent = as.numeric(unname(c(s0[[x]]["count"], sapply(gates[[x]], "[", "count"))))/as.numeric(stats::lag(unname(c(s0[[x]]["count"], sapply(gates[[x]], "[", "count"))), k=1)),
               sampleID = unname(s0[[x]]["sampleID"]))
  }))
  table <- merge(table, groups, by = "sampleID")
  table <- merge(table, samples, by = c("sampleID", "FileName"))
  table <- merge(table, d, by = c("FileName", "PopulationFullPath"))
  table <- table[,which(names(table) != "sampleID")]
  table[,"ws"] <- basename(ws)

  return(table)
}


wsp_xml_get_roots <- function(x) {
  if (is.character(x)) {
    x <- xml2::read_xml(x)
  }
  if (!any(class(x) == "xml_document")) {
    stop("x must be a xml-document or a character path to its location on disk")
  }
  lapply(xml_children(xml2::xml_child(x, "SampleList")), function(y) {
    xml2::xml_attrs(xml2::xml_child(y, "SampleNode"))
  })
}


wsp_xml_get_gates <- function(x) {
  if (is.character(x)) {
    x <- xml2::read_xml(x)
  }
  if (!any(class(x) == "xml_document")) {
    stop("x must be a xml-document or a character path to its location on disk")
  }

  s0 <- lapply(xml2::xml_children(xml2::xml_child(x, "SampleList")), function(y) {
    xml_attrs(xml_child(y, "SampleNode"))
  })

  s <- lapply(xml2::xml_children(xml2::xml_child(x, "SampleList")), function(y) {
    b <- xml2::xml_attrs(xml2::xml_find_all(y, ".//Population"))
    fp <- unlist(lapply(xml2::xml_find_all(y, ".//Population"), function(z) {
      par <- unlist(lapply(xml2::xml_parents(z), function(a) {
        return(xml2::xml_attr(a, "name"))
      }))
      return(paste(rev(par[which(!is.na(par))]), collapse = "/"))
    }))
    fp <- gsub(fp[1], "root", fp)

    for (i in seq_along(b)) {
      b[[i]] <- c(b[[i]], "parents" = fp[i])
      b[[i]] <- c(b[[i]], "PopulationFullPath" = gsub("^/", "", paste(c(b[[i]]["parents"], b[[i]]["name"]), collapse = "/")))
    }
    return(b)
  })


  id <- lapply(xml2::xml_children(xml2::xml_child(x, "SampleList")), function(y) {
    xml2::xml_attrs(xml2::xml_find_all(y, ".//Gate"))
  })

  dims <- lapply(xml2::xml_children(xml2::xml_child(x, "SampleList")), function(z) {
    list(xDim = sapply(xml2::xml_find_all(z, ".//Gate"), function(y) {xml2::xml_attrs(xml2::xml_child(xml2::xml_child(xml2::xml_child(y), 1), 1))}),
         yDim = sapply(xml2::xml_find_all(z, ".//Gate"), function(y) {xml2::xml_attrs(xml2::xml_child(xml2::xml_child(xml2::xml_child(y), 2), 1))}))
  })

  out <- lapply(seq_along(s), function(y) {
    r <- data.frame(PopulationFullPath = unname(c("root", gsub("root/", "", sapply(s[[y]], "[", "PopulationFullPath")))),
                    Population = unname(c("root", sapply(s[[y]], "[", "name"))),
                    Count = as.numeric(unname(c(s0[[y]][["count"]], sapply(s[[y]], "[", "count")))),
                    gate_id = c("root", sapply(id[[y]], "[", "id")),
                    parentgate_id = c(NA, "root", sapply(id[[y]], "[", "parent_id")[-1]),
                    xDim = c(NA, dims[[y]][["xDim"]]),
                    yDim = c(NA, dims[[y]][["yDim"]]))
    rr <- dplyr::left_join(r, setNames(r[,c(3,4)], c("ParentCount", "gate_id")), by = c("parentgate_id"="gate_id"))
    rr <- dplyr::left_join(r, setNames(r[,c(1,4)], c("Parent", "gate_id")), by = c("parentgate_id"="gate_id"))
    return(rr)
  })

  return(out)
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

  g <- sapply(xml2::xml_children(xml2::xml_child(x, "Groups")), function(y){
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

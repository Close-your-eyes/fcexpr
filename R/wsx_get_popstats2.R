#' Import data from a wsp file
#'
#' Flowjo wsp files contain many information like gated event counts, statistics and keywords
#' from FCS files. These may be accessed without a dongle and can be read completely independent of
#' the respective FCS files once the gating has been conducted.
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#' @param return_stats logical (T,F) whether to return statistics next to cells counts
#' @param groups vector of flowjo group names to consider
#' @param invert_groups logical whether to exclude the selected groups
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
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
wsx_get_popstats2 <- function(ws,
                              groups = NULL,
                              invert_groups = F,
                              return_stats = T,
                              lapply_fun = lapply,
                              strip_data = T,
                              ...) {

  # ws <- "/Volumes/CMS_SSD_2TB/Experiment_data/20230605_Blood_cohort2/FJ_workspaces/ExpPart_6_for_pub.wsp"

  ## allow to pass mclapply
  lapply_fun <- match.fun(lapply_fun)
  ws <- check_ws(ws)
  group_df <- wsx_get_groups(ws, collapse_groups = F)
  if (is.null(groups)) {
    groups <- unique(group_df[,"FlowJoGroup", drop=T])
  }
  if (invert_groups) {
    group_df <- group_df[which(!group_df[,"FlowJoGroup", drop=T] %in% groups),]
  } else {
    group_df <- group_df[which(group_df[,"FlowJoGroup", drop=T] %in% groups),]
  }

  if (nrow(group_df) == 0) {
    stop("Non of provided groups found.")
  }

  # each sample which may be in multiple groups only considered once
  ids <- unique(group_df[,"sampleID",drop=T])
  rel_nodes <- xml2::xml_children(xml2::xml_child(ws, "SampleList"))
  rel_nodes <- rel_nodes[which(purrr::map(rel_nodes, function(x) xml2::xml_attrs(xml2::xml_child(x, "DataSet"))[["sampleID"]]) %in% ids)]

  # in gg each sample is an own nodeset with gates that only belong to that sample
  # go sample wise
  library(xml2)
  library(ggraph)
  gg <- xml2::xml_find_all(rel_nodes, ".//Gate|.//Dependents", flatten = FALSE)
  gg_meta <- purrr::flatten(xml2::xml_find_all(rel_nodes, "SampleNode", flatten = FALSE))
  gg_meta_names <- do.call(rbind, purrr::map(gg_meta, xml2::xml_attrs))
  names(gg) <- gg_meta_names[,"name"]

  # create data with gate ids
  gate_ids <- purrr::map(gg, xml2::xml_attrs)
  gate_ids_df <- purrr::map_dfr(gate_ids, function(x) as.data.frame(do.call(rbind, x)), .id = "FileName")
  gate_ids_df$parent_id <- ifelse(gate_ids_df$parent_id == gate_ids_df$id, paste0("root_", gate_ids_df$FileName), gate_ids_df$parent_id)

  # pull all counts with and associated ids
  gate_details <- purrr::map(gg, xml2::xml_parents)
  counts <- purrr::map_dfr(gate_details, get_count_and_ids)
  gate_ids_df <- dplyr::left_join(gate_ids_df, counts)

  ## missing yet: root counts
  ## next: follow graph to derive population full paths

  # make graph
  gate_graph <- igraph::graph_from_data_frame(data.frame(from = gate_ids_df$parent_id,
                                                         to = gate_ids_df$id),
                                              directed = F)

  ggraph::ggraph(ggraph::create_layout(gate_graph, layout = "tree")) +
    geom_node_point() +
    geom_edge_link()






  gates <- lapply_fun(gg, function(n) {

    prnts <- xml2::xml_parents(n)

    s_node <- prnts[which(xml2::xml_name(prnts) == "Sample")]
    sampleID <- xml2::xml_attr(xml2::xml_child(s_node, "DataSet")[[1]], "sampleID")

    FilePath <- gsub("^file:", "", xml2::xml_attr(xml2::xml_child(s_node, "DataSet")[[1]], "uri"))
    FileName <- basename(FilePath)

    p_nodes <- prnts[which(xml2::xml_name(prnts) %in% c("AndNode", "OrNode", "NotNode", "Population"))]
    PopulationFullPath <- paste(rev(xml2::xml_attr(p_nodes, "name")), collapse = "/")
    Parent <- if (PopulationFullPath == basename(PopulationFullPath)) {"root"} else {dirname(PopulationFullPath)}
    Population <- basename(PopulationFullPath)

    Count <- xml2::xml_attr(p_nodes[1], "count")
    if (Count == -1) {
      stop("Count = -1 detected. One or more nodes a boolean gate (Or/And) depends may not have been found.
            Cannot derive correct Count.
            Did you rename the nodes an Or- or And-Gate depends on?
                 If so, please re-define the respective boolean gate.")
    }

    ParentCount <- if (length(p_nodes) > 1) {xml2::xml_attr(p_nodes[2], "count")} else {xml2::xml_attr(xml2::xml_child(s_node, "SampleNode"), "count")}
    gate_level <- length(p_nodes)

    xDim <- tryCatch({
      xml2::xml_attr(xml2::xml_child(xml2::xml_child(xml2::xml_child(n), 1)), "name")
    }, error = function(e) {
      NA
    })

    yDim <- tryCatch({
      xml2::xml_attr(xml2::xml_child(xml2::xml_child(xml2::xml_child(n), 2)), "name")
    }, error = function(e) {
      NA
    })

    if (xml2::xml_name(n) == "Dependents") {
      origin <- "Dependents"
    } else {
      origin <- "Gate"
    }

    # deps for Or, And or NotNodes; correct counts afterwards
    #deps <- xml2::xml_child(p_nodes[1], "Dependents")
    #deps <- list(xml2::xml_attr(xml2::xml_children(deps), "name"))

    gate_id <- xml2::xml_attr(n, "id")
    parentgate_id <- xml2::xml_attr(n, "parent_id")
    eventsInside <- xml2::xml_attr(xml2::xml_child(n), "eventsInside")

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
                      #deps = I(deps),
                      #n = n,
                      stringsAsFactors = F)
    )
  }) #, ...

  roots <- do.call(rbind, lapply_fun(rel_nodes, function(y) {
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
               #deps = I(list(character(0))),
               #n = 0,
               stringsAsFactors = F)
  })) #, ...

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

  # if any sample is in at least two groups, the group column becomes a list
  if (any(table(group_df$sampleID) > 1)) {
    group_df <- dplyr::summarise(dplyr::group_by(group_df, sampleID), FlowJoGroup = list(FlowJoGroup))
  }
  gates_out <- do.call(rbind, gates_list)
  gates_out <- dplyr::left_join(gates_out, group_df, by = "sampleID") # ...
  gates_out[,"ws"] <- basename(xml2::xml_attr(ws, "nonAutoSaveFileName"))
  gates_out <- gates_out[order(gates_out$FileName, gates_out$gate_level, factor(gates_out$origin, levels = c("root", "Gate", "Dependents"))),]
  rownames(gates_out) = seq(1,nrow(gates_out),1)

  if (strip_data) {
    gates_out <- gates_out[,which(!names(gates_out) %in% c("gate_id", "parentgate_id", "sampleID", "origin", "n", "gate_level"))]
  }


  if (return_stats) {
    stats_out <- do.call(rbind, lapply_fun(seq_along(rel_nodes), function(n) {
      node <- rel_nodes[n]
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
    }, ...))
    return(list(counts = gates_out, stats = stats_out))
  }
  return(gates_out)
}

get_count_and_ids <- function(nodeset) {
  ## how to do this more elegantly?
  attr_list <- purrr::map(nodeset, xml2::xml_attrs)
  inds <- which(sapply(purrr::map(attr_list, names), "[", 1) == "name")
  attr_list <- attr_list[inds]
  nodeset <- nodeset[inds]
  attr_list[[which(lengths(attr_list) == 7)]] <- attr_list[[which(lengths(attr_list) == 7)]][-7]
  df <- as.data.frame(do.call(rbind, attr_list))
  df$id <- purrr::map_chr(purrr::map(nodeset, xml2::xml_child, search = "Gate"), xml2::xml_attr, attr = "id")
  df$id <- ifelse(is.na(df$id), paste0("root_", df$name), df$id)

  return(df)
}


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
                              strip_data = T,
                              more_gate_data = F,
                              show_progress = F,
                              ...) {


  # check stats

  groups = NULL
  invert_groups = F
  more_gate_data = T
  show_progress = T
  ws <- check_ws(ws)
  group_df <- wsx_get_groups(ws, collapse = NULL)
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


  # jannis wsp: multiple samples from same FCS file

  # each sample which may be in multiple groups is only considered once
  ids <- unique(group_df[,"sampleID",drop=T])
  samples <- xml2::xml_children(xml2::xml_child(ws, "SampleList"))
  sample_ids <- which(sapply(xml2::xml_attrs(xml2::xml_child(samples, "DataSet")), "[[", "sampleID") %in% ids)
  if (anyDuplicated(basename(sapply(xml2::xml_attrs(xml2::xml_child(samples, "DataSet")), "[[", "uri"))) != 0) {
    message("Duplicate filenames detected.")
  }
  samples <- samples[sample_ids]

  file_paths <- as.data.frame(do.call(dplyr::bind_rows, xml2::xml_attrs(xml2::xml_child(samples, "DataSet"))))
  file_paths$uri <- gsub("^file:", "", file_paths$uri)

  file_names <- purrr::map_chr(xml2::xml_children(xml2::xml_child(ws, "SampleList")), function(x) xml2::xml_attrs(xml2::xml_child(x, "SampleNode"))[["name"]])
  keys_list <- purrr::map(setNames(xml2::xml_children(xml2::xml_child(ws, "SampleList")), ids), function(x) {
    keys <- xml2::xml_attrs(xml2::xml_contents(xml2::xml_child(x, "Keywords")))
    keys <- stats::setNames(sapply(keys, "[", 2), sapply(keys, "[", 1))
    keys <- utils::stack(keys)
    names(keys) <- c("value", "name")
    keys <- keys[,c(2,1)]
    keys$name <- as.character(keys$name)
    return(keys)
  })

  return_dup_files <- T
  dup_files <- purrr::map_dfr(keys_list, function(x) x[which(grepl("\\$FIL|\\$DATE|\\$TOT|\\$BTIM|\\$ETIM", x$name)),], .id = "sampleID")
  dup_files <- tidyr::pivot_wider(dup_files, names_from = "name", values_from = "value")
  dup_files <- tidyr::drop_na(dup_files)
  if (nrow(dup_files) > 1) {
    dup_files <- dplyr::group_by(dup_files, !!!rlang::syms(c("$FIL","$DATE","$TOT","$BTIM","$ETIM")))
    dup_files <- dplyr::mutate(dup_files, group_id = dplyr::cur_group_id())
    dup_files <- dplyr::ungroup(dup_files)
    dup_files <- dplyr::left_join(dup_files, file_paths, by = "sampleID")
    dup_table <- table(dup_files$group_id)
    if (any(dup_table) > 1) {
      return_dup_files <- T
      warning("One or multiple FCS files seem to be duplicates. A data.frame with details will be returned.")
      print(dup_files)
    }
  }
  names(keys_list) <- file_names

  # in samples each sample is an own nodeset with gates that only belong to that sample
  gates <- xml2::xml_find_all(samples, ".//Gate", flatten = F) # .//Dependents
  names(gates) <- do.call(dplyr::bind_rows, xml2::xml_attrs(xml2::xml_find_all(samples, "SampleNode", flatten = T)))[,"name",drop=T]
  #depends <- xml2::xml_find_all(samples, ".//Dependents", flatten = F)
  #names(depends) <- do.call(dplyr::bind_rows, xml2::xml_attrs(xml2::xml_find_all(samples, "SampleNode", flatten = T)))[,"name"]

  # !!!!!
  # do.call(rbind, x) fills missing values (e.g. missing parent_id) with value id, this is unexpected and a mistake!
  # use do.call(dplyr::bind_rows, x) which does it corectly and insert NA for missing values
  # !!!!!
  # children of OrNodes or AndNodes gates have no parent_id as the OrNodes or AndNodes has no id; this feature is (no parent_id) is shared with the root
  # create data with gate ids
  pop_df <- purrr::map_dfr(purrr::map(gates, xml2::xml_attrs), function(x) as.data.frame(do.call(dplyr::bind_rows, x)), .id = "FileName")
  pop_df$GateType <- unlist(purrr::map(gates, purrr::compose(xml2::xml_name, xml2::xml_children)))

  if (more_gate_data) {
    GateDef <- purrr::map(gates, xml2::xml_children)
    GateDef <- purrr::map(GateDef, xml2::as_list)
    pop_df$GateDef <- purrr::list_flatten(GateDef)
  }

  # make sure root counts are joined below: manually add the root 'gate' for each sample
  pop_df <- add_root_node(pop_df)

  # pull all counts with and associated ids
  # nodeset <- purrr::map(gates, xml2::xml_parents)[[1]]
  # nodeset <- purrr::map(gates, xml2::xml_parents)[["Compensation Controls_CD19 APC Stained Control.fcs"]]
  node_details_list <- purrr::map(.x = purrr::map(gates, xml2::xml_parents),
                                  .f = get_node_details2,
                                  .progress = show_progress,
                                  more_gate_data = more_gate_data)
  node_details_df <- bind_rows_chunked(df_list = sapply(node_details_list, "[", "df"))
  if (anyDuplicated(node_details_df[which(!is.na(node_details_df$id)), "id"]) != 0) {
    stop("Duplicate gate id detected. FlowJo wsp needs fixing?!")
  }

  # join in this way; pop_df has no And-/Or-Depend Gates but node_details_df has
  pop_df <- dplyr::left_join(node_details_df, pop_df, by = c("id", "FileName"))

  # from OrNodes and AndNodes in node_details_list: find all children of these nodes. the remaining row in pop_df must have the samples root as parent

  # assign root as parent_id
  #(i) name not in ornodenames and andnodenames
  #(ii) !is.na(id) - these the OrNodes and AndNodes
  #(iii) name_root != "root" - this is the root itself and has not parent_id
  # could be made simpler if there are no OrNodes and AndNodes (then there is only two types nodes without parent_id: the root itself and its direct children)

  if (any(!sapply(sapply(node_details_list, "[", "OrNodes"), is.null)) || any(!sapply(sapply(node_details_list, "[", "AndNodes"), is.null))) {
    #  || any(!sapply(sapply(node_details_list, "[", "NotNodes"), is.null))
    # filter for !is.null first
    for (i in names(node_details_list)) {
      # one entry for each OrNode or AndNode
      temp_nodes <- node_details_list[[i]][["OrNodes"]]
      temp_nodes <- temp_nodes[which(!sapply(temp_nodes, is.null))]
      ornodenames <- NULL
      if (length(temp_nodes) > 0) {
        temp <- xml2::xml_child(temp_nodes, "Subpopulations")
        temp <- temp[which(!is.na(temp))]
        ornodenames <- unlist(lapply(temp, function(x) {
          xml2::xml_attr(xml2::xml_children(x), "name")
        }))
      }

      temp_nodes <- node_details_list[[i]][["AndNodes"]]
      temp_nodes <- temp_nodes[which(!sapply(temp_nodes, is.null))]
      andnodenames <- NULL
      if (length(temp_nodes) > 0) {
        temp <- xml2::xml_child(temp_nodes, "Subpopulations")
        temp <- temp[which(!is.na(temp))]
        andnodenames <- unlist(lapply(temp, function(x) {
          xml2::xml_attr(xml2::xml_children(x), "name")
        }))
      }

      # needed?
      '      temp_nodes <- node_details_list[[i]][["NotNodes"]]
      temp_nodes <- temp_nodes[which(!sapply(temp_nodes, is.null))]
      notnodenames <- NULL
      if (length(temp_nodes) > 0) {
        temp <- xml2::xml_child(temp_nodes, "Subpopulations")
        temp <- temp[which(!is.na(temp))]
        notnodenames <- unlist(lapply(temp, function(x) {
          xml2::xml_attr(xml2::xml_children(x), "name")
        }))
      }'

      inds <- Reduce(intersect, list(which(!pop_df$name %in% c(ornodenames, andnodenames)),
                                     which(!is.na(pop_df$id)),
                                     which(pop_df$name_root != "root"),
                                     which(is.na(pop_df$parent_id)),
                                     which(pop_df$FileName == i)))
      pop_df[inds,"parent_id"] <- paste0("root_", pop_df[inds,"FileName"])
    }
  } else {
    inds <- Reduce(intersect, list(which(pop_df$name_root != "root"),
                                   which(is.na(pop_df$parent_id))))
    pop_df[inds,"parent_id"] <- paste0("root_", pop_df[inds,"FileName"])
  }

  ## error below here with
  # ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Complicated_OrAndGates_OrGate_at_diff_hierachies_sameGatingTree.wsp"

  # to make a complete graph, all id and parent_id have to be assigned (also of OrNodes and AndNodes)
  # but how to match unambiguously without PopulationFullPath? - use name and count (see add_boolean_gate_data)
  ## only run this if OrNodes and/or AndNodes exist
  if (any(!sapply(sapply(node_details_list, "[", "OrNodes"), is.null)) || any(!sapply(sapply(node_details_list, "[", "AndNodes"), is.null)) || any(!sapply(sapply(node_details_list, "[", "NotNodes"), is.null))) {
    # check for if one or the other is null #######
    if (any(!sapply(sapply(node_details_list, "[", "OrNodes"), is.null))) {
      pop_df <- add_boolean_gate_data(df = pop_df,
                                      node_details_list = node_details_list,
                                      nodes_name = "OrNodes",
                                      more_gate_data = more_gate_data)
    }
    if (any(!sapply(sapply(node_details_list, "[", "AndNodes"), is.null))) {
      pop_df <- add_boolean_gate_data(df = pop_df,
                                      node_details_list = node_details_list,
                                      nodes_name = "AndNodes",
                                      more_gate_data = more_gate_data)
    }
    # do not nodes at last, because notnodes stemming from OrNodes/AndNodes - those OrNodes/AndNodes need an id themselves first (before they are NA)
    # why can the be nodesets of length zero and not NULL?
    if (any(!sapply(sapply(node_details_list, "[", "NotNodes"), is.null)) & any(sapply(sapply(node_details_list, "[", "NotNodes"), length) > 0)) {
      pop_df <- add_boolean_gate_data(df = pop_df,
                                      node_details_list = node_details_list,
                                      nodes_name = "NotNodes",
                                      more_gate_data = more_gate_data)
    }

    # add "not" to id for duplicated ids from not nodes
    dup_id <- unique(pop_df$id[duplicated(pop_df$id)])
    inds <- intersect(which(pop_df$id %in% dup_id), which(pop_df$NodeType == "NotNode"))
    pop_df[inds,"id"] <- paste0(pop_df[inds,"id"], "_Not")

    # when NotNodes, somehow from same OrNodes/AndNodes still have the id - make it unique; all correct?
    # check somewhen. but this only occurs in kind of exotic gating strategies. will not cause problems, regularly
    #"/Volumes/CMS_SSD_2TB/example_workspaces/Multiple_OrNodes_AndNodes_NotNode_on_OrAndNodes_sameDims_sameGatingTrees.wsp"
    pop_df$id <- make.unique(pop_df$id)
  }

  # add channel desc here, after boolean gates have been enriched by data, including xChannel, yChannel
  if (more_gate_data) {
    pop_df <- add_channel_desc(df = pop_df, ws = ws, keys_list = keys_list)
  }

  ## next: follow graph to derive population full paths
  # make graph
  # find end nodes (vertices) by checking degree (number of outgoing edges); graph has to be directed
  gate_graph <- igraph::graph_from_data_frame(data.frame(from = pop_df[which(!is.na(pop_df$parent_id)), "parent_id"],
                                                         to = pop_df[which(!is.na(pop_df$parent_id)), "id"]), directed = T)

  ## maybe provide separate graphs with respective end edges? speed?
  edge_degrees <- igraph::degree(gate_graph, mode = "out")
  end_edges <- edge_degrees[which(edge_degrees == 0)]


  # adds grandparent_id, PopulationFullPath, PopulationFullPathID, GateDepth
  pop_df <- add_full_paths(df = pop_df, graph = gate_graph, end_edges = names(end_edges), show_progress = show_progress) # takes a bit with many samples

  # TODO: add more info to vertices?!
  igraph::V(gate_graph)$label <- pop_df[match(names(igraph::V(gate_graph)), unique(pop_df$id)),"name_root"]
  igraph::V(gate_graph)$parent_id <- pop_df[match(names(igraph::V(gate_graph)), unique(pop_df$id)),"parent_id"]
  # get one graph for each sample
  graph_subgroups <- igraph::components(gate_graph)
  subgroup_vertex_id_list <- split(names(graph_subgroups$membership), graph_subgroups$membership)
  gate_graph_samples <- lapply(stats::setNames(subgroup_vertex_id_list, gsub("root_", "", sapply(subgroup_vertex_id_list, "[", 1))), function(x) igraph::subgraph(graph = gate_graph, vids = x))


  pop_df <- add_parent_count(pop_df) # add parent count
  pop_df <- add_grandparent_count(pop_df) # add grandparent count
  pop_df <- add_total_count(pop_df)

  # sample wise of for whole pop_df at once?
  auto_paths <- shortest_unique_path(unique(pop_df$PopulationFullPath)[which(unique(pop_df$PopulationFullPath) != "")])
  names(auto_paths) <- unique(pop_df$PopulationFullPath)[which(unique(pop_df$PopulationFullPath) != "")]
  pop_df$Population <- auto_paths[pop_df$PopulationFullPath]
  pop_df$Population[which(is.na(pop_df$Population))] <- pop_df$name_root[which(is.na(pop_df$Population))]
  pop_df$PopulationFullPath[which(pop_df$Population == "root")] <- "root"

  pop_df$FlowJoWsp <- basename(xml2::xml_attr(ws, "nonAutoSaveFileName"))

  # harmonize with previous? compare? - with dplyr::anti_join
  # notify of number of different gating trees - how to do with igraph?

  # sampleID - done before
  # FilePath
  #file_paths <- as.data.frame(do.call(dplyr::bind_rows, xml2::xml_attrs(xml2::xml_child(samples, "DataSet"))))
  #file_paths$uri <- gsub("^file:", "", file_paths$uri)
  names(file_paths)[1] <- "FilePath"
  pop_df <- dplyr::left_join(pop_df, group_df, by = c("sampleID"))
  pop_df <- dplyr::left_join(pop_df, file_paths, by = c("sampleID"))

  # fill whitespaces with NA

  cols <- if (more_gate_data) {
    c("FileName",
      "PopulationFullPath",
      "PopulationFullPathID",
      "Population",
      "Count",
      "ParentCount",
      "GrandparentCount",
      "TotalCount",
      "FractionOfParent",
      "FractionOfGrandparent",
      "FractionOfTotal",
      "xChannel",
      "yChannel",
      "xDesc",
      "yDesc",
      "GateDef",
      "GateType",
      "NodeType",
      "eventsInside",
      "GateDepth",
      "id",
      "parent_id",
      "grandparent_id",
      "sampleID",
      "FlowJoGroup",
      "FlowJoWsp",
      "FilePath")
  } else {
    c("FileName",
      "PopulationFullPath",
      "PopulationFullPathID",
      "Population",
      "Count",
      "ParentCount",
      "GrandparentCount",
      "TotalCount",
      "FractionOfParent",
      "FractionOfGrandparent",
      "FractionOfTotal",
      "GateType",
      "NodeType",
      "GateDepth",
      "id",
      "parent_id",
      "grandparent_id",
      "sampleID",
      "FlowJoGroup",
      "FlowJoWsp",
      "FilePath")
  }
  # order rows - how?
  pop_df <- dplyr::group_by(pop_df, FileName)
  pop_df <- dplyr::arrange(pop_df, GateDepth, .by_group = T)
  pop_df <- dplyr::ungroup(pop_df)

  if (return_dup_files) {
    return(list(counts = pop_df[,cols], graph = gate_graph, graph_sample = gate_graph_samples, duplicate_FCS_files = dup_files))
  } else {
    return(list(counts = pop_df[,cols], graph = gate_graph, graph_sample = gate_graph_samples))
  }
  # stats

  # matrix/df with all parent gates as ref and FractionOfXXX


  stats <- xml2::xml_find_all(ws, ".//Statistic")

  '  if (strip_data) {
    gates_out <- gates_out[,which(!names(gates_out) %in% c("gate_id", "parentgate_id", "sampleID", "origin", "n", "gate_level"))]
  }'

  'if (return_stats) {
    stats_out <- do.call(dplyr::bind_rows, lapply_fun(seq_along(rel_nodes), function(n) {
      node <- rel_nodes[n]
      stats <- xml2::xml_find_all(node, ".//Statistic")

      stats_df <- do.call(dplyr::bind_rows, lapply(stats, function(x) {
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
  return(gates_out)'

}


'  ggraph::ggraph(ggraph::create_layout(gate_graph_samples[[2]], layout = "tree")) +
    ggraph::geom_edge_link() +
    ggraph::geom_node_point(ggplot2::aes(color = label), size = 2) +
    ggraph::geom_node_text(ggplot2::aes(label = name)) +
    ggplot2::theme_void() +
    ggplot2::scale_color_manual(values = fcexpr::col_pal("custom"))'

'  gate_graphs <- lapply(unique(pop_df$FileName), function(x) {
    pop_df2 <- dplyr::filter(pop_df, FileName == x)
    igraph::graph_from_data_frame(data.frame(from = pop_df2[which(!is.na(pop_df2$parent_id)), "parent_id"],
                                             to = pop_df2[which(!is.na(pop_df2$parent_id)), "id"]), directed = T)
  })

  igraph::subgraph_isomorphic(gate_graphs[[2]], gate_graphs[[1]])
  igraph::isomorphic(gate_graphs[[2]], gate_graphs[[3]])'

'  igraph::V(gate_graph)$label <- pop_df[match(names(igraph::V(gate_graph)), unique(pop_df$id)),"name_root"]
  ggraph::ggraph(ggraph::create_layout(gate_graph, layout = "tree")) +
    ggraph::geom_edge_link() +
    ggraph::geom_node_point(ggplot2::aes(color = label), size = 2) +
    ggraph::geom_node_text(ggplot2::aes(label = name)) +
    ggplot2::theme_void() +
    ggplot2::scale_color_manual(values = fcexpr::col_pal("custom"))'

'ggraph::ggraph(ggraph::create_layout(gate_graph_samples[[2]], layout = "tree")) +
  ggraph::geom_edge_link() +
  ggraph::geom_node_point() +
  ggraph::geom_node_text(ggplot2::aes(label = name))'


get_node_details2 <- function(nodeset, more_gate_data = F) {

  temp_attr_list <- purrr::map(nodeset, xml2::xml_attrs)
  temp_attr_list_names <- sapply(temp_attr_list, names)
  temp_attr_list <- temp_attr_list[which(grepl("sampleID", temp_attr_list_names))]
  sampleID <- temp_attr_list[[1]][7]
  temp_attr_list <- list(temp_attr_list[[1]][-7])
  # sampleNode becomes last row in df

  pops <- xml2::xml_find_all(nodeset, "Population|AndNode|OrNode|NotNode") # do not omit NotNode here
  gate_list <- xml2::xml_find_all(pops, ".//Gate|.//Dependents", flatten = T)

  # get parents to maintain correct order of id and associated gates in df
  gate_list_par <- xml2::xml_parent(gate_list)
  attr_list <- c(purrr::map(gate_list_par, xml2::xml_attrs), temp_attr_list)
  df <- as.data.frame(do.call(dplyr::bind_rows, attr_list))

  id <- c(xml2::xml_attr(gate_list, attr = "id"), NA) # add NA for SampleNode at the end
  node_types <- c(xml2::xml_name(gate_list_par), "SampleNode")

  # sometimes, boolean gates are returned duplicated, but only of them has an ID
  # in this case ID and node_types are of different lengths which causes error
  # filter for indices (rows), that are duplicates but have an ID
  if (length(id) != length(node_types)) {
    gate_list_par2 <- lapply(gate_list, xml2::xml_parent)
    gate_list_par2_attrs <- dplyr::bind_rows(lapply(gate_list_par2, xml2::xml_attrs))
    gate_list_par2_attrs$id <- xml2::xml_attr(gate_list, attr = "id")
    gate_list_par2_attrs <- dplyr::group_by(gate_list_par2_attrs, !!!rlang::syms(names(gate_list_par2_attrs)[1:6]))
    gate_list_par2_attrs <- dplyr::mutate(gate_list_par2_attrs, n = dplyr::n())
    gate_list_par2_attrs <- dplyr::filter(gate_list_par2_attrs, !(is.na(id) & n > 1))
    gate_list_par2_attrs <- dplyr::ungroup(gate_list_par2_attrs)
    #gate_list_par2_attrs$node_type <- node_types[-length(node_types)] # just for checking
    id <- c(gate_list_par2_attrs$id, NA) # add NA for SampleNode at the end
  }

  ### NotNode - yes no ?? if not node is from OrNode or AndNode then it should be added but if NotNode is normal, do not add?
  node_select <- node_types %in% c("SampleNode", "OrNode", "AndNode")
  node_select[which(node_types == "NotNode" & is.na(id))] <- T
  id <- id[which(!is.na(id))]
  # "NotNode"
  for (i in which(node_select)-1) {
    id <- append(id, NA, after = i)
  }
  df[,"id"] <- id

  df$count <- as.numeric(df$count)
  df$name_root <- ifelse(df$count == max(df$count), "root", df$name)
  df$id <- ifelse(df$count == max(df$count), paste0("root_", df$name), df$id)
  df$sampleID <- sampleID
  df$FileName <- df[df$count == max(df$count), "name"]
  df$NodeType <- node_types
  names(df)[which(names(df) == "count")] <- "Count"
  depend_list <- xml2::xml_find_all(pops, ".//Dependents", flatten = T) #|.//Dependents
  depend_list_parents <- xml2::xml_parent(depend_list)

  # handle assignment of parent nodes for OrNodes and AndNodes outside this function
  parent_node_types <- xml2::xml_name(depend_list_parents)
  ## add dependent notnodes?
  if ("OrNode" %in% parent_node_types) {
    OrNodes <- depend_list_parents[which(parent_node_types == "OrNode")]
  } else {
    OrNodes <- NULL
  }
  if ("AndNode" %in% parent_node_types) {
    AndNodes <- depend_list_parents[which(parent_node_types == "AndNode")]
  } else {
    AndNodes <- NULL
  }
  ## NotNodes only relevant if they stem from OrNodes/AndNodes, because then they have no id
  if ("NotNode" %in% parent_node_types) {
    NotNodes <- depend_list_parents[which(parent_node_types == "NotNode")]
    # check if NotNodes from node_details_list already have id in df - then they could be removed from node_details_list
    # they have an id if they do not stem from OrNodes/AndNodes; if they do they have no id
    # match returns only the first match !! similar to duplicated
    NotNodes <- NotNodes[which(is.na(df[match(xml2::xml_attr(NotNodes, "name"), df$name),"id"]))]
  } else {
    NotNodes <- NULL
  }


  if (more_gate_data) {
    # OrNodes and AndNodes: add xChannel and yChannel and eventInside later from originating gates in add_boolean_gate_data
    #gate_list2 <- gate_list[which(node_types == "Population")]

    gate_list_children <- xml2::xml_children(gate_list)
    # for id: filter for relevant children first, then obtain parents
    gate_list_children <- gate_list_children[which(lengths(xml2::xml_attrs(gate_list_children)) > 1)]
    id <- do.call(dplyr::bind_rows, xml2::xml_attrs(xml2::xml_parent(gate_list_children)))[,"id",drop=T]
    eventsInside <- do.call(dplyr::bind_rows, xml2::xml_attrs(gate_list_children))[,"eventsInside",drop=T]

    gate_list_children2 <- xml2::xml_children(gate_list_children)
    gate_list_children2 <- gate_list_children2[which(xml2::xml_name(gate_list_children2) == "dimension")]
    gate_list_children2 <- xml2::xml_children(gate_list_children2)
    dims <- xml2::xml_attrs(gate_list_children2)
    if (length(dims)/2 != length(id)) {
      ndims <- sapply(lapply(gate_list_children, function(x) xml2::xml_name(xml2::xml_children(x))), function(y) sum(y == "dimension"))
      xChannel <- character(length(id))
      yChannel <- character(length(id))
      n = 1
      for (i in 1:length(ndims)) {
        xChannel[i] <- dims[n]
        if (ndims[i] == 2) {
          yChannel[i] <- dims[n+1]
        } else {
          yChannel[i] <- NA
        }
        n <- n + i
        xChannel <- unlist(xChannel)
        yChannel <- unlist(yChannel)
      }
      # count dimension for each gate and split dims accordingly
    } else {
      # does not work with 1D-gates which only have one dimension and hence one channel
      xChannel <- unlist(dims[seq(1, length(dims), 2)])
      yChannel <- unlist(dims[seq(2, length(dims), 2)])
    }

    df2 <- data.frame(id = id,
                      eventsInside = eventsInside,
                      xChannel = xChannel,
                      yChannel = yChannel)
    df <- dplyr::left_join(df, df2[which(!is.na(df2$id)),], by = "id")
  }
  return(list(df = df, OrNodes = OrNodes, AndNodes = AndNodes, NotNodes = NotNodes))
}



add_full_paths <- function(df, graph, end_edges = NULL, show_progress = F) {
  # providing end_edges speeds up the process
  # starting from end edges should catch all gates (nodes) at least once (logic, maybe)
  if (is.null(end_edges)) {
    edges <- df$id
  } else {
    edges <- end_edges
  }

  # different OrNodes / AndNodes with same id (multiple id assigned by same originating gates?!)
  full_paths_df <- purrr::map(edges, function(x) {
    path_to_root <- igraph::shortest_paths(graph,
                                           mode = "all",
                                           from = x,
                                           to = paste0("root_", df[which(df$id == x), "FileName"]),
                                           algorithm = "unweighted")
    ## derive full paths
    pops_to_root <- rev(df[match(names(path_to_root[[1]][[1]]), df$id),"name"])
    pops_to_root[1] <- ""
    full_paths <- rev(sapply(1:length(pops_to_root), function(x) paste(pops_to_root[1:x], collapse = "/")))
    full_id_paths <- rev(sapply(1:length(names(path_to_root[[1]][[1]])), function(x) paste(rev(names(path_to_root[[1]][[1]]))[1:x], collapse = "/")))
    full_paths <- gsub("^/", "", full_paths)
    full_paths <- data.frame(name = ifelse(full_paths == "", rev(df[match(names(path_to_root[[1]][[1]]), df$id),"name"])[1], basename(full_paths)),
                             id = names(path_to_root[[1]][[1]]),
                             parent_id = c(names(path_to_root[[1]][[1]]), NA)[-1],
                             grandparent_id = c(names(path_to_root[[1]][[1]]), NA, NA)[-c(1:2)],
                             PopulationFullPath = full_paths,
                             PopulationFullPathID = full_id_paths)
    return(full_paths)
  }, .progress = show_progress)

  while (length(full_paths_df) > 20) {
    full_paths_df <- purrr::map(split(c(1:length(full_paths_df)), ceiling(seq_along(c(1:length(full_paths_df)))/10)), function(x) purrr::reduce(full_paths_df[x], dplyr::bind_rows))
  }
  full_paths_df <- unique(purrr::reduce(full_paths_df, dplyr::bind_rows))
  full_paths_df$GateDepth <- nchar(full_paths_df$PopulationFullPath) - nchar(gsub("/", "", full_paths_df$PopulationFullPath)) + 1
  full_paths_df$GateDepth <- ifelse(full_paths_df$PopulationFullPath == "", 0, full_paths_df$GateDepth)
  df <- dplyr::left_join(df, full_paths_df, by = c("id" = "id", "parent_id" = "parent_id", "name" = "name"))



  ## too complicated and slow:
  #grandparents <- purrr::map(igraph::incident_edges(graph, df$parent_id[which(!is.na(df$parent_id))]), igraph::ends, graph = graph)
  #grandparents <- purrr::map(grandparents, as.data.frame)
  #grandparents <- dplyr::bind_rows(grandparents)
  #grandparents <- dplyr::right_join(grandparents, grandparents, by = c("V2" = "V1"), relationship = "many-to-many")
  #grandparents <- grandparents[,-2]
  #names(grandparents) <- c("grandparent_id", "id")
  #grandparents <- unique(grandparents)
  #df <- dplyr::left_join(df, grandparents, by = "id")
  return(df)
}

add_root_node <- function(df) {
  dplyr::bind_rows(df, data.frame(FileName = unique(df$FileName),
                                  id = paste0("root_", unique(df$FileName)),
                                  parent_id = NA))
}

add_parent_count <- function(df) {
  df_parent <- df[,match(c("id", "Count"), names(df))]
  df_parent <- df_parent[which(!is.na(df_parent$id)),] # depends
  names(df_parent) <- c("parent_id", "ParentCount")
  df <- dplyr::left_join(df, df_parent, by = "parent_id")
  df$FractionOfParent <- df$Count/df$ParentCount
  return(df)
}

add_grandparent_count <- function(df) {
  df_grandparent <- df[,match(c("id", "Count"), names(df))]
  df_grandparent <- df_grandparent[which(!is.na(df_grandparent$id)),] # depends
  names(df_grandparent) <- c("grandparent_id", "GrandparentCount")
  df <- dplyr::left_join(df, df_grandparent, by = "grandparent_id")
  df$FractionOfGrandparent <- df$Count/df$GrandparentCount
  return(df)
}

add_total_count <- function(df) {
  df_total <- df[which(df$name == df$FileName),c("FileName", "Count")]
  names(df_total)[2] <- "TotalCount"
  df <- dplyr::left_join(df, df_total, by = "FileName")
  df$FractionOfTotal <- df$Count/df$TotalCount
  return(df)
}


add_boolean_gate_data <- function(df,
                                  node_details_list,
                                  nodes_name = c("OrNodes", "AndNodes", "NotNodes"),
                                  more_gate_data = F) {

  nodes_name <- match.arg(nodes_name, c("OrNodes", "AndNodes", "NotNodes"))

  # identify parents of OrNodes/AndNodes by name and Count but without PopulationFullPath
  ## add id and parent_id to OrNodes/AndNodes themselves
  ## something is weird here
  temp_df <- purrr::map_dfr(sapply(node_details_list, "[", nodes_name), function(x) {
    purrr::map_dfr(x, function(y) {
      xx <- do.call(dplyr::bind_rows, xml2::xml_attrs(xml2::xml_children(xml2::xml_parent(y))))
      yy <- basename(xml2::xml_attr(xml2::xml_children(xml2::xml_child(x, "Dependents")[[1]]), "name"))
      xx <- xx[which(xx$name %in% yy), c("name", "count")]
      xx$name2 <- xml2::xml_attr(y, "name")
      xx$Count2 <- xml2::xml_attr(y, "count")
      return(xx)
    })
  }, .id = "FileName")

  if (nrow(temp_df) == 0) {
    # this happened with NotNodes; but should now be caught outside, see above
    return(df)
  }

  #temp_df <- unique(temp_df)
  names(temp_df)[which(names(temp_df) == "count")] <- "Count"
  temp_df$Count <- as.numeric(temp_df$Count)
  temp_df$Count2 <- as.numeric(temp_df$Count2)
  temp_df$FileName <- gsub(paste0("\\.", nodes_name, "$"), "", temp_df$FileName)
  if (more_gate_data) {
    cols <- c("FileName", "Count", "name", "id", "parent_id", "xChannel", "yChannel", "eventsInside")
  } else {
    cols <- c("FileName", "Count", "name", "id", "parent_id")
  }

  temp_df <- dplyr::left_join(temp_df,
                              df[,cols],
                              by = c("FileName", "Count", "name"))
  # or gates in different dimension will cause error
  if (more_gate_data) {
    temp_df <- dplyr::group_by(temp_df, FileName, name2, Count2, parent_id, eventsInside)
  } else {
    temp_df <- dplyr::group_by(temp_df, FileName, name2, Count2, parent_id)
  }
  # do not group by xChannel and yChannel, but summarise them as below
  # if OrNodes/AndNodes are from different dimension, then summarizing would not work correctly if xChannel, yChannel were used for grouping
  # if done like all dimension from different feeding are gates saved

  temp_df <- dplyr::summarise(temp_df,
                              id = paste(sort(id), collapse = ","),
                              xChannel = paste(sort(unique(xChannel)), collapse = ","),
                              yChannel = paste(sort(unique(yChannel)), collapse = ","),
                              .groups = "drop")

  names(temp_df)[which(names(temp_df) == "name2")] <- "name"
  names(temp_df)[which(names(temp_df) == "Count2")] <- "Count"

  df <- coalesce_join(df, temp_df, by = c("FileName", "name", "Count")) # join via name and Count: only in a super rare case when there are two OrNodes with same name and same count, this will give a conflict


  ## add parent_id to children of OrNodes/AndNodes themselves
  temp_df2 <- purrr::map_dfr(sapply(node_details_list, "[", nodes_name), function(x) {
    purrr::map_dfr(x, function(y) {
      if (!is.na(xml2::xml_child(y, "Subpopulations"))) {
        # is na when there are no children to OrNodes or AndNodes
        xx <- do.call(dplyr::bind_rows, xml2::xml_attrs(xml2::xml_children(xml2::xml_child(y, "Subpopulations"))))[, c("name", "count"),drop=T]
        names(xx)[which(names(xx) == "name")] <- "name2"
        names(xx)[which(names(xx) == "count")] <- "count2"
        # is x[[1]] correct? or should is be y?
        xx$name <- xml2::xml_attrs(x[[1]])[["name"]]
        xx$Count <- xml2::xml_attrs(x[[1]])[["count"]]
        return(xx)
      } else {
        return(NULL)
      }
    })
  }, .id = "FileName")

  if (nrow(temp_df2) > 0) {
    temp_df2$Count <- as.numeric(temp_df2$Count)
    temp_df2$count2 <- as.numeric(temp_df2$count2)
    temp_df2$FileName <- gsub(paste0("\\.", nodes_name, "$"), "", temp_df2$FileName)
    # join temp_df
    temp_df2 <- dplyr::left_join(temp_df2, temp_df, by = c("FileName", "name", "Count"))
    temp_df2 <- temp_df2[,which(names(temp_df2) %in% c("FileName", "name2", "count2", "id"))]
    names(temp_df2)[2:4] <- c("name", "Count", "parent_id")
    # check for duplicate rows (very rare case)
    df <- coalesce_join(df, temp_df2, by = c("FileName", "name", "Count")) # join via name and Count: only in a super rare case when there are two OrNodes with same name and same count, this will give a conflict
  }

  return(df)
}


add_OrNode_AndNode_data <- function(df, node_details_list, nodes_name = c("OrNodes", "AndNodes")) {

  #grandparent_id
  #GateDepth
  #PopulationFullPathID
  nodes_name <- match.arg(nodes_name, c("OrNodes", "AndNodes"))

  temp_df <- purrr::map_dfr(sapply(node_details_list, "[", nodes_name), function(x) {
    data.frame(name = xml2::xml_attr(x, "name"),
               ParentFullPath = dirname(attr(xml2::as_list(xml2::xml_child(x, "Dependents"))[[1]][[1]], "name")))
  }, .id = "FileName")
  temp_df$FileName <- gsub(paste0("\\.", nodes_name, "$"), "", temp_df$FileName)
  temp_df <- dplyr::left_join(temp_df,
                              df[,c("FileName", "PopulationFullPath", "id")],
                              by = c("FileName"= "FileName", "ParentFullPath" = "PopulationFullPath"))
  names(temp_df)[which(names(temp_df) == "ParentFullPath")] <- "PopulationFullPath"
  names(temp_df)[which(names(temp_df) == "id")] <- "parent_id"
  temp_df$PopulationFullPath <- paste0(temp_df$PopulationFullPath, "/", temp_df$name)
  df <- coalesce_join(df, temp_df, by = c("FileName", "name", "Count")) # join via name and Count: only in a super rare case when there are two OrNodes with same name and same count, this will give a conflict

  ## add ids to OrNode
  temp_df <- purrr::map_dfr(sapply(node_details_list, "[", nodes_name), function(x) {
    temp <- stack(unlist(sapply(xml2::as_list(xml2::xml_find_all(x, "Dependents"))[[1]], attributes)))[,-2,drop=F]
    temp$name <- xml2::xml_attr(x, "name")
    temp$Count <- as.numeric(xml2::xml_attr(x, "count"))
    return(temp)
  }, .id = "FileName")
  temp_df$FileName <- gsub(paste0("\\.", nodes_name, "$"), "", temp_df$FileName)
  names(temp_df)[which(names(temp_df) == "values")] <- "PopulationFullPath"
  temp_df <- dplyr::left_join(temp_df, df[,c("FileName", "PopulationFullPath", "id")], by = c("FileName", "PopulationFullPath"))
  temp_df <- dplyr::group_by(temp_df, FileName, name, Count)
  temp_df <- dplyr::summarise(temp_df, id = paste(id, collapse = ","), .groups = "drop")
  df <- coalesce_join(df, temp_df, by = c("FileName", "name", "Count")) # join via name and Count: only in a super rare case when there are two OrNodes with same name and same count, this will give a conflict

  return(df)
}

bind_rows_chunked <- function(df_list, chunk_size = 10) {
  while (length(df_list) > chunk_size*2) {
    df_list <- purrr::map(.x = split(c(seq_along(df_list)), ceiling(seq_along(c(seq_along(df_list)))/chunk_size)),
                          .f = function(x) purrr::reduce(df_list[x], dplyr::bind_rows))
  }
  return(purrr::reduce(df_list, dplyr::bind_rows))
}

add_channel_desc <- function(df, ws, keys_list = NULL) {

  'keys <- purrr::map_dfr(xml2::xml_children(xml2::xml_child(ws, "SampleList")), function(x) {
    keys <- xml2::xml_attrs(xml2::xml_contents(xml2::xml_child(x, "Keywords")))
    keys <- stats::setNames(sapply(keys, "[", 2), sapply(keys, "[", 1))
    keys <- utils::stack(keys)
    names(keys) <- c("value", "name")
    keys <- keys[,c(2,1)]
    keys$name <- as.character(keys$name)
    keys <- keys[which(grepl("^\\$P[[:digit:]]{1,2}[NS]", keys$name)),]
    keys$channel_digit <- gsub("[NS]", "", gsub("\\$P", "", keys$name))
    keys$name <- gsub("\\$P[[:digit:]]{1,2}", "", keys$name)
    keys <- as.data.frame(tidyr::pivot_wider(keys, names_from = name, values_from = value))
    keys$FileName <- xml2::xml_attrs(xml2::xml_child(x, "SampleNode"))[["name"]]
    keys <- keys[,-which(names(keys) == "channel_digit")]
    return(keys)
  })'

  #get keys from above
  if (is.null(keys_list)) {
    file_names <- purrr::map_chr(xml2::xml_children(xml2::xml_child(ws, "SampleList")), function(x) xml2::xml_attrs(xml2::xml_child(x, "SampleNode"))[["name"]])
    keys_list <- purrr::map(setNames(xml2::xml_children(xml2::xml_child(ws, "SampleList")), file_names), function(x) {
      keys <- xml2::xml_attrs(xml2::xml_contents(xml2::xml_child(x, "Keywords")))
      keys <- stats::setNames(sapply(keys, "[", 2), sapply(keys, "[", 1))
      keys <- utils::stack(keys)
      names(keys) <- c("value", "name")
      keys <- keys[,c(2,1)]
      keys$name <- as.character(keys$name)
      return(keys)
    })
  }
  keys <- purrr::map_dfr(keys_list, function(x) {
    x <- x[which(grepl("^\\$P[[:digit:]]{1,2}[NS]", x$name)),]
    x$channel_digit <- gsub("[NS]", "", gsub("\\$P", "", x$name))
    x$name <- gsub("\\$P[[:digit:]]{1,2}", "", x$name)
    x <- as.data.frame(tidyr::pivot_wider(x, names_from = name, values_from = value))
    x <- x[,-which(names(x) == "channel_digit")]
    return(x)
  }, .id = "FileName")

  ## for OrNodes/AndNodes from different dimension: expand rows first, then join, then collapse them again
  ## check for comma in xChannel, yChannel to determine if this special procedure is needed

  if (any(grepl(",", df$xChannel)) || any(grepl(",", df$yChannel))) {
    if (any(grepl(",", df$xChannel))) {
      df$xChannel <- strsplit(df$xChannel, ",")
      df <- tidyr::unnest(df, xChannel)
      names(keys)[which(names(keys) == "N")] <- "xChannel"
      names(keys)[which(names(keys) == "S")] <- "xDesc"
      df <- dplyr::left_join(df, keys, by = c("FileName", "xChannel"))
      df <- dplyr::group_by(df, !!!rlang::syms(names(df)[which(!names(df) %in% c("xChannel", "xDesc"))]))
      df <- dplyr::summarise(df, xDesc = paste0(xDesc, collapse = ","), xChannel = paste0(xChannel, collapse = ","), .groups = "drop")
    }
    if (any(grepl(",", df$yChannel))) {
      df$yChannel <- strsplit(df$yChannel, ",")
      df <- tidyr::unnest(df, yChannel)
      if (any(grepl(",", df$xChannel))) {
        names(keys)[which(names(keys) == "xChannel")] <- "yChannel"
        names(keys)[which(names(keys) == "xDesc")] <- "yDesc"
      } else {
        names(keys)[which(names(keys) == "N")] <- "yChannel"
        names(keys)[which(names(keys) == "S")] <- "yDesc"
      }
      df <- dplyr::left_join(df, keys, by = c("FileName", "yChannel"))
      df <- dplyr::group_by(df, !!!rlang::syms(names(df)[which(!names(df) %in% c("yChannel", "yDesc"))]))
      df <- dplyr::summarise(df, yDesc = paste0(yDesc, collapse = ","), yChannel = paste0(yChannel, collapse = ","), .groups = "drop")
    }
  } else {
    names(keys)[which(names(keys) == "N")] <- "xChannel"
    names(keys)[which(names(keys) == "S")] <- "xDesc"
    df <- dplyr::left_join(df, keys, by = c("FileName", "xChannel"))
    names(keys)[which(names(keys) == "xChannel")] <- "yChannel"
    names(keys)[which(names(keys) == "xDesc")] <- "yDesc"
    df <- dplyr::left_join(df, keys, by = c("FileName", "yChannel"))
  }
  df$xDesc[which(df$xDesc == "")] <- NA
  df$yDesc[which(df$yDesc == "")] <- NA
  df$xDesc[which(df$xDesc == "NA")] <- NA
  df$yDesc[which(df$yDesc == "NA")] <- NA
  df$xChannel[which(df$xChannel == "NA")] <- NA
  df$yChannel[which(df$yChannel == "NA")] <- NA

  return(as.data.frame(df))
}

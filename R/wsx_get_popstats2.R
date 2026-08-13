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
#' @param strip_data remove internal columns from count df and speed up the function
#'
#' @return list with counts and statistics and graphs
#' @export
#' @examples
#' \dontrun{
#' # When the script is saved to R_scripts in the experiment folder,
#' # get the absolute path to the folder
#' wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
#' # find workspaces
#' ws <- list.files(path = wd, pattern = '\\.wsp$', recursive = T, full.names = T)
#' # import the population counts:
#' lst <- wsx_get_popstats(ws = ws[[1]])
#' #' # plot graph
#' ggraph::ggraph(lst$graph_sample[[1]]) +
#'   ggraph::geom_node_point() +
#'   ggraph::geom_edge_link() +
#'   ggraph::geom_node_label(ggplot2::aes(label = name))
#' }
wsx_get_popstats2 <- function(ws,
                              groups = NULL,
                              invert_groups = F,
                              return_stats = F,
                              strip_data = F) {

  if (!requireNamespace("brathering", quietly = T)) {
    pak::pak("Close-your-eyes/brathering")
  }

  ws_raw <- ws
  ws <- fcexpr:::check_ws(ws) #fcexpr:::
  group_df <- fcexpr:::get_group_df(ws, groups, invert_groups) #fcexpr:::
  samples <- fcexpr:::get_sample_nodes(ws, group_df) #fcexpr:::
  samplenodenames <- do.call(dplyr::bind_rows, xml2::xml_attrs(xml2::xml_find_all(samples, "SampleNode", flatten = T)))
  fcsfilenames <- basename(xml2::xml_attr(xml2::xml_find_all(samples, "DataSet", flatten = T), "uri"))

  keys_list <- wsx_get_keywords(
    ws = ws,
    samples = samplenodenames$name,
    verbose = F)
  if (is.null(keys_list)) {
    # when samples are added to wsp before renaming on disk
    keys_list <- wsx_get_keywords(ws = ws,
                                  samples = fcsfilenames,
                                  verbose = F)
  }
  if (is.null(keys_list)) {
    keys_list <- wsx_get_keywords(ws = ws,
                                  fverbose = F)
  }

  fcs_idents <-
    utils::stack(fcexpr:::get_fcs_identities(kwl = keys_list[["vec"]])) |>
    dplyr::rename(identity = "values", "FileName" = ind) |>
    dplyr::mutate(FileName = as.character(FileName))

  conv <- check_fcs_namechange(filenames = fcs_idents$FileName,
                               nodenames = samplenodenames[["name"]])
  # samplenodenames[["name"]] <- unname(conv[samplenodenames[["name"]]])

  tt <- purrr::map_dfr(samples, ~recursive_walk_xml(x = xml2::xml_child(.x, "SampleNode")))
  if (!is.null(conv)) {
    fullp <- strsplit(tt$PopulationFullPath, "/")
    fullp <- purrr::map2_chr(fullp, unname(conv[sapply(fullp, "[", 1)]), function(x,y) {
      x[1] <- y
      x <- paste(x, collapse = "/")
      return(x)
    })
    tt$PopulationFullPath <- fullp
  }

  tt <- tt |>
    tidyr::fill(sampleID) |>
    dplyr::mutate(parent = dirname(PopulationFullPath)) |>
    dplyr::mutate(grandparent = dirname(parent)) |>
    dplyr::rename("Count" = count) |>
    dplyr::mutate(Count = as.numeric(Count)) |>
    dplyr::mutate(GateDepth = stringr::str_count(PopulationFullPath, "/")) |>
    dplyr::mutate(id = ifelse(NodeType == "SampleNode", paste0("ID", sampleID), id),
                  parent_id = ifelse(is.na(parent_id) & GateDepth == 1, paste0("ID", sampleID), parent_id)) |>
    dplyr::mutate(PopulationFullPath2 = rm_root(PopulationFullPath)) |>
    dplyr::mutate(FlowJoWsp = ws_raw)  |>
    dplyr::left_join(group_df, by = c("sampleID")) |>
    dplyr::mutate(FileName = get_root(PopulationFullPath)) |>
    dplyr::left_join(fcs_idents, by = "FileName")

  tt <- fix_missing_ids(df = tt)
  tt <- add_grandparent_id(df = tt)
  tt <- add_count(df = tt, type = "parent")
  tt <- add_count(df = tt, type = "grandparent")
  tt <- add_total_count2(df = tt)
  tt <- add_population(df = tt, full_path_col = "PopulationFullPath2")
  tt <- add_filpath(df = tt, sample_nodeset = samples)

  if (!strip_data) {
    tt <- add_channel_desc(df = tt, ws = ws, keys_list = keys_list[["df"]])
  }

  graphs <- make_graphs(tt)

  cols <- c(
    "FileName",
    "PopulationFullPath",
    "Population",
    "Count",
    "ParentCount",
    "GrandparentCount",
    "TotalCount",
    "FractionOfParent",
    "FractionOfGrandparent",
    "FractionOfTotal",
    "FlowJoGroup",
    "FlowJoWsp",
    "FilePath",
    "identity"
  )
  if (!strip_data) {
    cols <- c(
      cols,
      "parent",
      "grandparent",
      "GateType",
      "NodeType",
      "id",
      "parent_id",
      "grandparent_id",
      "sampleID",
      "GateDepth",
      "xChannel",
      "yChannel",
      "xDesc",
      "yDesc",
      "eventsInside"
    )
  }

  tt <- tt |>
    dplyr::select(-PopulationFullPath) |>
    dplyr::rename("PopulationFullPath" = PopulationFullPath2) |>
    dplyr::mutate(PopulationFullPath = unname(PopulationFullPath))

  if (!strip_data) {
    tt <- dplyr::rename(tt[,cols], "ID" = id,
                        "ParentID" = parent_id,
                        "GrandparentID" = grandparent_id,
                        "ParentPath" = parent,
                        "GrandparentPath" = grandparent) |>
      dplyr::mutate(ParentPath = rm_root(ParentPath),
                    GrandparentPath = rm_root(GrandparentPath))
    # dplyr::mutate(ParentPath = ifelse(ParentPath == ".", NA, ParentPath)) |>
    # dplyr::mutate(GrandparentPath = ifelse(GrandparentPath == ".", NA, GrandparentPath))
    tt$ParentPath[which(tt$ParentPath == ".")] <- NA
    tt$GrandparentPath <- NA
  } else {
    tt <- tt[,cols]
  }


  stats_out <- NULL
  if (return_stats) {
    stats_out <- do.call(rbind, lapply(seq_along(samples), function(n) {
      node <- samples[n]
      stats <- xml2::xml_find_all(node, ".//Statistic")
      stats_df <- do.call(rbind, lapply(stats, function(x) {
        prnts <- xml2::xml_parents(x)
        p_nodes <- prnts[which(xml2::xml_name(prnts) %in% c("AndNode", "OrNode", "NotNode", "Population"))]

        sampleID <- xml2::xml_attr(xml2::xml_child(x, "DataSet"), "sampleID")
        FilePath <- gsub("^file:", "", xml2::xml_attr(xml2::xml_child(node, "DataSet"), "uri"))
        FileName <- basename(FilePath)
        PopulationFullPath <- if (length(p_nodes) == 0) {"root"} else {paste(rev(xml2::xml_attr(p_nodes, "name")), collapse = "/")}

        df <- data.frame(FileName = FileName,
                         PopulationFullPath = PopulationFullPath,
                         statistic = xml2::xml_attr(x, "name"),
                         channel = xml2::xml_attr(x, "id"),
                         value = suppressWarnings(as.numeric(xml2::xml_attr(x, "value"))),
                         FilePath = FilePath,
                         stringsAsFactors = F)
        if (is.na(df$value)) {
          message("stats: statistic not a number.")
          print(df)
        }
        return(df)
      }))

      return(stats_df)
    }))
    if (!is.null(stats_out)) {
      stats_out <- tibble::as_tibble(stats_out)
      if (anyNA(stats_out$value)) {
        print(dplyr::filter(stats_out, is.na(value)), n = Inf)
      }
    }
  }

  message("FileName is as in FlowJo wsp file. Join counts to sampledescription (sd) via FileName and identity.")
  message("dplyr::left_join(counts, sd, by = c('FileName', 'identity'))")
  # check for equal file names on hdd
  fcexpr:::check_filenames_hdd(pop_df = tt)

  ## check duplicated files
  dup_files <- fcs_idents[duplicated(fcs_idents$identity),]
  if (nrow(dup_files) > 1) {
    message("One or multiple FCS files seem to be duplicates. A data.frame with details will be returned.")
    print(dup_files)
  } else {
    dup_files <- NULL
  }

  return(list(
    counts = tibble::as_tibble(tt),
    stats = stats_out,
    graph = graphs[["graph"]],
    graph_sample = graphs[["graph_samples"]],
    duplicate_FCS_files = dup_files
  ))

}



recursive_walk_list <- function(x, parent_path = NULL) {

  path <- paste(c(parent_path, attr(x, "name")), collapse = "/")

  df <- data.frame(full_path = path,
                   count = attr(x, "count"),
                   sampleID = if (is.null(attr(x, "sampleID"))) NA else attr(x, "sampleID"))


  if ("Subpopulations" %in% names(x)) {
    for (i in 1:length(x[["Subpopulations"]])) {
      df <- dplyr::bind_rows(df, recursive_walk_list(x[["Subpopulations"]][[i]], parent_path = path))
    }
  }
  return(df)
}

recursive_walk_xml <- function(x, parent_path = NULL) {

  path <- paste(c(parent_path, xml2::xml_attr(x, "name")), collapse = "/")
  # print(path)

  gate <- xml2::xml_child(x, "Gate")
  node_type <- xml2::xml_name(x)
  gate_type <- if (!length(gate)) NA else xml2::xml_name(xml2::xml_child(gate))

  xChannel <- NA
  yChannel <- NA
  eventsInside <- NA
  if (length(gate)) {
    channels <- get_channels_dims(gate)
    xChannel <- channels[1]
    if (length(channels) == 2) {
      yChannel <- channels[2]
    }
    eventsInside <- xml2::xml_attr(xml2::xml_child(gate), "eventsInside")
  }

  df <- data.frame(PopulationFullPath = path,
                   count = xml2::xml_attr(x, "count"),
                   sampleID = if (is.null(xml2::xml_attr(x, "sampleID"))) NA else xml2::xml_attr(x, "sampleID"),
                   id = xml2::xml_attr(gate, "id"),
                   parent_id = xml2::xml_attr(gate, "parent_id"),
                   NodeType = node_type,
                   GateType = gate_type,
                   xChannel = xChannel,
                   yChannel = yChannel,
                   eventsInside = eventsInside,
                   depends = if (node_type %in% c("Population", "SampleNode")) NA else I(list(paste0(strsplit(path, "/")[[1]][1], "/", xml2::xml_attr(xml2::xml_children(xml2::xml_child(x, "Dependents")), "name")))))


  if ("Subpopulations" %in% xml2::xml_name(xml2::xml_children(x))) {
    nodeset <- xml2::xml_contents(xml2::xml_child(x, "Subpopulations"))
    for (i in which(xml2::xml_name(nodeset) != "Statistic")) {
      df <- dplyr::bind_rows(df, recursive_walk_xml(x = xml2::xml_contents(xml2::xml_child(x, "Subpopulations"))[[i]],
                                                    parent_path = path))
    }
  }

  return(df)
}

fix_missing_ids <- function(df) {

  repeat {
    df$old_id <- df$id

    df <- df |>
      dplyr::mutate(
        id = ifelse(
          is.na(id),
          purrr::map_chr(depends, function(dep) {
            # handle NULL or empty
            if (is.null(dep) || length(dep) == 0) return(NA_character_)

            # match paths → ids
            ids <- df$id[match(dep, df$PopulationFullPath)]

            # clean
            ids <- ids[!is.na(ids)]
            if (length(ids) == 0) return(NA_character_)

            # unique + sorted + collapsed
            paste(sort(unique(ids)), collapse = "_")
          }),
          id
        )
      )

    if (identical(df$old_id, df$id)) break
  }

  dup_counts <- ave(df$old_id, df$old_id, FUN = length)
  dup_index  <- ave(df$old_id, df$old_id, FUN = seq_along)

  df$id <- ifelse(
    dup_counts > 1,
    paste0(df$old_id, "_", dup_index),
    df$old_id
  )

  df_parent <- df |>
    dplyr::select(PopulationFullPath, id, sampleID) |>
    dplyr::rename("parent_id" = id, "parent" = PopulationFullPath)

  df <- brathering::coalesce_join(df, df_parent, by = c("parent" = "parent", "sampleID" = "sampleID"))
  df <- dplyr::select(df, -old_id)
  return(df)
}


add_total_count2 <- function(df) {
  df <- dplyr::left_join(df, df |>
                           dplyr::filter(GateDepth == 0) |>
                           dplyr::select(Count, sampleID) |>
                           dplyr::rename("TotalCount" = Count),
                         by = "sampleID") |>
    dplyr::mutate(FractionOfTotal = Count/TotalCount)
  return(df)
}

add_grandparent_id <- function(df) {
  df <- df |>
    dplyr::left_join(df |>
                       dplyr::select(PopulationFullPath, sampleID, id) |>
                       dplyr::rename("grandparent_id" = id),
                     by = c("grandparent" = "PopulationFullPath", "sampleID" = "sampleID"))
  return(df)
}

get_root <- function(x) {
  sapply(strsplit(x, "/"), "[", 1)
}

rm_root <- function(x) {
  sapply(brathering::strsplit2(x, "/"), "[", 2)
}

add_population <- function(df, full_path_col = "PopulationFullPath2") {
  auto_paths <- fcexpr:::shortest_unique_path(df[[full_path_col]]) #fcexpr:::
  df$Population <- unname(auto_paths[df[[full_path_col]]])
  return(df)
}

add_filpath <- function(df, sample_nodeset) {
  file_paths <-
    as.data.frame(do.call(dplyr::bind_rows, xml2::xml_attrs(xml2::xml_child(sample_nodeset, "DataSet")))) |>
    dplyr::rename("FilePath" = uri) |>
    dplyr::mutate(FilePath = gsub("^file:", "", FilePath))
  df <- dplyr::left_join(df, file_paths, by = "sampleID")
  return(df)
}

get_channels_dims <- function(gate) {
  y <- xml2::xml_children(xml2::xml_child(gate))
  y <- y[which(xml2::xml_name(y) == "dimension")]
  y <- xml2::xml_children(y)
  channels <- unname(unlist(xml2::xml_attrs(y)))
  return(channels)
}

make_graphs <- function(df) {
  df2 <- df |>
    dplyr::filter(parent != ".") |>
    dplyr::rename("from" = PopulationFullPath, "to" = parent) |>
    dplyr::select(from, to)

  graph <- igraph::graph_from_data_frame(df2, directed = T)
  edge_degrees <- igraph::degree(graph, mode = "out")
  # end_edges <- names(edge_degrees[which(edge_degrees == 0)])

  igraph::V(graph)$full_path <- igraph::V(graph)$name
  igraph::V(graph)$full_path2 <- unname(rm_root(igraph::V(graph)$full_path))
  igraph::V(graph)$final_node <- basename(igraph::V(graph)$name)
  igraph::V(graph)$population <- unname(df[match(names(igraph::V(graph)), unique(df$PopulationFullPath)), "Population"])
  igraph::V(graph)$PopulationFullPath <- igraph::V(graph)$full_path
  igraph::V(graph)$PopulationFullPath2 <- igraph::V(graph)$full_path2

  # get one graph for each sample
  graph_subgroups <- igraph::components(graph)
  subgroup_vertex_id_list <- split(names(graph_subgroups$membership), graph_subgroups$membership)
  subgroup_vertex_id_list <- purrr::map(subgroup_vertex_id_list, rev)
  graph_samples <- purrr::map(stats::setNames(subgroup_vertex_id_list, sapply(subgroup_vertex_id_list, "[", 1)),
                              ~igraph::subgraph(graph = graph, vids = .x))


  return(list(graph = graph,
              graph_samples = graph_samples))
}

check_fcs_namechange <- function(filenames, nodenames) {
  conv <- NULL
  if (any(!filenames %in% nodenames)) {
    message("did you change filenames after loading fcs files into flowjo?")
    match_inds <- apply(stringdist::stringdistmatrix(unique(nodenames), filenames),
                        1,
                        which.min)
    if (anyDuplicated(match_inds)) {
      stop("could not match old and new FileName unambigously. Make new wsp and re-import renamed FCS files.")
    }
    message("change old to new FileName based on best match. You may re-create the wsp with renamed FCS files from scratch.")
    # iteration over rows makes sense!
    conv <- stats::setNames(filenames[match_inds], unique(nodenames))
  }
  return(conv)
}

add_count <- function(df, type = c("parent", "grandparent")) {
  type <- rlang::arg_match(type)
  id <- paste0(type, "_id")
  name <- paste0(stringr::str_to_title(type), "Count")
  name2 <- paste0("FractionOf", stringr::str_to_title(type))
  df2 <- df[,match(c("id", "Count"), names(df))]
  df2 <- df2[which(!is.na(df2$id)),] # depends
  names(df2) <- c(id, name)
  # browser()
  df <- dplyr::left_join(df, df2, by = id)
  df[[name2]] <- df[["Count"]]/df[[name]]
  return(df)
}

add_channel_desc <- function(df, ws, keys_list = NULL) {

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

  ## why is this needed - check somewhen why, above
  keys <- dplyr::distinct(keys)

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

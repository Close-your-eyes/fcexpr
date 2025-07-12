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
wsx_get_popstats_legacy <- function(ws,
                                    groups = NULL,
                                    invert_groups = F,
                                    return_stats = T,
                                    lapply_fun = lapply,
                                    strip_data = T,
                                    ...) {

  ## allow to pass mclapply
  lapply_fun <- match.fun(lapply_fun)

  ws_raw <- ws
  ws <- fcexpr:::check_ws(ws) #fcexpr:::
  group_df <- get_group_df(ws, groups) #fcexpr:::
  samples <- get_sample_nodes(ws, group_df) #fcexpr:::


  gates_list <- do.call(rbind, lapply_fun(xml2::xml_find_all(samples, ".//Gate|.//Dependents"), function(n) {

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
    GateDepth <- length(p_nodes)

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

    ID <- xml2::xml_attr(n, "id")
    ParentID <- xml2::xml_attr(n, "parent_id")
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
                      ID = ID,
                      ParentID = ParentID,
                      eventsInside = eventsInside,
                      sampleID = sampleID,
                      FilePath = FilePath,
                      GateDepth = GateDepth,
                      origin = origin,
                      #deps = I(deps),
                      #n = n,
                      stringsAsFactors = F)
    )
  }, ...)) |> #, ...
    dplyr::bind_rows(do.call(rbind, lapply_fun(samples, function(y) {
      # roots
      data.frame(FileName = basename(xml2::xml_attr(xml2::xml_child(y, "DataSet"), "uri")),
                 PopulationFullPath = "root",
                 Parent = NA,
                 Population = "root",
                 Count = as.numeric(xml2::xml_attr(xml2::xml_child(y, "SampleNode"), "count")),
                 ParentCount = NA,
                 FractionOfParent = NA,
                 xDim = NA,
                 yDim = NA,
                 ID = NA,
                 ParentID = NA,
                 eventsInside = NA,
                 sampleID = xml2::xml_attr(xml2::xml_child(y, "DataSet"), "sampleID"),
                 FilePath = gsub("^file:", "", xml2::xml_attr(xml2::xml_child(y, "DataSet"), "uri")),
                 GateDepth = 0,
                 origin = "root",
                 #deps = I(list(character(0))),
                 #n = 0,
                 stringsAsFactors = F)
    }, ...))) |> #, ...
    dplyr::group_split(sampleID)

  # remove duplicate rows from gate+dependents
  pop_df <- purrr::map_dfr(gates_list, ~dplyr::filter(.x, !(dplyr::n() > 1 & origin == "Dependents"), .by = PopulationFullPath))

  if (any(dplyr::count(pop_df, FileName, PopulationFullPath)$n >1)) {
    stop("PopulationFullPaths not unique which cannot or should not be. Check.")
  }
  auto_paths <- shortest_unique_path(pop_df$PopulationFullPath) #fcexpr:::
  pop_df$Population <- unname(auto_paths[pop_df$PopulationFullPath])

  pop_df <-
    pop_df |>
    dplyr::left_join(group_df, by = c("sampleID")) |>
    dplyr::mutate(ws = ws_raw)
  pop_df <- pop_df[order(pop_df$FileName, pop_df$GateDepth, factor(pop_df$origin, levels = c("root", "Gate", "Dependents"))),]

  if (strip_data) {
    pop_df <- pop_df[,which(!names(pop_df) %in% c("ID", "ParentID", "sampleID", "origin", "n", "GateDepth"))]
  }

  ## get keywords to derive identity of fcs files
  keys_list <- fcexpr::wsx_get_keywords(ws = ws) # return = "data.frame"
  fcs_idents <-
    stack(fcexpr:::get_fcs_identities(keys_list[["vec"]])) |>
    dplyr::rename(identity = "values", "FileName" = ind) |>
    dplyr::mutate(FileName = as.character(FileName))
  pop_df <- dplyr::left_join(pop_df, fcs_idents, by = "FileName")

  message("FileName is as in FlowJo wsp file. Join counts to sampledescription (sd) via FileName and identity.")
  message("dplyr::left_join(counts, sd, by = c('FileName', 'identity'))")
  # check for equal file names on hdd
  check_filenames_hdd(pop_df = pop_df)

  # dup_files <- fcs_idents[duplicated(fcs_idents$identity),]
  # if (nrow(dup_files) > 1) {
  #   message("One or multiple FCS files seem to be duplicates. A data.frame with details will be returned.")
  #   print(dup_files)
  # } else {
  #   dup_files <- NULL
  # }


  if (return_stats) {
    stats_out <- do.call(rbind, lapply_fun(seq_along(samples), function(n) {
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
    if (anyNA(stats_out$value)) {
      print(tibble::as_tibble(stats_out) |> dplyr::filter(is.na(value)), n = Inf)
    }
    return(list(counts = tibble::as_tibble(pop_df), stats = tibble::as_tibble(stats_out)))
  }

  return(list(counts = tibble::as_tibble(pop_df), stats = NULL))
}

#' Get subsetted flowFrames from one or many flowjo workspaces
#'
#'
#'
#' @param wsp vector of paths to flowjo workspaces
#' @param FCS.file.folder path to folder of FCS files; has to be one common folder for all wsp
#' @param groups vector or list of groups in flowjo to consider; if a list, each index corresponds to the index in wsp
#' @param population which population (=node, =gate) to subset flowFrames one
#' @param invert_groups logical whether to invert group selection
#' @param samples vector of samples to select (names of FCS files); not-unique FCS file names across wsp cannot be selected individually yet
#' @param invert_samples logical whether to invert sample selection
#' @param inverse_transform return inverse- (T) or logicle- (F) transform or both (c(T,F))
#' @param downsample numeric, if < 0 then a fraction of each flowFrame is sampled, if > 0 an absolute number of each flowFrame is subsetted
#' @param remove_redundant_channels remove channels that are not part of the gating tree
#' @param lapply_fun lapply function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#'
#' @return
#' @export
#'
#' @examples
wsp_get_ff <- function(wsp,
                       FCS.file.folder,
                       groups = NULL,
                       population,
                       invert_groups = F,
                       samples = NULL,
                       invert_samples = F,
                       inverse_transform = c(T,F),
                       downsample = 1,
                       remove_redundant_channels = F,
                       lapply_fun = lapply,
                       ...) {

  if (!"CytoML" %in% rownames(installed.packages())) {BiocManager::install("CytoML")}
  if (!"flowWorkspace" %in% rownames(installed.packages())) {BiocManager::install("flowWorkspace")}

  lapply_fun <- match.fun(lapply_fun)

  if (missing(wsp) || class(wsp) != "character") {stop("Please provide a vector of paths to flowjo workspaces.")}
  if (!is.null(groups)) {
    if (class(groups) == "list" && length(groups) != length(wsp)) {stop("list of groups has to have the same length as wsp. Alternatively pass a vector groups to use for all workspace.")}
    if (class(groups) != "list") {groups <- rep(list(groups), length(wsp))}
  }
  if (!dir.exists(FCS.file.folder)) {stop(paste0(FCS.file.folder, " not found."))}


  smpl <- do.call(rbind, lapply(seq_along(wsp), function(x) {
    y <- wsx_get_fcs_paths(wsp[x], basename = T, split = F)
    y$wsp <- wsp[x]

    key <- sapply(wsx_get_keywords(ws), function(z) {
      z[which(z$name == "$FIL"),"value"]
    })
    y$FIL <- key[y$FilePath]

    if (!is.null(groups)) {
      if (invert_groups) {
        y <- y[which(!y$group %in% groups[[x]]),]
      } else {
        y <- y[which(y$group %in% groups[[x]]),]
      }
    }
    return(y)
  }))

  if (!is.null(samples)) {
    if (invert_samples) {
      smpl <- smpl[which(!smpl$FilePath %in% samples),]
    } else {
      smpl <- smpl[which(smpl$FilePath %in% samples),]
    }
  }

  smpl <- dplyr::distinct(smpl, FilePath, wsp, .keep_all = T)


  ff.list <- lapply_fun(split(smpl, 1:nrow(smpl)),
                        get_ff,
                        FCS.file.folder = FCS.file.folder,
                        inverse_transform = inverse_transform,
                        downsample = downsample,
                        remove_redundant_channels = remove_redundant_channels,
                        population = population)

  ffs <- sapply(ff.list, "[", 1)
  names(ffs) <- smpl$FilePath
  ffs <- lapply(seq_along(ffs[[1]]), function(x) {
    sapply(ffs, "[", x, simplify = T)
  })
  names(ffs) <- stats::setNames(c("inverse", "logicle"), c(T,F))[as.character(inverse_transform)]

  inds <- sapply(ff.list, "[", 2)
  names(inds) <- smpl$FilePath

  return(list(ffs, inds))
}


get_ff <- function (i, FCS.file.folder, inverse_transform, downsample, remove_redundant_channels, population) {
  gs <- .flowjo_to_gatingset2(ws = CytoML::open_flowjo_xml(i$wsp), name = i$group, path = FCS.file.folder, subset = i$FIL, truncate_max_range = F)

  if (remove_redundant_channels) {
    gs <- suppressMessages(flowWorkspace::gs_remove_redundant_channels(gs))
  }

  ex <- lapply(inverse_transform, function (x) {
    flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gh_pop_get_data(gs[[1]], inverse.transform = x))
  })

  inds <- flowWorkspace::gh_pop_get_indices(gs[[1]], y = population)
  inds_non_subset <- inds

  s <- if (downsample < 1) {
    sort(sample(which(inds), ceiling(length(which(inds))*downsample)))
  } else if (downsample > 1) {
    sort(sample(which(inds), min(length(which(inds)),downsample)))
  }
  inds[which(inds)[!which(inds) %in% s]] <- F

  if (length(which(inds)) > 1) {
    for (x in seq_along(ex)) {
      ex[[x]] <- subset(ex[[x]], inds)
    }
  }

  flowWorkspace::gs_cleanup_temp(gs)
  return(list(ex, inds_non_subset))
}

.flowjo_to_gatingset2 <- function (ws, name = NULL, subset = list(), execute = TRUE,
                                   path = "", cytoset = NULL, backend_dir = tempdir(), backend = flowWorkspace::get_default_backend(),
                                   includeGates = TRUE, additional.keys = "$TOT", additional.sampleID = FALSE,
                                   keywords = character(), keywords.source = "XML", keyword.ignore.case = FALSE,
                                   extend_val = 0, extend_to = -4000, channel.ignore.case = FALSE,
                                   leaf.bool = TRUE, include_empty_tree = FALSE, skip_faulty_gate = FALSE,
                                   compensation = NULL, transform = TRUE, fcs_file_extension = ".fcs",
                                   greedy_match = FALSE, mc.cores = 1, ...) {
  if (is.null(cytoset))
    cytoset <- flowWorkspace::cytoset()
  backend <- match.arg(backend, c("h5", "tile"))
  g <- CytoML::fj_ws_get_sample_groups(ws)
  groups <- g[!duplicated(g$groupName), ]
  groups <- groups[order(groups$groupID), "groupName"]
  if (is.null(name)) {
    groupInd <- menu(groups, graphics = FALSE, "Choose which group of samples to import:")
  } else if (is.numeric(name)) {
    if (length(groups) < name)
      stop("Invalid sample group index.")
    groupInd <- name
  } else if (is.character(name)) {
    if (is.na(match(name, groups)))
      stop("Invalid sample group name.")
    groupInd <- match(name, groups)
  }

  if (is(subset, "character")) {
    subset <- list(name = subset)
  }
  if (!is(subset, "list")) {
    stop("invalid 'subset' argument!")
  }
  if (is.null(additional.keys)) {additional.keys <- character(0)}
  if (is.null(path)) {path <- ""}
  args <- list(...)
  if (!is.null(args[["isNcdf"]])) {
    warning("'isNcdf' argument is deprecated!Data is always stored in h5 format by default!")
    args[["isNcdf"]] <- NULL
  }
  if (is.null(compensation)) {
    compensation <- list()
  } else {
    if (is.list(compensation) && !is.data.frame(compensation)) {
      compensation <- sapply(compensation, CytoML:::check_comp,
                             simplify = FALSE)
    } else compensation <- CytoML:::check_comp(compensation)
  }
  args <- list(ws = ws@doc, group_id = groupInd - 1, subset = subset,
               execute = execute, path = suppressWarnings(normalizePath(path)),
               cytoset = cytoset@pointer, backend_dir = suppressWarnings(normalizePath(backend_dir)),
               backend = backend, includeGates = includeGates, additional_keys = additional.keys,
               additional_sampleID = additional.sampleID, keywords = keywords,
               is_pheno_data_from_FCS = keywords.source == "FCS", keyword_ignore_case = keyword.ignore.case,
               extend_val = extend_val, extend_to = extend_to, channel_ignore_case = channel.ignore.case,
               leaf_bool = leaf.bool, include_empty_tree = include_empty_tree,
               skip_faulty_gate = skip_faulty_gate, comps = compensation,
               transform = transform, fcs_file_extension = fcs_file_extension,
               greedy_match = greedy_match, fcs_parse_arg = args, num_threads = mc.cores)
  p <- do.call(CytoML:::parse_workspace, args)
  gs <- new("GatingSet", pointer = p)
  gslist <- suppressMessages(flowWorkspace:::gs_split_by_tree(gs))
  if (length(gslist) > 1) {
    msg <- "GatingSet contains different gating tree structures and must be cleaned before using it!\n "
    if (grepl("all samples", groups[groupInd], ignore.case = TRUE)) {
      msg <- c(msg, "It seems that you selected the 'All Samples' group,",
               " which is a generic group and typically contains samples with different gating schemes attached.",
               "Please choose a different sample group and try again.")
    }
    warning(msg)
  }
  gs
}

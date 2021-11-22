#' Get subsetted flowFrames from FCS files in one or many flowjo workspaces (wsp)
#'
#'
#'
#' @param wsp vector of paths to flowjo workspaces
#' @param FCS.file.folder path to folder(s) of FCS files; may be one path for all wsp or a vector of paths, one for each wsp;
#' if not provided fcs file paths are derived individually from the wsp (xml)
#' @param groups vector or list of groups in flowjo to consider; if a list, each index corresponds to the index in wsp;
#' if NULL samples from all groups are read
#' @param population which population (=node, =gate) to subset flowFrames one; use wsx_get_poppaths to get paths
#' @param invert_groups logical whether to invert group selection
#' @param samples vector or list of samples to select (names of FCS files), each index corresponds to the index in wsp;
#' if NULL all samples (from selected groups) are read
#' @param invert_samples logical whether to invert sample selection
#' @param inverse_transform return inverse- (T) or logicle- (F) transform or both (c(T,F))
#' @param downsample numeric, if < 0 then a fraction of each flowFrame is sampled, if > 0 an absolute number of each flowFrame is subsetted
#' @param remove_redundant_channels remove channels that are not part of the gating tree, mainly to reduce memory load
#' @param lapply_fun lapply function name, unquoted; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#'
#' @return
#' @export
#'
#' @examples
wsp_get_ff <- function(wsp,
                       FCS.file.folder = NULL,
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

  if (!requireNamespace("CytoML", quietly = T)){
    BiocManager::install("CytoML")
  }
  if (!requireNamespace("flowWorkspace", quietly = T)){
    BiocManager::install("flowWorkspace")
  }
  lapply_fun <- match.fun(lapply_fun)


  checked_in <- check_in(wsp = wsp, groups = groups, samples = samples, FCS.file.folder = FCS.file.folder)
  groups <- checked_in[["groups"]]
  samples <- checked_in[["samples"]]
  FCS.file.folder <- checked_in[["FCS.file.folder"]]

  smpl <- get_smpl_df(wsp = wsp, groups = groups, invert_groups = invert_groups, samples = samples, invert_samples = invert_samples, FCS.file.folder = FCS.file.folder)
  if (is.null(smpl)) {
    return(NULL)
  }
  # remove doublets due to "All Samples" association
  smpl <- dplyr::distinct(smpl, FilePath, wsp, .keep_all = T)

  if (any(table(smpl$FilePath) > 1)) {
    print("Same FCS files found in multiple workspaces. This cannot be handled. Please provide the samples and/or groups argument or fix manually.")
    stop(print(smpl$FilePath[which(table(smpl$FilePath) > 1)]))
  }

  # check if population exists for each sample
  pp <- do.call(rbind, lapply(wsp, function(x) {
    wsx_get_poppaths(x, collapse = F)
  }))
  pp <- pp[which(pp$FileName %in% smpl$FileName),]
  pp <- pp %>%
    dplyr::group_by(PopulationFullPath, Population, ws) %>%
    dplyr::summarise(FileName = list(FileName), .groups = "drop")
  pp <- as.data.frame(pp)
  pp_is <- unique(unlist(pp[unique(c(which(pp$PopulationFullPath == population), which(pp$Population == population))), "FileName"]))
  if (length(smpl$FileName[which(!smpl$FileName %in% pp_is)]) > 0) {
    print(paste0("For ", paste(smpl$FileName[which(!smpl$FileName %in% pp_is)], collapse = ", "), " population was not found."))
  }
  smpl <- smpl[which(smpl$FileName %in% pp_is),]

  # writing h5 files to disk may be the limiting factor, so doing this by multicore may not improve speed
  ff.list <- lapply_fun(split(smpl, 1:nrow(smpl)),
                        get_ff,
                        inverse_transform = inverse_transform,
                        downsample = downsample,
                        remove_redundant_channels = remove_redundant_channels,
                        population = population,
                        ...)

  ffs <- sapply(ff.list, "[", 1)
  names(ffs) <- smpl$FileName
  ffs <- lapply(seq_along(ffs[[1]]), function(x) {
    sapply(ffs, "[", x, simplify = T)
  })
  names(ffs) <- stats::setNames(c("inverse", "logicle"), c(T,F))[as.character(inverse_transform)]

  inds <- sapply(ff.list, "[", 2)
  names(inds) <- smpl$FileName

  return(list(ffs, inds))
}

get_ff <- function(x, inverse_transform, downsample, remove_redundant_channels, population) {

  # one file at a time avoids problems due to different gating trees, but this may leave unintentional different gating trees undetected
  # pass full path as attr and check consistency later?

  if (nrow(x) > 1) {
    stop("Only one fcs file at a time.")
  }

  if (is.na(x$FCS.file.folder)) {
    #path <- x[,which(names(x) %in% c("sampleID", "FilePath")),drop=F]
    #names(path)[which(names(path) == "FilePath")] <- "file"
    path <- dirname(x$FilePath)
    if (!file.exists(path)) {
      stop(paste0(path, " not found. Was the workspace saved on another computer? If so, reconnect FCS files in flowjo or provide the FCS.file.folder(s) on the current computer."))
    }
  } else {
    path <- x$FCS.file.folder
  }


  gs <- CytoML::flowjo_to_gatingset(ws = CytoML::open_flowjo_xml(x$wsp), name = x$group, path = path, subset = `$FIL` == x$FIL, truncate_max_range = F, keywords = "$FIL")

  if (remove_redundant_channels) {
    gs <- suppressMessages(flowWorkspace::gs_remove_redundant_channels(gs))
  }

  ex <- lapply(inverse_transform, function (y) {
    flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gh_pop_get_data(gs[[1]], inverse.transform = y))
  })

  inds <- flowWorkspace::gh_pop_get_indices(gs[[1]], y = population)
  s <- if (downsample < 1) {
    sort(sample(which(inds), ceiling(length(which(inds))*downsample)))
  } else if (downsample > 1) {
    sort(sample(which(inds), min(length(which(inds)),downsample)))
  } else {
    which(inds)
  }
  inds[which(inds)[!which(inds) %in% s]] <- F


  for (i in seq_along(ex)) {
    ex[[i]] <- subset(ex[[i]], inds)
  }

  flowWorkspace::gs_cleanup_temp(gs)
  return(list(ex, inds))
}


check_in <- function(wsp, samples, groups, FCS.file.folder) {
  if (missing(wsp) || class(wsp) != "character") {stop("Please provide a vector of paths to flowjo workspaces.")}
  if (!is.null(groups)) {
    if (class(groups) == "list" && length(groups) != length(wsp)) {stop("list of groups has to have the same length as wsp. Alternatively pass a vector groups to use for all workspace.")}
    if (class(groups) != "list") {groups <- rep(list(groups), length(wsp))}
  }
  if (!is.null(samples)) {
    if (class(samples) == "list" && length(samples) != length(wsp)) {stop("list of samples has to have the same length as wsp. Alternatively pass a vector samples to use for all workspace.")}
    if (class(samples) != "list") {samples <- rep(list(samples), length(wsp))}
  }
  if (!is.null(FCS.file.folder)) {
    if (any(!dir.exists(FCS.file.folder))) {stop(paste0(FCS.file.folder[which(!dir.exists(FCS.file.folder))], " not found."))}
    if (length(FCS.file.folder) != length(wsp)) {stop("FCS.file.folder has to have the same length as wsp or 1.")}
    if (length(FCS.file.folder) == 1) {FCS.file.folder <- rep(FCS.file.folder, length(wsp))}
  }
  return(list(groups = groups, samples = samples, FCS.file.folder = FCS.file.folder))
}


get_smpl_df <- function(wsp, groups, invert_groups, samples, invert_samples, FCS.file.folder) {
  smpl <- do.call(rbind, lapply(seq_along(wsp), function(x) {
    y <- wsx_get_fcs_paths(wsp[x], split = F)
    y$wsp <- wsp[x]
    y$FileName <- basename(y$FilePath)

    key <- sapply(wsx_get_keywords(wsp[x]), function(z) {
      z[which(z$name == "$FIL"),"value"]
    })
    y$FIL <- key[y$FileName]

    if (!is.null(groups)) {
      if (invert_groups) {
        y <- y[which(!y$group %in% groups[[x]]),]
      } else {
        y <- y[which(y$group %in% groups[[x]]),]
      }
    }
    if(nrow(y) == 0) {
      return(NULL)
    }

    if (!is.null(samples)) {
      if (invert_samples) {
        y <- y[which(!y$FileName %in% samples[[x]]),]
      } else {
        y <- y[which(y$FileName %in% samples[[x]]),]
      }
    }
    if(nrow(y) == 0) {
      return(NULL)
    }

    if (is.null(FCS.file.folder)) {
      y$FCS.file.folder <- NA
    } else {
      y$FCS.file.folder <- FCS.file.folder[x]
    }

    return(y)
  }))
  return(smpl)
}

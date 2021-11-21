#' Get subsetted flowFrames from FCS files in one or many flowjo workspaces (wsp)
#'
#'
#'
#' @param wsp vector of paths to flowjo workspaces
#' @param FCS.file.folder path to folder(s) of FCS files; may be one path for all wsp or a vector of paths, one for each wsp;
#' if not provided fcs file paths are derived individually from the wsp (xml)
#' @param groups vector or list of groups in flowjo to consider; if a list, each index corresponds to the index in wsp
#' @param population which population (=node, =gate) to subset flowFrames one; use wsx_get_poppaths (make fun!) to get poppaths
#' @param invert_groups logical whether to invert group selection
#' @param samples vector or list of samples to select (names of FCS files)
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

  if (!"CytoML" %in% rownames(installed.packages())) {BiocManager::install("CytoML")}
  if (!"flowWorkspace" %in% rownames(installed.packages())) {BiocManager::install("flowWorkspace")}

  lapply_fun <- match.fun(lapply_fun)

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

    if (!is.null(samples)) {
      if (invert_samples) {
        y <- y[which(!y$FileName %in% samples),]
      } else {
        y <- y[which(y$FileName %in% samples),]
      }
    }

    if (is.null(FCS.file.folder)) {
      y$FCS.file.folder <- NA
    }

    return(y)
  }))

  # remove doublets due to "All Samples" association
  smpl <- dplyr::distinct(smpl, FilePath, wsp, .keep_all = T)

  if (any(table(smpl$FilePath) > 1)) {
    print("Same FCS files found in multiple workspaces. This cannot be handled. Please provide the samples argument or fix manually.")
    stop(print(smpl$FilePath[which(table(smpl$FilePath) > 1)]))
  }

  # writing h5 files to disk may be the limiting factor, so doing this by multicore may not imrove speed
  ff.list <- lapply_fun(split(smpl, 1:nrow(smpl)),
                        get_ff,
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

get_ff <- function (x, inverse_transform, downsample, remove_redundant_channels, population) {

  # one file at a time avoids problems due to different gating trees, but this may leave unintentional different gating trees undetected
  # pass full path as attr and check consistency later

  if (nrow(x) > 1) {
    stop("Only one fcs file at a time.")
  }

  if (is.na(x$FCS.file.folder)) {
    path <- x[,which(names(x) %in% c("sampleID", "FilePath")),drop=F]
    names(path)[which(names(path) == "FilePath")] <- "file"
    if (!file.exists(path$file)) {
      stop(paste0(path$file, " not found. Was the workspace saved on another computer? If so, reconnect FCS files in flowjo or provdide the FCS.file.folder(s) on the current computer."))
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
  }
  inds[which(inds)[!which(inds) %in% s]] <- F

  if (length(which(inds)) > 1) {
    for (i in seq_along(ex)) {
      ex[[i]] <- subset(ex[[i]], inds)
    }
  }

  flowWorkspace::gs_cleanup_temp(gs)
  return(list(ex, inds))
}


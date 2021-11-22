#' Get (subsetted) flowFrames from FCS files in one or many flowjo workspaces
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
#' @return a list of flowframes
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

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("CytoML", quietly = T)){
    BiocManager::install("CytoML")
  }
  if (!requireNamespace("flowWorkspace", quietly = T)){
    BiocManager::install("flowWorkspace")
  }
  lapply_fun <- match.fun(lapply_fun)


  checked_in <- check_in(wsp = wsp, groups = groups, samples = samples, FCS.file.folder = FCS.file.folder, inverse_transform = inverse_transform)
  groups <- checked_in[["groups"]]
  samples <- checked_in[["samples"]]
  FCS.file.folder <- checked_in[["FCS.file.folder"]]

  smpl <- get_smpl_df(wsp = wsp, groups = groups, invert_groups = invert_groups, samples = samples, invert_samples = invert_samples, FCS.file.folder = FCS.file.folder)
  if (is.null(smpl)) {
    return(NULL)
  }

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

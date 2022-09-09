#' Get (subsetted) flowFrames from one or many flowjo workspaces
#'
#' From a flowjo workspace with gated populations and the respective fcs files flowframes are generated.
#' Under the hood CytoML::flowjo_to_gatingset applies the geometric gate definitions and filters
#' relevant rows (indices) of fcs files. The compensation matrix as defined in flowjo will be used
#' to compensate fluorescence intensities.
#'
#' If it is intended to pass flowframes to fcexpr::dr_to_fcs, it is recommended to have transformed
#' and untransformed expression values returned.
#'
#' @param wsp vector of paths to flowjo workspaces
#' @param FCS.file.folder path to folder(s) of FCS files; may be one path for all wsp or a vector of paths, one for each wsp;
#' if not provided fcs file paths are derived individually from the wsp (xml). If the workspace was generated and saved on
#' another computer you will have to provide the path to FCS files on the current computer.
#' @param groups vector or list of groups in flowjo to consider; if a list, each index corresponds to the same index in wsp;
#' if NULL samples from all groups are read
#' @param population which population (=node, =gate) to subset flowFrames on; use fcexpr::wsx_get_poppaths to get paths
#' @param invert_groups logical whether to invert group selection
#' @param samples vector or list of samples to select (names of FCS files), each index corresponds to the index in wsp;
#' if NULL all samples (from selected groups) are read
#' @param invert_samples logical whether to invert sample selection
#' @param downsample numeric, if < 0 then a fraction of events is sampled, if > 0 an absolute number of events is sampled; or set to "min"
#' which will lead to downsampling each flowframe to the number of events in the flowframe with lowest number of events; can be a single value to treat all
#' FCS files equally or can be a vector of same length as FCS files
#' @param remove_redundant_channels remove channels that are not part of the gating tree, mainly to reduce memory load
#' @param lapply_fun lapply function name, unquoted; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#' @param return_untransformed logical; do return untransformed (inverse) data
#' @param return_logicle_transformed logical; do return logicle-transformed data
#' @param seed set a seed to reproduce downsampling
#'
#' @return a list of (subsetted) flowframes with events that are within the gated population only
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#'\dontrun{
#' ff_list <- fcexpr::wsp_get_ff(wsp = "mypath/my.wsp", population = "CD8+")
#' # ff.list[[1]] may be passed to fcexpr::dr_to_fcs for instance
#'}
wsp_get_ff <- function(wsp,
                       FCS.file.folder = NULL,
                       groups = NULL,
                       population,
                       invert_groups = F,
                       samples = NULL,
                       invert_samples = F,
                       return_untransformed = T,
                       return_logicle_transformed = T,
                       downsample = 1,
                       remove_redundant_channels = F,
                       lapply_fun = lapply,
                       seed = 42,
                       ...) {

  if (!requireNamespace("BiocManager", quietly = T)) {
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("CytoML", quietly = T)) {
    BiocManager::install("CytoML")
  }
  if (!requireNamespace("flowWorkspace", quietly = T)) {
    BiocManager::install("flowWorkspace")
  }
  lapply_fun <- match.fun(lapply_fun)

  checked_in <- check_in(wsp = wsp, groups = groups, samples = samples, FCS.file.folder = FCS.file.folder,
                         return_untransformed = return_untransformed,
                         return_logicle_transformed = return_logicle_transformed)
  groups <- checked_in[["groups"]]
  samples <- checked_in[["samples"]]
  FCS.file.folder <- checked_in[["FCS.file.folder"]]

  smpl <- get_smpl_df(wsp = wsp,
                      groups = groups,
                      invert_groups = invert_groups,
                      samples = samples,
                      invert_samples = invert_samples,
                      FCS.file.folder = FCS.file.folder)
  if (is.null(smpl)) {
    return(NULL)
  }

  if (is.numeric(downsample)) {
    ds <- downsample
  } else if (all(downsample == "min")) {
    ds <- 1
  } else {
    stop("downsample has to be numeric or 'min'. With min all flowframes will be downsampled to that flowframe with the lowest number of events.")
  }

  # check length of downsample equal to length of ind_mat or equal to 1
  if (length(ds) != 1 && length(ds) != nrow(smpl)) {
    stop("downsample has to have length 1 or length of ind_mat (one value for each FCS file).")
  }

  if (length(ds) != 1) {
    smpl$downsample <- 1
    for (x in 1:nrow(smpl)) {
      smpl$downsample[x] <- ds[x]
    }
  }

  if (any(table(smpl$FilePath) > 1)) {
    warning("Same FCS files found in multiple workspaces. This cannot be handled. Please provide the samples and/or groups argument or fix manually.")
    stop(smpl$FilePath[which(table(smpl$FilePath) > 1)])
  }

  pp <- do.call(rbind, lapply(wsp, function(x) {
    wsx_get_poppaths(x, collapse = F)
  }))
  pp <- pp[which(pp$FileName %in% smpl$FileName),]
  pp <-
    pp %>%
    dplyr::group_by(PopulationFullPath, Population, ws) %>%
    dplyr::summarise(FileName = list(FileName), .groups = "drop")
  pp <- as.data.frame(pp)
  pp_is <- unique(unlist(pp[unique(c(which(pp$PopulationFullPath == population), which(pp$Population == population))), "FileName"]))
  if (is.null(pp_is)) {
    stop("Population was not found for any sample.")
  }
  if (length(smpl$FileName[which(!smpl$FileName %in% pp_is)]) > 0) {
    message("For ", paste(smpl$FileName[which(!smpl$FileName %in% pp_is)], collapse = ", "), " population was not found.")
  }
  smpl <- smpl[which(smpl$FileName %in% pp_is),]

  ff.list <- lapply_fun(split(smpl, 1:nrow(smpl)),
                        get_ff,
                        return_untransformed = return_untransformed,
                        return_logicle_transformed = return_logicle_transformed,
                        downsample = ds,
                        remove_redundant_channels = remove_redundant_channels,
                        population = population,
                        seed = seed,
                        ...)

  if (all(downsample == "min")) {
    min <- min(unlist(lapply(sapply(sapply(ff.list, "[", 1), "[", 1), nrow)))
    for (x in seq_along(ff.list)) {
      inds <- sample(c(rep(T, min), rep(F, nrow(ff.list[[x]][[1]][[1]])-min)))
      for (y in seq_along(ff.list[[x]][[1]])) {
        ff.list[[x]][[1]][[y]] <- subset(ff.list[[x]][[1]][[y]], inds)
      }
    }
  }

  ffs <- sapply(ff.list, "[", 1)
  names(ffs) <- smpl$FileName

  ffs <- lapply(seq_along(ffs[[1]]), function(x) sapply(ffs, "[", x, simplify = T))

  if (return_untransformed && !return_logicle_transformed) {
    names(ffs) <- "untransformed"
  } else if (!return_untransformed && return_logicle_transformed) {
    names(ffs) <- "transformed"
  } else if (return_untransformed && return_logicle_transformed) {
    names(ffs) <- c("transformed", "untransformed")
  }

  #names(ffs) <- stats::setNames(c("inverse", "logicle"), c(T,F))[as.character(inverse_transform)]

  inds <- sapply(ff.list, "[", 2)
  names(inds) <- smpl$FileName

  return(list(flowframes = ffs, indices = inds))
}

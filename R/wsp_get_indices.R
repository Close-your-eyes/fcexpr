#' Get indices for all gated populations in a flowjo workspace
#'
#' Apply geometric definitions with CytoML::flowjo_to_gatingset and save the respective indices. This may take a little
#' bit of time as .h5 files are written to disk for every fcs file. Obtained indices for selected populations may then be applied to fcs files
#' by ...
#'
#'
#' @param wsp vector of paths to flowjo workspaces
#' @param FCS.file.folder path to folder(s) of FCS files; may be one path for all wsp or a vector of paths, one for each wsp;
#' if not provided fcs file paths are derived individually from the wsp (xml)
#' @param groups vector or list of groups in flowjo to consider; if a list, each index corresponds to the index in wsp;
#' if NULL samples from all groups are read
#' @param invert_groups logical whether to invert group selection
#' @param samples vector or list of samples to select (names of FCS files), each index corresponds to the index in wsp;
#' if NULL all samples (from selected groups) are read
#' @param lapply_fun lapply function name, unquoted; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#'
#' @return list of of matrices, one entry for each selected sample
#' @export
#'
#' @examples
wsp_get_indices <- function(wsp,
                            FCS.file.folder = NULL,
                            groups = NULL,
                            invert_groups = F,
                            samples = NULL,
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
  # remove duplicates due to "All Samples" association, which group to keep does not matter
  smpl <- dplyr::distinct(smpl, FilePath, wsp, .keep_all = T)

  if (any(table(smpl$FilePath) > 1)) {
    print("Same FCS files found in multiple workspaces. This cannot be handled. Please provide the samples and/or groups argument or fix manually.")
    stop(print(smpl$FilePath[which(table(smpl$FilePath) > 1)]))
  }

  ind.list <- lapply_fun(split(smpl, 1:nrow(smpl)),
                         get_inds,
                         ...)
  names(ind.list) <- smpl$FileName
  return(ind.list)
}

get_inds <- function(x) {
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
  inds <- flowWorkspace::gh_pop_get_indices_mat(gs[[1]], y = gh_get_pop_paths(gs[[1]]))
  attr(inds, "short_names") <- stats::setNames(shortest_unique_path(colnames(inds)), nm = colnames(inds))
  attr(inds, "ws") <- x$wsp
  attr(inds, "FilePath") <- x$FilePath

  return(inds)
}





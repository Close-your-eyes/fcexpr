#' Get indices for all gated populations in a flowjo workspace
#'
#' Every row in a fcs file represents an event. Every column is a parameter (channel). Gates select a subset (rows) of events by applying limits to usually one or two channels.
#' A sub-population of gated events hence may be defined by a vector of (row-) indices (TRUE - events is included in gate; FALSE - event is not included). A whole gating tree
#' may be represented by a matrix with n columns for n gates and m rows for m events. The output of this function may be saved to disk and applied to fcs files with fcexpr::inds_get_ff
#' in order to obtain subsetted flowfframes representing gated populations in flowjo.
#'
#' Geometric gate definitions from flowjo are applied with CytoML::flowjo_to_gatingset and indices matrices are obtained with flowWorkspace::gh_pop_get_indices_mat.
#' This process may take a while depend upon size of fcs files as a .h5 file is written to disk for every fcs file before indices can be derived. Hence, it is recommended
#' to save the indices-matrices in case of large FCS files.
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
#' @param invert_samples logical whether to invert sample selection
#'
#' @return list of of matrices, one for each selected sample
#' @export
#'
#' @examples
wsp_get_indices <- function(wsp,
                            FCS.file.folder = NULL,
                            groups = NULL,
                            invert_groups = F,
                            samples = NULL,
                            invert_samples = F,
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

  checked_in <- check_in(wsp = wsp, groups = groups, samples = samples, FCS.file.folder = FCS.file.folder, inverse_transform = c(T,F))
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

  ind.list <- lapply_fun(split(smpl, 1:nrow(smpl)),
                         get_inds,
                         ...)
  names(ind.list) <- smpl$FileName
  return(ind.list)
}




#' Get (subsetted) gatingsets from flowjo workspaces
#'
#' @param wsp vector of paths to flowjo workspaces
#' @param groups vector or list of groups in flowjo to consider; if a list, each index corresponds to the index in wsp;
#' if NULL samples from all groups are read
#' @param FCS.file.folder path to folder(s) of FCS files; may be one path for all wsp or a vector of paths, one for each wsp;
#' if not provided (NULL) fcs file paths are derived individually from wsps (xml)
#' @param invert_groups logical whether to invert group selection
#' @param samples vector or list of samples to select (names of FCS files), each index corresponds to the index in wsp;
#' if NULL all samples (from selected groups) are read
#' @param invert_samples logical whether to invert sample selection
#' @param remove_redundant_channels remove channels that are not part of the gating tree, mainly to reduce memory load
#'
#' @return list of gatingsets
#' @export
#'
#' @examples
#'\dontrun{
#' gs_list <- fcexpr::wsp_get_gs(wsp = "mypath/my.wsp")
#'}
wsp_get_gs <- function(wsp,
                       FCS.file.folder = NULL,
                       groups = NULL,
                       invert_groups = F,
                       samples = NULL,
                       invert_samples = F,
                       remove_redundant_channels = F) {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("CytoML", quietly = T)){
    BiocManager::install("CytoML")
  }
  if (!requireNamespace("flowWorkspace", quietly = T)){
    BiocManager::install("flowWorkspace")
  }


  checked_in <- check_in(wsp = wsp, samples = samples, groups = groups, FCS.file.folder = FCS.file.folder)

  groups <- checked_in[["groups"]]
  samples <- checked_in[["samples"]]
  FCS.file.folder <- checked_in[["FCS.file.folder"]]

  smpl <- get_smpl_df(wsp = wsp, groups = groups, invert_groups = invert_groups, samples = samples, invert_samples = invert_samples, FCS.file.folder = FCS.file.folder)
  if (is.null(smpl)) {
    return(NULL)
  }
  if (any(table(smpl$FilePath) > 1)) {
    print(smpl$FilePath[which(table(smpl$FilePath) > 1)])
    stop("Same FCS files found in multiple workspaces or groups. This cannot be handled. Please provide the samples and/or groups argument or fix manually.")
  }

  smpl_list <- split(smpl, paste0(basename(smpl$wsp), "_-_", smpl$group))
  smpl_list <- lapply(smpl_list, function(x) {
    if (any(is.na(x$FCS.file.folder))) {
      path <- x$FilePath
      n <- length(unique(path))
      while (n != 1) {
        path <- dirname(path)
        n <- length(unique(path))
      }
      message("Common FCS.file.folder determined: ", unique(path), ".")
      x$FCS.file.folder <- path
    }
    return(x)
  })

  gs_list <- lapply(smpl_list,
                    get_gs,
                    remove_redundant_channels = remove_redundant_channels)
  names(gs_list) <- names(smpl_list)

  return(gs_list)
}



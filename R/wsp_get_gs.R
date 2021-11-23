#'
#'
#' @param wsp
#' @param groups
#' @param FCS.file.folder
#' @param invert_groups
#' @param samples
#' @param invert_samples
#' @param remove_redundant_channels
#'
#' @return
#' @export
#'
#' @examples
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

  checked_in <- check_in(wsp = wsp, samples = samples, groups = groups, FCS.file.folder = FCS.file.folder, inverse_transform = c(T,F))
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

  smpl_list <- split(smpl, paste0(basename(smpl$wsp), "_-_", smpl$group))
  smpl_list <- lapply(smpl_list, function(x) {
    if (any(is.na(x$FCS.file.folder))) {
      path <- x$FilePath
      n <- length(unique(path))
      while (n != 1) {
        path <- dirname(path)
        n <- length(unique(path))
      }
    }
    print(paste0("Common FCS.file.folder determined: ", unique(path), "."))
    x$FCS.file.folder <- path
    return(x)
  })

  gs_list <- lapply(smpl_list,
                    get_gs,
                    remove_redundant_channels = remove_redundant_channels)
  names(gs_list) <- names(smpl_list)

  return(gs_list)
}



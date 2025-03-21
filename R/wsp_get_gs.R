#' Get (subsetted) gatingsets from flowjo workspaces
#'
#' Wsp files are checked for corresponding gating hierarchies and FCS files are
#' read into separate GatingSets (gs), respectively. That means FCS files from different
#' wsp are combined if their gatings are equal. But simple import from one wsp is
#' also possible. If file paths to FCS files in flowjo are not valid on the current computer
#' (e.g. because flowjo analysis was done on another computer), you have to provide the
#' local top folder which contains all of the necessary FCS files: FCS.files.folder.
#' Therein, the files are searched for recursively. Samples from wsp(s) with missing FCS files are skipped.
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
#' @return list of gatingsets, FCS files used, gating hierarchies
#' @export
#'
#' @examples
#'\dontrun{
#' # check fcs files in wsp first
#' samples <- fcexpr::wsx_get_fcs_paths(ws = "mypath/my.wsp", split = F, filter_AllSamples = T)
#' # read-in Gating set with selected samples
#' gs_list <- fcexpr::wsp_get_gs(wsp = "mypath/my.wsp", FCS.files.folder = "myLocalTopFolder/subfolder",
#' groups = "Compensation", invert_groups = T, samples = samples$FileName[1:5])
#'}
wsp_get_gs <- function(wsp,
                       FCS.file.folder = NULL,
                       groups = NULL,
                       invert_groups = F,
                       samples = NULL,
                       invert_samples = F,
                       remove_redundant_channels = F
) {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("CytoML", quietly = T)){
    BiocManager::install("CytoML")
  }
  if (!requireNamespace("flowWorkspace", quietly = T)){
    BiocManager::install("flowWorkspace")
  }

  checked_in <- check_in(
    wsp = wsp,
    samples = samples,
    groups = groups,
    FCS.file.folder = FCS.file.folder
  )

  smpl <- get_smpl_df(
    wsp = wsp,
    groups = checked_in[["groups"]],
    invert_groups = invert_groups,
    samples = checked_in[["samples"]],
    invert_samples = invert_samples,
    FCS.file.folder = checked_in[["FCS.file.folder"]]
  )
  if (is.null(smpl)) {
    return(NULL)
  }

  # check gating trees for equality
  gatings_list <- purrr::list_flatten(purrr::map(wsp, function(ws) {
    ps <-
      wsx_get_popstats(ws = ws, return_stats = F) |>
      dplyr::filter(identity %in% smpl$identity)
    gatings_list <- purrr::map(setNames(unique(ps$identity), unique(ps$identity)), function(x) {
      fcexpr::gating_tree_plot(ps |> dplyr::filter(identity == x) |> dplyr::pull("PopulationFullPath"))
    })
    return(gatings_list)
  }))

  gatings_df <-
    dplyr::bind_rows(purrr::map(gatings_list, `[[`, 3), .id = "identity") |>
    dplyr::left_join(ps |>
                       dplyr::select(PopulationFullPath, identity, xDim, yDim) |>
                       dplyr::rename("xDimfrom" = xDim, "yDimfrom" = yDim),
                     by = c("identity" = "identity", "from" = "PopulationFullPath")) |>
    dplyr::left_join(ps |>
                       dplyr::select(PopulationFullPath, identity, xDim, yDim) |>
                       dplyr::rename("xDimto" = xDim, "yDimto" = yDim),
                     by = c("identity" = "identity", "to" = "PopulationFullPath"))
  gatings_df <- split(gatings_df, gatings_df$identity)
  gatings_df <- lapply(gatings_df, function(x) {
    x <-
      x |>
      dplyr::select(-identity) |>
      dplyr::arrange(from, to)
    #x <- dplyr::select(x, from, to) # ignore channels for splitting?!
    rownames(x) <- seq(1, nrow(x))
    return(x)
  })
  dfhashes <-
    stack(purrr::map_chr(gatings_df, digest::digest, algo = "sha256")) |>
    dplyr::group_by(values) |>
    dplyr::mutate(gatinggroup = dplyr::cur_group_id()) |>
    dplyr::rename("hash" = values, "identity" = ind) |>
    dplyr::mutate(identity = as.character(identity))
  smpl <- dplyr::left_join(smpl, dfhashes, by = "identity")

  #waldo::compare(gatings_list[[1]], gatings_list[[2]])

  smpl_list <- split(smpl, smpl$gatinggroup)
  names(smpl_list) <- paste(names(smpl_list),
                            purrr::map_chr(smpl_list, function(x) {
                              paste(paste(unique(basename(x$wsp)), sep = "_"), paste(unique(x$FlowJoGroup), sep = "_"), sep = "_")
                            }),
                            sep = "_")
  if (length(length(smpl_list)) > 1) {
    message("Splitting fcs files into ", length(smpl_list), " groups of different gating hierarchies.")
    for (i in smpl_list) {
      message(paste(i$FlowJoFileName, collapse = ", "), "\n")
    }
  }

  gs_list <- lapply(smpl_list,
                    get_gs,
                    remove_redundant_channels = remove_redundant_channels)
  names(gs_list) <- names(smpl_list)

  return(list(
    gs_list = gs_list,
    FCS_files = smpl,
    gating_hierachies = gatings_list,
    gatings_compare = gatings_df
  ))
}



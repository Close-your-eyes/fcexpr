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
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param split_size chunk size of splits of fcs files parsed to CytoML::flowjo_to_gatingset;
#' splitting allows to have a progress bar from pblapply or to use multicore processing by parallel::mclapply
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
                       remove_redundant_channels = F,
                       lapply_fun = lapply
                       #split_size = 2,
                       #additional.sampleID = F,
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

  lapply_fun <- match.fun(lapply_fun)

  checked_in <- check_in(wsp = wsp, samples = samples, groups = groups, FCS.file.folder = FCS.file.folder)

  smpl <- get_smpl_df(
    wsp = wsp,
    groups = checked_in[["groups"]],
    invert_groups = invert_groups,
    samples = checked_in[["samples"]],
    invert_samples = invert_samples,
    FCS.file.folder = checked_in[["FCS.file.folder"]],
    lapply_fun = lapply_fun
  )
  if (is.null(smpl)) {
    return(NULL)
  }
  # check gating trees for equality
  ps <-
    wsx_get_popstats(ws = wsp, return_stats = F) |>
    dplyr::filter(identity %in% smpl$identity)
  gatings_list <- purrr::map(setNames(unique(ps$identity), unique(ps$identity)), function(x) {
    fcexpr::gating_tree_plot(ps |> dplyr::filter(identity == x) |> dplyr::pull("PopulationFullPath"))
  })
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
  if (any(table(smpl$identity) > 1)) {
    message("Same FCS files found in multiple workspaces. This may cause problems. Provide samples argument as a list of vectors, one for each wsp and filter the duplicate files.")
    message(paste(names(table(smpl$identity))[which(table(smpl$identity) > 1)], collapse = ", "))
  }

  smpl_list <- split(smpl, smpl$gatinggroup)
  message("Splitting fcs files into ", length(smpl_list), " groups of different gating hierarchies.")
  for (i in smpl_list) {
    message(paste(i$FlowJoFileName, collapse = ", "), "\n")
  }
  #x<- smpl
  gs_list <- lapply(smpl_list,
                    get_gs,
                    remove_redundant_channels = remove_redundant_channels,
                    lapply_fun = lapply_fun)
  # split_size = split_size,
  # additional.sampleID = additional.sampleID,
  # ...)
  names(gs_list) <- names(smpl_list)

  return(
    gs_list = gs_list,
    FCS_files = smpl,
    gating_hierachies = gatings_list,
    gatings_compare = gatings_df
  )
}



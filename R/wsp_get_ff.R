#' Get (subsetted) flowFrames from one or many flowjo workspaces
#'
#' From a flowjo workspace with gated populations and the respective fcs files flowframes are generated.
#' Under the hood CytoML::flowjo_to_gatingset applies the geometric gate definitions and filters
#' relevant rows (indices) of fcs files. The compensation matrix as defined in flowjo will be used
#' to compensate fluorescence intensities. The index matrix (ind_mat) returned cannot be used with
#' the flow frames as these contain subsetted rows. Rather ind_mat can be used to read defined
#' rows from the fcs files without the need to re-read gate definitions from the workspace again.
#' ind_mats only, without any flow frames of gated populations, are returned by fcexpr::wsp_get_indices.
#' When reading data with ind_mats from fcs files, compensation needs to be applied manually. One
#' has to extract the compensation matrix either from the fcs file or from the flowjo workspace
#' where it is usually defined.
#'
#' If it is intended to pass flowframes to fcexpr::dr_to_fcs, it is recommended to have both,
#' transformed and untransformed, expression values returned.
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
#' value for all FCS files or a vector of same length as FCS files
#' @param remove_redundant_channels remove channels that are not part of the gating tree, mainly to reduce memory load
#' @param lapply_fun lapply function name, unquoted; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#' @param return_untransformed logical; do return untransformed (inverse) data
#' @param return_transformed logical; do return transformed data (transformation as in flowjo?!)
#' @param seed set a seed to reproduce downsampling
#' @param channels channels to use for leverage score calculation; use wsx_get_keywords to retrieve channel names/descriptions
#' @param leverage_score_for_sampling logical whether to use leverage scores for downsampling
#'
#' @return a list of (subsetted) flowframes with events that are within the gated population only
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#'\dontrun{
#' ff_list <- fcexpr::wsp_get_ff(wsp = "mypath/my.wsp", population = c("CD8+", "CD4+"), groups = "Compensation", invert_groups = TRUE)
#' # ff.list[[1]] may be passed to fcexpr::dr_to_fcs for instance
#' ## how to to work with the returned nested list
#' # ignore ind_mats and create one big data.frame w/o any filtering from flow frames
#' ff1 <- purrr::list_flatten(x = purrr::list_flatten(x = ff.list[["flowframes"]],
#'                                                    name_spec = "{outer}___{inner}"),
#'                            name_spec = "{outer}___{inner}")
#' ff1 <-
#'   purrr::map_dfr(ff1, purrr::compose(tibble::as_tibble, flowCore::exprs), .id = "type") |>
#'   tidyr::separate(type, into = c("sample", "pop", "trafo"), sep = "___")
#'
#' # pluck a subset of flowframe only
#' ff2 <- purrr::map(ff.list[["flowframes"]], purrr::pluck, "CD4", "transformed")
#' ff2 <- purrr::map_dfr(ff2, purrr::compose(tibble::as_tibble, flowCore::exprs), .id = "sample")
#'
#' # discard a transformation, then extract expression values
#' ff3 <- purrr::map(ff.list[["flowframes"]], ~ purrr::map(.x, ~ purrr::discard(.x, names(.x) == "untransformed")))
#' ff3 <- purrr::map(ff3, ~ purrr::map(.x, ~ purrr::pluck(.x, "transformed")))
#' ff3 <- purrr::map(ff3, ~ purrr::map(.x, ~ tibble::as_tibble(flowCore::exprs(.x))))
#'
#' ## how to work with ind_mats
#' # see ?fcexpr::inds_get_ff
#' # or ?fcexpr::fcs_files_df_get_ff
#'}
wsp_get_ff <- function(wsp,
                       FCS.file.folder = NULL,
                       groups = NULL,
                       population,
                       invert_groups = F,
                       samples = NULL,
                       invert_samples = F,
                       return_untransformed = T,
                       return_transformed = T,
                       downsample = 1,
                       remove_redundant_channels = F,
                       lapply_fun = lapply,
                       seed = 42,
                       channels = NULL,
                       leverage_score_for_sampling = F,
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
  if (missing(population)) {
    stop("population missing. please provide one or more.")
  }
  if ((!is.null(return_untransformed) && !return_untransformed) && (!is.null(return_transformed) && !return_transformed)) {
    stop("At least one of return_transformed or return_untransformed has to be TRUE.")
  }
  lapply_fun <- match.fun(lapply_fun)

  all_smpl <- wsx_get_sample_df_light(wsp)
  summary <-
    all_smpl[["summary"]] |>
    dplyr::select(-sample_info) |>
    dplyr::filter(n_samples > 1)
  if (nrow(summary) > 0) {
    message("Some samples have equal $FIL keywords and are in the same group.")
    message("Even when FCS files are in different local folders (on disk), a warning may be thrown that one of the FCS does not match.")
    message("This is save to ignore.")
    print(summary, n = Inf)
  }

  smpl <- wsx_get_sample_df(
    ws = wsp,
    groups = groups,
    samples = samples,
    FCS.file.folder = FCS.file.folder,
    invert_groups = invert_groups,
    invert_samples = invert_samples
  )

  if (is.null(smpl)) {
    return(NULL)
  }

  groups <- split(smpl$wsp, smpl$FlowJoGroup)
  groups <- groups[wsp] # order correctly

  if (is.numeric(downsample)) {
    ds <- downsample
  } else if (all(downsample == "min")) {
    ds <- 1
  } else {
    stop("downsample has to be numeric or 'min'. With min, all flowframes will be downsampled to that flowframe with the lowest number of events.")
  }
  if (length(ds) != 1 && length(ds) != nrow(smpl)) {
    stop("downsample has to have length 1 or length of FCS files.")
  }
  smpl$downsample <- ds


  # check for filenames and populations
  pp <-
    do.call(rbind, purrr::map2(wsp, groups, function(x,y) {
      wsx_get_poppaths(ws = x, groups = y, collapse = F)
    })) %>%
    dplyr::filter(FileName %in% smpl$FlowJoFileName) |>
    dplyr::filter(PopulationFullPath  %in% population | Population %in% population)
  if (nrow(pp) == 0) {
    stop("Population(s) not found for any sample.")
  }
  ## message which population was not found by samples
  pp_split <- split(pp, pp$Population)
  for (i in pp_split) {
    if (any(!smpl$FlowJoFileName %in% i$FileName)) {
      message("Population ", unique(i$Population), " not found for ", paste(smpl$FlowJoFileName[which(!smpl$FlowJoFileName %in% i$FileName)], collapse = ", "), ".")
    }
  }

  gs_list <- lapply(list(smpl),
                    get_gs,
                    remove_redundant_channels = remove_redundant_channels,
                    merge_to_gs = F)[[1]]

  ff.list <- lapply_fun(gs_list,
                        get_ff,
                        return_untransformed = return_untransformed,
                        return_transformed = return_transformed,
                        population = pp,
                        seed = seed,
                        channels = channels,
                        leverage_score_for_sampling = leverage_score_for_sampling
  )
  message("temporary .h5 were removed.")

  ff.list <- list(flowframes = purrr::map(ff.list, purrr::pluck, "ff"),
                  ind_mats = purrr::map(ff.list, purrr::pluck, "ind_mat"))

  if (all(downsample == "min") && length(ff.list[["flowframes"]]) > 1) {
    pops <- unique(unlist(purrr::map(ff.list[["flowframes"]], names)))
    names(pops) <- pops
    minrows <- purrr::map_int(pops, function(pop) {
      min(purrr::map_int(purrr::map(ff.list[["flowframes"]], ~ purrr::pluck(.x, pop, 1)), purrr::compose(nrow, flowCore::exprs)))
    })
    for (pop in pops) {
      for (z in seq_along(ff.list[["flowframes"]])) {
        if (pop %in% names(ff.list[["flowframes"]][[z]])) {
          set.seed(seed)
          inds <- sample(c(rep(T, minrows[pop]), rep(F, nrow(ff.list[["flowframes"]][[z]][[pop]][[1]])-minrows[pop])))
          for (y in seq_along(ff.list[["flowframes"]][[z]][[pop]])) {
            ff.list[["flowframes"]][[z]][[pop]][[y]] <- flowCore::Subset(ff.list[["flowframes"]][[z]][[pop]][[y]], inds)
          }
        }
      }
    }
  }

  return(c(ff.list, list(FCS_files = smpl)))

  ## how to to work with the returned nested list
  # ignore ind_mats and create one big data.frame w/o any filtering from flow frames
  # ff1 <- purrr::list_flatten(x = purrr::list_flatten(x = ff.list[["flowframes"]],
  #                                                    name_spec = "{outer}___{inner}"),
  #                            name_spec = "{outer}___{inner}")
  # ff1 <-
  #   purrr::map_dfr(ff1, purrr::compose(tibble::as_tibble, flowCore::exprs), .id = "type") |>
  #   tidyr::separate(type, into = c("sample", "pop", "trafo"), sep = "___")
  #
  # # pluck a subset of flowframe only
  # ff2 <- purrr::map(ff.list[["flowframes"]], purrr::pluck, "CD4", "transformed")
  # ff2 <- purrr::map_dfr(ff2, purrr::compose(tibble::as_tibble, flowCore::exprs), .id = "sample")
  #
  # # discard a transformation, then extract expression values
  # ff3 <- purrr::map(ff.list[["flowframes"]], ~ purrr::map(.x, ~ purrr::discard(.x, names(.x) == "untransformed")))
  # ff3 <- purrr::map(ff3, ~ purrr::map(.x, ~ purrr::pluck(.x, "transformed")))
  # ff3 <- purrr::map(ff3, ~ purrr::map(.x, ~ tibble::as_tibble(flowCore::exprs(.x))))

  ## how to work with ind_mat

  ## how to rename exprs columns by stained markers
  # make a function that takes flowframes
}

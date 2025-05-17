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
#' @param pData data frame to join via rownames of pData(gs), this will make meta
#' data columns available for plotting with fcexpr::plot_gates
#' @param pData_join_col column name in pData for joining; joinable columns
#' are found with fcexpr:::find_unique_level_columns(pData, pData_join_col)
#' @param get_gates run fcexpr::gs_get_gates and include results in the return list
#' @param force_gs_merge try to eliminate non-congruent populations from gating
#' hierarchies to allow merging of all samples into one gs
#' @param get_gates_args argument list for get_gates
#'
#' @return list of gatingsets, FCS files used, gating hierarchies
#' @export
#'
#' @examples
#'\dontrun{
#' # check fcs files in wsp first
#' samples <- fcexpr::wsx_get_fcs_paths(ws = "mypath/my.wsp", split = F, filter_AllSamples = T)
#' # read-in Gating set with selected samples
#' gs_list <- fcexpr::wsp_get_gs(wsp = "mypath/my.wsp", FCS.files.folder = "myLocalTopFolder/subfolder", groups = "Compensation", invert_groups = T, samples = samples$FileName[1:5])
#'}
wsp_get_gs <- function(wsp,
                       FCS.file.folder = NULL,
                       groups = NULL,
                       invert_groups = F,
                       samples = NULL,
                       invert_samples = F,
                       remove_redundant_channels = F,
                       pData = NULL,
                       pData_join_col = "identity",
                       get_gates = T,
                       get_gates_args = list(n_bins = 200^2),
                       force_gs_merge = F

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

  if (!is.null(pData)) {
    if (!pData_join_col %in% names(pData)) {
      stop("pData_join_col not in pData.")
    }
  }
  if (!is.logical(get_gates)) {
    stop("get_gates must be T or F (logical).")
  }

  for (i in wsp) {
    all_smpl <- wsx_get_sample_df_light(i)
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

  # check gating trees for equality
  gh_comparison <- compare_gating_hierarchies(wsp = wsp, sample_df = smpl)
  smpl <- dplyr::left_join(smpl, gh_comparison[["dfhashes"]], by = c("identity" = "identity", "wsp" = "wsp"))
  smpl_list <- make_smpl_list(smpl)



  gs_list <- lapply(smpl_list,
                    get_gs,
                    remove_redundant_channels = remove_redundant_channels)
  names(gs_list) <- names(smpl_list)

  if (length(gs_list) > 1 && force_gs_merge) {
    '
  tt <- purrr::map_dfr(gh_comparison[["gatings_list"]], purrr::pluck, "ps") |>
    dplyr::left_join(gh_comparison[["dfhashes"]])
  tt2 <- tt |>
    dplyr::distinct(gatinggroup, PopulationFullPath, Population)
  tt3 <- split(tt2$PopulationFullPath, tt2$gatinggroup)
  common <- Reduce(intersect, tt3)
  rm_pop <- unique(unlist(purrr::map(tt3, setdiff, common)))'
    message("Trying to merge ", length(gs_list), " gs by removing non-common populations.\n")

    gs_pop <- lapply(gs_list, function(x) flowWorkspace::gh_get_pop_paths(x[[1]]))
    common_pop <- Reduce(intersect, gs_pop)
    rm_pop <- purrr::map(gs_pop, setdiff, common_pop)

    message("common pops:\n", paste(common_pop, collapse = "\n"))
    message("\n")
    message("Non-common pops to be removed: \n", paste(names(rm_pop), rm_pop, collapse = "\n\n", sep = "\n"))
    message("\n")
    for (i in seq_along(gs_list)) {
      for (j in seq_along(rm_pop[[i]])) {
        flowWorkspace::gs_pop_remove(gs_list[[i]], rm_pop[[i]][[j]])
      }
    }

    gs_list <- list(flowWorkspace::merge_list_to_gs(gs_list))
    ## writ messages
    ## collapse other vars?
    # gatings_list remain what they are
    smpl$gatinggroup <- 1
    gh_comparison[["dfhashes"]]$gatinggroup <- 1
    smpl_list <- make_smpl_list(smpl)
  }

  ## join pData
  if (!is.null(pData)) {
    # which colukmns only have one level per pData_join_col level
    joincols <- find_unique_level_columns(df = pData, refcol = pData_join_col)
    message("pData joinable columns: ", paste(joincols[-1], collapse = ", "))
    gs_list <- lapply(gs_list, function(x) {
      pd <- flowCore::pData(x) |>
        tibble::rownames_to_column(pData_join_col) |>
        dplyr::left_join(pData |> dplyr::distinct(dplyr::pick(dplyr::all_of(joincols))), by = pData_join_col) |>
        tibble::column_to_rownames(pData_join_col)
      flowCore::pData(x) <- pd
      return(x)
    })
  }

  ## get gate data frames to avoid extra variables on global env (naming is hard)
  gate_dfs = NULL
  if (get_gates) {
    message("running gs_get_gates")
    gate_dfs <- lapply(gs_list, function(x) {
      tryCatch(
        expr = {
          Gmisc::fastDoCall(gs_get_gates, args = c(list(x), get_gates_args))
        },
        error = function(e) {
          print(x)
          NULL
        }
      )
    })
  }

  return(list(
    gs_list = gs_list,
    gate_dfs = gate_dfs,
    FCS_files = smpl,
    gating_hierachies = purrr::map(gh_comparison[["gatings_list"]], purrr::pluck, "gatings_list"),
    gatings_compare = gh_comparison[["gatings_df"]]
  ))
}



find_unique_level_columns <- function(df, refcol) {
  if (!refcol %in% names(df)) {
    stop("refcol not in df")
  }
  col_level_df <-
    df |>
    dplyr::group_by(!!rlang::sym(refcol)) |>
    dplyr::summarise(dplyr::across(.cols = dplyr::everything(), ~dplyr::n_distinct(.x)), .groups = "drop")
  uniquelevelcols <- c(refcol, names(which(apply(col_level_df == 1, 2, all))))
  return(uniquelevelcols)
}



make_smpl_list <- function(smpl) {
  smpl_list <- split(smpl, smpl$gatinggroup)
  names(smpl_list) <- paste(names(smpl_list),
                            purrr::map_chr(smpl_list, function(x) {
                              paste(paste(unique(basename(x$wsp)), collapse = "_"), paste(unique(x$FlowJoGroup), collapse = "_"), sep = "_")
                            }),
                            sep = "_")

  if (length(smpl_list) > 1) {
    message("Splitting fcs files into ", length(smpl_list), " groups of different gating hierarchies.\n\n")
    for (i in smpl_list) {
      message(paste(i$FlowJoFileName, collapse = ", "), "\n")
    }
  }

  return(smpl_list)
}



compare_gating_hierarchies <- function(wsp, sample_df = NULL) {
  gatings_list <- purrr::map(stats::setNames(wsp, wsp), function(ws) {
    ps <- wsx_get_popstats(ws = ws, return_stats = F)[["counts"]]
    if (!is.null(sample_df)) {
      ps <- dplyr::filter(ps, identity %in% sample_df$identity)
    }
    gatings_list <- purrr::map(stats::setNames(unique(ps$identity), unique(ps$identity)), function(x) {
      pssub <- dplyr::filter(ps, identity == x)
      # no, wsp is iterated over
      # if (length(unique(pssub$ws)) > 1) {
      #   # one identity can only exist once in a workspace; so if multiple ws: same fcs file with potentially different gating in multiple ws
      #   message("gating tree comparison: duplicate identity in sample df found. if they dont have different gatings, this will cause a problem.")
      # }
      #
      fcexpr::gating_tree_plot(PopulationFullPath = dplyr::pull(pssub, "PopulationFullPath"))
    })
    return(list(ps = ps, gatings_list = gatings_list))
  })

  gatings_df <-
    dplyr::bind_rows(purrr::map(gatings_list, function(x) {
      purrr::map_dfr(x[["gatings_list"]], `[[`, 3, .id = "identity") |>
        dplyr::left_join(x[["ps"]] |>
                           dplyr::select(PopulationFullPath, identity, xDim, yDim) |>
                           dplyr::rename("xDimfrom" = xDim, "yDimfrom" = yDim),
                         by = c("identity" = "identity", "from" = "PopulationFullPath")) |>
        dplyr::left_join(x[["ps"]] |>
                           dplyr::select(PopulationFullPath, identity, xDim, yDim) |>
                           dplyr::rename("xDimto" = xDim, "yDimto" = yDim),
                         by = c("identity" = "identity", "to" = "PopulationFullPath"))
    }), .id = "wsp")
  seperator <- "_"
  while (any(c(grepl(seperator, gatings_df$identity), grepl(seperator, gatings_df$wsp)))) {
    seperator <- paste0("_", seperator)
  }
  gatings_df <- split(gatings_df, paste(gatings_df$identity, gatings_df$wsp, sep = seperator)) # propagate identity and wsp for joining with smpl below
  gatings_df <- lapply(gatings_df, function(x) {
    x <-
      x |>
      dplyr::select(-c(identity, wsp)) |> # remove to only compare df for from, to and dims
      dplyr::arrange(from, to)
    #x <- dplyr::select(x, from, to) # ignore channels for splitting?!
    rownames(x) <- seq(1, nrow(x))
    return(x)
  })
  dfhashes <-
    utils::stack(purrr::map_chr(gatings_df, digest::digest, algo = "sha256")) |>
    dplyr::group_by(values) |>
    dplyr::mutate(gatinggroup = dplyr::cur_group_id()) |>
    dplyr::ungroup() |>
    tidyr::separate(ind, into = c("identity", "wsp"), sep = seperator) |>
    dplyr::rename("hash" = values) |>
    dplyr::select(-hash)
  return(list(gatings_list = gatings_list, gatings_df = gatings_df, dfhashes = dfhashes))
}

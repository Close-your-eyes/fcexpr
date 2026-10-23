#' Add dimension reductions and cluster annotations to merged flow data
#'
#' Combine flow frames, optionally scale the selected expression values and
#' correct sample effects with Harmony, calculate dimension reductions and
#' Louvain clusters, and append the results as channels. Optionally calculate
#' cluster markers and save FCS and RDS files.
#'
#' @param ff_list A nonempty list of `flowCore::flowFrame` objects,
#' with matching channel names and order. Their expression values must match
#' `ff_transformation`. For transformed data from [wsp_get_ff()], select
#' the transformed frames obtained with `return_transformed = TRUE` and
#' `return_untransformed = FALSE`. Inverse transformation requires the
#' transformation attribute specified by `trafoname` on every frame.
#' @param ff_transformation Transformation state of the input expression values:
#' `"transformed"` (default) or `"untransformed"`. This declares the input state;
#' it does not transform values before analysis. With untransformed input,
#' `trafo_to_fcs` is forced to `"untransformed"`.
#' @param channels Channels to use for analysis, passed to [ff_get_channels()].
#' The default, `NULL`, uses that helper's channel selection rules. The selected
#' channels must be present in every input frame.
#' @param grouplist A named list of sample annotations to append as channels;
#' see [ff_add_groups()]. Each element must contain one value per input frame,
#' in the same order as `ff_list`, with no missing values. Default: `NULL`.
#' @param scale_whole Scaling applied after sample matrices are concatenated:
#' `"none"` (default), `"z.score"`, or `"min.max"`. Z-score scaling uses [scale()]
#' on each channel; min-max scaling uses [scales::rescale()] on the entire
#' selected matrix, using a common range across channels.
#' @param scale_samples Scaling applied separately to each sample before
#' concatenation. Accepts the same options as `scale_whole`; default: `"none"`.
#' For `"min.max"`, each sample uses a common range across its selected channels.
#' @param downsample One finite nonnegative number controlling random sampling
#' within each input frame before analysis and merging. Default: `1`, which
#' keeps all events. Values below `1` retain that fraction of each sample,
#' rounded down, with a minimum target of one event. Values above `1` specify
#' an event count per sample; use whole numbers. Samples no larger than the
#' target are kept intact. Sampling uses `seed`, without replacement, and
#' preserves the original order of retained events.
#' @param run_pca Number of principal components to retain for downstream UMAP,
#' SOM, and clustering, and append to the output. Default: `0`, which skips PCA.
#' Positive values are capped at one fewer than the number of selected channels.
#' The returned PCA object contains the full result from [ff_calc_pca()].
#' @param run_harmony Whether to correct sample effects with
#' `harmony::RunHarmony`. Default: `FALSE`. Runs only with more than one input
#' frame, treating each frame as a batch. Uses the selected expression matrix
#' after scaling and optional PCA; requires the harmony package.
#' @param run_umap Whether to calculate UMAP coordinates with
#' [ff_calc_umap_tsne()] using `uwot::umap`. Default: `TRUE`.
#' @param run_som Whether to calculate SOM coordinates with [ff_calc_som()],
#' using `EmbedSOM::SOM` followed by `EmbedSOM::EmbedSOM`. Default: `FALSE`.
#' @param run_find_clusters Whether to calculate Louvain clusters with
#' [get_louvain_cluster()] and append numeric cluster labels. Default: `TRUE`.
#' @param calc_cluster_marker Whether to calculate global cluster markers with
#' [ff_calc_marker()]. Requires cluster assignments from `run_find_clusters`.
#' Marker calculation uses expression values from the retained events before
#' scaling, PCA, and Harmony, in the transformation state supplied by the caller.
#' Default: `FALSE`.
#' @param UMAP_args A list passed as `fun_args` to [ff_calc_umap_tsne()].
#' Default: `list(metric = "cosine", verbose = TRUE, scale = TRUE)`.
#' An explicit list replaces these defaults; it is not merged with them.
#' Other UMAP settings use the defaults of `uwot::umap`.
#' @param SOM_args A list of arguments passed to `EmbedSOM::SOM` through
#' [ff_calc_som()]. Default: `list()`, which uses the backend defaults.
#' @param EmbedSOM_args A list of arguments passed to `EmbedSOM::EmbedSOM`
#' through [ff_calc_som()]. Default: `list()`, which uses the backend defaults.
#' @param FindNeighbors_args A list of arguments passed to
#' `Seurat::FindNeighbors` through [get_louvain_cluster()]. Default: `list()`.
#' @param FindClusters_args A list of arguments passed to `Seurat::FindClusters`
#' through [get_louvain_cluster()]. Default: `list(resolution = c(0.1))`.
#' Supply multiple resolutions to calculate multiple clusterings.
#' @param mc_cores Number of cores used by [get_louvain_cluster()] for multiple
#' clustering resolutions. Default: `1`; capped at one fewer than the detected
#' number of cores.
#' This argument is not forwarded to marker calculation.
#' @param save_to_disk Output formats: `"fcs"`, `"rds"`, or both (the default).
#' Set to `NULL` to skip file writing. RDS output contains the same list that is
#' returned by this function.
#' @param save_path Output directory. By default, a timestamped directory ending
#' in `_FCS_dr` is created in the working directory. Set to `NULL` or `NA` to
#' disable file writing. A non-missing path may be created even when
#' `save_to_disk = NULL`.
#' @param save_name Output file stem. Default: `NULL`, which uses timestamped
#' names. A supplied trailing `.fcs` extension is removed, ignoring case;
#' `.fcs` or `.rds` is then appended for the respective output format.
#' @param trafo_to_fcs Expression channel states to include in the returned
#' flow frame and any saved FCS file: `"transformed"`, `"untransformed"`, or
#' both (the default). For transformed input, requesting both uses [ff_merge()]
#' to retain transformed fluorescence channels with the `_trans` suffix and
#' inverse-transformed channels with their original names; non-fluorescence
#' channels follow that helper's exclusion rules. Requesting either state alone
#' keeps the original channel names. Untransformed input is always returned
#' as untransformed channels.
#' @param trafoname Name of the inverse-transformation attribute on each input
#' frame. Default: `"trafolistinv"`. Used to inverse-transform transformed input
#' when untransformed output channels are requested.
#' @param timeChannel Candidate time channel names passed to [ff_get_channels()].
#' Default: `c("Time", "HDR-T")`. That helper controls time-channel exclusion and
#' automatic detection, including when this argument is `NULL`.
#' @param seed Random seed used for downsampling, Harmony, and the UMAP and
#' SOM wrappers. Default: `42`.
#' Clustering seeds are controlled separately through `FindClusters_args`.
#' @param ... Additional arguments passed to [ff_get_channels()], such as
#' `part_match` and `rm_scatter`.
#'
#' @details Downsampling is applied first, followed by sample scaling,
#' whole-data scaling, optional PCA, and optional Harmony correction. UMAP, SOM,
#' and clustering use the resulting analysis matrix. Scaling, PCA, and Harmony
#' do not replace expression channels in the output. Appended PCA coordinates
#' and the returned PCA object contain the result before Harmony correction;
#' Harmony coordinates are not appended separately or returned.
#' With all analyses disabled, the retained input events are still merged,
#' with an `ident` channel identifying their input frame. `trafo_to_fcs` controls
#' the expression channel states in this merged output.
#'
#' Errors from the SOM and marker helper calls are reported as messages, and
#' their results remain `NULL`. Other errors propagate. Incompatible channel
#' names cause an early return of the diagnostics from [ff_compare_channels()].
#'
#' @return Normally, a list with `flowframe` (the merged flow frame with analysis
#' channels), `marker` (the result of [ff_calc_marker()], or `NULL`), and `pca`
#' (the result of [ff_calc_pca()], or `NULL` when PCA is disabled).
#' If channel names are incompatible, the diagnostics from
#' [ff_compare_channels()] are returned instead.
#' @export
#'
#'
#' @examples
#' \dontrun{
#' # Simulated input has no inverse-transformation attribute.
#' frames <- ff_simulate2(npop = 5, dims = 10, m = 3)
#' channels <- setdiff(colnames(flowCore::exprs(frames[[1]])), "Time")
#' result <- dimred_to_fcs(
#'   ff_list = frames,
#'   ff_transformation = "untransformed",
#'   channels = channels,
#'   downsample = 0.1,
#'   run_pca = 5,
#'   run_umap = FALSE,
#'   run_find_clusters = FALSE,
#'   save_path = NULL
#' )
#' result$flowframe
#' result$pca
#'
#' # With transformed frames from wsp_get_ff(), keep transformed channels only.
#' # This avoids requiring inverse-transformation attributes on the frames.
#' transformed <- wsp_get_ff(
#'   wsp = "workspace.wsp",
#'   FCS.file.folder = "FCS_files",
#'   groups = "All Samples",
#'   population = "live",
#'   return_transformed = TRUE,
#'   return_untransformed = FALSE
#' )
#' frames <- lapply(transformed, function(x) x[["transformed"]][[1]])
#' out <- dimred_to_fcs(
#'   ff_list = frames,
#'   trafo_to_fcs = "transformed",
#'   run_pca = 5,
#'   run_harmony = TRUE,
#'   save_path = NULL
#' )
#' marker <- ff_calc_marker(
#'   out[["flowframe"]],
#'   cluster_cols = "res.0.1",
#'   channels = paste0(channels, "_trans")
#' )
#' df <- marker[["res.0.1"]][["global"]]
#' df <-
#'   df |>
#'   dplyr::mutate(pvalue2 = -log10(pvalue)) |>
#'   dplyr::mutate(pvalue2 = ifelse(is.infinite(pvalue2), 304, pvalue2))
#' heatmap_long_df(df = df,
#'                 groups = "cluster",
#'                 features = "channel",
#'                 values = "mean_cluster_scale",
#'                 dotsizes = "pvalue2")
#' }
dimred_to_fcs <- function(ff_list,
                          ff_transformation = c("transformed", "untransformed"),
                          channels = NULL,
                          grouplist = NULL,
                          scale_whole = c("none", "z.score", "min.max"),
                          scale_samples = c("none", "z.score", "min.max"),
                          downsample = 1,
                          run_pca = 0,
                          run_harmony = F,
                          run_umap = T,
                          run_som = F,
                          run_find_clusters = T,
                          calc_cluster_marker = F,
                          UMAP_args = list(metric = "cosine", verbose = T, scale = T),
                          SOM_args = list(),
                          EmbedSOM_args = list(),
                          FindNeighbors_args = list(),
                          FindClusters_args = list(resolution = c(0.1)),

                          mc_cores = 1,
                          save_to_disk = c("fcs", "rds"),
                          save_path = file.path(getwd(), paste0(substr(gsub("\\.", "", make.names(as.character(Sys.time()))), 2, 15), "_FCS_dr")),
                          save_name = NULL,

                          trafo_to_fcs = c("transformed", "untransformed"),
                          trafoname = "trafolistinv",
                          timeChannel = c("Time", "HDR-T"),
                          seed = 42,
                          ...) {
  fcexpr:::.ensure_packages(c("flowCore", "parallel"))


  ## ---- checks --------

  if (!is.null(save_to_disk)) {
    save_to_disk <- match.arg(save_to_disk, c("fcs", "rds"), several.ok = T)
  }

  mc_cores <- min(mc_cores, parallel::detectCores() - 1)
  trafo_to_fcs <- rlang::arg_match(trafo_to_fcs, multiple = TRUE)
  ff_transformation <- rlang::arg_match(ff_transformation)

  if (ff_transformation == "untransformed" && (length(trafo_to_fcs) > 1 || trafo_to_fcs != "untransformed")) {
    message("ff_transformation is untransformed, setting trafo_to_fcs to untransformed.")
    trafo_to_fcs <- "untransformed"
  }

  # check grouplist
  if (!is.null(grouplist)) {
    if (!is.list(grouplist)) {
      stop("grouplist has to be a list.")
    }
    if (is.null(names(grouplist))) {
      stop("grouplist has to have names. These names will become channel names in the FCS file.")
    }
    # if (!all(unlist(lapply(grouplist, function(x) is.numeric(x))))) {
    #   stop("Please provide numeric values only in grouplist: E.g. use as.numeric(as.factor(x)), where x is numeric.")
    # }
    if (any(unlist(lapply(grouplist, function(x) is.na(x))))) {
      stop("NA found in sample infos.")
    }
    if (!all(unlist(lapply(grouplist, function(x) length(x) == length(ff_list))))) {
      stop("Length of each additional sample information has to match the length of selected samples, which is: ", length(ff_list),".")
    }
  }

  ## check trafoname

  # set scaling funs
  scale_samples <- switch(match.arg(scale_samples, c("none", "z.score", "min.max")),
                          z.score = scale,
                          min.max = scales::rescale,
                          none = function(x) return(x))

  scale_whole <- switch(match.arg(scale_whole, c("none", "z.score", "min.max")),
                        z.score = scale,
                        min.max = scales::rescale,
                        none = function(x) return(x))


  ## channel names from first ff
  channels <- ff_get_channels(ff = ff_list[[1]], timeChannel = timeChannel, channels = channels, ...)
  message("Channels: ", paste(channels, collapse = ", "), "\n")

  # check if channel names and desc are equal
  if (!is.null(fcs_check <- ff_compare_channels(ff_list = ff_list, channels = channels, strict = T))) {
    return(fcs_check)
  }

  if (downsample != 1) {
    if (!is.numeric(downsample) ||
        length(downsample) != 1L ||
        !is.finite(downsample)) {
      stop("downsample must be NULL or one finite number.")
    }

    set.seed(seed)
    ff_list <- purrr::map(ff_list, function(ff) {
      n <- nrow(ff)

      if (downsample<1) {
        downsample <- floor(downsample*n)
      }

      if (downsample == 0) {
        downsample <- 1
      }

      if (n <= downsample) {
        return(ff)
      }
      keep <- sort(sample.int(n, size = downsample, replace = FALSE))
      ff[keep, ]
    })
  }

  ## apply scaling which was selected above and select channels to use for dimension reduction.
  expr_select <- scale_whole(do.call(rbind, purrr::map(ff_list, ~scale_samples(flowCore::exprs(.x)[, channels, drop = FALSE])))) # names to desc here?!

  message("using ", nrow(expr_select), " cells from ", length(ff_list), " FCS files.\n")

  pca_result <- NULL
  pca_dims <- NULL
  if (run_pca > 0) {
    run_pca <- min(ncol(expr_select) - 1, run_pca)
    message("Calculating PCA.\nStart: ", Sys.time())
    pca_result <- ff_calc_pca(exprs = expr_select)
    pca_dims <- pca_result[["x"]][, seq_len(run_pca), drop = FALSE]
    expr_select <- pca_dims # overwrite original data with PCA embedding
    message("Done. ", Sys.time())
  }


  ## ---- harmonize ------
  if (length(ff_list) > 1L && run_harmony) {
    fcexpr:::.ensure_packages("harmony")

    harmony_meta <- data.frame(
      file = factor(rep(
        seq_along(ff_list),
        times = vapply(ff_list, nrow, integer(1))
      ))
    )

    set.seed(seed)
    expr_select <- harmony::RunHarmony(
      data_mat = expr_select,
      meta_data = harmony_meta,
      vars_use = "file",
      return_object = FALSE,
      verbose = TRUE
    )
  }

  ## ---- dim reds -------
  # umap
  umap_dims <- NULL
  if (run_umap) {
    message("Calculating UMAP.\nStart: ", Sys.time())

    umap_dims <- ff_calc_umap_tsne(exprs = expr_select,
                                   seed = seed,
                                   fun_args = UMAP_args)

    message("End: ", Sys.time())
  }

  ## SOMs do not like the input to be scaled
  ## whereas UMAP and tSNE do like it

  # SOM
  som_dims <- NULL
  if (run_som) {
    tryCatch({
      message("Calculating SOM.\nStart: ", Sys.time())

      som_dims <- ff_calc_som(exprs = expr_select,
                              seed = seed,
                              som_args = SOM_args,
                              embedsom_args = EmbedSOM_args)

      message("End: ", Sys.time())
    }, error = function(err) {
      message("SOM with error: ", err)
      som_dims <- NULL
    })
  }

  clust_idents <- NULL
  if (run_find_clusters) {
    clust_idents <- get_louvain_cluster(exprs = expr_select,
                                        FindNeighbors_args = FindNeighbors_args,
                                        FindClusters_args = FindClusters_args,
                                        mc.cores = mc_cores,
                                        return_as = "numeric")
  }

  ## ---- prepare final fcs file -------
  if (ff_transformation == "transformed") {
    if (all(c("transformed", "untransformed") %in% trafo_to_fcs)) {
      ff <- ff_merge(ff_list = ff_list,
                     grouplist = grouplist,
                     add_transformed_channels = T,
                     trafoname = trafoname)
    } else if ("transformed" %in% trafo_to_fcs) {
      ff <- ff_merge(ff_list = ff_list,
                     grouplist = grouplist,
                     add_transformed_channels = F)
    } else if ("untransformed" %in% trafo_to_fcs) {
      ff <- ff_merge(ff_list = ff_transform(ff = ff_list, trafolist = trafoname),
                     grouplist = grouplist,
                     add_transformed_channels = F)
    }
  } else if (ff_transformation == "untransformed") {
    # transformed are not relevant then
    # untransformed is what people are used to
    ff <- ff_merge(ff_list = ff_list,
                   grouplist = grouplist,
                   add_transformed_channels = F)
  }


  dimred_dims <- do.call(cbind, list(pca_dims, umap_dims, som_dims, clust_idents))
  if (!is.null(dimred_dims) && ncol(dimred_dims) > 0) {
    ff <- ff_add_columns(ff, mat = dimred_dims)
  }

  # no empty desc. may overload fcs, but easier for desc transfer to colnames of exprs
  desc <- ff@parameters@data$desc
  desc_na <- which(is.na(desc))
  desc[desc_na] <- ff@parameters@data$name[desc_na]
  ff@parameters@data$desc <- desc

  # prepare channel desc
  # name_desc <- stats::setNames(ff_list[[1]]@parameters@data[["desc"]], ff_list[[1]]@parameters@data[["name"]])
  # name_desc <- name_desc[which(!is.na(name_desc))]
  # channel_desc <- rep("", ncol(dim_red_data))
  #
  # for (i in seq_along(name_desc)) {
  #   channel_desc[grep(names(name_desc)[i], colnames(dim_red_data))] <- name_desc[i]
  # }
  #
  # channel_desc_augment <- channel_desc
  # #channel_desc_augment[intersect(which(channel_desc_augment != ""), which(grepl("scaled$", colnames(dim_red_data))))] <- paste0(channel_desc_augment[intersect(which(channel_desc_augment != ""), which(grepl("scaled$", colnames(dim_red_data))))], "_scaled")
  # channel_desc_augment[intersect(which(channel_desc_augment != ""), which(grepl(paste0(trafosuffix, "$"), colnames(dim_red_data))))] <- paste0(channel_desc_augment[intersect(which(channel_desc_augment != ""), which(grepl(paste0(trafosuffix, "$"), colnames(dim_red_data))))], paste0("_", trafosuffix))
  #
  # channel_desc_augment[which(channel_desc_augment == "")] <- colnames(dim_red_data)[which(channel_desc_augment == "")]
  # channel_desc_augment <- make.names(channel_desc_augment)
  # names(channel_desc_augment) <- colnames(dim_red_data)

  # https://github.com/RGLab/flowCore/issues/201
  #flowCore::keyword(ff) <- flowCore:::updateTransformKeywords(ff)

  ## ---- calc cluster markers -------
  marker <- NULL
  if (calc_cluster_marker) {
    message("Calculating markers.")
    if (all(c("transformed", "untransformed")) %in% trafo_to_fcs) {
      channels2 <- paste0(channels, formals(ff_merge)$suffix)
      exprs <- NULL
    } else if ("transformed" %in% trafo_to_fcs) {
      channels2 <- channels
      exprs <- NULL
    } else if ("untransformed" %in% trafo_to_fcs) {
      channels2 <- channels
      exprs <- do.call(rbind, lapply(ff_list, function(x) {
        flowCore::exprs(x)[, channels, drop = FALSE]
      }))
      exprs <- cbind(exprs, clust_idents)
    }
    tryCatch({
      marker <- ff_calc_marker(ff = ff,
                               exprs = exprs,
                               channels = channels2,
                               cluster_cols = colnames(clust_idents))
    }, error = function(err) {
      message("cluster marker calculation caused an error: ", err)
    })
  }


  ## ---- write to disk -------
  if (!is.null(save_path) && !is.na(save_path)) {
    message("Writing files to disk: ", save_path)
    t <- format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%Y%m%d_%H%M%S")
    dir.create(save_path, showWarnings = F)
    if (!is.null(save_name)) {
      save_name <- gsub("\\.fcs$", "", save_name, ignore.case = T)
      save_name <- save_name[1]
    }
    if (is.null(save_name)) {
      sv_pth <- file.path(save_path, paste0(t, "_dr_ff_list.rds"))
    } else {
      sv_pth <- file.path(save_path, paste0(save_name, ".rds"))
    }
    if ("rds" %in% save_to_disk) {
      brathering::saverds2(list(flowframe = ff, marker = marker, pca = pca_result), file = sv_pth, compress = F)
    }
    if (is.null(save_name)) {
      sv_pth <- file.path(save_path, paste0(t, "_dr.fcs"))
    } else {
      sv_pth <- file.path(save_path, paste0(save_name, ".fcs"))
    }
    if ("fcs" %in% save_to_disk) {
      flowCore::write.FCS(ff, sv_pth)
    }
  }
  return(list(flowframe = ff, marker = marker, pca = pca_result))
}



# tsne
'  if (run_tsne) {
    message("Calculating tSNE.\nStart: ", Sys.time())

    if ("perplexity" %in% names(tSNE_args) && length(tSNE_args[["perplexity"]]) > 1) {
      tsne_dims <- do.call(cbind, parallel::mclapply(tSNE_args[["perplexity"]], function(z) {
        set.seed(seed)
        out <- Gmisc::fastDoCall(Rtsne::Rtsne, args = c(list(X = expr_select, verbose = F, perplexity = z), tSNE_args[which(!names(tSNE_args) %in% c("X", "verbose", "perplexity"))]))$Y
        colnames(out) <- c(paste0("tSNE_1_", z), paste0("tSNE_2_", z))
        return(out)
      }, mc.cores = mc_cores))
    } else {
      set.seed(seed)
      tsne_dims <- Gmisc::fastDoCall(Rtsne::Rtsne, args = c(list(X = expr_select), tSNE_args[which(names(tSNE_args) != "X")]))$Y
      colnames(tsne_dims) <- c("tSNE_1", "tSNE_2")
    }

    message("End: ", Sys.time())
  }'


## ---- clusterings -------
# find communities (clusters)
# louvain is run on the original data, not on any metaclusters
' if (run_louvain) {
    tryCatch({
      message("Calculating snn for louvain.\nStart: ", Sys.time())
      rownames(expr_select) <- 1:nrow(expr_select)
      snn <- Gmisc::fastDoCall(Seurat::FindNeighbors, args = c(list(object = Gmisc::fastDoCall(louvain_scale, expr_select)), FindNeighbors_args[which(names(FindNeighbors_args) != "object")]))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("error in snn calculation: ", err)
      run_louvain <- F
    })
  }

  if (run_louvain) {
    tryCatch({
      message("Finding clusters with Seurats implementation of the Louvain algorithm and parallel::mclapply using ", mc_cores," cores.\nStart: ",Sys.time())
      if ("resolution" %in% names(FindClusters_args) && length(FindClusters_args[["resolution"]]) > 1) {
        clust_idents <- do.call(cbind, parallel::mclapply(FindClusters_args[["resolution"]], function(x) {
          apply(Gmisc::fastDoCall(Seurat::FindClusters, args = c(list(object = snn$snn, resolution = x, verbose = F), FindClusters_args[which(!names(FindClusters_args) %in% c("object", "resolution", "verbose"))])), 2, as.numeric)
        }, mc.cores = mc_cores))
      } else {
        clust_idents <- apply(Gmisc::fastDoCall(Seurat::FindClusters, args = c(list(object = snn$snn), FindClusters_args[which(!names(FindClusters_args) %in% c("object"))])), 2, as.numeric)
      }


      colnames(clust_idents) <- paste0("louvain_", FindClusters_args[["resolution"]])
      dim_red_data <- do.call(cbind, list(dim_red_data, clust_idents))
      print(apply(clust_idents, 2, function(x) length(unique(x))))
      message("End: ", Sys.time())

    }, error = function(err) {
      message("run_louvain with error: ", err)
    })
  }'


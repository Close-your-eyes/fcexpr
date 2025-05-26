#' Calculate dimension reduction and cluster annotation with data from one or more flow frames and add these parameters to a (concatenated) fcs file
#'
#' It is expected that ff.list was made with wsp_get_ff with return_transformed
#' TRUE and return_untransformed FALSE.
#' After optional scaling, this runs ff_calc_pca, ff_calc_umap, ff_calc_som,
#' get_louvain_cluster. flow frames from ff.list are joined with ff_merge,
#' data from pca, umap etc added with ff_add_columns.
#' If write.untransformed.channels.to.FCS is TRUE, ff_transform is called.
#'
#' @param ff.list a list of flowFrames as received from fcexpr::wsp_get_ff
#' (compensated with Compensation Matrix as defined in FlowJo by default) or
#' as received from fcexpr::inds_get_ff (directly from FCS files, not
#' compensated by default)
#' @param channels channel selection, passed to ff_get_channels
#' @param grouplist named list of additional channels to identify samples or
#' group them in flowjo; see ff_add_groups
#' @param scale.whole if and how to scale channels after concatenation of
#' flowframes in ff.list
#' @param scale.samples if and how to scale channels of flowframes in ff.list
#' individually before concatenation
#' @param run.pca select number of pcs; if 0 no pca is computed
#' @param run.umap logical, whether to calculate UMAP dimension reduction
#' with uwot::umap
#' @param run.som logical, whether to calculate SOM dimension reduction
#' EmbedSOM::SOM followed by EmbedSOM::EmbedSOM
#' @param run.find.clusters logical, whether to detect clusters with
#' get_louvain_clusters
#' @param calc_cluster_marker calculate global cluster markers;
#' see ff_calc_markers
#' @param mc.cores number of cores to use for parallel computing of
#' cluster markers and multiple cluster resolutions;
#' mc.cores is passed to parallel::mclapply or parallel::mcmapply;
#' limited to parallel::detectCores()-1
#' @param save.to.disk what to save to disk: (concatenated) and appended FCS
#' file and/or rds file with several elements in a list
#' @param save.path where to save elements specified in save.to.disk;
#' set to NULL to have nothing written to disk
#' @param write.transformed.channels.to.FCS do save transformed channels
#' (e.g. logicle or arcsinh) to FCS file
#' @param write.untransformed.channels.to.FCS do save untransformed channels
#' (inverse) to FCS file; requires trafoname in attributes of ff
#' @param timeChannel name of the Time channel to exclude from all
#' analyses and calculation; if NULL will be attempted to be detected
#' automatically
#' @param trafoname name of transformatio in attributes of ff
#' @param seed set a seed for reproduction of dimension reductions
#' @param UMAP_args arguments to uwot::umap
#' @param SOM_args args to EmbedSOM::SOM
#' @param EmbedSOM_args args to EmbedSOM::EmbedSOM
#' @param save.name name to use for save rds and or fcs file
#' @param FindNeighbors_args arugments to Seurat function fof louvain clustering
#' @param FindClusters_args arugments to Seurat function fof louvain clustering
#' @param ... arugments to ff_get_channels
#'
#' @return
#' A list with 3 elements: (i) The matrix of fluorescence intensities and appended information (dim red, clustering). This is the table which is written into a newly generated fcs file.
#' (ii) A character vector of meaningful column names which may be used for the table in R (rather for convenience). (iii) Tables of marker features (each cluster vs all other events and all clusters pairwise).
#' @export
#'
#'
#' @examples
#'\dontrun{
#'############################
#'### Plot cluster markers ###
#'############################
#' dr <- fcexpr::dr_to_fcs(ff.list = ffs,
#' channels = channels,
#' louvain__resolution = 0.5,
#' run.lda = "louvain_0.5",
#' clustering.for.marker.calc = c("louvain_0.5"),
#' save.path = NULL)
#' marker <- dr[[3]][[1]][[1]]
#'marker$channel_desc2 <- sapply(strsplit(marker$channel_desc, "_"), "[", 1)
#'marker <-
#'  marker %>%
#'  dplyr::mutate(pvalue = ifelse(pvalue == 0, 1e-300, marker$pvalue)) %>%
#'  dplyr::group_by(channel_desc2) %>%
#'  dplyr::mutate(mean_scaled = fcexpr:::min.max.normalization(mean))
#'
#'ggplot(marker, aes(x = as.factor(cluster), y = channel_desc2, fill = -log10(pvalue))) +
#'  geom_tile(color = "black") +
#'  theme_bw() +
#'  geom_text(aes(label = diff_sign)) +
#'  scale_fill_viridis_c()
#'
#'ggplot(marker, aes(x = as.factor(cluster), y = channel_desc2, fill = mean_diff)) +
#'  geom_tile(color = "black") +
#'  theme_bw() +
#'  geom_text(aes(label = diff_sign)) +
#'  scale_fill_viridis_c()
#'
#'ggplot(marker, aes(x = as.factor(cluster), y = channel_desc2, fill = mean_scaled)) +
#'  geom_tile(color = "black") +
#'  theme_bw() +
#'  scale_fill_viridis_c()
#'
#'ggplot(marker, aes(x = mean_diff, y = -log10(pvalue), label = channel_desc2)) + #color = channel_desc2,
#'  geom_point() +
#'  theme_bw() +
#'  labs(title = "cluster markers (vs all other cells each)") +
#'  ggrepel::geom_text_repel(max.overlaps = 20, show.legend = F) +
#'  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), strip.background = element_rect(fill = "hotpink2")) +
#'  geom_vline(xintercept = 0, col = "tomato2", linetype = "dashed") +
#'  geom_hline(yintercept = 100, col = "tomato2", linetype = "dashed") +
#'  facet_wrap(vars(cluster))
#'
#'  ##############################################
#'  ### find clusters which may be subject #######
#'  ### to bi- or multimodality in any channel ###
#'  ##############################################
#'  # make one data frame
#'  marker_all <- purrr::map_dfr(setNames(names(out[["marker"]]), names(out[["marker"]])),
#'  function(x) out[["marker"]][[x]][["marker_table"]], .id = "clustering")
#'
#'  # sort by diptest p value; or low p indicates bi- or multimodality
#'  marker_all <- dplyr::arrange(marker_all, diptest_p)
#'  # see ?diptest::dip.test
#'  # a low p-value indicates bi- or multimodality (multiple peaks)
#'  # a high p-value (close to 1) indicates unimodality
#'  # multimodality indicates heterogeneity within in the cluster
#'  # and may justify to separate that cluster further into sub-clusters
#'  # this depends on the interpretation of the scientist though
#'
#'
#'  ##############################################
#'  ### overlay gated populations from flowjo ####
#'  ### on dimension reduction plot (SOM/UMAP) ###
#'  ##############################################
#'
#' common_cols <- Reduce(intersect, purrr::map(ff[["indices"]], colnames))
#' ff[["indices"]] <- purrr::map(ff[["indices"]], function(x) x[,which(colnames(x) %in% common_cols)])
#' ind_mat <- do.call(rbind, ff[["indices"]])
#'
#' som <- ggplot(out[["df"]], aes(SOM_1, SOM_2, color = as.factor(cutree_30))) +
#' geom_point(size = 0.5) +
#' geom_density2d(data = . %>% dplyr::filter(ind_mat[,which(basename(colnames(ind_mat)) == "NK cells")]), color = "black",
#' contour_var = "ndensity", breaks = c(0.1, 0.2, 0.4, 0.6, 0.8)) +
#' guides(color = guide_legend(override.aes = list(size = 2)))
#'
#'}
dr_to_fcs3 <- function(ff.list,
                       channels = NULL,
                       grouplist = NULL,
                       scale.whole = c("none", "z.score", "min.max"),
                       scale.samples = c("none", "z.score", "min.max"),
                       run.pca = 0,
                       run.umap = T,
                       run.som = T,
                       run.find.clusters = T,
                       calc_cluster_marker = F,
                       UMAP_args = list(metric = "cosine", verbose = T, scale = T),
                       SOM_args = list(),
                       EmbedSOM_args = list(),
                       FindNeighbors_args = list(),
                       FindClusters_args = list(resolution = c(0.1)),

                       mc.cores = 1,
                       save.to.disk = c("fcs", "rds"),
                       save.path = file.path(getwd(), paste0(substr(gsub("\\.", "", make.names(as.character(Sys.time()))), 2, 15), "_FCS_dr")),
                       save.name = NULL,

                       write.transformed.channels.to.FCS = T,
                       write.untransformed.channels.to.FCS = T,
                       trafoname = "trafolistinv",
                       timeChannel = c("Time", "HDR-T"),
                       seed = 42,
                       ...) {


  ## ---- checks --------
  if (!requireNamespace("diptest", quietly = T)) {
    utils::install.packages("diptest")
  }
  if (!requireNamespace("matrixStats", quietly = T)) {
    utils::install.packages("matrixStats")
  }
  if (run.umap && !requireNamespace("uwot", quietly = T)) {
    utils::install.packages("uwot")
  }
  if (!requireNamespace("devtools", quietly = T)) {
    utils::install.packages("devtools")
  }
  if (calc_cluster_marker && !requireNamespace("presto", quietly = T)) {
    devtools::install_github('immunogenomics/presto')
  }
  if (calc_cluster_marker && !requireNamespace("matrixTests", quietly = T)) {
    utils::install.packages("matrixTests")
  }
  if (run.som && !requireNamespace("EmbedSOM", quietly = T)) {
    devtools::install_github("exaexa/EmbedSOM")
  }

  if (!is.null(save.to.disk)) {
    save.to.disk <- match.arg(save.to.disk, c("fcs", "rds"), several.ok = T)
  }

  mc.cores <- min(mc.cores, parallel::detectCores() - 1)

  # check grouplist
  if (!is.null(grouplist)) {
    if (!is.list(grouplist)) {
      stop("grouplist has to be a list.")
    }
    if (is.null(names(grouplist))) {
      stop("grouplist has to have names. These names will become channel names in the FCS file.")
    }
    if (!all(unlist(lapply(grouplist, function(x) is.numeric(x))))) {
      stop("Please provide numeric values only in grouplist: E.g. use as.numeric(as.factor(x)), where x is numeric.")
    }
    if (any(unlist(lapply(grouplist, function(x) is.na(x))))) {
      stop("NA found in sample infos.")
    }
    if (!all(unlist(lapply(grouplist, function(x) length(x) == length(ff.list[[1]]))))) {
      stop("Length of each additional sample information has to match the length of selected samples, which is: ", length(ff.list[[1]]),".")
    }
  }

  ## check trafoname

  # set scaling funs
  scale.samples <- switch(match.arg(scale.samples, c("none", "z.score", "min.max")),
                          z.score = scale,
                          min.max = min.max.normalization,
                          none = function(x) return(x))

  scale.whole <- switch(match.arg(scale.whole, c("none", "z.score", "min.max")),
                        z.score = scale,
                        min.max = min.max.normalization,
                        none = function(x) return(x))


  ## channel names from first ff
  channels <- ff_get_channels(ff = ff.list[[1]], timeChannel = timeChannel, channels = channels, ...)
  message("Channels: ", paste(channels, collapse = ", "))

  # check if channel names and desc are equal
  if (!is.null(fcs_check <- ff_compare_channels(ff_list = ff.list, channels = channels, strict = T))) {
    return(fcs_check)
  }

  ## apply scaling which was selected above and select channels to use for dimension reduction.
  expr.select <- scale.whole(do.call(rbind, lapply(ff.list, function(x) scale.samples(flowCore::exprs(x)[, channels])))) # names to desc here?!

  pca.result <- NULL
  if (run.pca > 0) {
    run.pca <- min(ncol(expr.select) - 1, run.pca)
    message("Calculating PCA.\nStart: ", Sys.time())
    pca.result <- ff_calc_pca(exprs = expr.select)
    expr.select <- pca.result[["x"]][, 1:run.pca] # overwrite original data with PCA embedding
    message("Done. ", Sys.time())
  }

  ## ---- dim reds -------
  # umap
  umap.dims <- NULL
  if (run.umap) {
    message("Calculating UMAP.\nStart: ", Sys.time())

    umap.dims <- ff_calc_umap(exprs = expr.select,
                              seed = seed,
                              args = UMAP_args)

    message("End: ", Sys.time())
  }

  ## SOMs do not like the input to be scaled
  ## whereas UMAP and tSNE do like it

  # SOM
  som.dims <- NULL
  if (run.som) {
    tryCatch({
      message("Calculating SOM.\nStart: ", Sys.time())

      som.dims <- ff_calc_som(exprs = expr.select,
                              seed = seed,
                              som_args = SOM_args,
                              embedsom_args = EmbedSOM_args)

      message("End: ", Sys.time())
    }, error = function(err) {
      message("SOM with error: ", err)
      som.dims <- NULL
    })
  }

  clust_idents <- NULL
  if (run.find.clusters) {
    clust_idents <- get_louvain_cluster(exprs = expr.select,
                                        FindNeighbors_args = FindNeighbors_args,
                                        FindClusters_args = FindClusters_args,
                                        mc.cores = mc.cores)
  }


  ## ---- prepare final fcs file -------

  if (write.transformed.channels.to.FCS && write.untransformed.channels.to.FCS) {
    ff <- ff_merge(ff.list,
                   grouplist = grouplist,
                   add_transformed_channels = T)
  } else if (write.transformed.channels.to.FCS) {
    ff <- ff_merge(ff.list,
                   grouplist = grouplist,
                   add_transformed_channels = F)
  } else if (write.untransformed.channels.to.FCS) {
    ff <- ff_merge(ff_transform(ff.list, trafolist = trafoname),
                   grouplist = grouplist,
                   add_transformed_channels = F)
  }

  ff <- ff_add_columns(ff, mat = do.call(cbind, list(pca.result[["x"]][, 1:run.pca],
                                                     umap.dims,
                                                     som.dims,
                                                     clust_idents)))


  # prepare channel desc
  # name.desc <- stats::setNames(ff.list[[1]]@parameters@data[["desc"]], ff.list[[1]]@parameters@data[["name"]])
  # name.desc <- name.desc[which(!is.na(name.desc))]
  # channel.desc <- rep("", ncol(dim.red.data))
  #
  # for (i in seq_along(name.desc)) {
  #   channel.desc[grep(names(name.desc)[i], colnames(dim.red.data))] <- name.desc[i]
  # }
  #
  # channel.desc_augment <- channel.desc
  # #channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("scaled$", colnames(dim.red.data))))] <- paste0(channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("scaled$", colnames(dim.red.data))))], "_scaled")
  # channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl(paste0(trafosuffix, "$"), colnames(dim.red.data))))] <- paste0(channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl(paste0(trafosuffix, "$"), colnames(dim.red.data))))], paste0("_", trafosuffix))
  #
  # channel.desc_augment[which(channel.desc_augment == "")] <- colnames(dim.red.data)[which(channel.desc_augment == "")]
  # channel.desc_augment <- make.names(channel.desc_augment)
  # names(channel.desc_augment) <- colnames(dim.red.data)

  # https://github.com/RGLab/flowCore/issues/201
  #flowCore::keyword(ff) <- flowCore:::updateTransformKeywords(ff)

  ## ---- calc cluster markers -------
  marker <- NULL
  if (calc_cluster_marker) {
    message("Calculating markers.")
    if (write.transformed.channels.to.FCS && write.untransformed.channels.to.FCS) {
      channels2 <- paste0(channels, formals(ff_merge)$suffix)
      exprs <- NULL
    } else if (write.transformed.channels.to.FCS) {
      channels2 <- channels
      exprs <- NULL
    } else if (write.untransformed.channels.to.FCS) {
      channels2 <- channels
      exprs <- do.call(rbind, lapply(ff.list, flowCore::exprs))
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
  if (!is.null(save.path) && !is.na(save.path)) {
    message("Writing files to disk.")
    t <- format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%Y%m%d_%H%M%S")
    dir.create(save.path, showWarnings = F)
    if (!is.null(save.name)) {
      save.name <- gsub("\\.fcs$", "", save.name, ignore.case = T)
      save.name <- save.name[1]
    }
    if (is.null(save.name)) {
      sv_pth <- file.path(save.path, paste0(t, "_dr_ff_list.rds"))
    } else {
      sv_pth <- file.path(save.path, paste0(save.name, ".rds"))
    }
    if ("rds" %in% save.to.disk) {
      saveRDS(list(flowframe = ff, marker = marker, pca = pca.result), file = sv_pth, compress = F)
    }
    if (is.null(save.name)) {
      sv_pth <- file.path(save.path, paste0(t, "_dr.fcs"))
    } else {
      sv_pth <- file.path(save.path, paste0(save.name, ".fcs"))
    }
    if ("fcs" %in% save.to.disk) {
      flowCore::write.FCS(ff, sv_pth)
    }
  }
  return(list(flowframe = ff, marker = marker, pca = pca.result))
}



# tsne
'  if (run.tsne) {
    message("Calculating tSNE.\nStart: ", Sys.time())

    if ("perplexity" %in% names(tSNE_args) && length(tSNE_args[["perplexity"]]) > 1) {
      tsne.dims <- do.call(cbind, parallel::mclapply(tSNE_args[["perplexity"]], function(z) {
        set.seed(seed)
        out <- Gmisc::fastDoCall(Rtsne::Rtsne, args = c(list(X = expr.select, verbose = F, perplexity = z), tSNE_args[which(!names(tSNE_args) %in% c("X", "verbose", "perplexity"))]))$Y
        colnames(out) <- c(paste0("tSNE_1_", z), paste0("tSNE_2_", z))
        return(out)
      }, mc.cores = mc.cores))
    } else {
      set.seed(seed)
      tsne.dims <- Gmisc::fastDoCall(Rtsne::Rtsne, args = c(list(X = expr.select), tSNE_args[which(names(tSNE_args) != "X")]))$Y
      colnames(tsne.dims) <- c("tSNE_1", "tSNE_2")
    }

    message("End: ", Sys.time())
  }'


## ---- clusterings -------
# find communities (clusters)
# louvain is run on the original data, not on any metaclusters
' if (run.louvain) {
    tryCatch({
      message("Calculating snn for louvain.\nStart: ", Sys.time())
      rownames(expr.select) <- 1:nrow(expr.select)
      snn <- Gmisc::fastDoCall(Seurat::FindNeighbors, args = c(list(object = Gmisc::fastDoCall(louvain.scale, expr.select)), FindNeighbors_args[which(names(FindNeighbors_args) != "object")]))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("error in snn calculation: ", err)
      run.louvain <- F
    })
  }

  if (run.louvain) {
    tryCatch({
      message("Finding clusters with Seurats implementation of the Louvain algorithm and parallel::mclapply using ", mc.cores," cores.\nStart: ",Sys.time())
      if ("resolution" %in% names(FindClusters_args) && length(FindClusters_args[["resolution"]]) > 1) {
        clust_idents <- do.call(cbind, parallel::mclapply(FindClusters_args[["resolution"]], function(x) {
          apply(Gmisc::fastDoCall(Seurat::FindClusters, args = c(list(object = snn$snn, resolution = x, verbose = F), FindClusters_args[which(!names(FindClusters_args) %in% c("object", "resolution", "verbose"))])), 2, as.numeric)
        }, mc.cores = mc.cores))
      } else {
        clust_idents <- apply(Gmisc::fastDoCall(Seurat::FindClusters, args = c(list(object = snn$snn), FindClusters_args[which(!names(FindClusters_args) %in% c("object"))])), 2, as.numeric)
      }


      colnames(clust_idents) <- paste0("louvain_", FindClusters_args[["resolution"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, clust_idents))
      print(apply(clust_idents, 2, function(x) length(unique(x))))
      message("End: ", Sys.time())

    }, error = function(err) {
      message("run.louvain with error: ", err)
    })
  }'



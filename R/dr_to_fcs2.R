#' Calculate dimension reduction and cluster annotation with data from one or more flow frames and add these parameters to a (concatenated) fcs file
#'
#'
#' @param ff.list a list of flowFrames as received from fcexpr::wsp_get_ff (compensated with Compensation Matrix as defined in FlowJo by default) or
#' as received from fcexpr::inds_get_ff (directly from FCS files, not compensated by default)
#' @param channels a named vector of channels to use for dimension reduction. values can be channel names (v-450LP..., b-520..., or so), or channel descriptions (e.g. CD3 or LAG3-PeCy7 for example)
#' names will be used as new description in the fcs file to be written; if no names provided, names of the very first FCS file will be used
#' @param add.sample.info named list of additional channels to identify samples or group them in flowjo;
#' e.g. for 9 fcs files to be concatenated: add.sample.info = list(condition = c(1,2,3,1,2,3,1,2,3), donor = c(1,1,1,2,2,2,3,3,3))
#' @param scale.whole if and how to scale channels after concatenation of flowframes in ff.list
#' @param scale.samples if and how to scale channels of flowframes in ff.list individually before concatenation
#' @param run.pca select number of pcs; if 0 no pca is computed
#' @param run.umap logical, whether to calculate UMAP dimension reduction with uwot::umap
#' @param run.som logical, whether to calculate SOM dimension reduction EmbedSOM::SOM followed by EmbedSOM::EmbedSOM
#' @param run.gqtsom logical, whether to calculate GQTSOM dimension reduction EmbedSOM::GQTSOM followed by EmbedSOM::EmbedSOM
#' @param find.clusters logical, whether to detect clusters with stats::dist, stats::hclust and stats::cutree
#' @param extra.cols numeric vector of an extra column (or matrix of multiple columns) with arbitraty information to add to the final fcs file;
#' has to be equal to the number of rows in all flowframes provided; colnames will be the channel names in the FCS file;
#' could be a previously calculated dimension reduction or cluster annotation.
#' @param clustering.for.marker.calc if NULL nothing is calculated; otherwise marker features (stained markers) are determined by wilcox test
#' using \href{https://github.com/immunogenomics/presto}{presto::wilcoxauc}
#' @param mc.cores number of cores to use for parallel computing of cluster markers and multiple cluster resolutions;
#' mc.cores is passed to parallel::mclapply or parallel::mcmapply; limited to parallel::detectCores()-1
#' @param save.to.disk what to save to disk: (concatenated) and appended FCS file and/or rds file with several elements in a list
#' @param save.path where to save elements specified in save.to.disk; set to NULL to have nothing written to disk
#' @param exclude.extra.channels when scaled and transform channels are written to FCS file, some channels may be redundant
#' and will only occupy disk space, those are specified here; matched with grepl
#' @param write.scaled.channels.to.FCS do save scaled channels (scale.whole, scale.samples) to FCS file
#' @param write.transformed.channels.to.FCS do save transformed channels (e.g. logicle or arcsinh) to FCS file
#' @param write.untransformed.channels.to.FCS do save untransformed channels (inverse) to FCS file
#' @param timeChannel name of the Time channel to exclude from all analyses and calculation; if NULL will be attempted
#' to be detected automatically
#' @param transformation_name name of the applied transformation (will appear in FCS file as channel desc)
#' @param return_processed_raw_data_only do not calculate dimension reduction etc but only return the preprocessed
#' data for external calculations or tryouts
#' @param seed set a seed for reproduction of dimension reductions
#' @param calc.global.markers logical whether to calculate global cluster markers: so each cluster vs. all other cells
#' @param calc.pairwise.markers logical whether to calculate pairwise cluster markers: so each cluster vs. each cluster
#' @param metacluster_map which map to use for metaclustering; SOM or GQTSOM
#' @param metacluster_map_source which source to use for metaclustering; actual expression valueS (fluorescence intensites)
#' or UMAP dimension (the latter of which was found to work really well)
#' @param UMAP_args arguments to uwot::umap; use pca and scale arguments to have data only modified for UMAP but not for SOM
#' @param SOM_args args to EmbedSOM::SOM
#' @param GQTSOM_args args to EmbedSOM::GQTSOM
#' @param EmbedSOM_args args to EmbedSOM::EmbedSOM
#' @param dist_args arguments to stats::dist
#' @param hclust_args arguments to stats::hclust
#' @param cutree_args arguments to stats::cutree; only k is relevant
#'
#' @return
#' A list with 3 elements: (i) The matrix of fluorescence intensities and appended information (dim red, clustering). This is the table which is written into a newly generated fcs file.
#' (ii) A character vector of meaningful column names which may be used for the table in R (rather for convenience). (iii) Tables of marker features (each cluster vs all other events and all clusters pairwise).
#' @export
#'
#' @importFrom magrittr "%>%"
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
dr_to_fcs2 <- function(ff.list,
                       channels = NULL,
                       add.sample.info = NULL,
                       scale.whole = c("none", "z.score", "min.max"),
                       scale.samples = c("none", "z.score", "min.max"),
                       run.pca = 0,
                       run.umap = T,
                       run.som = T,
                       run.gqtsom = F,
                       find.clusters = T,
                       metacluster_map = c("SOM", "GQTSOM"),
                       metacluster_map_source = c("UMAP", "expr"),
                       clustering.for.marker.calc = c("cutree_5", "cutree_10", "cutree_15"),
                       calc.global.markers = F,
                       calc.pairwise.markers = F,
                       UMAP_args = list(metric = "cosine", verbose = T, scale = T),
                       SOM_args = list(),
                       GQTSOM_args = list(),
                       EmbedSOM_args = list(),

                       dist_args = list(),
                       hclust_args = list(method = "average"),
                       cutree_args = list(k = c(5,10,15)),

                       extra.cols = NULL,
                       mc.cores = 1,
                       save.to.disk = c("fcs", "rds"),
                       save.path = file.path(getwd(), paste0(substr(gsub("\\.", "", make.names(as.character(Sys.time()))), 2, 15), "_FCS_dr")),
                       save.name = NULL,
                       exclude.extra.channels = ifelse(length(ff.list) == 1 && names(ff.list) == "transformed", "cluster.id", "FSC|SSC|Time|cluster.id"),
                       write.transformed.channels.to.FCS = T,
                       write.untransformed.channels.to.FCS = T,
                       write.scaled.channels.to.FCS = F,
                       timeChannel = c("Time", "HDR-T"),
                       transformation_name = "trans",
                       return_processed_raw_data_only = F,
                       seed = 42) {


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
  if (!is.null(clustering.for.marker.calc) && !requireNamespace("presto", quietly = T)) {
    devtools::install_github('immunogenomics/presto')
  }
  if (!is.null(clustering.for.marker.calc) && !requireNamespace("matrixTests", quietly = T)) {
    utils::install.packages("matrixTests")
  }
  if (run.som && !requireNamespace("EmbedSOM", quietly = T)) {
    devtools::install_github("exaexa/EmbedSOM")
  }

  # check names of ff.list entries
  if (any(!names(ff.list) %in% c("untransformed", "transformed"))) {
    stop("ff.list has to contain a list of flowframes named 'transformed' and optionally an additional list named 'untransformed' (inverse transformed, original as in flowjo.). Check the output of fcexpr::wsp_get_ff for valid input.")
  }

  # check if names in ff.list are the same
  if (length(unique(lengths(ff.list))) != 1) {
    stop("number of flowframes in untransformed and transformed has to be equal.")
  }

  if (!is.null(save.to.disk)) {
    save.to.disk <- match.arg(save.to.disk, c("fcs", "rds"), several.ok = T)
  }

  mc.cores <- min(mc.cores, parallel::detectCores() - 1)

  metacluster_map_source <- match.arg(metacluster_map_source, c("UMAP", "expr"))
  metacluster_map <- match.arg(metacluster_map, c("SOM", "GQTSOM"))

  if (!is.null(extra.cols)) {
    if (!is.matrix(extra.cols)) {
      extra.cols <- as.matrix(extra.cols)
    }
    if (any(!apply(extra.cols, 2, is.numeric))) {
      stop("extra.cols has to be a numeric matrix.")
    }
  }

  # check add.sample.info
  if (!is.null(add.sample.info)) {
    if (!is.list(add.sample.info)) {
      stop("add.sample.info has to be a list.")
    }
    if (is.null(names(add.sample.info))) {
      stop("add.sample.info has to have names. These names will become channel names in the FCS file.")
    }
    if (!all(unlist(lapply(add.sample.info, function(x) is.numeric(x))))) {
      stop("Please provide numeric values only in add.sample.info: E.g. use as.numeric(as.factor(x)), where x is numeric.")
    }
    if (any(unlist(lapply(add.sample.info, function(x) is.na(x))))) {
      stop("NA found in sample infos.")
    }
    if (!all(unlist(lapply(add.sample.info, function(x) length(x) == length(ff.list[[1]]))))) {
      stop("Length of each additional sample information has to match the length of selected samples, which is: ", length(ff.list[[1]]),".")
    }
  }

  # set scaling funs
  scale.samples <- switch(match.arg(scale.samples, c("none", "z.score", "min.max")),
                          z.score = scale,
                          min.max = min.max.normalization,
                          none = function(x) return(x))

  scale.whole <- switch(match.arg(scale.whole, c("none", "z.score", "min.max")),
                        z.score = scale,
                        min.max = min.max.normalization,
                        none = function(x) return(x))


  # check if channel names and desc are equal
  if (!is.null(fcs_check <- .check.ff.list(ff.list = ff.list, channels = channels, strict = T))) {
    return(fcs_check)
  }

  ## channel names from first ff
  channels <- .get.channels(ff = ff.list[[1]][[1]], timeChannel = timeChannel, channels = channels)


  # first step to create matrix for fcs file; put here to allow early cluster calculation which then
  # allows for lda calculation before umap, som, etc.
  # prepare matrix for FCS file

  ## option to save inverse (~untransformed) data only to FCS
  if (write.untransformed.channels.to.FCS && !"untransformed" %in% names(ff.list)) {
    message("Untransformed data not provided. Hence, cannot be saved to FCS.")
    write.untransformed.channels.to.FCS <- F
  }

  ## ---- prepare data -------
  ## initiate empty dim.red.data
  dim.red.data <- data.frame()
  ## write original data (transformed and/or untransformed) to fcs
  if (write.untransformed.channels.to.FCS && "untransformed" %in% names(ff.list)) {
    dim.red.data <- do.call(rbind, lapply(ff.list[["untransformed"]], flowCore::exprs))
  }

  if (write.transformed.channels.to.FCS && "transformed" %in% names(ff.list)) {
    expr_trans <- do.call(rbind, lapply(ff.list[["transformed"]], flowCore::exprs))
    expr_trans <- expr_trans[, which(!grepl(exclude.extra.channels, colnames(expr_trans)))]
    colnames(expr_trans) <- paste0(colnames(expr_trans), "_", transformation_name)
    if (nrow(dim.red.data) > 0) {
      dim.red.data <- cbind(dim.red.data, expr_trans)
    } else {
      dim.red.data <- expr_trans
    }
  }

  ## if transformed data is provided, these are used, if not then untransformed
  ff.list_index <- ifelse("transformed" %in% names(ff[["flowframes"]]), "transformed", "untransformed")

  ## write ident to fcs; any transformed or untransformed is fine
  dim.red.data <- cbind(dim.red.data, ident = rep(1:length(ff.list[[ff.list_index]]), sapply(ff.list[[ff.list_index]], nrow)))



  ## if data is scaled and pca applied this affects all algorithms (also SOM and GQTSOM which do not like it I think)
  ## if only for umap: pass parameters to umap

  ## apply scaling which was selected above and select channels to use for dimension reduction.
  expr.select <- scale.whole(do.call(rbind, lapply(ff.list[[ff.list_index]], function(x) scale.samples(flowCore::exprs(x)[, channels])))) # names to desc here?!


  ### allow to pass expr.select here.
  ## requires a lot of checking though

  pca.result <- NULL
  if (run.pca > 0) {
    run.pca <- min(ncol(expr.select) - 1, run.pca)
    message("Calculating PCA.\nStart: ", Sys.time())
    pca.result <- stats::prcomp(expr.select, scale. = F, center = F)
    pca.dims <- pca.result[["x"]]
    expr.select <- pca.dims[, 1:run.pca] # overwrite original data with PCA embedding
    message("Done. ", Sys.time())
  }

  # first exit here
  if (return_processed_raw_data_only) {
    return(expr.select)
  }

  ## ---- dim reds -------
  # umap
  if (run.umap) {
    message("Calculating UMAP.\nStart: ", Sys.time())

    if ("n_neighbors" %in% names(UMAP_args) && length(UMAP_args[["n_neighbors"]]) > 1) {
      umap.dims <- do.call(cbind, parallel::mclapply(UMAP_args[["n_neighbors"]], function(z) {
        set.seed(seed)
        out <- Gmisc::fastDoCall(uwot::umap, args = c(list(X = expr.select, verbose = F, n_neighbors = z), UMAP_args[which(!names(UMAP_args) %in% c("X", "verbose", "n_neighbors"))]))
        colnames(out) <- c(paste0("UMAP_1_", z), paste0("UMAP_2_", z))
        return(out)
      }, mc.cores = mc.cores))
    } else {
      set.seed(seed)
      umap.dims <- Gmisc::fastDoCall(uwot::umap, args = c(list(X = expr.select), UMAP_args[which(names(UMAP_args) != "X")]))
      colnames(umap.dims) <- c("UMAP_1", "UMAP_2")
    }

    message("End: ", Sys.time())
  }

  ## SOMs do not like the input to be scaled
  ## whereas UMAP and tSNE do like it

  # SOM
  if (run.som) {
    message("Calculating SOM.\nStart: ", Sys.time())
    set.seed(seed)
    som.map.dr <- Gmisc::fastDoCall(EmbedSOM::SOM, args = c(list(data = expr.select), SOM_args[which(!names(SOM_args) %in% c("data"))]))
    som.dims <- Gmisc::fastDoCall(EmbedSOM::EmbedSOM, args = c(list(data = expr.select, map = som.map.dr), EmbedSOM_args[which(!names(EmbedSOM_args) %in% c("data", "map"))]))
    colnames(som.dims) <- c("SOM_1", "SOM_2")
    message("End: ", Sys.time())
  }

  # GQTSOM
  if (run.gqtsom) {
    message("Calculating GQTSOM.\nStart: ", Sys.time())
    set.seed(seed)
    gqtsom.map.dr <- Gmisc::fastDoCall(EmbedSOM::GQTSOM, args = c(list(data = expr.select), GQTSOM_args[which(!names(GQTSOM_args) %in% c("data"))]))
    gqtsom.dims <- Gmisc::fastDoCall(EmbedSOM::EmbedSOM, args = c(list(data = expr.select, map = gqtsom.map.dr), EmbedSOM_args[which(!names(EmbedSOM_args) %in% c("data", "map"))]))
    colnames(som.dims) <- c("GQTSOM_1", "GQTSOM_2")
    message("End: ", Sys.time())
  }

  ## ---- metaclustering -------
  ## decide here how to run metaclustering
  ## and whether to run it on umap or tsne dimensions
  # metacluster_map = SOM or GQTSOM
  # metacluster_map_source = expression or UMAP

  ## always run metaclustering as other procedures are too slow
  if (metacluster_map_source == "expr") {
    if (metacluster_map == "SOM") {
      som.map.clust <- som.map.dr
    }
    if (metacluster_map == "GQTSOM") {
      som.map.clust <- gqtsom.map.dr
    }

  } else if (metacluster_map_source == "UMAP") {
    if (metacluster_map == "SOM") {
      som.map.clust <- Gmisc::fastDoCall(EmbedSOM::SOM, args = c(list(data = umap.dims), SOM_args[which(!names(SOM_args) %in% c("data"))]))
    }
    if (metacluster_map == "GQTSOM") {
      som.map.clust <- Gmisc::fastDoCall(EmbedSOM::GQTSOM, args = c(list(data = umap.dims), GQTSOM_args[which(!names(GQTSOM_args) %in% c("data"))]))
    }
  }

  if (find.clusters) {
    tryCatch({
      message("Finding clusters with hclust.\nStart: ", Sys.time())

      d <- Gmisc::fastDoCall(stats::dist, args = c(list(x = som.map.clust$codes), dist_args[which(names(dist_args) != "x")]))
      h <- Gmisc::fastDoCall(stats::hclust, args = c(list(d = d), hclust_args[which(names(hclust_args) != "d")]))
      ks <- Gmisc::fastDoCall(cbind, lapply(cutree_args[which(names(cutree_args) == "k")], function(x) stats::cutree(tree = h, k = x)))

      ks <- apply(ks, 2, function (x) x[som.map.clust[["mapping"]][,1]])

      # make sure that cluster 1 is the largest and so on
      ks <- .cluster_ordering(ks = ks)

      colnames(ks) <- paste0("cutree_", cutree_args[which(names(cutree_args) == "k")])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("find.clusters with error: ", err)
    })
  }


  ## ---- prepare final fcs file -------
  if (write.scaled.channels.to.FCS) {
    scaled.expr <- scale.whole(do.call(rbind, lapply(ff.list[[ff.list_index]], function(x) scale.samples(flowCore::exprs(x)[, channels]))))
    scaled.expr <- scaled.expr[, which(!grepl(exclude.extra.channels, colnames(scaled.expr)))]
    colnames(scaled.expr) <- paste0(colnames(scaled.expr), "_scaled")
    dim.red.data <- do.call(cbind, list(dim.red.data, scaled.expr))
  }

  if (run.pca) {
    dim.red.data <- do.call(cbind, list(dim.red.data, pca.dims[,1:n.pca.dims]))
  }
  if (run.umap) {
    dim.red.data <- do.call(cbind, list(dim.red.data, umap.dims))
  }
  #if (run.tsne) {dim.red.data <- do.call(cbind, list(dim.red.data, tsne.dims))}

  if (run.som) {
    dim.red.data <- do.call(cbind, list(dim.red.data, som.dims))
  }
  if (run.gqtsom) {
    dim.red.data <- do.call(cbind, list(dim.red.data, gqtsom.dims))
  }

  if (!is.null(extra.cols)) {
    if (nrow(extra.cols) == nrow(dim.red.data)) {
      dim.red.data <- cbind(dim.red.data, extra.cols)
    } else {
      message("extra.cols not added due to wrong row number.")
    }
  }

  # 2021 06 17 necessary
  dim.red.data <- as.data.frame(dim.red.data)

  if (!is.null(add.sample.info)) {
    tryCatch({
      for (i in names(add.sample.info)) {
        dim.red.data <- do.call(cbind, list(dim.red.data, rep(add.sample.info[[i]], times = as.numeric(table(rep(1:length(ff.list[[ff.list_index]]), sapply(ff.list[[ff.list_index]], nrow)))))))
        names(dim.red.data)[length(dim.red.data)] <- i
      }
    }, error = function(err) {
      message("Addition of sample info caused an error: ", err)
    })
  }

  # prepare channel desc
  name.desc <- stats::setNames(ff.list[[1]][[1]]@parameters@data[["desc"]], ff.list[[1]][[1]]@parameters@data[["name"]])
  name.desc <- name.desc[which(!is.na(name.desc))]
  channel.desc <- rep("", ncol(dim.red.data))

  for (i in seq_along(name.desc)) {
    channel.desc[grep(names(name.desc)[i], colnames(dim.red.data))] <- name.desc[i]
  }

  channel.desc_augment <- channel.desc
  channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("scaled$", colnames(dim.red.data))))] <- paste0(channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("scaled$", colnames(dim.red.data))))], "_scaled")
  channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl(paste0(transformation_name, "$"), colnames(dim.red.data))))] <- paste0(channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl(paste0(transformation_name, "$"), colnames(dim.red.data))))], paste0("_", transformation_name))

  channel.desc_augment[which(channel.desc_augment == "")] <- colnames(dim.red.data)[which(channel.desc_augment == "")]
  channel.desc_augment <- make.names(channel.desc_augment)
  names(channel.desc_augment) <- colnames(dim.red.data)


  ## ---- write flowframe file -------

  # get common (intersecting keywords)
  # a bit unhandy but vectorized version (below) did not work
  all <- lapply(ff.list[[1]], flowCore::keyword)
  new_kw <- as.list(unlist(lapply(names(all[[1]]), function(x) {
    if (length(unique(unlist(sapply(all, "[", x)))) == 1) {
      return(stats::setNames(unique(unlist(sapply(all, "[", x))), nm = x))
    }
  })))

  #new_kw <- Reduce(intersect, lapply(ff.list[[1]], flowCore::keyword))
  #names(new_kw) <- names(flowCore::keyword(ff.list[[1]][[1]]))[match(new_kw, flowCore::keyword(ff.list[[1]][[1]]))]
  #new_kw <- flowCore::keyword(ff.list[[1]][[1]])
  new_pars <- flowCore::parameters(ff.list[[1]][[1]])

  new_kw[["$PAR"]] <- as.character(ncol(dim.red.data))
  n <- nrow(flowCore::parameters(ff.list[[1]][[1]])) + 1
  for (z in n:ncol(dim.red.data)) {
    new_p <- flowCore::parameters(ff.list[[1]][[1]])[1,]
    new_p_number <- z
    rownames(new_p) <- c(paste0("$P", new_p_number))

    new_pars <- BiocGenerics::combine(new_pars, new_p)

    new_p_name <- names(dim.red.data)[z]
    new_p_desc <- channel.desc[z]
    flowCore::pData(new_pars)$name[new_p_number] <- new_p_name
    flowCore::pData(new_pars)$desc[new_p_number] <- new_p_desc
    flowCore::pData(new_pars)$minRange[new_p_number] <- as.integer(round(min(dim.red.data[, z])))
    flowCore::pData(new_pars)$maxRange[new_p_number] <- as.integer(round(max(dim.red.data[, z])))
    flowCore::pData(new_pars)$range[new_p_number] <- as.integer(round(max(dim.red.data[, z])))

    #### WRITE KEYWORD WITH 2 BRACKETS!! OTHERWISE FLOWJO DOES NOT READ IT #### ????????????
    new_kw[[paste0("$P", new_p_number, "N")]] <- new_p_name
    new_kw[[paste0("$P", new_p_number, "S")]] <- new_p_desc
    new_kw[[paste0("$P", new_p_number, "E")]] <- "0,0"
    new_kw[[paste0("$P", new_p_number, "G")]] <- "1"
    new_kw[[paste0("$P", new_p_number, "B")]] <- new_kw[["$P1B"]]
    new_kw[[paste0("$P", new_p_number, "R")]] <- as.integer(round(max(dim.red.data[, z])))
    new_kw[[paste0("$P", new_p_number, "V")]] <- "1"
    new_kw[[paste0("$P", new_p_number, "DISPLAY")]] <- "LIN"
    new_kw[[paste0("flowCore_$P", new_p_number, "Rmin")]] <- as.integer(round(min(dim.red.data[, z])))
    new_kw[[paste0("flowCore_$P", new_p_number, "Rmax")]] <- as.integer(round(max(dim.red.data[, z])))
  }

  for (z in 1:n) {
    flowCore::pData(new_pars)$minRange[z] <- as.integer(round(min(dim.red.data[, z])))
    flowCore::pData(new_pars)$maxRange[z] <- as.integer(round(max(dim.red.data[, z])))
    flowCore::pData(new_pars)$range[z] <- as.integer(round(max(dim.red.data[, z])))
    if (!grepl("FSC|SSC", new_kw[paste0("$P", as.character(z), "N")])) {
      new_kw[[paste0("$P", as.character(z), "R")]] <- as.integer(round(max(dim.red.data[, z])))
    }
    new_kw[[paste0("flowCore_$P", as.character(z), "Rmin")]] <- as.integer(round(min(dim.red.data[, z])))
    new_kw[[paste0("flowCore_$P", as.character(z), "Rmax")]] <- as.integer(round(max(dim.red.data[, z])))
  }
  ff <- methods::new("flowFrame", exprs = as.matrix(dim.red.data), parameters = new_pars, description = new_kw)
  # https://github.com/RGLab/flowCore/issues/201
  #flowCore::keyword(ff) <- flowCore:::updateTransformKeywords(ff)

  ## ---- calc cluster markers -------
  ## always use transformed data if provided
  marker <- NULL
  if (!is.null(clustering.for.marker.calc)) {
    message("Calculating markers.")
    tryCatch({
      marker <- lapply(clustering.for.marker.calc, function (clust_col) {
        message("Clustering: ", clust_col)
        if (ff.list_index == "transformed") {
          # when transformed data is provided always use them for marker calc
          channels <- paste0(channels, "_", transformation_name)
        }
        # do not use expr.select which may have become dimensions of LDA
        dat <- as.matrix(dim.red.data[,which(colnames(dim.red.data) %in% channels),drop=F])
        cluster <- dim.red.data[,which(colnames(dim.red.data) == clust_col),drop=T]

        message("Cluster sizes: ")
        print(table(cluster))
        #cluster <- cluster[which(table(cluster) > cluster.marker.min.cells)]


        ## alternatives to dip.test:
        # silverman test: https://github.com/jenzopr/silvermantest
        # multimode: https://cran.r-project.org/web/packages/multimode/multimode.pdf

        # global markers
        if (calc.global.markers) {
          message("Global markers.\nStart: ", Sys.time())
          marker_table <- .calc.global.cluster.marker(dat = dat, cluster = cluster, mc.cores = mc.cores)
          marker_table[,"channel_desc"] <- channel.desc_augment[marker_table[,"channel"]]
          message("End: ", Sys.time())
        } else {
          marker_table <- NULL
        }

        ## pairwise markers
        if (calc.pairwise.markers) {
          message("Pairwise markers.\nStart: ", Sys.time())
          pair_marker_table <- .calc.pairwise.cluster.marker(dat = dat, cluster = cluster, mc.cores = mc.cores)
          pair_marker_table[,"channel_desc"] <- channel.desc_augment[pair_marker_table[,"channel"]]
          message("End: ", Sys.time())
        } else {
          pair_marker_table <- NULL
        }

        return(list(marker_table = marker_table, pairwise_marker_table = pair_marker_table, cluster_sizes = table(cluster)))
      })
      names(marker) <- clustering.for.marker.calc
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
      sv_pth <- file.path(save.path, paste0(t, "_dr_ff_list.rds")
    } else {
      sv_pth <- file.path(save.path, paste0(save.name, ".rds")
    }
    if ("rds" %in% save.to.disk) {
      saveRDS(list(df = dim.red.data, col_names = channel.desc_augment, marker = marker, flowframe = ff, pca = pca.result), file = sv_pth, compress = F)
    }
    if (is.null(save.name)) {
      sv_pth <- file.path(save.path, paste0(t, "_dr.fcs"))
    } else {
      sv_pth <- file.path(save.path, paste0(save.name, ".fcs")
    }
    if ("fcs" %in% save.to.disk) {
      flowCore::write.FCS(ff, sv_pth)
    }
  }
  return(list(df = dim.red.data, col_names = channel.desc_augment, marker = marker, flowframe = ff, pca = pca.result))
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



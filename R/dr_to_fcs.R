#' Calculate dimension reductions and cluster annotations with data from one or more flow frames and add these parameters to a (concanetated fcs file)
#'
#'
#'
#'
#'  logicle trans: 2006 - Parks_A New Logicle Display Method Avoids Deceptive Effects of Logarithmic Scaling for Low Signals and Compensated Data
#'
#' @param ff.list a list of flowFrames as received from fcexpr::wsp_get_ff (compensated with Compensation Matrix as defined in FlowJo by default) or
#' as received from fcexpr::inds_get_ff (directly from FCS files, not compensated by default)
#' @param channels a named vector of channels to use for dimension reduction. values can be channel names (v-450LP..., b-520..., or so), or channel descriptions (e.g. CD3 or LAG3-PeCy7 for example)
#' names will be used as new description in the fcs file to be written; if no names provided, names of the very first FCS file will be used
#' @param add.sample.info named list of additional channels to identify samples or group them in flowjo;
#' e.g.: add.sample.info = list(condition = c(1,2,3,1,2,3,1,2,3), donor = c(1,1,1,2,2,2,3,3,3))
#' @param scale.whole if and how to scale channels after concatenation of flowframes in ff.list
#' @param scale.samples if and how to scale channels of flowframes in ff.list individually before concatenation
#' @param run.harmony attempt batch correction using harmony::HarmonyMatrix; if TRUE, then harmony__meta_data has to be provided in ... indicating the groups to be batch corrected;
#' harmony is conducted before run.pca; to calculate a pca before harmony, pass harmony__do_pca = T and optional a number of pcs with harmony__npcs. Set run.pca = F
#' when a pca is calculated in harmony.
#' @param run.pca run principle component analysis before dimension reduction. may not be necessary when less than 10 or so are provided in channels.
#' on the other hand if you decide to pass all channels (including those without any actual staining), then pca will focus on the ones with highest variances.
#' @param n.pca.dims number of principle components to calculate; generally, with 10 channels or less PCA may not be necessary to calculate;
#' e.g. if you choose only channel which do have an amount a variation (different populations) then PCA is not required; if you are lazy
#' and simply select all available channels, PCA will select the most relevant number of dimensions for you
#' @param run.lda run Linear Discriminant Analysis before dimension reduction;
#' should be F (FALSE) or a clustering calculated before, e.g. louvain_0.8 or leiden_1.1, kmeans_12 etc.; respective clustering calculation
#' has to be provided in ...
#' @param run.umap calculate UMAP dimension reduction with uwot::umap
#' @param run.tsne calculate tsne dimension reduction with Rtsne::Rtsne
#' @param run.som calculate SOM dimension reduction EmbedSOM::SOM followed by EmbedSOM::EmbedSOM
#' @param run.gqtsom calculate GQTSOM dimension reduction EmbedSOM::GQTSOM followed by EmbedSOM::EmbedSOM
#' @param run.louvain detect clusters (communities, groups) of cells with the louvain algorithm, implemented in Seurat::FindClusters (subsequent to snn detection by Seurat::FindNeighbors)
#' @param run.leiden detect clusters (communities, groups) of cells with the leiden algorithm, with leiden::leiden (subsequent to snn detection by Seurat::FindNeighbors)
#' @param run.kmeans detect clusters with stats::kmeans; will be the quickest way to get cluster annotation!!!
#' @param run.hclust detect clusters with stats::dist, stats::hclust and stats::cutree
#' @param run.flowClust detect clusters with flowClust::flowClust
#' @param run.MUDAN detect clusters with MUDAN::getComMembership (k = as.numeric(names(run.MUDAN))); e.g. run.MUDAN = setNames(TRUE, 8)
#' @param extra.cols vector of one extra column (or matrix of multiple columns) to add to the final fcs file;
#' has to be numeric; has to be equal to the number of rows of all flowframes provided; colnames of matrix dictate
#' channel names in the FCS file
#' @param calc.cluster.markers if NULL nothing is calculated; otherwise provide the clustering(s) for which cluster markers are to be determined,
#' using matrixStats::col_wilcoxon_twosample every cluster is compared to all other cells as well as all clusters pairwise.
#' respective clustering calculation has to be provided in ...; e.g. if louvain__resolution = 0.5 is provided set calc.cluster.markers = louvain_0.5;
#' and if in addition leiden__resolution_parameter = 0.7 then set calc.cluster.markers = c(louvain_0.5, leiden_0.7).
#' @param mc.cores mc.cores to calculate clusterings, limited to parallel::detectCores()-1
#' @param save.to.disk what to save to disk: (concatenated) and appended FCS file and/or rds file with several elements in a list
#' @param save.path where to save elements specified in save.to.disk; set to NULL to have nothing written to disk
#' @param exclude.extra.channels when scaled and transform channels are written to FCS file, some channels may be redundant
#' and will only occupy disk space, those are specified here; matched with grepl
#' @param write.scaled.channels.to.FCS do save scaled channels (scale.whole, scale.samples) to FCS file
#' @param timeChannel name of the Time channel to exclude from all analyses and calculation; if NULL will be attempted
#' to be detected automatically
#' @param ... additional parameters to calculations of UMAP, tSNE, som, gqtsom, EmbedSOM, louvain, leiden, harmony, hclust, flowClust, MUDAN, kmeans;
#' provide arguments as follows: UMAP__n_neighbors = c(15,20,25), or tsne__theta = 0.3, etc.
#' see respected help files to get to know which arguments can be passed:
#' uwot::umap, Rtsne::Rtsne, EmbedSOM::SOM, EmbedSOM::GQTSOM, EmbedSOM::EmbedSOM, harmony::HarmonyMatrix, flowClust::flowClust,
#' louvain: Seurat::FindNeighbors and Seurat::FindCluster, leiden: Seurat::FindNeighbors and leiden::leiden.
#' hclust: stats::dist and stats::hclust, MUDAN: MUDAN::getComMembership, stats::kmeans
#'
#' @return
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#'\dontrun{
#'########################################
#'### Plot cluster markers with ggplot ###
#'########################################
#' dr <- fcexpr::dr_to_fcs(ff.list = ffs,
#' channels = channels,
#' louvain__resolution = 0.5,
#' run.lda = "louvain_0.5",
#' calc.cluster.markers = c("louvain_0.5"),
#' save.path = NULL)
#' marker <- dr[[3]][[1]][[1]]
#'marker$channel_desc2 <- sapply(strsplit(marker$channel_desc, "_"), "[", 1)
#'marker <-
#'  marker %>%
#'  dplyr::mutate(pvalue = ifelse(pvalue == 0, 1e-300, marker$pvalue))
#'dplyr::group_by(channel_desc2) %>%
#'  dplyr::mutate(mean_scaled = fcexpr:::min.max.normalization(mean))

#'ggplot(marker, aes(x = as.factor(cluster), y = channel_desc2, fill = -log10(pvalue))) +
#'  geom_tile(color = "black") +
#'  theme_bw() +
#'  geom_text(aes(label = diff_sign)) +
#'  scale_fill_viridis_c()

#'ggplot(marker, aes(x = as.factor(cluster), y = channel_desc2, fill = mean_diff)) +
#'  geom_tile(color = "black") +
#'  theme_bw() +
#'  geom_text(aes(label = diff_sign)) +
#'  scale_fill_viridis_c()

#'ggplot(marker, aes(x = as.factor(cluster), y = channel_desc2, fill = mean_scaled)) +
#'  geom_tile(color = "black") +
#'  theme_bw() +
#'  scale_fill_viridis_c()

#'ggplot(marker, aes(x = mean_diff, y = -log10(pvalue), label = channel_desc2)) + #color = channel_desc2,
#'  geom_point() +
#'  theme_bw() +
#'  labs(title = "cluster markers (vs all other cells each)") +
#'  ggrepel::geom_text_repel(max.overlaps = 20, show.legend = F) +
#'  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), strip.background = element_rect(fill = "hotpink2")) +
#'  geom_vline(xintercept = 0, col = "tomato2", linetype = "dashed") +
#'  geom_hline(yintercept = 100, col = "tomato2", linetype = "dashed") +
#'  facet_wrap(vars(cluster))
#'}
dr_to_fcs <- function(ff.list,
                      channels = NULL,
                      add.sample.info,
                      scale.whole = c("z.score", "min.max", "none"),
                      scale.samples = c("none", "z.score", "min.max"),
                      run.harmony = F,
                      run.pca = F,
                      run.lda = F,
                      run.umap = T,
                      run.tsne = F,
                      run.som = T,
                      run.gqtsom = T,
                      run.louvain = T,
                      run.kmeans = F,
                      run.leiden = F,
                      run.hclust = F,
                      run.flowClust = F,
                      run.MUDAN = F,
                      n.pca.dims = 0,
                      calc.cluster.markers = NULL,
                      extra.cols,
                      mc.cores = 1,
                      save.to.disk = c("fcs", "rds"),
                      save.path = file.path(getwd(), paste0(substr(gsub("\\.", "", make.names(as.character(Sys.time()))), 2, 15), "_FCS_dr")),
                      exclude.extra.channels = ifelse(length(ff.list) == 1 && names(ff.list) == "logicle", "cluster.id", "FSC|SSC|Time|cluster.id"),
                      write.scaled.channels.to.FCS = T,
                      timeChannel = "Time",
                      ...) {
  # batch effect correction: https://cytekbio.com/blogs/blog/how-to-identify-and-prevent-batch-effects-in-longitudinal-flow-cytometry-research-studies
  # cytonorm (https://github.com/saeyslab/CytoNorm) requires reference sample for every batch - not always available
  # check somewhen: https://github.com/casanova-lab/iMUBAC
  # harmony: https://github.com/immunogenomics/harmony
  # MUDAN: https://github.com/JEFworks/MUDAN


  ### to do:
  # add: MUDAN::clusterBasedBatchCorrect
  # harmony correction in clusters - then recompute everything?! - kind of circular
  # harmony could be done at beginning (e.g unbiased or based on detected clusters)

  ## preprocessCore::normalize.quantiles() - allow to normalize channel of ffs within defined groups
  # but this can put a really strong bias on the data:
  'df <- data.frame(x1 = c(rnorm(1e5,0,1), rnorm(1e4,15,1)),
                   x2 = c(rnorm(1e5,1,1), rnorm(1e4,9,1)),
                   x3 = c(rnorm(1e4,-1,4), rnorm(1e5,30,5)))
  df <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(df)))
  names(df) <- paste0("x", 1:ncol(df))
  df <-
    df %>%
    tidyr::pivot_longer(names_to = "name", values_to = "value", cols = c(x1,x2,x3))


  ggplot(df, aes(x = value, y = name))+
    ggridges::geom_density_ridges() +
    ggridges::theme_ridges()'

  ## allow to provide expr.select directly instead of ff.list

  if (!requireNamespace("Rtsne", quietly = T)) utils::install.packages("Rtsne")
  if (run.umap && !requireNamespace("uwot", quietly = T)) {
    utils::install.packages("uwot")
  }
  if (run.leiden && !requireNamespace("Seurat", quietly = T)) {
    utils::install.packages("Seurat")
  }
  if (!requireNamespace("parallel", quietly = T)) {
    utils::install.packages("parallel")
  }
  if (run.leiden &&!requireNamespace("leiden", quietly = T)) {
    utils::install.packages("leiden")
  }
  if (!requireNamespace("Biobase", quietly = T)) {
    BiocManager::install("Biobase")
  }
  if (run.flowClust && !requireNamespace("flowClust", quietly = T)) {
    BiocManager::install("flowClust")
  }
  if (!requireNamespace("devtools", quietly = T)) {
    utils::install.packages("devtools")
  }
  if (run.harmony && !requireNamespace("run.harmony", quietly = T)) {
    devtools::install_github("immunogenomics/harmony")
  }
  if (run.MUDAN && !requireNamespace("MUDAN", quietly = T)) {
    devtools::install_github("JEFworks/MUDAN")
  }
  if (!"logicle" %in% names(ff.list)) {
    stop("logicle transformed has to be in ff.list.")
  }

  dots <- list(...)

  if (any(!names(ff.list) %in% c("inverse", "logicle"))) {
    stop("ff.list has to contain a list of flowframes named 'logicle' (logicle transformed)
         and optionally an additional list named 'inverse' (inverse transformed, original as in flowjo.).")
  }

  if (length(unique(lengths(ff.list))) != 1) {
    stop("number of flowframes in inverse and logicle has to be equal.")
  }
  # check if names in ff.list inverse and logicle are the same

  for (par in c("louvain", "leiden", "umap", "tsne", "som", "gqtsom", "harmony")) {
    if (any(grepl(paste0("^", par, "__"), names(dots), ignore.case = T)) &&!eval(rlang::sym(paste0("run.", par)))) {
      message(paste0(par, " parameters provided in ... but ", "'run.", par, " = F'."))
    }
  }

  if (run.harmony && !any(grepl("^harmony__meta_data", names(dots)))) {
    stop("When 'run.harmony = T' harmony__meta_data has to be provided in ..., see ?harmony::HarmonyMatrix.")
  }

  if (run.hclust && !any(grepl("^hclust__k", names(dots)))) {
    stop("When 'run.hclust = T' hclust__k has to be provided in ..., see ?stats::cutree.")
  }

  if (run.flowClust && !any(grepl("^flowClust__K", names(dots)))) {
    stop("When 'run.flowClust = T' flowClust__K has to be provided in ..., see ?flowClust::flowClust.")
  }

  if (run.MUDAN && !any(grepl("^MUDAN__k", names(dots)))) {
    stop("When 'run.MUDAN = T' MUDAN__k (e.g. MUDAN__k = 30), has to be provided in ..., see ?MUDAN::getComMembership.")
  }

  if (run.kmeans && !any(grepl("^kmeans__centers", names(dots)))) {
    stop("When 'run.kmeans = T' kmeans__centers, has to be provided in ..., see ?stats::kmeans")
  }

  if (run.louvain && !any(grepl("^louvain__resolution", names(dots)))) {
    stop("When 'run.louvain = T' louvain__resolution, has to be provided in ..., see ?Seurat::FindClusters")
  }

  if (run.leiden && !any(grepl("^leiden__resolution_parameter", names(dots)))) {
    stop("When 'run.leiden = T' leiden__resolution_parameter, has to be provided in ..., see ?leiden::leiden")
  }

  if (n.pca.dims < 0) {
    run.pca <- F
  }

  if (run.pca && run.lda) {
    stop("run.pca = T and run.lda = T at the same time is not possible.")
  }

  if (is.logical(run.lda) && run.lda) {
    stop("Do not set 'run.lda = T' but provide a clustering that should be used to calculate it, e.g. a pattern like louvain_0.4.")
  }

  if (!is_logical(run.lda)) {
    # clustering columns need to follow the pattern name_resolution.
    if (!run.lda %in% paste0(sapply(strsplit(names(dots), "__"), "[", 1), "_", dots)) {
      stop("Value for run.lda not found in ... . They respective clustering to use for lda has to be computed! E.g. if run.lda = 'louvain_0.9' then louvain__resolution = 0.9 has to be passed.")
    }
  }

  if (run.harmony && any(grepl("^harmony__do_pca", names(dots))) && dots[["harmony__do_pca"]] == T && run.pca) {
    warning("harmony is calculated with do_pca = T, hence a subsequnt pca (run.pca = T) is not required. Consider setting run.pca to FALSE.")
  }

  if (!is.null(calc.cluster.markers)) {
    if (!any(calc.cluster.markers %in% paste0(sapply(strsplit(names(dots), "__"), "[", 1), "_", dots))) {
      stop("calc.cluster.markers: ", calc.cluster.markers[which(!calc.cluster.markers %in% paste0(sapply(strsplit(names(dots), "__"), "[", 1), "_", dots))], " not found in ... .")
    }
  }

  if (!is.null(save.to.disk)) {
    save.to.disk <- match.arg(save.to.disk, c("fcs", "rds"), several.ok = T)
  }

  mc.cores <- min(mc.cores, parallel::detectCores() - 1)

  if (!missing(extra.cols)) {
    if (!is.matrix(extra.cols)) {
      extra.cols <- as.matrix(extra.cols)
    }
    if (any(!apply(extra.cols, 2, is.numeric))) {
      stop("extra.cols has to be a numeric matrix.")
    }
  }

  # check add.sample.info
  if (!missing(add.sample.info)) {
    if (!is.list(add.sample.info)) {
      stop("add.sample.info has to be a list.")
    }
    if (is.null(names(add.sample.info))) {
      stop("add.sample.info has to have names. These names become channel names in the FCS file.")
    }
    if (!all(unlist(lapply(add.sample.info, function(x) is.numeric(x))))) {
      stop("Please only provide numeric values as additional sample infos.")
    }
    if (any(unlist(lapply(add.sample.info, function(x) is.na(x))))) {
      stop("NA found in sample infos.")
    }
  }

  if (!missing(add.sample.info)) {
    if (!all(unlist(lapply(add.sample.info, function(x) length(x) == length(ff.list[[1]]))))) {
      stop(paste0("Length of each additional sample information has to match the length of selected samples, which is: ", length(ff.list[[1]]),".")
      )
    }
  }

  scale.samples <-
    switch(match.arg(scale.samples, c("none", "z.score", "min.max")),
           z.score = scale,
           min.max = min.max.normalization,
           none = function(x) {
             return(x)
           }
    )

  scale.whole <-
    switch(match.arg(scale.whole, c("z.score", "min.max", "none")),
           z.score = scale,
           min.max = min.max.normalization,
           none = function(x) {
             return(x)
           }
    )

  # check if channel names and desc are equal
  .check.ff.list(ff.list = ff.list)

  channels <- .get.channels(
    ff = ff.list[["logicle"]][[1]],## channel names from first ff
    timeChannel = timeChannel,
    channels = channels
  )

  # overwrite channel desc in ffs
  # correct order is important, as provided by .get.channels
  for (i in seq_along(ff.list[["logicle"]])) {
    ff.list[["logicle"]][[i]]@parameters@data[["desc"]][which(ff.list[["logicle"]][[i]]@parameters@data[["name"]] %in% channels)] <- names(channels)
  }
  if ("inverse" %in% names(ff.list)) {
    for (i in seq_along(ff.list[["inverse"]])) {
      ff.list[["inverse"]][[i]]@parameters@data[["desc"]][which(ff.list[["inverse"]][[i]]@parameters@data[["name"]] %in% channels)] <- names(channels)
    }
  }

  # first step to create matrix for fcs file; put here to allow early cluster calculation which then
  # allows for lda calculation before umap, som, etc.
  # prepare matrix for FCS file
  expr.logicle <- do.call(rbind, lapply(ff.list[["logicle"]], function(x) flowCore::exprs(x)))
  expr.logicle <- expr.logicle[, which(!grepl(exclude.extra.channels, colnames(expr.logicle)))]
  colnames(expr.logicle) <- paste0(colnames(expr.logicle), "_logicle")

  if ("inverse" %in% names(ff.list)) {
    dim.red.data <- do.call(cbind, list(do.call(rbind, lapply(ff.list[["inverse"]], function(x) flowCore::exprs(x))), expr.logicle, ident = rep(1:length(ff.list[["logicle"]]), sapply(ff.list[["logicle"]], nrow))))
  } else {
    dim.red.data <- do.call(cbind, list(expr.logicle, ident = rep(1:length(ff.list[["logicle"]]), sapply(ff.list[["logicle"]], nrow))))
  }

  ## apply scaling which was selected above and select channels to use for dimension reduction.
  expr.select <- scale.whole(do.call(rbind, lapply(ff.list[["logicle"]], function(x) scale.samples(flowCore::exprs(x)[, channels]))))


  ### allow to pass expr.select here.


  # run harmony
  # if harmony is run with do_pca = T a subsequent pca should not be computed
  if (run.harmony) {
    # https://portals.broadinstitute.org/harmony/articles/quickstart.html
    # https://portals.broadinstitute.org/harmony/advanced.html
    temp_dots <- dots[which(grepl("^harmony__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^harmony__", "", names(temp_dots), ignore.case = T)
    if (!any(grepl("do_pca", names(temp_dots)), ignore.case = T)) {
      print("'do_pca set to FALSE' in harmony::HarmonyMatrix.")
      temp_dots <- c(temp_dots, do_pca = F)
    }
    expr.select <- do.call(harmony::HarmonyMatrix, args = c(list(data_mat = expr.select), temp_dots))
  }

  pca.result <- NULL
  if (run.pca) {
    if (n.pca.dims == 0) {
      n.pca.dims <- ncol(expr.select) - 1
    } else if (n.pca.dims > 0) {
      n.pca.dims <- min(ncol(expr.select) - 1, n.pca.dims)
    }
    print(paste0("Calculating PCA. Start: ", Sys.time()))
    # https://slowkow.com/notes/pca-benchmark/
    # mat_irlba2 <- irlba::irlba(A = expr.select, nv = n.pca.dims)
    # mat_irlba2$x <- mat_irlba2$u %*% diag(mat_irlba2$d)
    pca.result <- stats::prcomp(expr.select, scale. = F, center = F)
    pca.matrix <- pca.result[["x"]]
    expr.select <- pca.matrix[, 1:n.pca.dims]
    print(paste0("Done. ", Sys.time()))
  }

  # find communities (clusters)
  tryCatch(
    if (run.louvain || run.leiden) {
      temp_dots <- dots[which(grepl("^louvain__|^leiden__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^louvain__|^leiden__", "", names(temp_dots), ignore.case = T)

      if (!any(grepl("annoy.metric", names(temp_dots), ignore.case = T))) {
        print("snn calculation with annoy.metric = 'cosine' by default.")
        temp_dots <- c(temp_dots, annoy.metric = "cosine")
      }
      print(paste0("Calculating snn for louvain and/or leiden. Start: ", Sys.time()))
      rownames(expr.select) <- 1:nrow(expr.select)
      # Seurat::FindNeighbors ignores all 'wrong' arguments; suppress the warnings though
      snn <- suppressMessages(suppressWarnings(do.call(Seurat::FindNeighbors, args = c(list(object = expr.select), temp_dots))))
      print(paste0("End: ", Sys.time()))
    },
    error = function(e) {
      print("louvain / leiden: error in snn calculation.")
      run.louvain <- F
      run.leiden <- F
    }
  )

  tryCatch(
    if (run.louvain) {
      temp_dots <- dots[which(grepl("^louvain__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^louvain__", "", names(temp_dots), ignore.case = T)
      print(paste0("Finding clusters with Seurats implementation of the Louvain algorithm and parallel::mclapply using ", mc.cores," cores. Start: ",Sys.time()))
      if (any(grepl("resolution", names(temp_dots), ignore.case = T))) {
        temp_dots[["resolution"]] <- as.numeric(temp_dots[["resolution"]])
        if (any(is.na(temp_dots[["resolution"]]))) {
          warning("Provide numeric values for louvain__resolution! Non numeric elements are ignored.")
          temp_dots[["resolution"]] <-temp_dots[["resolution"]][which(!is.na(temp_dots[["resolution"]]))]
        }
        clust_idents <- do.call(cbind,parallel::mclapply(temp_dots[["resolution"]], function(x) {
          apply(do.call(Seurat::FindClusters, args = c(list(object = snn$snn, resolution = x, verbose = F, algorithm = 1), temp_dots[which(names(temp_dots) != "resolution")])), 2, as.numeric)
        }, mc.cores = mc.cores))
      } else {
        clust_idents <- apply(do.call(Seurat::FindClusters, args = c(list(object = snn$snn, resolution = x, verbose = F, algorithm = 1), temp_dots)), 2, as.numeric)
      }
      colnames(clust_idents) <- paste0("louvain_", temp_dots[["resolution"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, clust_idents))
      print(apply(clust_idents, 2, function(x) length(unique(x))))
      print(paste0("End: ", Sys.time()))
    },
    error = function(e) {
      print("run.louvain with error")
    }
  )

  tryCatch(
    if (run.leiden) {
      temp_dots <- dots[which(grepl("^leiden__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^leiden__", "", names(temp_dots), ignore.case = T)

      print(paste0("Finding clusters with leiden algorithm and parallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time()))

      if (any(grepl("resolution_parameter", names(temp_dots), ignore.case = T))) {
        temp_dots[["resolution_parameter"]] <- as.numeric(temp_dots[["resolution_parameter"]])
        if (any(is.na(temp_dots[["resolution_parameter"]]))) {
          warning("Provide numeric values for leiden__resolution! Non numeric elements are ignored.")
          temp_dots[["resolution_parameter"]] <- temp_dots[["resolution_parameter"]][which(!is.na(temp_dots[["resolution_parameter"]]))]
        }
        clust_idents <- do.call(cbind, parallel::mclapply(temp_dots[["resolution_parameter"]], function(x) {
          do.call(leiden::leiden, args = c( list(object = snn$snn, resolution_parameter = x), temp_dots[which(names(temp_dots) != "resolution_parameter")]))
        }, mc.cores = mc.cores))
      } else {
        clust_idents <- do.call(leiden::leiden, args = c(list(object = snn$snn), temp_dots))
      }
      colnames(clust_idents) <- paste0("leiden_", temp_dots[["resolution_parameter"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, clust_idents))
      print(apply(clust_idents, 2, function(x) length(unique(x))))
      print(paste0("End: ", Sys.time()))
    },
    error = function(e) {
      print("run.leiden with error")
    }
  )

  tryCatch(
    if (run.kmeans) {
      temp_dots <- dots[which(grepl("^kmeans__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^kmeans__", "", names(temp_dots), ignore.case = T)
      print(paste0("Finding clusters with kmeans andparallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time()))

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["centers"]], function(x) {
        do.call(stats::kmeans, args = c(list(x = expr.select, centers = x), temp_dots[which(names(temp_dots) != "centers")]))$cluster
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("kmeans_", temp_dots[["centers"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      print(paste0("End: ", Sys.time()))
    },
    error = function(e) {
      print("run.kmeans with error")
    }
  )

  tryCatch(
    if (run.hclust) {
      temp_dots <- dots[which(grepl("^hclust__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^hclust__", "", names(temp_dots), ignore.case = T)
      print(paste0("Finding clusters with hclust. Start: ", Sys.time()))

      d <- do.call(stats::dist, args = c(list(x = expr.select), temp_dots[which(names(temp_dots) %in% names(formals(stats::dist))[-1])]))
      h <- do.call(stats::hclust, args = c(list(d = d), temp_dots[which(names(temp_dots) %in% names(formals(stats::hclust))[-1])]))
      ks <- do.call(cbind, parallel::mclapply(temp_dots[["k"]], function(x) {
        stats::cutree(tree = h, k = x)
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("hclust_", temp_dots[["k"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      print(paste0("End: ", Sys.time()))
    },
    error = function(e) {
      print("run.hclust with error")
    }
  )

  tryCatch(
    if (run.flowClust) {

      temp_dots <- dots[which(grepl("^flowClust__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^flowClust__", "", names(temp_dots), ignore.case = T)

      print(paste0("Finding clusters with flowClust andparallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time()))
      ks <- do.call(cbind, parallel::mclapply(temp_dots[["K"]], function(x) {
        suppressMessages(do.call(flowClust::flowClust, args = c(list(x = expr.select, K = x), temp_dots[which(names(temp_dots) != "K")]))@label)
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("flowClust_", temp_dots[["K"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      print(paste0("End: ", Sys.time()))
    },
    error = function(e) {
      print("run.flowClust with error")
    }
  )

  tryCatch(
    if (run.MUDAN) {

      temp_dots <- dots[which(grepl("^MUDAN__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^MUDAN__", "", names(temp_dots), ignore.case = T)

      print(paste0("Finding clusters with MUDAN.", Sys.time()))

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["k"]], function(x) {
        ## documentation is wrong (mat: cells as rows and features as cols!)
        do.call(MUDAN::getComMembership, args = c(list(mat = expr.select, k = x), temp_dots[which(names(temp_dots) != "k")]))
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("MUDAN_", temp_dots[["k"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      print(paste0("End: ", Sys.time()))
    },
    error = function(e) {
      print("run.MUDAN with error")
    }
  )

  ### optionally run MUDAN::clusterBasedBatchCorrect here
  ## actually though harmony performs something similar with multiple iterations: https://portals.broadinstitute.org/harmony/articles/quickstart.html
  ## MUDAN advertises harmony: http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/mudan.html

  if (!is.logical(run.lda)) {
    # optimized embedding with LDA, so embedding biased by detected clusters. is that cheating or just the best dimensions to separate clusters on a map?!
    # https://highdemandskills.com/lda-clustering/
    ldam <- MASS::lda(com ~ ., data = as.data.frame(cbind(expr.select, com = dim.red.data[, run.lda])))
    expr.select <- stats::predict(ldam, as.data.frame(cbind(expr.select, com = dim.red.data[, run.lda])))$x
  }


  if (run.umap) {
    temp_dots <- dots[which(grepl("^UMAP__", names(dots), ignore.case = T))]
    names(temp_dots) <-gsub("^UMAP__", "", names(temp_dots), ignore.case = T)

    if (!any(grepl("metric", names(temp_dots)), ignore.case = T)) {
      print("UMAP metric set to 'cosine' by default.")
      temp_dots <- c(temp_dots, metric = "cosine")
    }
    print(paste0("Calculating UMAP. Start: ", Sys.time()))
    if (any(grepl("n_neighbors", names(temp_dots), ignore.case = T))) {
      umap.dims <- do.call(cbind,parallel::mclapply(temp_dots[["n_neighbors"]], function(z) {
        out <- do.call(uwot::umap, args = c(list(X = expr.select, verbose = F, n_neighbors = z),temp_dots[which(names(temp_dots) != "n_neighbors")]))
        colnames(out) <- c(paste0("UMAP_1_", z), paste0("UMAP_2_", z))
        return(out)
      }, mc.cores = mc.cores))
    } else {
      umap.dims <- do.call(uwot::umap, args = c(list(X = expr.select), temp_dots))
    }
    if (!any(grepl("n_neighbors", names(temp_dots), ignore.case = T)) || length(temp_dots[["n_neighbors"]]) == 1) {
      colnames(umap.dims) <- c("UMAP_1", "UMAP_2")
    }
    print(paste0("End: ", Sys.time()))
  }

  if (run.tsne) {
    temp_dots <- dots[which(grepl("^tSNE__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^tSNE__", "", names(temp_dots), ignore.case = T)

    if (run.pca) {
      print("As 'run.pca=T' inital pca in tSNE is not performed.")
      temp_dots <- c(temp_dots, pca = F)
    }
    if (!any(grepl("normalize", names(temp_dots), ignore.case = T))) {
      print("Rtsne normalize set to 'FALSE' by default. (scale.whole, scale.samples?!)")
      temp_dots <- c(temp_dots, normalize = F)
    }

    print(paste0("Calculating tSNE. Start: ", Sys.time()))
    if (any(grepl("perplexity", names(temp_dots), ignore.case = T))) {
      tsne.dims <- do.call(cbind, parallel::mclapply(temp_dots[["perplexity"]], function(z) {
        out <- do.call(Rtsne::Rtsne, args = c(list(X = expr.select, verbose = F, perplexity = z), temp_dots[which(names(temp_dots) != "perplexity")]))$Y
        colnames(out) <- c(paste0("tSNE_1_", z), paste0("tSNE_2_", z))
        return(out)
      }, mc.cores = mc.cores))
    } else {
      tsne.dims <- do.call(Rtsne::Rtsne, args = c(list(X = expr.select, verbose = T), temp_dots))$Y
    }
    if (!any(grepl("perplexity", names(temp_dots), ignore.case = T)) || length(temp_dots[["perplexity"]]) == 1) {
      colnames(tsne.dims) <- c("tSNE_1", "tSNE_2")
    }
    print(paste0("End: ", Sys.time()))
  }


  if (run.som) {
    temp_dots <- dots[which(grepl("^SOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^SOM__", "", names(temp_dots), ignore.case = T)
    map <-do.call(EmbedSOM::SOM, args = c(list(data = expr.select), temp_dots))

    temp_dots <- dots[which(grepl("^EmbedSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^EmbedSOM__", "", names(temp_dots), ignore.case = T)
    som.dims <- do.call(EmbedSOM::EmbedSOM, args = c(list(data = expr.select, map = map), temp_dots))
    colnames(som.dims) <- c("SOM_1", "SOM_2")
  }

  if (run.gqtsom) {
    temp_dots <- dots[which(grepl("^GQTSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^GQTSOM__", "", names(temp_dots), ignore.case = T)
    map <- do.call(EmbedSOM::GQTSOM, args = c(list(data = expr.select), temp_dots))

    temp_dots <- dots[which(grepl("^EmbedSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^EmbedSOM__", "", names(temp_dots), ignore.case = T)
    gqtsom.dims <- do.call(EmbedSOM::EmbedSOM, args = c(list(data = expr.select, map = map), temp_dots))
    colnames(gqtsom.dims) <- c("GQTSOM_1", "GQTSOM_2")
  }

  if (write.scaled.channels.to.FCS) {
    scaled.expr <- scale.whole(do.call(rbind, lapply(ff.list[["logicle"]], function(x) {
      scale.samples(flowCore::exprs(x)[, channels])
    })))
    scaled.expr <- scaled.expr[, which(!grepl(exclude.extra.channels, colnames(scaled.expr)))]
    colnames(scaled.expr) <- paste0(colnames(scaled.expr), "_scaled")
    dim.red.data <- do.call(cbind, list(dim.red.data, scaled.expr))
  }

  if (run.pca) {
    dim.red.data <- do.call(cbind, list(dim.red.data, pca.matrix))
  }
  if (run.umap) {
    dim.red.data <- do.call(cbind, list(dim.red.data, umap.dims))
  }
  if (run.tsne) {
    dim.red.data <- do.call(cbind, list(dim.red.data, tsne.dims))
  }
  if (run.som) {
    dim.red.data <- do.call(cbind, list(dim.red.data, som.dims))
  }
  if (run.gqtsom) {
    dim.red.data <- do.call(cbind, list(dim.red.data, gqtsom.dims))
  }

  if (!missing(extra.cols)) {
    if (nrow(extra.cols) == nrow(dim.red.data)) {
      dim.red.data <- cbind(dim.red.data, extra.cols)
    } else {
      print("extra.cols not added due to wrong row number.")
    }
  }

  # 2021 06 17 necessary
  dim.red.data <- as.data.frame(dim.red.data)

  tryCatch(
    if (!missing(add.sample.info)) {
      for (i in names(add.sample.info)) {
        dim.red.data <- do.call(cbind, list(dim.red.data, rep(add.sample.info[[i]], times = as.numeric(table(rep(1:length(ff.list[["logicle"]]), sapply(ff.list[["logicle"]], nrow)))))))
        names(dim.red.data)[length(dim.red.data)] <- i
      }
    },
    error = function(e) {
      print("Addition of sample info failed due to an error. Check add.sample.info.")
    }
  )

  # prepare channel desc
  name.desc <- setNames(ff.list[[1]][[1]]@parameters@data[["desc"]], ff.list[[1]][[1]]@parameters@data[["name"]])
  name.desc <- name.desc[which(!is.na(name.desc))]
  channel.desc <- rep("", ncol(dim.red.data))
  for (i in seq_along(name.desc)) {
    channel.desc[grep(names(name.desc)[i], colnames(dim.red.data))] <-
      name.desc[i]
  }

  channel.desc_augment <- channel.desc
  channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("scaled", colnames(dim.red.data))))] <-
    paste0(channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("scaled", colnames(dim.red.data))))], "_scaled")
  channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("logicle", colnames(dim.red.data))))] <-
    paste0(channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("logicle", colnames(dim.red.data))))], "_logicle")
  channel.desc_augment[which(channel.desc_augment == "")] <-
    colnames(dim.red.data)[which(channel.desc_augment == "")]
  channel.desc_augment <- make.names(channel.desc_augment)
  names(channel.desc_augment) <- colnames(dim.red.data)

  # write FCS file
  new_p <- flowCore::parameters(ff.list[[1]][[1]])[1, ]
  new_kw <- flowCore::keyword(ff.list[[1]][[1]])
  new_pars <- flowCore::parameters(ff.list[[1]][[1]])

  for (z in (nrow(new_pars) + 1):ncol(dim.red.data)) {
    new_p_number <- z
    rownames(new_p) <- c(paste0("$P", new_p_number))
    new_pars <- BiocGenerics::combine(new_pars, new_p)

    new_p_name <- names(dim.red.data)[z]
    new_p_desc <- channel.desc[z]
    flowCore::pData(new_pars)$name[new_p_number] <- new_p_name
    flowCore::pData(new_pars)$desc[new_p_number] <- new_p_desc

    new_kw["$PAR"] <- as.character(new_p_number)
    new_kw[paste0("$P", as.character(new_p_number), "N")] <- new_p_name
    new_kw[paste0("$P", as.character(new_p_number), "S")] <- new_p_desc
    new_kw[paste0("$P", as.character(new_p_number), "E")] <- "0,0"
    new_kw[paste0("$P", as.character(new_p_number), "G")] <- "1"
    new_kw[paste0("$P", as.character(new_p_number), "B")] <- new_kw["$P1B"]
    new_kw[paste0("$P", as.character(new_p_number), "R")] <- max(dim.red.data[, z])
    new_kw[paste0("$P", as.character(new_p_number), "DISPLAY")] <- "LIN"
    new_kw[paste0("flowCore_$P", as.character(new_p_number), "Rmin")] <- min(dim.red.data[, z])
    new_kw[paste0("flowCore_$P", as.character(new_p_number), "Rmax")] <- max(dim.red.data[, z])
  }
  ff <- methods::new("flowFrame", exprs = as.matrix(dim.red.data), parameters = new_pars, description = new_kw)

  # get cluster markers
  ## always used logicle transformed data?!?!
  marker <- lapply(calc.cluster.markers, function (clust_col) {
    # do not use expr.select which may have become dimenions of LDA
    dat <- dim.red.data[,c(which(colnames(dim.red.data) %in% paste0(channels, "_logicle")), which(colnames(dim.red.data) == clust_col))]
    split_var <- dat[,clust_col]
    split_var_levels <- sort(unique(split_var))
    ## keep dat a data frame until here to allow split (works only on data.frame); after that convert to matrix
    dat_split <- split(dat, split_var)
    dat <- as.matrix(dat)
    dat_split <- lapply(dat_split, function(x) as.matrix(x[,-which(names(x) == clust_col)]))
    all_pairs <- combn(split_var_levels, 2, simplify = F)

    ## all pairwise
    pair_marker_table <- dplyr::bind_rows(lapply(all_pairs, function(x) {
      out <- matrixTests::col_wilcoxon_twosample(dat_split[[as.character(x[1])]], dat_split[[as.character(x[2])]])
      out[,"mean_1"] <- matrixStats::colMeans2(dat_split[[as.character(x[1])]])
      out[,"mean_2"] <- matrixStats::colMeans2(dat_split[[as.character(x[2])]])
      out[,"mean_diff"] <- out[,"mean_1"] - out[,"mean_2"]
      out[,"diptest_p_1"] <- apply(dat_split[[as.character(x[1])]], 2, function(x) diptest::dip.test(x)[["p.value"]])
      out[,"diptest_p_2"] <- apply(dat_split[[as.character(x[2])]], 2, function(x) diptest::dip.test(x)[["p.value"]])
      out <- tibble::rownames_to_column(out, "channel")
      out[,"cluster_1"] <- x[1]
      out[,"cluster_2"] <- x[2]
      out[,"diff_sign"] <- ifelse(out[,"mean_diff"] == 0, "+/-", ifelse(out[,"mean_diff"] > 0, "+", "-"))
      out <- dplyr::select(out, channel, cluster_1, cluster_2, pvalue, mean_1, mean_2, mean_diff, diff_sign, diptest_p_1, diptest_p_2)
      out <- dplyr::arrange(out, pvalue)
      return(out)
    }))
    pair_marker_table[,"channel_desc"] <- channel.desc_augment[pair_marker_table[,"channel"]]

    marker_table <- dplyr::bind_rows(lapply(split_var_levels, function(x) {
      y <- dat[which(dat[,clust_col] == x),which(colnames(dat) != clust_col)]
      z <- dat[which(dat[,clust_col] != x),which(colnames(dat) != clust_col)]
      out <- matrixTests::col_wilcoxon_twosample(y, z)
      out[,"mean"] <- round(matrixStats::colMeans2(y), 2)
      out[,"mean_not"] <- round(matrixStats::colMeans2(z), 2)
      out[,"mean_diff"] <- round(out[,"mean"] - out[,"mean_not"], 2)
      out[,"diptest_p"] <- round(apply(y, 2, function(x) diptest::dip.test(x)[["p.value"]]), 2)
      out[,"diptest_not_p"] <- round(apply(z, 2, function(x) diptest::dip.test(x)[["p.value"]]), 2)
      out <- tibble::rownames_to_column(out, "channel")
      out[,"cluster"] <- x
      out[,"diff_sign"] <- ifelse(out[,"mean_diff"] == 0, "+/-", ifelse(out[,"mean_diff"] > 0, "+", "-"))
      out <- dplyr::select(out, channel, cluster, pvalue, mean, mean_not, mean_diff, diff_sign, diptest_p, diptest_not_p)
      out <- dplyr::arrange(out, pvalue)
      return(out)
    }))
    marker_table[,"channel_desc"] <- channel.desc_augment[marker_table[,"channel"]]
    return(list(marker_table = marker_table, pairwise_marker_table = pair_marker_table))
  })
  names(marker) <- calc.cluster.markers

  ## find out which contrasts are marked by a channel
  #xxx <- dplyr::filter(pair_marker_table, mean_diff > 0) %>% dplyr::group_by(channel, cluster) %>% dplyr::slice_min(pvalue, n = 3)
  ## group by cluster to find cluster markers
  #yyy <- dplyr::filter(marker_table, mean_diff > 0) %>% dplyr::group_by(cluster) %>% dplyr::slice_min(pvalue, n = 3)
  ## group by channel to find out which clusters they indicate
  #zzz <- dplyr::filter(marker_table, mean_diff > 0) %>% dplyr::group_by(channel) %>% dplyr::slice_min(pvalue, n = 3)


  # save results
  if (!is.null(save.path) && !is.na(save.path)) {
    print("Writing files to disk.")
    t <- format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%Y%m%d_%H%M%S")
    dir.create(save.path, showWarnings = F)
    if ("rds" %in% save.to.disk) {
      saveRDS(list(table = dim.red.data, flowframe = ff, pca = pca.result), file.path(save.path, paste0(t, "_dr_ff_list.rds")))
    }
    if ("fcs" %in% save.to.disk) {
      flowCore::write.FCS(ff, file.path(save.path, paste0(t, "_dr.fcs")))
    }
  }
  return(list(df = dim.red.data, col_names = channel.desc_augment, marker = marker))
}

#' Calculate dimension reductions and cluster annotations with data from one or more flow frames and add these parameters to a (concatenated) fcs file
#'
#' Prepare a list of flowframes (ff.list) with fcexpr::wsp_get_ff or fcexpr::inds_get_ff. The objects returned from these function will compatible with fcexpr::dr_to_fcs.
#' The list must contain logicle transformed fluorescence intensities (FI) and optionally may contain corresponding inverse transformed FI (which is the transformation you see in flowjo).
#' For dimension reduction logicle transformed FI are used. Currently, this is obligatory. This transformation, optionally additional scaling operations alogn the way (scale.whole, scale.samples)
#' as well as cluster annotation will be written to the resulting (concatenated) fcs file for manual inspection in flowjo. Certain dimension reduction algorithms and cluster detection algorithm
#' may become slow with a large number of events (e.g. > 1e6, see the details). In order to get a quick impression of what the algorithms can pull out for you, you may use the 'downsample'
#' argument in fcexpr::wsp_get_ff or fcexpr::inds_get_ff to conveniently sample a random subset of events from each fcs files (flowframe).
#'
#' Logicle transformation of FI has been described here: \href{https://pubmed.ncbi.nlm.nih.gov/16604519/}{Parks, et al., 2006, PMID: 16604519  DOI: 10.1002/cyto.a.20258}.
#' Different transformations have been compared for instance here: \href{https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR-Transformations.html}{Transformations}.
#' Dimension reduction (low dimension embedding) algorithms are: UMAP, tSNE, SOM, GQTSOM, PCA. tSNE is expected to be too slow for more then 1e4 or 1e5 cells (depending
#' precision set by tsne__theta). UMAP has been developed later and is well-know. SOM and GQTSOM are similar to \href{https://github.com/SofieVG/FlowSOM}{flowSOM}
#' but are much more convenient to use in programming as they accept matrices whereas flowSOM strictly requires flowFrames. Also SOM and GQTSOM may be superior with
#' respect to calculation speed. \href{"https://www.youtube.com/watch?v=FgakZw6K1QQ"}{PCA} will be orders of magnitude faster than UMAP especially on large numbers of events (1e6 and above).
#' On the other hand it will not produce as nicely separated clusters as it is a linear algorithm.
#' Cluster (community) detection algorithms are louvain, leiden, kmeans, hclust, flowClust and MUDAN. These assign events with similar properties in the high dimensional space a common discrete label. kmeans will be the quickest.
#' Start with this one and progress to using louvain which is slower but may yield better results. For a low number of events (e.g. below 1e6) louvain will perform in a reasonable amount of time
#' and could used immediately in parallel to kmeans. leiden will require the original python library and is approximately similar to louvain with respect to calculation time,
#' maybe a bit slower. leiden is an enhancement of louvain, though louvain does not produce "wrong" results per se. flowClust, MUDAN, and hclust have not been tested extensively.
#' For kmeans, hclust, flowClust the number of expected cluster can be provided. louvain and leiden take a resolution parameter (lower resolution --> less clusters). MUDAN takes
#' the a number of nearest neighbors (less neighbors --> more clusters).
#' An initial PCA may not be required if the selected channels are chosen thoughtfully (only those with stained markers). If you decide though to simply
#' use all channels (with and and without stained marker) for dimension reduction, PCA will automatically lead to focus on the ones with highest variation (so those channels which contain staining
#' on a subset of events). In this case you may set n.pca.dims roughly equal to the number of channels which contain a staining.
#'
#' @param ff.list a list of flowFrames as received from fcexpr::wsp_get_ff (compensated with Compensation Matrix as defined in FlowJo by default) or
#' as received from fcexpr::inds_get_ff (directly from FCS files, not compensated by default)
#' @param channels a named vector of channels to use for dimension reduction. values can be channel names (v-450LP..., b-520..., or so), or channel descriptions (e.g. CD3 or LAG3-PeCy7 for example)
#' names will be used as new description in the fcs file to be written; if no names provided, names of the very first FCS file will be used
#' @param add.sample.info named list of additional channels to identify samples or group them in flowjo;
#' e.g. for 9 fcs files to be concatenated: add.sample.info = list(condition = c(1,2,3,1,2,3,1,2,3), donor = c(1,1,1,2,2,2,3,3,3))
#' @param scale.whole if and how to scale channels after concatenation of flowframes in ff.list
#' @param scale.samples if and how to scale channels of flowframes in ff.list individually before concatenation
#' @param run.harmony attempt batch correction using harmony::HarmonyMatrix; if TRUE, then harmony__meta_data has to be provided in ... indicating the groups to be batch corrected;
#' harmony is conducted before run.pca; to calculate a pca before harmony, pass harmony__do_pca = T and optional a number of pcs with harmony__npcs. Set run.pca = F
#' when a pca is before in harmony.
#' @param run.pca run principle component analysis only, or prior to other dimension reduction (umap, SOM, tSNE, ...)
#' @param n.pca.dims number of principle components to calculate
#' @param run.lda run Linear Discriminant Analysis before dimension reduction;
#' should be F (FALSE) or a clustering calculated before, e.g. louvain_0.8 or leiden_1.1, kmeans_12 etc.; respective clustering calculation
#' has to be provided in ...
#' @param run.umap calculate UMAP dimension reduction with uwot::umap
#' @param run.tsne calculate tsne dimension reduction with Rtsne::Rtsne
#' @param run.som calculate SOM dimension reduction EmbedSOM::SOM followed by EmbedSOM::EmbedSOM
#' @param run.gqtsom calculate GQTSOM dimension reduction EmbedSOM::GQTSOM followed by EmbedSOM::EmbedSOM
#' @param run.louvain detect clusters (communities, groups) of cells with the louvain algorithm, implemented in Seurat::FindClusters (subsequent to snn detection by Seurat::FindNeighbors)
#' @param run.leiden detect clusters (communities, groups) of cells with the leiden algorithm, with leiden::leiden (subsequent to snn detection by Seurat::FindNeighbors)
#' @param run.kmeans detect clusters with stats::kmeans
#' @param run.minibatchkmeans detect clusters with ClusterR::MiniBatchKmeans
#' @param run.kmeans_arma detect clusters with ClusterR::KMeans_arma
#' @param run.kmeans_rcpp detect clusters with ClusterR::KMeans_rcpp
#' @param run.hclust detect clusters with stats::dist, stats::hclust and stats::cutree
#' @param run.flowClust detect clusters with flowClust::flowClust
#' @param run.MUDAN detect clusters with MUDAN::getComMembership
#' @param extra.cols vector of one extra column (or matrix of multiple columns) to add to the final fcs file;
#' has to be numeric; has to be equal to the number of rows in all flowframes provided; colnames will be the channel names in the FCS file;
#' could be a previously calculated dimension reduction or cluster annotation.
#' @param calc.cluster.markers if NULL nothing is calculated; otherwise marker features (stained markers) are determined by wilcox test
#' using \href{https://github.com/immunogenomics/presto}{presto::wilcoxauc} for the provided clustering(s). each cluster
#' is tested against other events and clusters are compaired pairwise. respective clustering calculation has to be provided in ...;
#' e.g. if louvain__resolution = 0.5 is provided set calc.cluster.markers = louvain_0.5;
#' and if in addition leiden__resolution_parameter = 0.7 then set calc.cluster.markers = c(louvain_0.5, leiden_0.7).
#' @param mc.cores mc.cores to calculate clusterings, limited to parallel::detectCores()-1
#' @param save.to.disk what to save to disk: (concatenated) and appended FCS file and/or rds file with several elements in a list
#' @param save.path where to save elements specified in save.to.disk; set to NULL to have nothing written to disk
#' @param exclude.extra.channels when scaled and transform channels are written to FCS file, some channels may be redundant
#' and will only occupy disk space, those are specified here; matched with grepl
#' @param write.scaled.channels.to.FCS do save scaled channels (scale.whole, scale.samples) to FCS file
#' @param timeChannel name of the Time channel to exclude from all analyses and calculation; if NULL will be attempted
#' to be detected automatically
#' @param ... additional parameters to calculations of \href{https://github.com/jlmelville/uwot}{UMAP}, \href{https://github.com/jkrijthe/Rtsne}{tSNE}, \href{https://github.com/exaexa/EmbedSOM}{SOM, GQTSOM, EmbedSOM}, louvain, \href{https://github.com/TomKellyGenetics/leiden}{leiden},
#' \href{https://github.com/immunogenomics/harmony}{harmony}, \href{https://www.bioconductor.org/packages/release/bioc/html/flowClust.html}{flowClust}, hclust, \href{https://github.com/JEFworks/MUDAN}{MUDAN}, kmeans;
#' provide arguments as follows: UMAP__n_neighbors = c(15,20,25), or tsne__theta = 0.3, etc.
#' see respected help files to get to know which arguments can be passed:
#' uwot::umap, Rtsne::Rtsne, EmbedSOM::SOM, EmbedSOM::GQTSOM, EmbedSOM::EmbedSOM, harmony::HarmonyMatrix, flowClust::flowClust,
#' louvain: Seurat::FindNeighbors and Seurat::FindCluster, leiden: Seurat::FindNeighbors and leiden::leiden.
#' hclust: stats::dist and stats::hclust, MUDAN: MUDAN::getComMembership, stats::kmeans,
#' ClusterR::MiniBatchKmeans, ClusterR::KMeans_rcpp, ClusterR::KMeans_arma
#'
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
#'
#'  ##############################################
#'  ### find clusters which may be subject #######
#'  ### to bi- or multimodality in any channel ###
#'  ##############################################
#'  # make one data frame
#'  marker_all <- dplyr::bind_rows(lapply(names(dr[["marker"]]), function(x) dplyr::mutate(dr[["marker"]][[x]][["marker_table"]], clustering = x)))
#'
#'  # sort by diptest p value; or low p indicates bi- or multimodality
#'  marker_all <- dplyr::arrange(marker_all, diptest_p)
#'  # see ?diptest::dip.test
#'  # multimodality indicates heterogeneity within in the cluster
#'  # and may justify to separate that cluster further into sub-clusters
#'  # this depends on the interpretation of the scientist though
#'
#'
#'}
dr_to_fcs <- function(ff.list,
                      channels = NULL,
                      add.sample.info = NULL,
                      scale.whole = c("z.score", "min.max", "none"),
                      scale.samples = c("none", "z.score", "min.max"),
                      run.harmony = F,
                      run.pca = F,
                      run.lda = F,
                      run.umap = T,
                      run.tsne = F,
                      run.som = T,
                      run.gqtsom = F,
                      run.louvain = F,
                      run.kmeans = F,
                      run.minibatchkmeans = F,
                      run.kmeans_arma = F,
                      run.kmeans_rcpp = F,
                      run.leiden = F,
                      run.hclust = F,
                      run.flowClust = F,
                      run.MUDAN = F,
                      n.pca.dims = 0,
                      calc.cluster.markers = NULL,
                      extra.cols = NULL,
                      mc.cores = 1,
                      save.to.disk = c("fcs", "rds"),
                      save.path = file.path(getwd(), paste0(substr(gsub("\\.", "", make.names(as.character(Sys.time()))), 2, 15), "_FCS_dr")),
                      exclude.extra.channels = ifelse(length(ff.list) == 1 && names(ff.list) == "logicle", "cluster.id", "FSC|SSC|Time|cluster.id"),
                      write.scaled.channels.to.FCS = F,
                      timeChannel = "Time",
                      ...) {
  # batch effect correction: https://cytekbio.com/blogs/blog/how-to-identify-and-prevent-batch-effects-in-longitudinal-flow-cytometry-research-studies
  # cytonorm (https://github.com/saeyslab/CytoNorm) requires reference sample for every batch - not always available
  # check somewhen: https://github.com/casanova-lab/iMUBAC
  # harmony: https://github.com/immunogenomics/harmony
  # MUDAN: https://github.com/JEFworks/MUDAN

  # harmony: return whole object??

  # MUDAN::getComMembership offers to pass one of many algorithms (seen here: https://slowkow.com/notes/harmony-animation/), though how to pass individual arguments? e.g. for igraph::cluster_leiden

  # optionally add: MUDAN::clusterBasedBatchCorrect
  ## allow to provide expr.select directly instead of ff.list
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

  if (!requireNamespace("diptest", quietly = T)) {
    utils::install.packages("diptest")
  }
  if (!requireNamespace("matrixStats", quietly = T)) {
    utils::install.packages("matrixStats")
  }
  if (!requireNamespace("Rtsne", quietly = T)) {
    utils::install.packages("Rtsne")
  }
  if (run.umap && !requireNamespace("uwot", quietly = T)) {
    utils::install.packages("uwot")
  }
  if ((run.leiden || run.louvain) && !requireNamespace("Seurat", quietly = T)) {
    utils::install.packages("Seurat")
  }
  if (run.leiden &&!requireNamespace("leiden", quietly = T)) {
    utils::install.packages("leiden")
  }
  if (run.flowClust && !requireNamespace("flowClust", quietly = T)) {
    BiocManager::install("flowClust")
  }
  if (!requireNamespace("devtools", quietly = T)) {
    utils::install.packages("devtools")
  }
  if (!requireNamespace("presto", quietly = T)) {
    devtools::install_github('immunogenomics/presto')
  }
  if (run.harmony && !requireNamespace("harmony", quietly = T)) {
    devtools::install_github("immunogenomics/harmony")
  }
  if (run.MUDAN && !requireNamespace("MUDAN", quietly = T)) {
    devtools::install_github("JEFworks/MUDAN")
  }
  if (run.som && !requireNamespace("EmbedSOM", quietly = T)) {
    devtools::install_github("devtools::install_github('exaexa/EmbedSOM')")
  }
  if ((run.minibatchkmeans || run.kmeans_rcpp || run.kmeans_arma) && !requireNamespace("ClusterR", quietly = T)) {
    utils::install.packages("ClusterR")
  }

  dots <- list(...)

  expect_dots <- "^harmony__|^hclust__|^flowClust_|^MUDAN__|^kmeans__|^louvain__|^leiden__|^som__|^gqtsom__|^tsne__|^umap__|^EmbedSOM|^kmeans_arma__|^kmeans_rcpp__|^minibatchkmeans__"
  if (any(!names(dots) %in% names(formals(dr_to_fcs)) & !grepl(expect_dots, names(dots)))) {
    message("These arguments are unknown: ", paste(names(dots)[which(!names(dots) %in% names(formals(dr_to_fcs)) & !grepl(expect_dots, names(dots)))], collapse = ", "))
  }


  if (length(dots) > 0) {
    dots_expanded <- unname(unlist(mapply(paste, sapply(strsplit(names(dots), "__"), "[", 1), dots, sep = "_")))
  } else {
    dots_expanded <- NULL
  }

  if (!"logicle" %in% names(ff.list)) {
    stop("'logicle' not found in ff.list: one list of flowframes in ff.list has to be named 'logicle' as has to contain logicle transformed fluorescence intensities.")
  }

  if (any(!names(ff.list) %in% c("inverse", "logicle"))) {
    stop("ff.list has to contain a list of flowframes named 'logicle' (logicle transformed) and optionally an additional list named 'inverse' (inverse transformed, original as in flowjo.).")
  }

  # check if names in ff.list inverse and logicle are the same
  if (length(unique(lengths(ff.list))) != 1) {
    stop("number of flowframes in inverse and logicle has to be equal.")
  }

  for (par in c("louvain", "leiden", "umap", "tsne", "som", "gqtsom", "harmony", "kmeans", "kmeans_rcpp", "kmeans_arma", "minibatchkmeans", "flowclust", "hclust", "harmony", "mudan")) {
    if (any(grepl(paste0("^", par, "__"), names(dots), ignore.case = T)) &&!eval(rlang::sym(paste0("run.", par)))) {
      message(par, " parameters provided in ... but ", "'run.", par, " = F'.")
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

  if (run.minibatchkmeans && !any(grepl("^minibatchkmeans__clusters", names(dots)))) {
    stop("When 'run.minibatchkmeans = T' minibatchkmeans__clusters, has to be provided in ..., see ?ClusterR::MiniBatchKmeans")
  }

  if (run.kmeans_rcpp && !any(grepl("^kmeans_rcpp__clusters", names(dots)))) {
    stop("When 'run.kmeans_rcpp = T' kmeans_rcpp__clusters, has to be provided in ..., see ?ClusterR::KMeans_rcpp")
  }

  if (run.kmeans_arma && !any(grepl("^kmeans_arma__clusters", names(dots)))) {
    stop("When 'run.kmeans_arma = T' kmeans_arma__clusters, has to be provided in ..., see ?ClusterR::KMeans_arma")
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

  if (is.logical(run.lda) && run.lda) {
    stop("Do not set 'run.lda = T' but provide a clustering that should be used to calculate it, e.g. a pattern like louvain_0.4.")
  }

  if (!is.logical(run.lda)) {
    # clustering columns need to follow the pattern name_resolution.
    if (!run.lda %in% dots_expanded) {
      stop("Value for run.lda not found in ... . They respective clustering to use for lda has to be computed! E.g. if run.lda = 'louvain_0.9' then louvain__resolution = 0.9 has to be passed.")
    }
    if (length(run.lda) > 1) {
      stop("run.lda must be of length 1. So, it should only indicate one clustering, e.g. kmeans_9.")
    }
  }

  if (run.harmony && any(grepl("^harmony__do_pca", names(dots))) && dots[["harmony__do_pca"]] == T && run.pca) {
    warning("harmony is calculated with do_pca = T, hence a subsequnt pca (run.pca = T) is not required. Consider setting run.pca to FALSE.")
  }

  if (!is.null(calc.cluster.markers)) {
    if (!any(calc.cluster.markers %in% dots_expanded)) {
      stop("calc.cluster.markers: ", calc.cluster.markers[which(!calc.cluster.markers %in% dots_expanded)], " not found in ... .")
    }
  }

  if (!is.null(save.to.disk)) {
    save.to.disk <- match.arg(save.to.disk, c("fcs", "rds"), several.ok = T)
  }

  mc.cores <- min(mc.cores, parallel::detectCores() - 1)

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
      stop("add.sample.info has to have names. These names become channel names in the FCS file.")
    }
    if (!all(unlist(lapply(add.sample.info, function(x) is.numeric(x))))) {
      stop("Please only provide numeric values as additional sample infos.")
    }
    if (any(unlist(lapply(add.sample.info, function(x) is.na(x))))) {
      stop("NA found in sample infos.")
    }
    if (!all(unlist(lapply(add.sample.info, function(x) length(x) == length(ff.list[[1]]))))) {
      stop("Length of each additional sample information has to match the length of selected samples, which is: ", length(ff.list[[1]]),".")
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
      message("'do_pca set to FALSE' in harmony::HarmonyMatrix.")
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
    message("Calculating PCA. Start: ", Sys.time())
    # https://slowkow.com/notes/pca-benchmark/
    ## version 1
    # mat_irlba2 <- irlba::irlba(A = expr.select, nv = n.pca.dims)
    # mat_irlba2$x <- mat_irlba2$u %*% diag(mat_irlba2$d)
    ## version 2
    #X_center <- rowMeans(X)
    #X_scale <- proxyC::rowSds(X)
    #suppressWarnings({
    # retval <- irlba::irlba(A = t(X), nv = 20, center = X_center, scale = X_scale)
    #})
    #retval$x <- retval$u %*% diag(retval$d)
    #retval

    pca.result <- stats::prcomp(expr.select, scale. = F, center = F)
    pca.dims <- pca.result[["x"]]
    expr.select <- pca.dims[, 1:n.pca.dims]
    message("Done. ", Sys.time())
  }


  # find communities (clusters)
  tryCatch(
    if (run.louvain || run.leiden) {
      temp_dots <- dots[which(grepl("^louvain__|^leiden__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^louvain__|^leiden__", "", names(temp_dots), ignore.case = T)

      if (!any(grepl("annoy.metric", names(temp_dots), ignore.case = T))) {
        message("shared nearest neighbor (snn) calculated with annoy.metric = 'cosine' by default.")
        temp_dots <- c(temp_dots, annoy.metric = "cosine")
      }
      message("Calculating snn for louvain and/or leiden. Start: ", Sys.time())
      rownames(expr.select) <- 1:nrow(expr.select)
      # Seurat::FindNeighbors ignores all 'wrong' arguments; suppress the warnings though
      snn <- suppressMessages(suppressWarnings(do.call(Seurat::FindNeighbors, args = c(list(object = expr.select), temp_dots))))
      message("End: ", Sys.time())
    },
    error = function(e) {
      message("louvain / leiden: error in snn calculation.")
      run.louvain <- F
      run.leiden <- F
    }
  )

  tryCatch(
    if (run.louvain) {

      temp_dots <- dots[which(grepl("^louvain__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^louvain__", "", names(temp_dots), ignore.case = T)

      message("Finding clusters with Seurats implementation of the Louvain algorithm and parallel::mclapply using ", mc.cores," cores. Start: ",Sys.time())

      if (any(grepl("resolution", names(temp_dots), ignore.case = T))) {
        temp_dots[["resolution"]] <- as.numeric(temp_dots[["resolution"]])
        if (any(is.na(temp_dots[["resolution"]]))) {
          warning("Provide numeric values for louvain__resolution! Non numeric elements are ignored.")
          temp_dots[["resolution"]] <- temp_dots[["resolution"]][which(!is.na(temp_dots[["resolution"]]))]
        }
        clust_idents <- do.call(cbind,parallel::mclapply(temp_dots[["resolution"]], function(x) {
          apply(do.call(Seurat::FindClusters, args = c(list(object = snn$snn, resolution = x, verbose = F, algorithm = 1), temp_dots[which(names(temp_dots) != "resolution")])), 2, as.numeric)
        }, mc.cores = mc.cores))
      } else {
        clust_idents <- apply(do.call(Seurat::FindClusters, args = c(list(object = snn$snn, resolution = x, algorithm = 1), temp_dots)), 2, as.numeric)
      }

      colnames(clust_idents) <- paste0("louvain_", temp_dots[["resolution"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, clust_idents))
      message(apply(clust_idents, 2, function(x) length(unique(x))))
      message("End: ", Sys.time())
    },
    error = function(e) {
      message("run.louvain with error")
    }
  )

  tryCatch(
    if (run.leiden) {
      temp_dots <- dots[which(grepl("^leiden__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^leiden__", "", names(temp_dots), ignore.case = T)

      message("Finding clusters with leiden algorithm and parallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time())

      if (any(grepl("resolution_parameter", names(temp_dots), ignore.case = T))) {
        temp_dots[["resolution_parameter"]] <- as.numeric(temp_dots[["resolution_parameter"]])
        if (any(is.na(temp_dots[["resolution_parameter"]]))) {
          warning("Provide numeric values for leiden__resolution! Non numeric elements are ignored.")
          temp_dots[["resolution_parameter"]] <- temp_dots[["resolution_parameter"]][which(!is.na(temp_dots[["resolution_parameter"]]))]
        }
        clust_idents <- do.call(cbind, parallel::mclapply(temp_dots[["resolution_parameter"]], function(x) {
          do.call(leiden::leiden, args = c(list(object = snn$snn, resolution_parameter = x), temp_dots[which(names(temp_dots) != "resolution_parameter")]))
        }, mc.cores = mc.cores))
      } else {
        clust_idents <- do.call(leiden::leiden, args = c(list(object = snn$snn), temp_dots))
      }
      colnames(clust_idents) <- paste0("leiden_", temp_dots[["resolution_parameter"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, clust_idents))
      message(apply(clust_idents, 2, function(x) length(unique(x))))
      message("End: ", Sys.time())
    },
    error = function(e) {
      message("run.leiden with error")
    }
  )

  tryCatch(
    if (run.kmeans_arma) {
      temp_dots <- dots[which(grepl("^kmeans_arma__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^kmeans_arma__", "", names(temp_dots), ignore.case = T)
      message("Finding clusters with kmeans_arma and parallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time())

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["clusters"]], function(x) {
        ClusterR::predict_KMeans(data = expr.select, CENTROIDS = do.call(ClusterR::KMeans_arma, args = c(list(data = expr.select, clusters = x), temp_dots[which(names(temp_dots) != "clusters")])))
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("kmeans_arma_", temp_dots[["clusters"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    },
    error = function(e) {
      message("run.kmeans_arma with error")
    }
  )

  tryCatch(
    if (run.kmeans_rcpp) {
      temp_dots <- dots[which(grepl("^kmeans_rcpp__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^kmeans_rcpp__", "", names(temp_dots), ignore.case = T)
      message("Finding clusters with kmeans_rcpp and parallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time())

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["clusters"]], function(x) {
        do.call(ClusterR::KMeans_rcpp, args = c(list(data = expr.select, clusters = x), temp_dots[which(names(temp_dots) != "clusters")]))[["clusters"]]
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("kmeans_rcpp_", temp_dots[["clusters"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    },
    error = function(e) {
      message("run.kmeans_rcpp with error")
    }
  )

  tryCatch(
    if (run.minibatchkmeans) {
      temp_dots <- dots[which(grepl("^minibatchkmeans__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^minibatchkmeans__", "", names(temp_dots), ignore.case = T)
      message("Finding clusters with minibatchkmeans and parallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time())

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["clusters"]], function(x) {
        ClusterR::predict_MBatchKMeans(data = expr.select, CENTROIDS = do.call(ClusterR::MiniBatchKmeans, args = c(list(data = expr.select, clusters = x), temp_dots[which(names(temp_dots) != "clusters")]))[["centroids"]])
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("minibatchkmeans_", temp_dots[["clusters"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    },
    error = function(e) {
      message("run.minibatchkmeans with error")
    }
  )

  tryCatch(
    if (run.kmeans) {
      temp_dots <- dots[which(grepl("^kmeans__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^kmeans__", "", names(temp_dots), ignore.case = T)
      message("Finding clusters with kmeans and parallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time())

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["centers"]], function(x) {
        do.call(stats::kmeans, args = c(list(x = expr.select, centers = x), temp_dots[which(names(temp_dots) != "centers")]))$cluster
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("kmeans_", temp_dots[["centers"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    },
    error = function(e) {
      message("run.kmeans with error")
    }
  )

  tryCatch(
    if (run.hclust) {
      temp_dots <- dots[which(grepl("^hclust__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^hclust__", "", names(temp_dots), ignore.case = T)
      message("Finding clusters with hclust. Start: ", Sys.time())

      d <- do.call(stats::dist, args = c(list(x = expr.select), temp_dots[which(names(temp_dots) %in% names(formals(stats::dist))[-1])]))
      h <- do.call(stats::hclust, args = c(list(d = d), temp_dots[which(names(temp_dots) %in% names(formals(stats::hclust))[-1])]))
      ks <- do.call(cbind, parallel::mclapply(temp_dots[["k"]], function(x) {
        stats::cutree(tree = h, k = x)
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("hclust_", temp_dots[["k"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    },
    error = function(e) {
      message("run.hclust with error")
    }
  )

  tryCatch(
    if (run.flowClust) {

      temp_dots <- dots[which(grepl("^flowClust__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^flowClust__", "", names(temp_dots), ignore.case = T)

      message("Finding clusters with flowClust and parallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time())
      ks <- do.call(cbind, parallel::mclapply(temp_dots[["K"]], function(x) {
        suppressMessages(do.call(flowClust::flowClust, args = c(list(x = expr.select, K = x), temp_dots[which(names(temp_dots) != "K")]))@label)
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("flowClust_", temp_dots[["K"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    },
    error = function(e) {
      message("run.flowClust with error")
    }
  )

  tryCatch(
    if (run.MUDAN) {

      temp_dots <- dots[which(grepl("^MUDAN__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^MUDAN__", "", names(temp_dots), ignore.case = T)

      message("Finding clusters with MUDAN.", Sys.time())

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["k"]], function(x) {
        ## documentation is wrong (mat: cells as rows and features as cols!)
        do.call(MUDAN::getComMembership, args = c(list(mat = expr.select, k = x), temp_dots[which(names(temp_dots) != "k")]))
      }, mc.cores = mc.cores))

      colnames(ks) <- paste0("MUDAN_", temp_dots[["k"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    },
    error = function(e) {
      message("run.MUDAN with error")
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
      message("UMAP metric set to 'cosine' by default.")
      temp_dots <- c(temp_dots, metric = "cosine")
    }

    message("Calculating UMAP. Start: ", Sys.time())
    if (any(grepl("n_neighbors", names(temp_dots), ignore.case = T))) {
      umap.dims <- do.call(cbind,parallel::mclapply(temp_dots[["n_neighbors"]], function(z) {
        out <- do.call(uwot::umap, args = c(list(X = expr.select, verbose = F, n_neighbors = z),temp_dots[which(names(temp_dots) != "n_neighbors")]))
        colnames(out) <- c(paste0("UMAP_1_", z), paste0("UMAP_2_", z))
        return(out)
      }, mc.cores = mc.cores))
    } else {
      if (!any(grepl("verbose", names(temp_dots)), ignore.case = T)) {
        temp_dots <- c(temp_dots, verbose = T)
      }
      umap.dims <- do.call(uwot::umap, args = c(list(X = expr.select), temp_dots))
    }
    if (!any(grepl("n_neighbors", names(temp_dots), ignore.case = T)) || length(temp_dots[["n_neighbors"]]) == 1) {
      colnames(umap.dims) <- c("UMAP_1", "UMAP_2")
    }
    message("End: ", Sys.time())
  }

  if (run.tsne) {
    temp_dots <- dots[which(grepl("^tSNE__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^tSNE__", "", names(temp_dots), ignore.case = T)

    if (run.pca) {
      message("As 'run.pca=T' inital pca in tSNE is not performed.")
      temp_dots <- c(temp_dots, pca = F)
    }
    if (!any(grepl("normalize", names(temp_dots), ignore.case = T))) {
      message("Rtsne normalize set to 'FALSE' by default. (scale.whole, scale.samples?!)")
      temp_dots <- c(temp_dots, normalize = F)
    }

    message("Calculating tSNE. Start: ", Sys.time())
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
    message("End: ", Sys.time())
  }

  if (run.som) {
    message("Calculating SOM. Start: ", Sys.time())
    temp_dots <- dots[which(grepl("^SOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^SOM__", "", names(temp_dots), ignore.case = T)
    map <- do.call(EmbedSOM::SOM, args = c(list(data = expr.select), temp_dots))

    temp_dots <- dots[which(grepl("^EmbedSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^EmbedSOM__", "", names(temp_dots), ignore.case = T)
    som.dims <- do.call(EmbedSOM::EmbedSOM, args = c(list(data = expr.select, map = map), temp_dots))
    colnames(som.dims) <- c("SOM_1", "SOM_2")
    message("End: ", Sys.time())
  }

  if (run.gqtsom) {
    message("Calculating GQTSOM. Start: ", Sys.time())
    temp_dots <- dots[which(grepl("^GQTSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^GQTSOM__", "", names(temp_dots), ignore.case = T)
    map <- do.call(EmbedSOM::GQTSOM, args = c(list(data = expr.select), temp_dots))

    temp_dots <- dots[which(grepl("^EmbedSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^EmbedSOM__", "", names(temp_dots), ignore.case = T)
    gqtsom.dims <- do.call(EmbedSOM::EmbedSOM, args = c(list(data = expr.select, map = map), temp_dots))
    colnames(gqtsom.dims) <- c("GQTSOM_1", "GQTSOM_2")
    message("End: ", Sys.time())
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
    dim.red.data <- do.call(cbind, list(dim.red.data, pca.dims[,1:n.pca.dims]))
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

  if (!is.null(extra.cols)) {
    if (nrow(extra.cols) == nrow(dim.red.data)) {
      dim.red.data <- cbind(dim.red.data, extra.cols)
    } else {
      message("extra.cols not added due to wrong row number.")
    }
  }

  # 2021 06 17 necessary
  dim.red.data <- as.data.frame(dim.red.data)

  tryCatch(
    if (!is.null(add.sample.info)) {
      for (i in names(add.sample.info)) {
        dim.red.data <- do.call(cbind, list(dim.red.data, rep(add.sample.info[[i]], times = as.numeric(table(rep(1:length(ff.list[["logicle"]]), sapply(ff.list[["logicle"]], nrow)))))))
        names(dim.red.data)[length(dim.red.data)] <- i
      }
    },
    error = function(e) {
      message("Addition of sample info failed due to an error. Check add.sample.info.")
    }
  )

  # prepare channel desc
  name.desc <- stats::setNames(ff.list[[1]][[1]]@parameters@data[["desc"]], ff.list[[1]][[1]]@parameters@data[["name"]])
  name.desc <- name.desc[which(!is.na(name.desc))]
  channel.desc <- rep("", ncol(dim.red.data))

  for (i in seq_along(name.desc)) {
    channel.desc[grep(names(name.desc)[i], colnames(dim.red.data))] <- name.desc[i]
  }


  channel.desc_augment <- channel.desc
  channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("scaled", colnames(dim.red.data))))] <- paste0(channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("scaled", colnames(dim.red.data))))], "_scaled")
  channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("logicle", colnames(dim.red.data))))] <- paste0(channel.desc_augment[intersect(which(channel.desc_augment != ""), which(grepl("logicle", colnames(dim.red.data))))], "_logicle")

  channel.desc_augment[which(channel.desc_augment == "")] <- colnames(dim.red.data)[which(channel.desc_augment == "")]
  channel.desc_augment <- make.names(channel.desc_augment)
  names(channel.desc_augment) <- colnames(dim.red.data)

  # write FCS file

  # get common (intersecting keywords)
  # a bit unhandy but belows vectorized version did not work
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

  # get cluster markers
  ## always used logicle transformed data?!?!
  if (!is.null(calc.cluster.markers)) {
    message("Calculating markers.")
  }
  tryCatch({
    marker <- lapply(calc.cluster.markers, function (clust_col) {
      # do not use expr.select which may have become dimenions of LDA
      dat <- dim.red.data[,c(which(colnames(dim.red.data) %in% paste0(channels, "_logicle")), which(colnames(dim.red.data) == clust_col))]
      split_var <- dat[,clust_col]
      split_var_levels <- sort(unique(split_var))
      ## keep dat a data frame until here to allow split (works only on data.frame); after that convert to matrix
      dat_split <- split(dat, split_var)
      dat <- as.matrix(dat)
      dat_split <- lapply(dat_split, function(x) as.matrix(x[,-which(names(x) == clust_col)]))
      all_pairs <- utils::combn(split_var_levels, 2, simplify = F)

      ## all pairwise
      pair_marker_table <- dplyr::bind_rows(parallel::mclapply(all_pairs, function(x) {
        #out <- matrixTests::col_wilcoxon_twosample(dat_split[[as.character(x[1])]], dat_split[[as.character(x[2])]])
        out <-
          presto::wilcoxauc(cbind(t(dat_split[[as.character(x[1])]]),t(dat_split[[as.character(x[2])]])), c(rep("y", nrow(dat_split[[as.character(x[1])]])), rep("z", nrow(dat_split[[as.character(x[2])]])))) %>%
          dplyr::filter(group == "y") %>%
          dplyr::select(feature, pval) %>%
          dplyr::rename("pvalue" = pval) %>%
          dplyr::rename("channel" = feature)

        ## wilcox.test procedure:
        '      p <- lapply(seq_along(colnames(dat_split[[as.character(x[1])]])), function(k) wilcox.test(dat_split[[as.character(x[1])]][,k], dat_split[[as.character(x[2])]][,k])[["p.value"]])
      names(p) <- colnames(dat_split[[as.character(x[1])]])
      out <- utils::stack(p)
      out[,2] <- as.character(out[,2])
      names(out) <- c("pvalue", "channel")'

        out[,"mean_1"] <- round(matrixStats::colMeans2(dat_split[[as.character(x[1])]]), 2)
        out[,"mean_2"] <- round(matrixStats::colMeans2(dat_split[[as.character(x[2])]]), 2)
        out[,"mean_diff"] <- round(out[,"mean_1"] - out[,"mean_2"], 2)
        out[,"diptest_p_1"] <- suppressMessages(round(apply(dat_split[[as.character(x[1])]], 2, function(x) diptest::dip.test(x)[["p.value"]]), 2))
        out[,"diptest_p_2"] <- suppressMessages(round(apply(dat_split[[as.character(x[2])]], 2, function(x) diptest::dip.test(x)[["p.value"]]), 2))
        #out <- tibble::rownames_to_column(out, "channel")
        out[,"cluster_1"] <- x[1]
        out[,"cluster_2"] <- x[2]
        out[,"diff_sign"] <- ifelse(out[,"mean_diff"] == 0, "+/-", ifelse(out[,"mean_diff"] > 0, "+", "-"))
        out <- dplyr::select(out, channel, cluster_1, cluster_2, pvalue, mean_1, mean_2, mean_diff, diff_sign, diptest_p_1, diptest_p_2)
        out <- dplyr::arrange(out, pvalue)
        return(out)
      }, mc.cores = mc.cores))
      pair_marker_table[,"channel_desc"] <- channel.desc_augment[pair_marker_table[,"channel"]]

      marker_table <- dplyr::bind_rows(parallel::mclapply(split_var_levels, function(x) {
        y <- t(dat[which(dat[,clust_col] == x),which(colnames(dat) != clust_col)])
        z <- t(dat[which(dat[,clust_col] != x),which(colnames(dat) != clust_col)])
        out <-
          presto::wilcoxauc(cbind(y,z), c(rep("y", ncol(y)), rep("z", ncol(z)))) %>%
          dplyr::filter(group == "y") %>%
          dplyr::select(feature, pval) %>%
          dplyr::rename("pvalue" = pval) %>%
          dplyr::rename("channel" = feature)

        #out <- matrixTests::col_wilcoxon_twosample(y, z) # produced error once
        # parallel::mclapply
        ## wilcox.test procedure:
        '
y <- dat[which(dat[,clust_col] == x),which(colnames(dat) != clust_col)]
      z <- dat[which(dat[,clust_col] != x),which(colnames(dat) != clust_col)]
p <- lapply(seq_along(colnames(y)), function(k) wilcox.test(y[,k], z[,k])[["p.value"]])
      names(p) <- colnames(y)
      out <- utils::stack(p)
      out[,2] <- as.character(out[,2])
      names(out) <- c("pvalue", "channel")'

        out[,"mean"] <- round(matrixStats::rowMeans2(y), 2)
        out[,"mean_not"] <- round(matrixStats::rowMeans2(z), 2)
        out[,"mean_diff"] <- round(out[,"mean"] - out[,"mean_not"], 2)
        out[,"diptest_p"] <- suppressMessages(round(apply(y, 1, function(x) diptest::dip.test(x)[["p.value"]]), 2))
        out[,"diptest_not_p"] <- suppressMessages(round(apply(z, 1, function(x) diptest::dip.test(x)[["p.value"]]), 2))
        #out <- tibble::rownames_to_column(out, "channel")
        out[,"cluster"] <- x
        out[,"diff_sign"] <- ifelse(out[,"mean_diff"] == 0, "+/-", ifelse(out[,"mean_diff"] > 0, "+", "-"))
        out <- dplyr::select(out, channel, cluster, pvalue, mean, mean_not, mean_diff, diff_sign, diptest_p, diptest_not_p)
        out <- dplyr::arrange(out, pvalue)
        return(out)
      }, mc.cores = mc.cores))
      marker_table[,"channel_desc"] <- channel.desc_augment[marker_table[,"channel"]]
      return(list(marker_table = marker_table, pairwise_marker_table = pair_marker_table))
    })
    names(marker) <- calc.cluster.markers
  }, error = function(e) {
    message("cluster marker calculation caused an error.")
    marker <- NULL
  }
  )




  # save results
  if (!is.null(save.path) && !is.na(save.path)) {
    message("Writing files to disk.")
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

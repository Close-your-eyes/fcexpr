#' Calculate dimension reduction and cluster annotation with data from one or more flow frames and add these parameters to a (concatenated) fcs file
#'
#' Prepare a list of flowframes (ff.list) with fcexpr::wsp_get_ff or fcexpr::inds_get_ff. The objects returned from these function will be compatible with fcexpr::dr_to_fcs.
#' The ff.list must contain a list of flowframes with untransformed fluorescence intensities (FI) and optionally may contain and additional list of flowframes with transformed FI. The transformation
#' is your choice and may be logicle, arcsinh or whatever. It may be individual to channels and/or to flowframes. You have to calculate the transformation outside this function (dr_to_fcs) though.
#' If a list of transformed flowframes is provided this one will be used for dimension reduction and clustering etc. If only a list untransformed flowframes is provided, this one will be used.
#' Both, the transformed and/or untransformed values will be written to the concatenated fcs file (see function arguments). Calculated dimension reductions and cluster annotations will
#' also be written as channels into that fcs file. Some dimension reduction and cluster detection algorithm may become slow with a large number of events (e.g. > 1e6, see the details).
#' In order to get a quick impression of what the algorithms can pull out for you, you may use the 'downsample'
#' argument in fcexpr::wsp_get_ff or fcexpr::inds_get_ff to conveniently sample a random subset of events from each fcs files (flowframe).
#'
#' Logicle transformation of FI has been described here: \href{https://pubmed.ncbi.nlm.nih.gov/16604519/}{Parks, et al., 2006, PMID: 16604519  DOI: 10.1002/cyto.a.20258}.
#' Different transformations have been compared for instance here: \href{https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR-Transformations.html}{Transformations}.
#' Dimension reduction (low dimension embedding) algorithms are: UMAP, tSNE, SOM, GQTSOM, PCA. tSNE is expected to be too slow for more then 1e4 or 1e5 cells (depending
#' precision set by tsne__theta). UMAP has been developed more recently, is well-know and quicker. SOM and GQTSOM are very quick and are similar to \href{https://github.com/SofieVG/FlowSOM}{flowSOM}
#' but are much more convenient to use in programming as they accept matrices whereas flowSOM strictly requires flowFrames. SOMs also enable the so-called meta-clustering
#' which massively speeds up cluster detection. The cost of preciseness has to be inspected manually.
#' \href{"https://www.youtube.com/watch?v=FgakZw6K1QQ"}{PCA} will be orders of magnitude faster than UMAP especially on large numbers of events (1e6 and above).
#' On the other hand it will not produce as nicely separated clusters as it is a linear algorithm.
#' Cluster (community) detection algorithms are louvain, leiden, kmeans, hclust, flowClust and MUDAN. They assign events with similar properties in the high dimensional space a common discrete label.
#' kmeans will be the quickest. Set metaclustering.on to 'SOM' or 'GQTSOM' to enable quick metaclustering. In this any of the clustering algorithms is performed on the SOM-map.
#' Louvain is slower but may yield better results on original data (not meta-clustering). For a low number of events (e.g. below 1e6) louvain will perform in a reasonable amount of time
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
#' @param run.umap logical, whether to calculate UMAP dimension reduction with uwot::umap
#' @param run.tsne logical, whether to calculate tsne dimension reduction with Rtsne::Rtsne
#' @param run.som logical, whether to calculate SOM dimension reduction EmbedSOM::SOM followed by EmbedSOM::EmbedSOM
#' @param run.gqtsom logical, whether to calculate GQTSOM dimension reduction EmbedSOM::GQTSOM followed by EmbedSOM::EmbedSOM
#' @param metaclustering.on SOM or GQTSOM; run clustering-algorithms not on original data but on SOM map (so called metaclustering); SOM or GQTSOM
#' may be selected for metaclustering; this is expected to yield a substantial improvement to calculation speed on large data sets;
#' if NULL, clustering is performed on original data; recommendation: when applying meta-clustering choose run.hclust = T and supply hclust__method = "average"
#' and cutree__k = c(xx) with xx being the number of desired/expected clusters
#' @param run.louvain logical, whether to detect clusters (communities, groups) of cells with the louvain algorithm, implemented in Seurat::FindClusters (subsequent to snn detection by Seurat::FindNeighbors)
#' @param run.leiden logical, whether to detect clusters (communities, groups) of cells with the leiden algorithm, with leiden::leiden (subsequent to snn detection by Seurat::FindNeighbors)
#' @param run.kmeans logical, whether to detect clusters with stats::kmeans
#' @param run.minibatchkmeans logical, whether to detect clusters with ClusterR::MiniBatchKmeans
#' @param run.kmeans_arma logical, whether to detect clusters with ClusterR::KMeans_arma
#' @param run.kmeans_rcpp logical, whether to detect clusters with ClusterR::KMeans_rcpp
#' @param run.hclust logical, whether to detect clusters with stats::dist, stats::hclust and stats::cutree
#' @param run.flowClust logical, whether to detect clusters with flowClust::flowClust
#' @param run.MUDAN logical, whether to detect clusters with MUDAN::getComMembership
#' @param extra.cols numeric vector of an extra column (or matrix of multiple columns) with arbitraty information to add to the final fcs file;
#' has to be equal to the number of rows in all flowframes provided; colnames will be the channel names in the FCS file;
#' could be a previously calculated dimension reduction or cluster annotation.
#' @param clustering.for.marker.calc if NULL nothing is calculated; otherwise marker features (stained markers) are determined by wilcox test
#' using \href{https://github.com/immunogenomics/presto}{presto::wilcoxauc} for the provided clustering(s). each cluster
#' is tested against other events and clusters are compared pairwise. respective clustering calculation has to be provided in ...;
#' e.g. if louvain__resolution = 0.5 is provided set clustering.for.marker.calc = louvain_0.5;
#' and if in addition leiden__resolution_parameter = 0.7 then set clustering.for.marker.calc = c(louvain_0.5, leiden_0.7).
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
#' @param ... additional parameters to calculations of \href{https://github.com/jlmelville/uwot}{UMAP}, \href{https://github.com/jkrijthe/Rtsne}{tSNE}, \href{https://github.com/exaexa/EmbedSOM}{SOM, GQTSOM, EmbedSOM}, louvain, \href{https://github.com/TomKellyGenetics/leiden}{leiden},
#' \href{https://github.com/immunogenomics/harmony}{harmony}, \href{https://www.bioconductor.org/packages/release/bioc/html/flowClust.html}{flowClust}, hclust, \href{https://github.com/JEFworks/MUDAN}{MUDAN}, kmeans;
#' provide arguments as follows: UMAP__n_neighbors = c(15,20,25), or tsne__theta = 0.3, etc.
#' see respected help files to get to know which arguments can be passed:
#' uwot::umap, Rtsne::Rtsne, EmbedSOM::SOM, EmbedSOM::GQTSOM, EmbedSOM::EmbedSOM, harmony::HarmonyMatrix, flowClust::flowClust,
#' louvain: Seurat::FindNeighbors and Seurat::FindCluster, leiden: Seurat::FindNeighbors and leiden::leiden.
#' hclust: stats::dist and stats::hclust, MUDAN: MUDAN::getComMembership, stats::kmeans,
#' ClusterR::MiniBatchKmeans, ClusterR::KMeans_rcpp, ClusterR::KMeans_arma
#' @param transformation_name name of the applied transformation (will appear in FCS file as channel desc)
#' @param return_processed_raw_data_only do not calculate dimension reduction etc but only return the preprocessed
#' data for external calculations or tryouts
#' @param seed set a seed for reproduction of dimension reductions
#' @param calc.global.markers logical whether to calculate global cluster markers: so each cluster vs. all other cells
#' @param calc.pairwise.markers logical whether to calculate pairwise cluster markers: so each cluster vs. each cluster
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
                      metaclustering.on = NULL,
                      run.louvain = F,
                      run.kmeans = F,
                      run.minibatchkmeans = F,
                      run.kmeans_arma = F,
                      run.kmeans_rcpp = F,
                      run.leiden = F,
                      run.hclust = F,
                      #run.mhclust = F,
                      run.flowClust = F,
                      run.MUDAN = F,
                      n.pca.dims = 0,
                      clustering.for.marker.calc = NULL,
                      calc.global.markers = F,
                      calc.pairwise.markers = F,
                      #cluster.marker.min.cells = 50,
                      extra.cols = NULL,
                      mc.cores = 1,
                      save.to.disk = c("fcs", "rds"),
                      save.path = file.path(getwd(), paste0(substr(gsub("\\.", "", make.names(as.character(Sys.time()))), 2, 15), "_FCS_dr")),
                      exclude.extra.channels = ifelse(length(ff.list) == 1 && names(ff.list) == "transformed", "cluster.id", "FSC|SSC|Time|cluster.id"),
                      write.transformed.channels.to.FCS = T,
                      write.untransformed.channels.to.FCS = T,
                      write.scaled.channels.to.FCS = F,
                      timeChannel = "Time",
                      transformation_name = "trans",
                      return_processed_raw_data_only = F,
                      seed = 42,
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
  #df <- data.frame(x1 = c(rnorm(1e5,0,1), rnorm(1e4,15,1)), x2 = c(rnorm(1e5,1,1), rnorm(1e4,9,1)),x3 = c(rnorm(1e4,-1,4), rnorm(1e5,30,5)))
  #df <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(df)))
  #names(df) <- paste0("x", 1:ncol(df))
  #df <- tidyr::pivot_longer(df, names_to = "name", values_to = "value", cols = c(x1,x2,x3))
  #ggplot(df, aes(x = value, y = name)) + ggridges::geom_density_ridges() + ggridges::theme_ridges()'

  ## make warning when calc.global.markers = T or calc.pairwise.markers = T but is.null(clustering.for.marker.calc) or vice versa


  ## ---- checks --------
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
  if (!is.null(clustering.for.marker.calc) && !requireNamespace("presto", quietly = T)) {
    devtools::install_github('immunogenomics/presto')
  }
  if (!is.null(clustering.for.marker.calc) && !requireNamespace("matrixTests", quietly = T)) {
    utils::install.packages("matrixTests")
  }
  if (run.harmony && !requireNamespace("harmony", quietly = T)) {
    devtools::install_github("immunogenomics/harmony")
  }
  if (run.MUDAN && !requireNamespace("MUDAN", quietly = T)) {
    devtools::install_github("JEFworks/MUDAN")
  }
  if (run.som && !requireNamespace("EmbedSOM", quietly = T)) {
    devtools::install_github("exaexa/EmbedSOM")
  }
  if ((run.minibatchkmeans || run.kmeans_rcpp || run.kmeans_arma) && !requireNamespace("ClusterR", quietly = T)) {
    utils::install.packages("ClusterR")
  }
  # 2022 09 07
  run.mhclust <- F
  if (run.mhclust && !is.null(metaclustering.on) && !requireNamespace("mhca", quietly = T)) {
    devtools::install_github("tsieger/mhca")
  }

  if (run.mhclust && is.null(metaclustering.on)) {
    run.mhclust <- F
    message("run.mhclust should only be ran with metaclustering.on = T. Otherwise it will be too slow. run.mhclust is set FALSE now.")
  }

  dots <- list(...)

  expect_dots <- "^harmony__|^hclust__|^dist__|^cutree__|^flowClust_|^MUDAN__|^kmeans__|^louvain__|^leiden__|^som__|^gqtsom__|^tsne__|^umap__|^EmbedSOM|^kmeans_arma__|^kmeans_rcpp__|^minibatchkmeans__"
  if (any(!names(dots) %in% names(formals(dr_to_fcs)) & !grepl(expect_dots, names(dots), ignore.case = T))) {
    message("These arguments are unknown: ", paste(names(dots)[which(!names(dots) %in% names(formals(dr_to_fcs)) & !grepl(expect_dots, names(dots), ignore.case = T))], collapse = ", "))
  }

  if (length(dots) > 0) {
    dots_sub <- dots[which(!sapply(dots, is.function))]
    dots_expanded <- unname(unlist(mapply(paste, sapply(strsplit(names(dots_sub), "__"), "[", 1), dots_sub, sep = "_")))
  } else {
    dots_expanded <- NULL
  }

  if (any(!names(ff.list) %in% c("untransformed", "transformed"))) {
    stop("ff.list has to contain a list of flowframes named 'transformed' and optionally an additional list named 'untransformed' (inverse transformed, original as in flowjo.). Check the output of fcexpr::wsp_get_ff for valid input.")
  }

  # check if names in ff.list are the same
  if (length(unique(lengths(ff.list))) != 1) {
    stop("number of flowframes in untransformed and transformed has to be equal.")
  }

  for (par in c("louvain", "leiden", "umap", "tsne", "som", "gqtsom", "harmony", "kmeans", "kmeans_rcpp", "kmeans_arma", "minibatchkmeans", "flowclust", "hclust", "mhclust", "harmony", "mudan")) {
    # cutree and dist not captured
    if (any(grepl(paste0("^", par, "__"), names(dots), ignore.case = T)) &&!eval(rlang::sym(paste0("run.", par)))) {
      message(par, " parameters provided in ... but ", "'run.", par, " = F'.")
    }
  }

  if (run.harmony && !any(grepl("^harmony__meta_data", names(dots)))) {
    stop("When 'run.harmony = T' harmony__meta_data has to be provided in ..., see ?harmony::HarmonyMatrix.")
  }

  if (run.hclust && !any(grepl("^cutree__k", names(dots)))) {
    stop("When 'run.hclust = T' cutree__k has to be provided in ..., see ?stats::cutree.")
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
    stop("Do not set 'run.lda = T' but provide a clustering that should be used to calculate it, e.g. a pattern like 'louvain_0.4' or 'leiden_1.1' or 'kmeans_20'.")
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

  if (!is.null(clustering.for.marker.calc)) {
    if (!any(clustering.for.marker.calc %in% dots_expanded)) {
      stop("clustering.for.marker.calc: ", clustering.for.marker.calc[which(!clustering.for.marker.calc %in% dots_expanded)], " not found in ... .")
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
  if (!is.null(metaclustering.on)) {
    metaclustering.on <- match.arg(arg = metaclustering.on, choices = c("SOM", "GQTSOM"))
  }

  if (!is.null(metaclustering.on) && metaclustering.on == "SOM" && !run.som) {
    run.som <- T
    message("run.som set to TRUE to allow metaclustering on it.")
  }
  if (!is.null(metaclustering.on) && metaclustering.on == "GQTSOM" && !run.gqtsom) {
    run.gqtsom <- T
    message("run.gqtsom set to TRUE to allow metaclustering on it.")
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

  scale.whole <- switch(match.arg(scale.whole, c("z.score", "min.max", "none")),
                        z.score = scale,
                        min.max = min.max.normalization,
                        none = function(x) return(x))


  # check if channel names and desc are equal
  if (!is.null(fcs_check <- .check.ff.list(ff.list = ff.list, channels = channels, strict = T))) {
    return(fcs_check)
  }

  ## channel names from first ff
  channels <- .get.channels(ff = ff.list[[1]][[1]], timeChannel = timeChannel, channels = channels)

  # overwrite channel desc in ffs
  # correct order is important, as provided by .get.channels
  #
  '  if ("transformed" %in% names(ff.list)) {
    for (i in seq_along(ff.list[["transformed"]])) {
      ff.list[["transformed"]][[i]]@parameters@data[["desc"]][which(ff.list[["transformed"]][[i]]@parameters@data[["name"]] %in% channels)] <- names(channels)
    }
  }
  if ("untransformed" %in% names(ff.list)) {
    for (i in seq_along(ff.list[["untransformed"]])) {
      ff.list[["untransformed"]][[i]]@parameters@data[["desc"]][which(ff.list[["untransformed"]][[i]]@parameters@data[["name"]] %in% channels)] <- names(channels)
    }
  }'

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

  ## apply scaling which was selected above and select channels to use for dimension reduction.
  expr.select <- scale.whole(do.call(rbind, lapply(ff.list[[ff.list_index]], function(x) scale.samples(flowCore::exprs(x)[, channels])))) # names to desc here?!


  ### allow to pass expr.select here.
  ## requires a lot of checking though

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

  ## pca in harmony has to be set to TRUE explicitly, then harmony is performed in pc-space
  ## in this case subsequent pca is not advisable
  pca.result <- NULL
  if (run.pca) {
    if (n.pca.dims == 0) {
      n.pca.dims <- ncol(expr.select) - 1
    } else if (n.pca.dims > 0) {
      n.pca.dims <- min(ncol(expr.select) - 1, n.pca.dims)
    }
    message("Calculating PCA.\nStart: ", Sys.time())
    pca.result <- stats::prcomp(expr.select, scale. = F, center = F)
    pca.dims <- pca.result[["x"]]
    expr.select <- pca.dims[, 1:n.pca.dims] # overwrite original data with PCA embedding
    message("Done. ", Sys.time())
  }

  # first exit here
  if (return_processed_raw_data_only) {
    return(expr.select)
  }

  ## ---- metaclustering -------

  if (!is.null(metaclustering.on) && metaclustering.on == "SOM") {
    if (run.lda) {
      warning("LDA not applied to SOM calculation.")
    }
    message("Calculating SOM.\nStart: ", Sys.time())

    temp_dots <- dots[which(grepl("^SOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^SOM__", "", names(temp_dots), ignore.case = T)
    set.seed(seed)
    map <- do.call(EmbedSOM::SOM, args = c(list(data = expr.select), temp_dots))

    temp_dots <- dots[which(grepl("^EmbedSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^EmbedSOM__", "", names(temp_dots), ignore.case = T)
    som.dims <- do.call(EmbedSOM::EmbedSOM, args = c(list(data = expr.select, map = map), temp_dots))
    colnames(som.dims) <- c("SOM_1", "SOM_2")

    message("End: ", Sys.time())
  }

  if (!is.null(metaclustering.on) && metaclustering.on == "GQTSOM") {
    if (run.lda) {
      warning("LDA not applied to GQTSOM calculation.")
    }
    message("Calculating GQTSOM.\nStart: ", Sys.time())

    temp_dots <- dots[which(grepl("^GQTSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^GQTSOM__", "", names(temp_dots), ignore.case = T)
    set.seed(seed)
    map <- do.call(EmbedSOM::GQTSOM, args = c(list(data = expr.select), temp_dots))

    temp_dots <- dots[which(grepl("^EmbedSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^EmbedSOM__", "", names(temp_dots), ignore.case = T)
    gqtsom.dims <- do.call(EmbedSOM::EmbedSOM, args = c(list(data = expr.select, map = map), temp_dots))
    colnames(gqtsom.dims) <- c("GQTSOM_1", "GQTSOM_2")
    message("End: ", Sys.time())
  }

  ## in case of metaclustering set expr.clust to map-codes
  if (!is.null(metaclustering.on)) {
    expr.clust <- map[["codes"]]
  } else {
    expr.clust <- expr.select
  }

  ## ---- clusterings -------
  # find communities (clusters)
  if (run.louvain || run.leiden) {
    tryCatch({
      temp_dots <- dots[which(grepl("^louvain__|^leiden__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^louvain__|^leiden__", "", names(temp_dots), ignore.case = T)
      message("Calculating snn for louvain and/or leiden.\nStart: ", Sys.time())

      if (!any(grepl("annoy.metric", names(temp_dots), ignore.case = T))) {
        message("shared nearest neighbor (snn) calculated with annoy.metric = 'cosine' by default.")
        temp_dots <- c(temp_dots, annoy.metric = "cosine")
      }
      rownames(expr.clust) <- 1:nrow(expr.clust)
      # Seurat::FindNeighbors ignores all 'wrong' arguments; suppress the warnings though
      snn <- suppressMessages(suppressWarnings(do.call(Seurat::FindNeighbors, args = c(list(object = expr.clust), temp_dots))))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("louvain / leiden: error in snn calculation: ", err)
      run.louvain <- F
      run.leiden <- F
    })
  }

  if (run.louvain) {
    tryCatch({
      temp_dots <- dots[which(grepl("^louvain__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^louvain__", "", names(temp_dots), ignore.case = T)
      message("Finding clusters with Seurats implementation of the Louvain algorithm and parallel::mclapply using ", mc.cores," cores.\nStart: ",Sys.time())

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

      if (!is.null(metaclustering.on)) {
        clust_idents <- apply(clust_idents, 2, function (x) x[map[["mapping"]][,1]])
      }

      colnames(clust_idents) <- paste0("louvain_", temp_dots[["resolution"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, clust_idents))
      print(apply(clust_idents, 2, function(x) length(unique(x))))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("run.louvain with error: ", err)
    })
  }


  if (run.leiden) {
    tryCatch({
      temp_dots <- dots[which(grepl("^leiden__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^leiden__", "", names(temp_dots), ignore.case = T)

      message("Finding clusters with leiden algorithm and parallel::mclapply using ", mc.cores, " cores.\nStart: ", Sys.time())

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

      if (!is.null(metaclustering.on)) {
        clust_idents <- apply(clust_idents, 2, function (x) x[map[["mapping"]][,1]])
      }

      colnames(clust_idents) <- paste0("leiden_", temp_dots[["resolution_parameter"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, clust_idents))
      print(apply(clust_idents, 2, function(x) length(unique(x))))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("run.leiden with error: ", err)
    })
  }

  if (run.kmeans_arma) {
    tryCatch({
      temp_dots <- dots[which(grepl("^kmeans_arma__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^kmeans_arma__", "", names(temp_dots), ignore.case = T)
      message("Finding clusters with kmeans_arma and parallel::mclapply using ", mc.cores, " cores.\nStart: ", Sys.time())

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["clusters"]], function(x) {
        ClusterR::predict_KMeans(data = expr.clust, CENTROIDS = do.call(ClusterR::KMeans_arma, args = c(list(data = expr.clust, clusters = x), temp_dots[which(names(temp_dots) != "clusters")])))
      }, mc.cores = mc.cores))

      if (!is.null(metaclustering.on)) {
        ks <- apply(ks, 2, function (x) x[map[["mapping"]][,1]])
      }

      # make sure that cluster 1 is the largest and so on
      ks <- .cluster_ordering(ks = ks)

      colnames(ks) <- paste0("kmeans_arma_", temp_dots[["clusters"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("run.kmeans_arma with error: ", err)
    })
  }

  if (run.kmeans_rcpp) {
    tryCatch({
      temp_dots <- dots[which(grepl("^kmeans_rcpp__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^kmeans_rcpp__", "", names(temp_dots), ignore.case = T)
      message("Finding clusters with kmeans_rcpp and parallel::mclapply using ", mc.cores, " cores.\nStart: ", Sys.time())

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["clusters"]], function(x) {
        do.call(ClusterR::KMeans_rcpp, args = c(list(data = expr.clust, clusters = x), temp_dots[which(names(temp_dots) != "clusters")]))[["clusters"]]
      }, mc.cores = mc.cores))

      if (!is.null(metaclustering.on)) {
        ks <- apply(ks, 2, function (x) x[map[["mapping"]][,1]])
      }

      # make sure that cluster 1 is the largest and so on
      ks <- .cluster_ordering(ks = ks)

      colnames(ks) <- paste0("kmeans_rcpp_", temp_dots[["clusters"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("run.kmeans_rcpp with error: ", err)
    })
  }

  if (run.minibatchkmeans) {
    tryCatch({
      temp_dots <- dots[which(grepl("^minibatchkmeans__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^minibatchkmeans__", "", names(temp_dots), ignore.case = T)
      message("Finding clusters with minibatchkmeans and parallel::mclapply using ", mc.cores, " cores.\nStart: ", Sys.time())

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["clusters"]], function(x) {
        ClusterR::predict_MBatchKMeans(data = expr.clust, CENTROIDS = do.call(ClusterR::MiniBatchKmeans, args = c(list(data = expr.clust, clusters = x), temp_dots[which(names(temp_dots) != "clusters")]))[["centroids"]])
      }, mc.cores = mc.cores))

      if (!is.null(metaclustering.on)) {
        ks <- apply(ks, 2, function (x) x[map[["mapping"]][,1]])
      }

      # make sure that cluster 1 is the largest and so on
      ks <- .cluster_ordering(ks = ks)

      colnames(ks) <- paste0("minibatchkmeans_", temp_dots[["clusters"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("run.minibatchkmeans with error: ", err)
    })
  }

  if (run.kmeans) {
    tryCatch({
      temp_dots <- dots[which(grepl("^kmeans__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^kmeans__", "", names(temp_dots), ignore.case = T)
      message("Finding clusters with kmeans and parallel::mclapply using ", mc.cores, " cores.\nStart: ", Sys.time())

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["centers"]], function(x) {
        do.call(stats::kmeans, args = c(list(x = expr.clust, centers = x), temp_dots[which(names(temp_dots) != "centers")]))$cluster
      }, mc.cores = mc.cores))

      if (!is.null(metaclustering.on)) {
        ks <- apply(ks, 2, function (x) x[map[["mapping"]][,1]])
      }

      # make sure that cluster 1 is the largest and so on
      ks <- .cluster_ordering(ks = ks)

      colnames(ks) <- paste0("kmeans_", temp_dots[["centers"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("run.kmeans with error: ", err)
    })
  }

  if (run.hclust) {
    tryCatch({
      message("Finding clusters with hclust.\nStart: ", Sys.time())

      temp_dots <- dots[which(grepl("^dist__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^dist__", "", names(temp_dots), ignore.case = T)
      d <- do.call(stats::dist, args = c(list(x = expr.clust), temp_dots[which(names(temp_dots) %in% names(formals(stats::dist))[-1])]))

      temp_dots <- dots[which(grepl("^hclust__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^hclust__", "", names(temp_dots), ignore.case = T)
      h <- do.call(stats::hclust, args = c(list(d = d), temp_dots[which(names(temp_dots) %in% names(formals(stats::hclust))[-1])]))

      temp_dots <- dots[which(grepl("^cutree__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^cutree__", "", names(temp_dots), ignore.case = T)
      ks <- do.call(cbind, parallel::mclapply(temp_dots[["k"]], function(x) stats::cutree(tree = h, k = x), mc.cores = mc.cores))

      if (!is.null(metaclustering.on)) {
        ks <- apply(ks, 2, function (x) x[map[["mapping"]][,1]])
      }
      # make sure that cluster 1 is the largest and so on
      ks <- .cluster_ordering(ks = ks)

      colnames(ks) <- paste0("cutree_", temp_dots[["k"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("run.hclust with error: ", err)
    })
  }

  if (run.mhclust) {
    tryCatch({
      #run.mhclust logical, whether to run Mahalanobis distance-based hierarchical cluster analysis \href{https://github.com/tsieger/mhca}{(mhca)} (similar to hclust)
      # 2022 09 07: how to supply apriori clusters? # see shinysom
      # slow
      message("Finding clusters with mhclust.\nStart: ", Sys.time())
      temp_dots <- dots[which(grepl("^mhclust__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^mhclust__", "", names(temp_dots), ignore.case = T)

      temp_dots <- c(list(g = ks[,1], quick = T, gintra = F), temp_dots) ## see shinysom for argument suggestion

      h <- do.call(mhca::mhclust, args = c(list(x = expr.select), temp_dots[which(names(temp_dots) %in% names(formals(mhca::mhclust))[-1])]))

      temp_dots <- dots[which(grepl("^cutree__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^cutree__", "", names(temp_dots), ignore.case = T)
      ks <- do.call(cbind, parallel::mclapply(temp_dots[["k"]], function(x) stats::cutree(tree = h, k = x), mc.cores = mc.cores))

      if (!is.null(metaclustering.on)) {
        ks <- apply(ks, 2, function (x) x[map[["mapping"]][,1]])
      }
      # make sure that cluster 1 is the largest and so on
      ks <- .cluster_ordering(ks = ks)

      colnames(ks) <- paste0("mhclust_", temp_dots[["k"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
    }, error = function(err) {
      message("run.mhclust with error: ", err)
    })
  }

  if (run.flowClust) {
    tryCatch({
      temp_dots <- dots[which(grepl("^flowClust__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^flowClust__", "", names(temp_dots), ignore.case = T)

      message("Finding clusters with flowClust and parallel::mclapply using ", mc.cores, " cores.\nStart: ", Sys.time())
      ks <- do.call(cbind, parallel::mclapply(temp_dots[["K"]], function(x) {
        suppressMessages(do.call(flowClust::flowClust, args = c(list(x = expr.clust, K = x), temp_dots[which(names(temp_dots) != "K")]))@label)
      }, mc.cores = mc.cores))

      if (!is.null(metaclustering.on)) {
        ks <- apply(ks, 2, function (x) x[map[["mapping"]][,1]])
      }

      # make sure that cluster 1 is the largest and so on
      ks <- .cluster_ordering(ks = ks)

      colnames(ks) <- paste0("flowClust_", temp_dots[["K"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("run.flowClust with error: ", err)
    })
  }

  if (run.MUDAN) {
    tryCatch({
      temp_dots <- dots[which(grepl("^MUDAN__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^MUDAN__", "", names(temp_dots), ignore.case = T)

      message("Finding clusters with MUDAN.\nStart: ", Sys.time())

      ks <- do.call(cbind, parallel::mclapply(temp_dots[["k"]], function(x) {
        ## documentation is wrong (mat: cells as rows and features as cols!)
        do.call(MUDAN::getComMembership, args = c(list(mat = expr.clust, k = x), temp_dots[which(names(temp_dots) != "k")]))
      }, mc.cores = mc.cores))

      if (!is.null(metaclustering.on)) {
        ks <- apply(ks, 2, function (x) x[map[["mapping"]][,1]])
      }

      # make sure that cluster 1 is the largest and so on
      ks <- .cluster_ordering(ks = ks)

      colnames(ks) <- paste0("MUDAN_", temp_dots[["k"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      message("End: ", Sys.time())
    }, error = function(err) {
      message("run.MUDAN with error: ", err)
    })
  }

  ### optionally run MUDAN::clusterBasedBatchCorrect here
  ## actually though harmony performs something similar with multiple iterations: https://portals.broadinstitute.org/harmony/articles/quickstart.html
  ## MUDAN advertises harmony: http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/mudan.html

  if (!is.logical(run.lda)) {
    # optimized embedding with LDA, so embedding biased by detected clusters. is that cheating or just the best dimensions to separate clusters on a map?!
    # https://highdemandskills.com/lda-clustering/
    ldam <- MASS::lda(com ~ ., data = as.data.frame(cbind(expr.select, com = dim.red.data[, run.lda])))
    expr.select <- stats::predict(ldam, as.data.frame(cbind(expr.select, com = dim.red.data[, run.lda])))$x
  }


  ## ---- dim reds -------
  if (run.umap) {
    temp_dots <- dots[which(grepl("^UMAP__", names(dots), ignore.case = T))]
    names(temp_dots) <-gsub("^UMAP__", "", names(temp_dots), ignore.case = T)

    if (!any(grepl("metric", names(temp_dots)), ignore.case = T)) {
      message("UMAP metric set to 'cosine' by default.")
      temp_dots <- c(temp_dots, metric = "cosine")
    }

    message("Calculating UMAP.\nStart: ", Sys.time())
    if (any(grepl("n_neighbors", names(temp_dots), ignore.case = T))) {
      umap.dims <- do.call(cbind, parallel::mclapply(temp_dots[["n_neighbors"]], function(z) {
        set.seed(seed)
        out <- do.call(uwot::umap, args = c(list(X = expr.select, verbose = F, n_neighbors = z),temp_dots[which(names(temp_dots) != "n_neighbors")]))
        colnames(out) <- c(paste0("UMAP_1_", z), paste0("UMAP_2_", z))
        return(out)
      }, mc.cores = mc.cores))
    } else {
      if (!any(grepl("verbose", names(temp_dots)), ignore.case = T)) {
        temp_dots <- c(temp_dots, verbose = T)
      }
      set.seed(seed)
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

    message("Calculating tSNE.\nStart: ", Sys.time())
    if (any(grepl("perplexity", names(temp_dots), ignore.case = T))) {
      tsne.dims <- do.call(cbind, parallel::mclapply(temp_dots[["perplexity"]], function(z) {
        set.seed(seed)
        out <- do.call(Rtsne::Rtsne, args = c(list(X = expr.select, verbose = F, perplexity = z), temp_dots[which(names(temp_dots) != "perplexity")]))$Y
        colnames(out) <- c(paste0("tSNE_1_", z), paste0("tSNE_2_", z))
        return(out)
      }, mc.cores = mc.cores))
    } else {
      set.seed(seed)
      tsne.dims <- do.call(Rtsne::Rtsne, args = c(list(X = expr.select, verbose = T), temp_dots))$Y
    }
    if (!any(grepl("perplexity", names(temp_dots), ignore.case = T)) || length(temp_dots[["perplexity"]]) == 1) {
      colnames(tsne.dims) <- c("tSNE_1", "tSNE_2")
    }
    message("End: ", Sys.time())
  }

  ## this part of code is a duplicate, for now
  ## only here lda applies if run.lda = TRUE
  if (run.som && (is.null(metaclustering.on) || metaclustering.on != "SOM")) {
    message("Calculating SOM.\nStart: ", Sys.time())
    temp_dots <- dots[which(grepl("^SOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^SOM__", "", names(temp_dots), ignore.case = T)
    set.seed(seed)
    map <- do.call(EmbedSOM::SOM, args = c(list(data = expr.select), temp_dots))

    temp_dots <- dots[which(grepl("^EmbedSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^EmbedSOM__", "", names(temp_dots), ignore.case = T)
    som.dims <- do.call(EmbedSOM::EmbedSOM, args = c(list(data = expr.select, map = map), temp_dots))
    colnames(som.dims) <- c("SOM_1", "SOM_2")
    message("End: ", Sys.time())
  }

  if (run.gqtsom && (is.null(metaclustering.on) || metaclustering.on != "GQTSOM")) {
    message("Calculating GQTSOM.\nStart: ", Sys.time())
    temp_dots <- dots[which(grepl("^GQTSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^GQTSOM__", "", names(temp_dots), ignore.case = T)
    set.seed(seed)
    map <- do.call(EmbedSOM::GQTSOM, args = c(list(data = expr.select), temp_dots))

    temp_dots <- dots[which(grepl("^EmbedSOM__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^EmbedSOM__", "", names(temp_dots), ignore.case = T)
    gqtsom.dims <- do.call(EmbedSOM::EmbedSOM, args = c(list(data = expr.select, map = map), temp_dots))
    colnames(gqtsom.dims) <- c("GQTSOM_1", "GQTSOM_2")
    message("End: ", Sys.time())
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
      marker <- NULL
    })
  } else {
    marker <- NULL
  }

  ## ---- write to disk -------
  if (!is.null(save.path) && !is.na(save.path)) {
    message("Writing files to disk.")
    t <- format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%Y%m%d_%H%M%S")
    dir.create(save.path, showWarnings = F)
    if ("rds" %in% save.to.disk) {
      saveRDS(list(df = dim.red.data, col_names = channel.desc_augment, marker = marker, flowframe = ff, pca = pca.result), file.path(save.path, paste0(t, "_dr_ff_list.rds")))
    }
    if ("fcs" %in% save.to.disk) {
      flowCore::write.FCS(ff, file.path(save.path, paste0(t, "_dr.fcs")))
    }
  }
  return(list(df = dim.red.data, col_names = channel.desc_augment, marker = marker, flowframe = ff, pca = pca.result))
}

.cluster_ordering <- function(ks) {
  if (methods::is(ks, "matrix")) {
    ks <- apply(ks, 2, function (x) {
      new_order <- stats::setNames(names(table(x)), nm = names(sort(table(x), decreasing = T)))
      return(as.numeric(new_order[as.character(x)]))
    })
  } else if (methods::is(ks, "integer")) {
    new_order <- stats::setNames(names(table(ks)), nm = names(sort(table(ks), decreasing = T)))
    ks <- as.numeric(new_order[as.character(ks)])
  } else {
    new_order <- stats::setNames(names(table(ks)), nm = names(sort(table(ks), decreasing = T)))
    ks <- unname(new_order[as.character(ks)])
  }

  return(ks)
}

.calc.pairwise.cluster.marker <- function(dat, cluster, levels = NULL, mc.cores = 1) {
  mc.cores <- min(mc.cores, parallel::detectCores() - 1)

  dat_split <- split_mat(x = dat, f = cluster)

  if (is.null(levels)) {
    levels <- sort(unique(cluster))
  }
  all_pairs <- utils::combn(levels, 2, simplify = T)

  # x = stats::setNames(as.character(all_pairs[1,]), nm = paste0(all_pairs[1,], "_____", all_pairs[2,]))
  dplyr::bind_rows(parallel::mcmapply(x = as.character(all_pairs[1,]), y = as.character(all_pairs[2,]), function(x,y) {
    'tryCatch({
      out <- suppressWarnings(matrixTests::col_wilcoxon_twosample(dat_split[[x]],
                                                                  dat_split[[y]])) %>%
        dplyr::select(pvalue) %>%
        tibble::rownames_to_column("channel")
    }, error=function(err) {
      message("Ran matrixTests::col_wilcoxon_twosample with error in level : ", x, " vs ", y, ": ")
      message(err)
      message("Trying presto::wilcoxauc.")
      out <-
        presto::wilcoxauc(cbind(t(dat_split[[x]]),t(dat_split[[y]])), c(rep("y", length(which(as.character(cluster) == x))),
                                                                        rep("z", length(which(as.character(cluster) == y))))) %>%
        dplyr::filter(group == "y") %>%
        dplyr::select(feature, pval) %>%
        dplyr::rename("pvalue" = pval, "channel" = feature)
    }, warning = function(war) {
      message("Ran matrixTests::col_wilcoxon_twosample with warning in level : ", x, " vs ", y, ": ")
      message(war)
      message("Trying presto::wilcoxauc.")
      out <-
        presto::wilcoxauc(cbind(t(dat_split[[x]]),t(dat_split[[y]])), c(rep("y", length(which(as.character(cluster) == x))),
                                                                        rep("z", length(which(as.character(cluster) == y))))) %>%
        dplyr::filter(group == "y") %>%
        dplyr::select(feature, pval) %>%
        dplyr::rename("pvalue" = pval, "channel" = feature)
    })'

    out <-
      presto::wilcoxauc(X = cbind(t(dat_split[[x]]),t(dat_split[[y]])), y = c(rep("y", length(which(as.character(cluster) == x))),
                                                                              rep("z", length(which(as.character(cluster) == y))))) %>%
      dplyr::filter(group == "y") %>%
      dplyr::select(feature, pval) %>%
      dplyr::rename("pvalue" = pval, "channel" = feature)

    out[,"mean_1"] <- round(matrixStats::colMeans2(dat_split[[x]]), 2)
    out[,"mean_2"] <- round(matrixStats::colMeans2(dat_split[[y]]), 2)
    out[,"mean_diff"] <- round(out[,"mean_1"] - out[,"mean_2"], 2)
    out[,"diptest_pvalue_1"] <- suppressWarnings(apply(dat_split[[x]], 2, function(z) diptest::dip.test(x = if(length(z) > 71999) {sample(z,71999)} else {z})[["p.value"]]))
    out[,"diptest_pvalue_2"] <- suppressWarnings(apply(dat_split[[y]], 2, function(z) diptest::dip.test(x = if(length(z) > 71999) {sample(z,71999)} else {z})[["p.value"]]))
    out[,"cluster_1"] <- as.character(x) #sapply(strsplit(out$cluster12, "_____"), "[", 1, simplify = T)
    out[,"cluster_2"] <- as.character(y) #sapply(strsplit(out$cluster12, "_____"), "[", 2, simplify = T)
    out[,"diff_sign"] <- ifelse(out[,"mean_diff"] == 0, "+/-", ifelse(out[,"mean_diff"] > 0, "+", "-"))

    #cluster_sizes <- utils::stack(table(cluster)) %>% dplyr::mutate(ind = as.character(ind))
    out <-
      out %>%
      #dplyr::left_join(cluster_sizes, by = c("cluster_1" = "ind")) %>%
      #dplyr::rename("cluster_1_events" = "values") %>%
      #dplyr::left_join(cluster_sizes, by = c("cluster_2" = "ind")) %>%
      #dplyr::rename("cluster_2_events" = "values") %>%
      #, cluster_1_events, cluster_2_events
      dplyr::select(channel, cluster_1, cluster_2, pvalue, mean_1, mean_2, mean_diff, diff_sign, diptest_pvalue_1, diptest_pvalue_2) %>%
      dplyr::arrange(pvalue)
    return(out)
  }, mc.cores = mc.cores, SIMPLIFY = F))

  out$cluster_1 <- factor(out$cluster_1, levels = levels)
  out$cluster_2 <- factor(out$cluster_2, levels = levels)
  return(out)
}

.calc.global.cluster.marker <- function(dat, cluster, levels = NULL, mc.cores = 1) {

  #method = c("presto", "matrixTests")
  #method <- match.arg(method, c("presto", "matrixTests"))

  mc.cores <- min(mc.cores, parallel::detectCores() - 1)

  # cluster is ident for each row in dat
  if (nrow(dat) != length(cluster)) {
    stop("nrow(dat) != length(cluster)")
  }
  if (is.null(levels)) {
    levels <- sort(unique(cluster))
  }

  #levels_out <<- levels
  #dat_out <<- dat
  #cluster_out <<- cluster

  ## try matrixStats first and on error run presto which requires transposation though
  ## matrixStats caused an error once (Integer Overflow with large matrices (approx. 1e6 cells as initial input))
  out <- dplyr::bind_rows(parallel::mclapply(levels, function(x) {
    ' tryCatch({
      out <- matrixTests::col_wilcoxon_twosample(dat[which(cluster == x),,drop = F],
                                                 dat[which(cluster != x),,drop = F]) %>%
        dplyr::select(pvalue) %>%
        tibble::rownames_to_column("channel")
    }, error=function(err) {
      message("Ran matrixTests::col_wilcoxon_twosample with error in level : ", x, ": ")
      message(err)
      message("Trying presto::wilcoxauc.")
      out <-
        presto::wilcoxauc(cbind(t(dat[which(cluster == x),,drop = F]),
                                t(dat[which(cluster != x),,drop = F])), c(rep("y", length(which(cluster == x))),
                                                                          rep("z", length(which(cluster != x))))) %>%
        dplyr::filter(group == "y") %>%
        dplyr::select(feature, pval) %>%
        dplyr::rename("pvalue" = pval, "channel" = feature)
    }, warning = function(war) {
      message("Ran matrixTests::col_wilcoxon_twosample with warning in level : ", x, ": ")
      message(war)
      message("Trying presto::wilcoxauc.")
      out <-
        presto::wilcoxauc(cbind(t(dat[which(cluster == x),,drop = F]),
                                t(dat[which(cluster != x),,drop = F])), c(rep("y", length(which(cluster == x))),
                                                                          rep("z", length(which(cluster != x))))) %>%
        dplyr::filter(group == "y") %>%
        dplyr::select(feature, pval) %>%
        dplyr::rename("pvalue" = pval, "channel" = feature)
    })'


    out <-
      presto::wilcoxauc(X = cbind(t(dat[which(cluster == x),,drop = F]),
                                  t(dat[which(cluster != x),,drop = F])), y = c(rep("y", length(which(cluster == x))),
                                                                                rep("z", length(which(cluster != x))))) %>%
      dplyr::filter(group == "y") %>%
      dplyr::select(feature, pval) %>%
      dplyr::rename("pvalue" = pval, "channel" = feature)

    out[,"mean_cluster"] <- round(matrixStats::colMeans2(dat[which(cluster == x),,drop = F]), 2)
    out[,"mean_not_cluster"] <- round(matrixStats::colMeans2(dat[which(cluster != x),,drop = F]), 2)
    out[,"mean_diff"] <- round(out[,"mean_cluster"] - out[,"mean_not_cluster"], 2)
    out[,"diptest_pvalue_cluster"] <- suppressWarnings(apply(dat[which(cluster == x),,drop = F], 2, function(z) diptest::dip.test(x = if(length(z) > 71999) {sample(z,71999)} else {z})[["p.value"]]))
    out[,"diptest_pvalue_notcluster"] <- suppressWarnings(apply(dat[which(cluster != x),,drop = F], 2, function(z) diptest::dip.test(x = if(length(z) > 71999) {sample(z,71999)} else {z})[["p.value"]]))
    out[,"cluster"] <- as.character(x)
    out[,"diff_sign"] <- ifelse(out[,"mean_diff"] == 0, "+/-", ifelse(out[,"mean_diff"] > 0, "+", "-"))
    out <-
      out %>%
      #dplyr::left_join(utils::stack(table(cluster)) %>% dplyr::mutate(ind = as.character(ind)), by = c("cluster" = "ind")) %>%
      #dplyr::rename("cluster_events" = "values") %>%
      #cluster_events
      dplyr::select(channel, cluster, pvalue, mean_cluster, mean_not_cluster, mean_diff, diff_sign, diptest_pvalue_cluster, diptest_pvalue_notcluster) %>%
      dplyr::arrange(pvalue)
    return(out)
  }, mc.cores = mc.cores))

  out$cluster <- factor(out$cluster, levels = levels)
  return(out)
}



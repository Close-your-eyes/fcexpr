#' Title
#'
#'
#'  logicle trans: 2006 - Parks_A New Logicle Display Method Avoids Deceptive Effects of Logarithmic Scaling for Low Signals and Compensated Data
#'
#' @param ff.list
#' @param channels
#' @param channel.desc
#' @param add.sample.info named list of additional channels to identify samples or group them in flowjo;
#' e.g.: add.sample.info = list(condition = c(1,2,3,1,2,3,1,2,3), donor = c(1,1,1,2,2,2,3,3,3))
#' @param scale.whole
#' @param scale.samples
#' @param run.harmony
#' @param harmony.subsets
#' @param run.pca
#' @param run.umap
#' @param run.tsne
#' @param run.louvain
#' @param run.kmeans
#' @param run.leiden
#' @param run.hclust
#' @param run.flowClust
#' @param run.MUDAN
#' @param MUDAN.k
#' @param n.pca.dims
#' @param umap.n_neighbors
#' @param umap.n_epochs
#' @param umap.metric
#' @param tsne.theta
#' @param louvain_cluster_resolutions
#' @param leiden_cluster_resolutions
#' @param n_cluster
#' @param extra.cols
#' @param mc.cores
#' @param save.RDS
#' @param wd
#' @param save.path
#' @param exclude.extra.channels
#' @param write.scaled.channels.to.FCS
#' @param check.channels
#' @param umap.ret_model
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
dr_to_fcs <- function(ff.list,
                      channels = NULL, # names vector
                      #channel.desc = NULL,
                      timeChannel = "Time",
                      add.sample.info,
                      scale.whole = "z.score",
                      scale.samples = "none",
                      run.harmony = F,
                      harmony.subsets,
                      run.pca = F,
                      run.umap = T,
                      run.tsne = F,
                      run.louvain = T,
                      run.kmeans = F,
                      run.leiden = F,
                      run.hclust = F,
                      run.flowClust = F,
                      run.MUDAN = F,
                      MUDAN.k = 30,
                      n.pca.dims = 0,
                      umap.n_neighbors = 15,
                      umap.n_epochs = 300,
                      umap.metric = "cosine",
                      tsne.theta = 0.1,
                      louvain_cluster_resolutions = c(seq(0.3,0.5,0.1)),
                      leiden_cluster_resolutions = c(seq(0.3,0.5,0.1)),
                      n_cluster = c(5:10),
                      extra.cols,
                      mc.cores = 1,
                      save.to.disk = c("fcs", "rds"),
                      save.path = file.path(getwd(), paste0(substr(gsub("\\.", "", make.names(as.character(Sys.time()))), 2, 15), "_FCS_dr")),
                      exclude.extra.channels = "FSC|SSC|Time|cluster.id",
                      write.scaled.channels.to.FCS = T,
                      check.channels = T,
                      umap.ret_model = F,
                      ...) {

  if (!requireNamespace("Rtsne", quietly = T)) utils::install.packages("Rtsne")
  if (run.umap && !requireNamespace("uwot", quietly = T)) utils::install.packages("uwot")
  if (run.leiden && !requireNamespace("Seurat", quietly = T)) utils::install.packages("Seurat")
  if (!requireNamespace("parallel", quietly = T)) utils::install.packages("parallel")
  if (run.leiden && !requireNamespace("leiden", quietly = T)) utils::install.packages("leiden")
  if (!requireNamespace("Biobase", quietly = T)) BiocManager::install("Biobase")
  if (run.flowClust && !requireNamespace("flowClust", quietly = T)) BiocManager::install("flowClust")
  if (!requireNamespace("devtools", quietly = T)) utils::install.packages("devtools")
  if (run.harmony && !requireNamespace("run.harmony", quietly = T)) devtools::install_github("immunogenomics/harmony")
  if (run.MUDAN && !requireNamespace("MUDAN", quietly = T)) devtools::install_github("JEFworks/MUDAN")

  if (!"logicle" %in% names(ff.list)) stop("logicle has to be in ff.list.")
  if (umap.ret_model && length(umap.n_neighbors) > 1) stop("UMAP can only be returned if length(umap.n_neighbors) == 1. Set umap.ret_model = F or provide only one value for umap.n_neighbors")
  if (n.pca.dims < 0) run.pca <- F
  save.to.disk <- match.arg(save.to.disk, c("fcs", "rds"), several.ok = T)

  if (length(ff.list) == 1 && names(ff.list) == "logicle" && check.channels) {
    print("ff.list only contains logicle. FSC and SSC and Time are removed from exclude.extra.channels. If not desired like this, set check.channels = F.")
    exclude.extra.channels <- gsub("FSC\\|", "", exclude.extra.channels)
    exclude.extra.channels <- gsub("SSC\\|", "", exclude.extra.channels)
    exclude.extra.channels <- gsub("Time\\|", "", exclude.extra.channels)
  }

  mc.cores <- min(mc.cores,  parallel::detectCores()-1)
  #if (!is.null(channel.desc) && length(channel.desc) != length(channels)) {stop("Please provide a channel description for every channel selected. Use a '-' for empty channels.")}

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
    if (!is.list(add.sample.info)) {stop("add.sample.info has to be a list.")}
    if (is.null(names(add.sample.info))) {stop("add.sample.info has to have names. These names become channel names in the FCS file.")}
    if (!all(unlist(lapply(add.sample.info, function(x) is.numeric(x))))) {stop("Please only provide numeric values as additional sample infos.")}
    if (any(unlist(lapply(add.sample.info, function(x) is.na(x))))) {stop("NA found in sample infos.")}
  }

  if (!missing(add.sample.info)) {
    if (!all(unlist(lapply(add.sample.info, function(x) length(x) == length(ff.list[[1]]))))) {stop(paste0("Length of each additional sample information has to match the length of selected samples, which is: ", length(ff.list[[1]]), "."))}
  }

  scale.samples <- switch(match.arg(scale.samples, c("z.score", "min.max", "none")),
                          z.score = scale,
                          min.max = min.max.normalization,
                          none = function(x) {return(x)})

  scale.whole <- switch(match.arg(scale.whole, c("z.score", "min.max", "none")),
                        z.score = scale,
                        min.max = min.max.normalization,
                        none = function(x) {return(x)})

  # check if channel names and desc are equal
  .check.ff.list(ff.list = ff.list)


  channels <- .get.channels(ff = ff.list[["logicle"]][[1]],
                            timeChannel = timeChannel,
                            channels = channels,
                            return = "inds")
  browser()
  # set channel.desc if provided; which is used later then
  if (!is.null(channel.desc)) {
    for (i in seq_along(ff.list[["logicle"]])) {
      ff.list[["logicle"]][[i]]@parameters@data[["desc"]][channels] <- channel.desc
    }
  }

  expr.select <- scale.whole(do.call(rbind, lapply(ff.list[["logicle"]], function(x) {scale.samples(flowCore::exprs(x)[,channels.index])})))

  # run harmony
  # if harmony is run with do_pca = T a subsequent pca should not be computed
  if (run.harmony) {
    row.ident <- rep(seq_along(ff.list[[1]]), sapply(ff.list[[1]], function(x) {nrow(flowCore::exprs(x))}))
    if (missing(harmony.subsets)) {
      if (!exists("meta_data")) {
        meta_data <- row.ident
      }
      expr.select <- harmony::HarmonyMatrix(expr.select, meta_data = meta_data, ...)
    } else {
      if (!exists("meta_data")) {
        meta_data <- lapply(harmony.subsets, function(x) {
          row.ident[which(row.ident %in% x)]
        })
      }
      for (i in seq_along(harmony.subsets)) {
        expr.select[which(row.ident %in% harmony.subsets[[i]]),] <- harmony::HarmonyMatrix(expr.select[which(row.ident %in% harmony.subsets[[i]]),], meta_data = meta_data[[i]], ...)
      }
    }
  }

  pca.result <- NULL
  if (run.pca) {
    if (n.pca.dims == 0) {n.pca.dims <- ncol(expr.select) - 1} else if (n.pca.dims > 0) {n.pca.dims <- min(ncol(expr.select), n.pca.dims)}
    print(paste0("Calculating PCA. Start: ", Sys.time()))
    pca.result <- prcomp(expr.select, scale. = F, center = F)
    pca.matrix <- pca.result[["x"]]
    expr.select <- pca.matrix[,1:n.pca.dims]
    print(paste0("Done. ", Sys.time()))
  }

  umap.model <- NULL
  if (run.umap) {
    print(paste0("Calculating UMAP. Start: ", Sys.time()))
    if (length(umap.n_neighbors) > 1) {
      print(paste0("Start: ", Sys.time()))
      umap.dims <- bind_cols(mclapply(umap.n_neighbors, function(z) {
        out <- uwot::umap(expr.select, verbose = F, metric = umap.metric, n_neighbors = z, n_epochs = umap.n_epochs, ret_model = F)
        colnames(out) <- c(paste0("UMAP_1_", z), paste0("UMAP_2_", z))
        return(as.data.frame(out))
      }, mc.cores = mc.cores))
      print(paste0("End: ", Sys.time()))
    } else {
      umap.model <- uwot::umap(expr.select, verbose = T, metric = umap.metric, n_neighbors = umap.n_neighbors, n_epochs = umap.n_epochs, ret_model = umap.ret_model)
      if (!umap.ret_model) {
        umap.dims <- as.data.frame(umap.model)
      } else {
        umap.dims <- as.data.frame(umap.model[["embedding"]])
      }
      colnames(umap.dims) <- c(paste0("UMAP_1_", umap.n_neighbors), paste0("UMAP_2_", umap.n_neighbors))
    }
    print(paste0("Done. ", Sys.time()))
  }

  if (run.tsne) {
    print(paste0("Calculating tSNE. Start: ", Sys.time()))
    tsne <- Rtsne(expr.select, verbose = T, pca = F, normalize = F, theta = tsne.theta)
    tsne.dims <- tsne$Y
    colnames(tsne.dims) <- c("tSNE_1", "tSNE_2")
    print(paste0("Done. ", Sys.time()))
  }

  # prepare matrix for FCS file
  expr.logicle <- do.call(rbind,lapply(ff.list[["logicle"]], function(x) {flowCore::exprs(x)}))
  expr.logicle <- expr.logicle[,which(!grepl(exclude.extra.channels, colnames(expr.logicle)))]
  colnames(expr.logicle) <- paste0(colnames(expr.logicle), "_logicle")

  if ("inverse" %in% names(ff.list)) {
    dim.red.data <- do.call(cbind, list(do.call(rbind,lapply(ff.list[["inverse"]], function(x) {flowCore::exprs(x)})), expr.logicle, ident = rep(1:length(ff.list[["logicle"]]), sapply(ff.list[["logicle"]], nrow))))
  } else {
    dim.red.data <- do.call(cbind, list(expr.logicle, ident = rep(1:length(ff.list[["logicle"]]), sapply(ff.list[["logicle"]], nrow))))
  }

  if (write.scaled.channels.to.FCS) {
    scaled.expr <- scale.whole(do.call(rbind, lapply(ff.list[["logicle"]], function(x) {scale.samples(flowCore::exprs(x)[,channels.index])})))
    scaled.expr <- scaled.expr[,which(!grepl(exclude.extra.channels, colnames(scaled.expr)))]
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

  # find communities (clusters)
  tryCatch(
    if (run.louvain) {
      print(paste0("Finding clusters with Seurats implementation of the Louvain algorithm and mclapply using ", mc.cores, " cores. Start: ", Sys.time()))
      rownames(expr.select) <- 1:nrow(expr.select)
      snn <- Seurat::FindNeighbors(expr.select, annoy.metric = "cosine")
      clust_idents <- as.matrix(do.call(cbind, mclapply(louvain_cluster_resolutions, function (x) {
        as.data.frame(apply(Seurat::FindClusters(snn$snn, resolution = x, verbose = T, algorithm = 1), 2, as.numeric))
      }, mc.cores = mc.cores)))
      colnames(clust_idents) <- paste0("louvain.", colnames(clust_idents))
      dim.red.data <- do.call(cbind, list(dim.red.data, clust_idents))
      print(apply(clust_idents, 2, function(x) length(unique(x))))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) print("run.louvain with error"))

  tryCatch(
    if (run.leiden) {
      print(paste0("Finding clusters leiden algorithm and mclapply using ", mc.cores, " cores. Start: ", Sys.time()))
      rownames(expr.select) <- 1:nrow(expr.select)
      if (!exists("snn")) {
        snn <- Seurat::FindNeighbors(expr.select, annoy.metric = "cosine")
      }
      clust_idents <- as.matrix(do.call(cbind, mclapply(leiden_cluster_resolutions, function (x) {
        leiden::leiden(snn$snn, resolution_parameter = x, n_iterations = 10, seed = 0)
      }, mc.cores = mc.cores)))
      colnames(clust_idents) <- paste0("leiden.res.", leiden_cluster_resolutions)
      dim.red.data <- do.call(cbind, list(dim.red.data, clust_idents))
      print(apply(clust_idents, 2, function(x) length(unique(x))))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) print("run.leiden with error"))

  tryCatch(
    if (run.kmeans) {
      print(paste0("Finding clusters with kmeans and mclapply using ", mc.cores, " cores. Start: ", Sys.time()))
      ks <- do.call(cbind, mclapply(n_cluster, function(x) {
        kmeans(expr.select, centers = x)$cluster
      }, mc.cores = mc.cores))
      colnames(ks) <- paste0("kmean.", n_cluster)
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) print("run.kmeans with error"))

  tryCatch(
    if (run.hclust) {
      print(paste0("Finding clusters with hclust. Start: ", Sys.time()))
      h <- stats::hclust(dist(expr.select))
      ks <- do.call(cbind, lapply(n_cluster, function(x) {
        stats::cutree(h, k = x)
      }))
      colnames(ks) <- paste0("hclust.", n_cluster)
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) print("run.hclust with error"))

  tryCatch(
    if (run.flowClust) {
      print(paste0("Finding clusters with flowClust and mclapply using ", mc.cores, " cores. Start: ", Sys.time()))
      ks <- do.call(cbind, mclapply(n_cluster, function(x) {
        flowClust::flowClust(expr.select, K = x)@label
      }, mc.cores = mc.cores))
      colnames(ks) <- paste0("flowClust.", n_cluster)
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) print("run.flowClust with error"))


  tryCatch(
    if (run.MUDAN) {
      print(paste0("Finding clusters with MUDAN.",  Sys.time()))
      dim.red.data <- do.call(cbind, list(dim.red.data, data.frame(MUDAN = as.numeric(MUDAN::getComMembership(expr.select, k = MUDAN.k)))))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) print("run.flowClust with error"))


  '  if (run.flowsom) {
    #fsom <- FlowSOM::ReadInput(flowSet(ff.list[["logicle"]]))
    #fsom <- FlowSOM::BuildSOM(fsom = fsom, colsToUse = channels.index)
    #fsom <- FlowSOM::BuildMST(fsom = fsom)
    #metacl <- FlowSOM::MetaClustering(fsom$data, "metaClustering_consensus", max = max(flowsom_n.cluster))
    browser()
    print(paste0("Finding clusters with flowSOMs hierachical clustering and mclapply using ", mc.cores, " cores. Start: ", Sys.time()))
    ks <- do.call(cbind, mclapply(n_cluster, function(x) {
      FlowSOM::metaClustering_consensus(expr.select, k = x)
    }, mc.cores = mc.cores))

    library(pbapply)
    ks <- do.call(cbind, pblapply(n_cluster, function(x) {
      FlowSOM::metaClustering_consensus(expr.select, k = x)
      metacl <- FlowSOM::MetaClustering(expr.select, "metaClustering_consensus", max = max(n_cluster))
    }))

    colnames(ks) <- paste0("flowSOM.", n_cluster)
    dim.red.data <- do.call(cbind, list(dim.red.data, ks))
    print(paste0("End: ", Sys.time()))
  }'

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
    error=function(e) print("Addition of sample info failed due to an error. Check add.sample.info."))

  '  new_p <- parameters(ff)[1,]
  new_kw <- flowCore::keyword(ff)
  new_pars <- parameters(ff)
  for (z in seq_along(ref.ffs)) {
    new_p_number <- as.integer(dim(ff)[2]+z)
    rownames(new_p) <- c(paste0("$P", new_p_number))
    new_pars <- BiocGenerics::combine(new_pars, new_p)
    new_p_name <- paste0("cluster.id_", names(ref.ffs)[z])
    new_p_desc <- "signal.correction"
    new_pars@data$name[new_p_number] <- new_p_name
    new_pars@data$desc[new_p_number] <- new_p_desc

    new_kw["$PAR"] <- as.character(new_p_number)
    new_kw[paste0("$P",as.character(new_p_number),"N")] <- new_p_name
    new_kw[paste0("$P",as.character(new_p_number),"S")] <- new_p_desc
    new_kw[paste0("$P",as.character(new_p_number),"E")] <- "0,0"
    new_kw[paste0("$P",as.character(new_p_number),"G")] <- "1"
    new_kw[paste0("$P",as.character(new_p_number),"B")] <- new_kw["$P1B"]
    new_kw[paste0("$P",as.character(new_p_number),"R")] <- max(tarExprs[,new_p_name])
    new_kw[paste0("$P",as.character(new_p_number),"DISPLAY")] <- "LIN"
    new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- min(tarExprs[,new_p_name])
    new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- max(tarExprs[,new_p_name])
  }
  ff <- new("flowFrame", exprs = tarExprs, parameters = new_pars, description = new_kw)'

  # prepare channel desc
  name.desc <- setNames(ff.list[[1]][[1]]@parameters@data[["desc"]], ff.list[[1]][[1]]@parameters@data[["name"]])
  name.desc <- name.desc[which(!is.na(name.desc))]
  channel.desc <- rep("", ncol(dim.red.data))
  for (i in seq_along(name.desc)) {
    channel.desc[grep(names(name.desc)[i], colnames(dim.red.data))] <- name.desc[i]
  }

  # write FCS file
  metadata <- data.frame(name = colnames(dim.red.data), desc = channel.desc, stringsAsFactors = F)
  metadata$minRange <- apply(dim.red.data,2,min)
  metadata$maxRange <- apply(dim.red.data,2,max)

  # add meta data to flowFrame (FCS)
  flowFrame.description <- list(channels = paste(channels, collapse = ", "), sample.idents = paste(as.character(1:length(ff.list[["logicle"]])), collapse = ","))
  ff <- new("flowFrame", exprs = as.matrix(dim.red.data), parameters = Biobase::AnnotatedDataFrame(metadata), description = flowFrame.description) # exprs must be a matrix

  # save results

  if (!is.null(save.path) && !is.na(save.path)) {
    print("Writing files to disk.")
    t <- format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%Y%m%d_%H%M%S")
    dir.create(save.path, showWarnings = F)
    if ("rds" %in% save.to.disk) {
      saveRDS(list(table = dim.red.data, flowframe = ff, pca = pca.result, umap.model = umap.model), file.path(save.path, paste0(t, "_dr_ff_list.rds")))
    }
    if ("fcs" %in% save.to.disk) {
      flowCore::write.FCS(ff, file.path(save.path, paste0(t, "_dr.fcs")))
    }
  }
  return(dim.red.data)
}


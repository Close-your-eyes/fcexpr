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
                      channels = NULL, # named vector
                      timeChannel = "Time",
                      add.sample.info,
                      scale.whole = "z.score",
                      scale.samples = "none",
                      run.harmony = F,
                      run.pca = F,
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
                      n_cluster = c(5:10),
                      extra.cols,
                      mc.cores = 1,
                      save.to.disk = c("fcs", "rds"),
                      save.path = file.path(getwd(), paste0(substr(gsub("\\.", "", make.names(as.character(Sys.time()))), 2, 15), "_FCS_dr")),
                      exclude.extra.channels = "FSC|SSC|Time|cluster.id",
                      write.scaled.channels.to.FCS = T,
                      check.channels = T,
                      ...) {

  ## louvain, leiden to Seurat::FindNeighbors
  # lovain to Seurat::FindClusters
  # leider to leiden::leiden

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
  if (!"logicle" %in% names(ff.list)) stop("logicle transformed has to be in ff.list.")

  dots <- list(...)

  for (par in c("louvain", "leiden", "umap", "tsne", "som", "gqtsom", "harmony")) {
    if (any(grepl(paste0("^", par, "__"), names(dots), ignore.case = T)) && !eval(rlang::sym(paste0("run.", par)))) {
      message(paste0(par, " parameters provided in ... but ", "'run.",par," = F'."))
    }
  }

  if (run.harmony && !any(grepl("^harmony__meta_data", names(dots)))) {
    stop("When 'run.harmony = T' harmony__meta_data has to be provided in ..., see ?harmony::HarmonyMatrix.")
  }

  if (n.pca.dims < 0) run.pca <- F

  if (!is.null(save.to.disk)) {
    save.to.disk <- match.arg(save.to.disk, c("fcs", "rds"), several.ok = T)
  }


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

  channels <- .get.channels(ff = ff.list[["logicle"]][[1]], ## channel names from first ff
                            timeChannel = timeChannel,
                            channels = channels)

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

  expr.select <- scale.whole(do.call(rbind, lapply(ff.list[["logicle"]], function(x) scale.samples(flowCore::exprs(x)[,channels]))))

  # run harmony
  # if harmony is run with do_pca = T a subsequent pca should not be computed
  if (run.harmony) {
    temp_dots <- dots[which(grepl("^harmony__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^harmony__", "", names(temp_dots), ignore.case = T)
    expr.select <- do.call(harmony::HarmonyMatrix, args = c(list(data_mat = expr.select), temp_dots))
  }

  pca.result <- NULL
  if (run.pca) {
    if (n.pca.dims == 0) {n.pca.dims <- ncol(expr.select) - 1} else if (n.pca.dims > 0) {n.pca.dims <- min(ncol(expr.select)-1, n.pca.dims)}
    print(paste0("Calculating PCA. Start: ", Sys.time()))
    # https://slowkow.com/notes/pca-benchmark/
    #mat_irlba2 <- irlba::irlba(A = expr.select, nv = n.pca.dims)
    #mat_irlba2$x <- mat_irlba2$u %*% diag(mat_irlba2$d)
    pca.result <- stats::prcomp(expr.select, scale. = F, center = F)
    pca.matrix <- pca.result[["x"]]
    expr.select <- pca.matrix[,1:n.pca.dims]
    print(paste0("Done. ", Sys.time()))
  }

  if (run.umap) {
    temp_dots <- dots[which(grepl("^UMAP__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^UMAP__", "", names(temp_dots), ignore.case = T)

    if (!any(grepl("metric", names(temp_dots)), ignore.case = T)) {
      print("UMAP metric set to 'cosine' by default.")
      temp_dots <- c(temp_dots, metric = "cosine")
    }
    print(paste0("Calculating UMAP. Start: ", Sys.time()))
    if (any(grepl("n_neighbors", names(temp_dots), ignore.case = T))) {
      umap.dims <- do.call(cbind, parallel::mclapply(temp_dots[["n_neighbors"]], function(z) {
        out <- do.call(uwot::umap, args = c(list(X = expr.select, verbose = F, n_neighbors = z), temp_dots[which(names(temp_dots) != "n_neighbors")]))
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
    map <- do.call(EmbedSOM::SOM, args = c(list(data = expr.select), temp_dots))

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
    scaled.expr <- scale.whole(do.call(rbind, lapply(ff.list[["logicle"]], function(x) {scale.samples(flowCore::exprs(x)[,channels])})))
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
  if (run.som) {
    dim.red.data <- do.call(cbind, list(dim.red.data, som.dims))
  }
  if (run.gqtsom) {
    dim.red.data <- do.call(cbind, list(dim.red.data, gqtsom.dims))
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
      # Seurat::FindNeighbors ignores all 'wrong' arguments; suppress the warning though
      snn <- suppressWarnings(do.call(Seurat::FindNeighbors, args = c(list(object = expr.select), temp_dots)))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) {
      print("louvain / leiden: error in snn calculation.")
      run.louvain <- F
      run.leiden <- F
    }
  )

  tryCatch(
    if (run.louvain) {
      temp_dots <- dots[which(grepl("^louvain__", names(dots), ignore.case = T))]
      names(temp_dots) <- gsub("^louvain__", "", names(temp_dots), ignore.case = T)

      print(paste0("Finding clusters with Seurats implementation of the Louvain algorithm and parallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time()))
      if (any(grepl("resolution", names(temp_dots), ignore.case = T))) {
        temp_dots[["resolution"]] <- as.numeric(temp_dots[["resolution"]])
        if (any(is.na(temp_dots[["resolution"]]))) {
          warning("Provide numeric values for louvain__resolution! Non numeric elements are ignored.")
          temp_dots[["resolution"]] <- temp_dots[["resolution"]][which(!is.na(temp_dots[["resolution"]]))]
        }
        clust_idents <- do.call(cbind, parallel::mclapply(temp_dots[["resolution"]], function (x) {
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
    error=function(e) print("run.louvain with error")
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
        clust_idents <- do.call(cbind, parallel::mclapply(temp_dots[["resolution_parameter"]], function (x) {
          do.call(leiden::leiden, args = c(list(object = snn$snn, resolution_parameter = x), temp_dots[which(names(temp_dots) != "resolution_parameter")]))
        }, mc.cores = mc.cores))
      } else {
        clust_idents <- do.call(leiden::leiden, args = c(list(object = snn$snn), temp_dots))
      }
      colnames(clust_idents) <- paste0("leiden_", temp_dots[["resolution_parameter"]])
      dim.red.data <- do.call(cbind, list(dim.red.data, clust_idents))
      print(apply(clust_idents, 2, function(x) length(unique(x))))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) print("run.leiden with error")
  )

  ### work on the 4 below!

  tryCatch(
    if (run.kmeans) {
      print(paste0("Finding clusters with kmeans andparallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time()))
      ks <- do.call(cbind, parallel::mclapply(n_cluster, function(x) {
        stats::kmeans(expr.select, centers = x)$cluster
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
      ks <- do.call(cbind, parallel::mclapply(n_cluster, function(x) {
        stats::cutree(h, k = x)
      }, mc.cores = mc.cores))
      colnames(ks) <- paste0("hclust_", n_cluster)
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) print("run.hclust with error"))

  tryCatch(
    if (run.flowClust) {
      print(paste0("Finding clusters with flowClust andparallel::mclapply using ", mc.cores, " cores. Start: ", Sys.time()))
      ks <- do.call(cbind, parallel::mclapply(n_cluster, function(x) {
        flowClust::flowClust(expr.select, K = x)@label
      }, mc.cores = mc.cores))
      colnames(ks) <- paste0("flowClust_", n_cluster)
      dim.red.data <- do.call(cbind, list(dim.red.data, ks))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) print("run.flowClust with error"))

  tryCatch(
    if (run.MUDAN) {
      print(paste0("Finding clusters with MUDAN.",  Sys.time()))
      dim.red.data <- do.call(cbind, list(dim.red.data, data.frame(MUDAN = as.numeric(MUDAN::getComMembership(t(expr.select), k = MUDAN.k)))))
      print(paste0("End: ", Sys.time()))
    },
    error=function(e) print("run.MUDAN with error"))

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
      saveRDS(list(table = dim.red.data, flowframe = ff, pca = pca.result), file.path(save.path, paste0(t, "_dr_ff_list.rds")))
    }
    if ("fcs" %in% save.to.disk) {
      flowCore::write.FCS(ff, file.path(save.path, paste0(t, "_dr.fcs")))
    }
  }
  return(list(df = dim.red.data, desc = channel.desc))
}


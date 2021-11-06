pca.to.FSC <- function(file_path,
                       channels = NULL,
                       compensate = T,
                       compMat,
                       timeChannel,
                       logicle_trans = T,
                       n_pca_dims = NULL,
                       remove_outliers = F,
                       outlier_limit = 10,
                       output_folder = NULL) {

  library(factoextra)
  library(Biobase)

  ff_orig <- flowCore::read.FCS(file_path, truncate_max_range = F, emptyValue = F)

  if (compensate) {
    if (missing(compMat)) {
      compMat <- flowCore::keyword(ff_orig)[["SPILL"]]
      print("SPILL keyword used for compensation.")
    }
    ff <- flowCore::compensate(ff_orig, compMat)
  } else {
    ff <- ff_orig
    print("No compensation applied.")
  }

  if(!missing(timeChannel)) {
    if (!timeChannel %in% colnames(flowCore::exprs(ff))) {
      stop("timeChannel not found in exprs of flowFrame.")
    }
  } else {
    timeChannel <- flowCore:::findTimeChannel(ff)
    print(paste0("time channel detected: ", timeChannel))
  }

  if (is.null(channels)) {
    channels <- colnames(flowCore::exprs(ff))
    channels <- channels[which(channels != timeChannel)]
  } else {
    # print warning if channels not found
    # allow desc for selection as well
    #flowCore::pData(flowCore::parameters(ff))
    channels <- channels[which(channels %in% colnames(flowCore::exprs(ff)))]
  }
  if (length(channels) == 0) {
    stop("no channels matched to those in the flowFrame.")
  }

  if (logicle_trans) {
    ff <- lgcl_trfm_ff(ff, channels = channels)
  }

  exprs <- flowCore::exprs(ff)
  exprs <- exprs[,which(colnames(exprs) %in% channels)]

  # remove outliers
  if (remove_outliers) {
    ind.outlier <- unique(unlist(apply(exprs, 2, function(x) which(abs((abs(x - median(x)) / mad(x))) > outlier_limit))))
    exprs <- exprs[-c(ind.outlier),]
    flowCore::exprs(ff_orig) <- flowCore::exprs(ff_orig)[-ind.outlier,]
    print(paste0(length(ind.outlier), " outliers are removed."))
  }

  pca <- stats::prcomp(exprs, center = T, scale. = T)
  if (is.null(n_pca_dims)) {
    n_pca_dims <- ncol(pca[["x"]])
  }
  n_pca_dims <- min(n_pca_dims, ncol(pca[["x"]]))

  ## make ff
  new_p <- flowCore::parameters(ff)[1,]
  new_kw <- flowCore::keyword(ff)
  new_pars <- flowCore::parameters(ff)
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
  ff <- new("flowFrame", exprs = tarExprs, parameters = new_pars, description = new_kw)


  flowCore::parameters(ff_orig)
  flowCore::exprs(ff_orig) <- cbind(flowCore::exprs(ff_orig), pca[["x"]][,1:n_pca_dims])


  desc <- c(flowCore::parameters(ff_orig)@data[["desc"]], rep(NA, n_pca_dims))
  names(desc)[which(names(desc) == "")] <- NA
  desc <- desc[which(is.na(desc))] <- ""

  metadata <- data.frame(name = colnames(flowCore::exprs(ff_orig)), desc = desc)
  metadata$minRange <- apply(flowCore::exprs(ff),2,min)
  metadata$maxRange <- apply(flowCore::exprs(ff),2,max)

  ff <- new("flowFrame", exprs = flowCore::exprs(ff), parameters = AnnotatedDataFrame(metadata), description = c(flowCore::keyword(ff), list(channels.for.pca = channels)))

  write.FCS(ff, paste0(export.folder.path, str_replace(basename(file_path), ".fcs", ""), "_pca.fcs"))

  if (is.null(output_folder)) {
    dir.create(output_folder, recursive = T, showWarnings = F)
  }

  return(list(pca.data = pca.data, pca.plot = pca.plot))
}

file_path <- "/Users/vonskopnik/Documents/20190225_CMS_Blut/Compensation_FCS_files/20190225_PBMC-comp_cd8-apc_003.fcs"


lgcl_trfm_ff <- function(ff, channels = NULL) {

  if (is.null(channels)) {
    channels <- colnames(flowCore::exprs(ff))
    channels <- channels[which(channels != flowCore:::findTimeChannel(ff))]
  }

  trfms <- lapply(channels, function(z) {
    m <- 4.5
    lgcl <- NULL
    while(is.null(lgcl)) {
      lgcl <- tryCatch(flowCore::estimateLogicle(ff, z, m = m),
                       error = function(e) {
                         #print(m)
                         return(NULL)
                       }
      )
      m <- m + 0.1
    }
    return(lgcl)
  })


  for (i in seq_along(trfms)) {
    ff <- flowCore::transform(ff, trfms[[i]])
  }
  return(ff)
}


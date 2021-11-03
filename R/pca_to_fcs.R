pca.to.FSC <- function(file_path,
                       channels = NA,
                       compensate = T,
                       time.channel.name = "time",
                       n.dims.to.fcs = NA,
                       remove.extreme.outliers = F,
                       export.folder.path = NA) {

  library(factoextra)
  library(Biobase)

  ff <- flowCore::read.FCS(fcs.path, truncate_max_range = F)
  if (compensate) {
    ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
  }

  if (missing(channels)) {
    print("No channels selected. Only excluding putative time-channel.")
    channels <- colnames(as.matrix(flowCore::exprs(ff)))[which(tolower(colnames(as.matrix(flowCore::exprs(ff)))) != time.channel.name)]
  }
  print(paste0("Channels selected for PCA: ", paste(colnames(as.matrix(flowCore::exprs(ff)))[grepl(paste(c(channels), collapse = "|"), colnames(as.matrix(flowCore::exprs(ff))))], collapse = ", ")))

  # calc pca
  expr.data <- log10(apply(as.matrix(flowCore::exprs(ff))[,which(colnames(as.matrix(flowCore::exprs(ff))) %in% colnames(as.matrix(flowCore::exprs(ff)))[grepl(paste(c(channels), collapse = "|"), colnames(as.matrix(flowCore::exprs(ff))))])], 2, shift.to.positive))

  # remove outliers
  if (remove.extreme.outliers) {
    ind.outlier <- unique(unlist(apply(expr.data, 2, function(x) which( abs((abs(x - median(x)) / mad(x))) > 10)))) # 10 is chosen arbitrarily to remove very extreme values only - maybe worth to check again, somewhen
    expr.data <- expr.data[-c(ind.outlier),]
    flowCore::exprs(ff) <- flowCore::exprs(ff)[-c(ind.outlier),]
    print(paste0(length(ind.outlier), " outliers are removed."))
  }

  pca.data <- prcomp(expr.data, center = T, scale. = T)

  if (missing(n.dims.to.fcs)) {
    n.dims.to.fcs <- ncol(pca.data[["x"]])
  } else if (n.dims.to.fcs <= 0 | n.dims.to.fcs %% 1 != 0 ) {
    print(paste0("n.dims.to.fcs has to be a positive integer. All pca dims are added to fcs file."))
    n.dims.to.fcs <- ncol(pca.data[["x"]])
  } else if (n.dims.to.fcs > ncol(pca.data[["x"]])) {
    print(paste0("n.dims.to.fcs has to be smaller or equal to the number of channels.All pca dims are added to fcs file."))
    n.dims.to.fcs <- ncol(pca.data[["x"]])
  }

  flowCore::exprs(ff) <- cbind(flowCore::exprs(ff), pca.data[["x"]][,1:n.dims.to.fcs])

  if (is.na(export.folder.path)) {
    export.folder.path <- dirname(fcs.path)
  }
  if (!str_detect(export.folder.path, "/$")) {
    export.folder.path <- paste0(export.folder.path, "/")
  }

  desc <- c(flowCore::parameters(ff)@data[["desc"]], rep(NA, n.dims.to.fcs))
  names(desc)[which(names(desc) == "")] <- NA

  metadata <- data.frame(name = colnames(flowCore::exprs(ff)), desc = desc) %>% dplyr::mutate(desc = ifelse(is.na(desc), "", desc))
  metadata$minRange <- apply(flowCore::exprs(ff),2,min)
  metadata$maxRange <- apply(flowCore::exprs(ff),2,max)

  ff <- new("flowFrame", exprs = flowCore::exprs(ff), parameters = AnnotatedDataFrame(metadata), description = c(flowCore::keyword(ff), list(channels.for.pca = channels)))

  write.FCS(ff, paste0(export.folder.path, str_replace(basename(fcs.path), ".fcs", ""), "_pca.fcs"))

  return(list(pca.data = pca.data, pca.plot = pca.plot))
}

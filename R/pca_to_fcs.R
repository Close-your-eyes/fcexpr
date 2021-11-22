#' Add principle components (PCs) to FCS file
#'
#' A PCA is calculated on selected events and selected channels with stats::prcomp.
#' Resulting PCs are added as channels to a newly generated FCS file.
#'
#' @param file_path character, path to a fcs file
#' @param which_lines numeric vector, passed to flowCore::read.FCS(..., which.lines = which_lines); every line in a fcs file is one acquired event.
#' If you know your events (~lines) of interest (e.g. a subpopulation to select or outliers to exclude) these events may be selected.
#' @param channels character vector, which channels to use for PC calculation. Either channels names (e.g. v-450/50-F-A) or descriptions (e.g. CD3 or CD4-PECy7) may be provided. Also a mixture is possible. If not provided, all channels (scatter and fluorescence) except for the Time channel are used.
#' @param compensate logical, should compensation be applied before PC calculation
#' @param compMat matrix, optional; a compensation matrix to use for compensation. If not provided the SPILL argument of the fcs file will be used. If you have generated a compensation matrix in FlowJo see ?fcexpr::wsx_compMats_to_fcs in order to have it copied to fcs files.
#' @param timeChannel character, optional; name of the time channel. If not provided flowCore:::findTimeChannel() is used to derive the time channel.
#' @param logicle_trans logical, should the logical transformation (Parks, 2006, https://pubmed.ncbi.nlm.nih.gov/16604519/) be applied before PCA calculation. Recommendation: yes.
#' @param processed_channels_to_FCS logical, should the processed fluorescence intensities (compensation and/or logicle transformation) be saved as extra channels to the newly generated fcs file?
#' Respective channels are suffixed by _comp, _lgcl, or _comp_lgcl depend upon the selections above.
#' @param n_pca_dims numeric, the number of PCs to add to the newly generated fcs file. Default: all.
#' @param output_folder character, optional, path to a folder where to save the newly generated fcs file. Default is dirname(file_path).
#' @param new_file_suffix character, the suffix to add to the the newly generated fcs file. Default is _pca.
#'
#' @return list of pca-object and appended flowFrame which is also saved as fcs file
#' @export
#'
#' @examples
pca_to_fcs <- function(file_path,
                       which_lines = NULL,
                       channels = NULL,
                       compensate = T,
                       compMat,
                       timeChannel,
                       logicle_trans = T,
                       processed_channels_to_FCS = T,
                       n_pca_dims = NULL,
                       output_folder = NULL,
                       new_file_suffix = "pca") {

  if (requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }
  if (requireNamespace("BiocGenerics", quietly = T)){
    BiocManager::install("BiocGenerics")
  }

  if (!file.exists(file_path)) {
    stop(paste0(file_path, " not found."))
  }

  ff_orig <- flowCore::read.FCS(file_path, which.lines = which_lines, truncate_max_range = F, emptyValue = F)

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
    inds <- unique(which(flowCore::pData(flowCore::parameters(ff))$name %in% channels),
                   which(flowCore::pData(flowCore::parameters(ff))$desc %in% channels))
    notfound <- channels[intersect(which(!channels %in% flowCore::pData(flowCore::parameters(ff))$name),
                                   which(!channels %in% flowCore::pData(flowCore::parameters(ff))$desc))]
    if (length(notfound) > 0) {
      print(paste0(paste(notfound, collapse = ", "), " channels not found in flowFrame."))
    }
    channels <- flowCore::pData(flowCore::parameters(ff))$name[inds]
  }
  if (length(channels) == 0) {
    stop("no channels matched to those in the flowFrame.")
  }

  if (logicle_trans) {
    ff <- lgcl_trsfrm_ff(ff, channels = channels)
  }

  exprs <- flowCore::exprs(ff)
  exprs <- exprs[,which(colnames(exprs) %in% channels)]

  pca <- stats::prcomp(exprs, center = T, scale. = T)
  if (is.null(n_pca_dims)) {
    n_pca_dims <- ncol(pca[["x"]])
  }
  n_pca_dims <- min(n_pca_dims, ncol(pca[["x"]]))

  new_p <- flowCore::parameters(ff_orig)[1,]
  new_kw <- flowCore::keyword(ff_orig)
  new_pars <- flowCore::parameters(ff_orig)
  for (z in 1:n_pca_dims) {
    new_p_number <- as.integer(dim(ff_orig)[2]+z)
    rownames(new_p) <- c(paste0("$P", new_p_number))
    new_pars <- BiocGenerics::combine(new_pars, new_p)
    new_p_name <- paste0("PC", z)
    new_p_desc <- paste0("PCA_dim", z)
    flowCore::pData(new_pars)$name[new_p_number] <- new_p_name
    flowCore::pData(new_pars)$desc[new_p_number] <- new_p_desc

    new_kw["$PAR"] <- as.character(new_p_number)
    new_kw[paste0("$P",as.character(new_p_number),"N")] <- new_p_name
    new_kw[paste0("$P",as.character(new_p_number),"S")] <- new_p_desc
    new_kw[paste0("$P",as.character(new_p_number),"E")] <- "0,0"
    new_kw[paste0("$P",as.character(new_p_number),"G")] <- "1"
    new_kw[paste0("$P",as.character(new_p_number),"B")] <- new_kw["$P1B"]
    new_kw[paste0("$P",as.character(new_p_number),"R")] <- max(pca[["x"]][,z])
    new_kw[paste0("$P",as.character(new_p_number),"DISPLAY")] <- "LIN"
    new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- min(pca[["x"]][,z])
    new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- max(pca[["x"]][,z])
  }
  ff_new <- methods::new("flowFrame", exprs = cbind(flowCore::exprs(ff_orig), pca[["x"]][,1:n_pca_dims]), parameters = new_pars, description = new_kw)

  if (processed_channels_to_FCS && (compensate || logicle_trans)) {
    if (compensate && logicle_trans) {
      ext <- "_comp_lgcl"
    } else if (logicle_trans) {
      ext <- "_lgcl"
    } else if (compensate) {
      ext <- "_comp"
    }
    new_p <- flowCore::parameters(ff_new)[1,]
    new_kw <- flowCore::keyword(ff_new)
    new_pars <- flowCore::parameters(ff_new)
    for (z in 1:ncol(exprs)) {
      new_p_number <- as.integer(dim(ff_new)[2]+z)
      rownames(new_p) <- c(paste0("$P", new_p_number))
      new_pars <- BiocGenerics::combine(new_pars, new_p)
      new_p_name <- paste0(colnames(exprs)[z], ext)
      ppp <- flowCore::pData(flowCore::parameters(ff))
      new_p_desc <- ppp[which(ppp$name == colnames(exprs)[z]), "desc"]
      if (is.na(new_p_desc)) {
        new_p_desc <- ""
      }
      flowCore::pData(new_pars)$name[new_p_number] <- new_p_name
      flowCore::pData(new_pars)$desc[new_p_number] <- new_p_desc

      new_kw["$PAR"] <- as.character(new_p_number)
      new_kw[paste0("$P",as.character(new_p_number),"N")] <- new_p_name
      new_kw[paste0("$P",as.character(new_p_number),"S")] <- new_p_desc
      new_kw[paste0("$P",as.character(new_p_number),"E")] <- "0,0"
      new_kw[paste0("$P",as.character(new_p_number),"G")] <- "1"
      new_kw[paste0("$P",as.character(new_p_number),"B")] <- new_kw["$P1B"]
      new_kw[paste0("$P",as.character(new_p_number),"R")] <- max(exprs[,z])
      new_kw[paste0("$P",as.character(new_p_number),"DISPLAY")] <- ifelse(logicle_trans, "LIN", "LOG")
      new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- min(exprs[,z])
      new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- max(exprs[,z])
    }
    colnames(exprs) <- paste0(colnames(exprs), ext)
    ff_new <- methods::new("flowFrame", exprs = cbind(flowCore::exprs(ff_new), exprs), parameters = new_pars, description = new_kw)
  }

  if (!is.null(output_folder)) {
    dir.create(output_folder, recursive = T, showWarnings = F)
  } else {
    output_folder <- dirname(file_path)
  }

  flowCore::write.FCS(ff_new, file.path(output_folder, paste0(gsub("\\.fcs$", "", basename(file_path)), "_", new_file_suffix, ".fcs")))
  return(list(pca = pca, ff = ff_new))
}

lgcl_trsfrm_ff <- function(ff, channels = NULL) {

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


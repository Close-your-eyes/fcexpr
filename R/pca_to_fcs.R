#' Add principle components (PCs) to FCS file
#'
#' A PCA is calculated on selected events and selected channels with stats::prcomp.
#' Resulting PCs are added as channels to a newly generated FCS file.
#'
#' @param file character, path to a fcs file or a flowFrame read with flowCore::read.FCS()
#' @param which_lines numeric vector, which events to keep; every line in a fcs file is one acquired event.
#' If you know your events (~lines) of interest (e.g. a subpopulation to select or outliers to exclude) these events may be selected.
#' @param channels character vector, which channels to use for PC calculation. Either channels names (e.g. v-450/50-F-A) or descriptions (e.g. CD3, CD4-PECy7) may be provided. Also a mixture is possible. If not provided, all channels (scatter and fluorescence) except for the Time channel are used.
#' @param compensate logical, should compensation be applied before PC calculation
#' @param compMat matrix, optional; a compensation matrix to use for compensation. If not provided the SPILL argument of the fcs file will be used. If you have generated a compensation matrix in FlowJo see ?fcexpr::wsx_compMats_to_fcs in order to have it copied to fcs files.
#' @param timeChannel character, optional; name of the time channel. If not provided flowCore:::findTimeChannel() is used to derive the time channel.
#' @param logicle_trans logical, should the logical transformation (Parks, 2006, https://pubmed.ncbi.nlm.nih.gov/16604519/) be applied before PCA calculation. Recommendation: yes.
#' @param processed_channels_to_FCS logical, should the processed fluorescence intensities (compensation and/or logicle transformation) be saved as extra channels to the newly generated fcs file?
#' Respective channels are suffixed by _comp, _lgcl, or _comp_lgcl depend upon the selections above.
#' @param n_pca_dims numeric, the number of PCs to add to the newly generated fcs file. Default: all.
#' @param output_folder character, optional, path to a folder where to save the newly generated fcs file. Default is dirname(file).
#' @param new_file_suffix character, the suffix to add to the the newly generated fcs file. Default is _pca.; Set this one to NULL as well as output_folder to overwrite the original file.
#'
#' @return list of pca-object and appended flowFrame which is also saved as fcs file
#' @export
#'
#' @examples
#' \dontrun{
#' pf <- pca_to_fcs(file = "mypath/file.fcs")
#' }
pca_to_fcs <- function(file,
                       which_lines = NULL,
                       channels = NULL,
                       compensate = T,
                       compMat,
                       timeChannel = NULL,
                       logicle_trans = T,
                       processed_channels_to_FCS = T,
                       n_pca_dims = NULL,
                       output_folder = NULL,
                       new_file_suffix = "pca") {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }
  if (!requireNamespace("BiocGenerics", quietly = T)){
    BiocManager::install("BiocGenerics")
  }
'  if (!requireNamespace("irlba", quietly = T)){
    utils::install.packages("irlba")
  }'

  if (is.character(file)) {
    if (length(file) > 1) {
      warning(paste0("Please only provide one file path. Only the first element of file will be used: ", file[1], "."))
    }
    if (!file.exists(file)) {
      stop(paste0(file, " not found."))
    }
    # which.lines is slow - filter afterwards
    ff_orig <- flowCore::read.FCS(file,
                                  truncate_max_range = F,
                                  emptyValue = F)
  } else if (methods::is(ff, "flowFrame")) {
    ff_orig <- file
  } else {
    stop("file must be a path to a fcs file on disk or a flowFrame read with flowCore::read.fcs().")
  }

  if (!is.null(which_lines)) {
    if (!is.numeric(which_lines)) {
      warning("which_lines has to be a numeric vector. It will be ignored as it is now.")
    }
    ff_orig <- ff_orig[which_lines,]
  }

  if (compensate) {
    if (missing(compMat)) {
      compMat <- flowCore::keyword(ff_orig)[["SPILL"]]
      if (is.null(compMat)) {
        compMat <- flowCore::keyword(ff_orig)[["$SPILLOVER"]]
        if (is.null(compMat)) {
          stop("compMat could not be determined.")
        }
        message("$SPILLOVER keyword used for compensation.")
      } else {
        message("SPILL keyword used for compensation.")
      }
    }
    ff <- flowCore::compensate(ff_orig, compMat)
  } else {
    ff <- ff_orig
    message("No compensation applied.")
  }

  channels <- .get.channels(ff = ff,
                            timeChannel = timeChannel,
                            channels = channels)

  if (logicle_trans) {
    ff <- lgcl_trsfrm_ff(ff, channels = channels)
  }
  exprs <- flowCore::exprs(ff)
  exprs <- exprs[,which(colnames(exprs) %in% channels)]

  #https://slowkow.com/notes/pca-benchmark/#irlbairlba
  '  mat_irlba2 <- irlba::irlba(
    A      = t(mat),
    nv     = n_pcs,
    center = Matrix::rowMeans(mat),
    scale  = proxyC::rowSds(mat)
  )
  mat_irlba2$x <- mat_irlba2$u %*% diag(mat_irlba2$d)
  '
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
    output_folder <- dirname(file)
  }

  if (!is.null(new_file_suffix)) {
    flowCore::write.FCS(ff_new, file.path(output_folder, paste0(gsub("\\.fcs$", "", basename(file)), "_", new_file_suffix, ".fcs")))
  } else {
    flowCore::write.FCS(ff_new, file.path(output_folder, paste0(gsub("\\.fcs$", "", basename(file)), ".fcs")))
  }

  return(list(pca = pca, ff = ff_new))
}


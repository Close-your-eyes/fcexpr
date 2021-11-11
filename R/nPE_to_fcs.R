#' Convert fluorescence intensity (FI) to number of Photoelectrons (nPE)
#'
#' Channels of selected flow cytometers have been analyzed with a QuantiFlash device. By generation of very defined pulses of light the
#' fluorescence intensity (FI) measured by a PMT (detector in a flow cytometer) can be correlated to an absolute number of detected photoelectrons (nPE).
#' Fluorescence signals acquired at different voltages and even by different flow cytometers become comparable on the nPE-scale.
#'
#' @param file_path character, path to a fcs file
#' @param compensate logical, should compensation be applied before PC calculation
#' @param compMat matrix, optional; a compensation matrix to use for compensation. If not provided the SPILL argument of the fcs file will be used. If you have generated a compensation matrix in FlowJo see ?fcexpr::wsx_compMats_to_fcs in order to have it copied to fcs files.
#' @param logicle_trans logical, should the logical transformation (Parks, 2006, https://pubmed.ncbi.nlm.nih.gov/16604519/) be applied before nPE calculation.
#' @param kfactor_df data.frame, table with k-factors for every channel and voltage per machine, defaults to system.file("extdata", "k_factors.rds", package = "fcexpr")
#' @param output_folder character, optional, path to a folder where to save the newly generated fcs file. Default is dirname(file_path).
#' @param new_file_suffix character, the suffix to add to the the newly generated fcs file. Default is _nPE.
#'
#' @return appended flowFrame which is also saved as fcs file
#' @export
#'
#' @examples
nPE_to_fcs <- function(file_path,
                       compensate = F,
                       compMat,
                       logicle_trans = F,
                       kfactor_df,
                       output_folder = NULL,
                       new_file_suffix = "nPE") {

  if (!"BiocManager" %in% rownames(utils::installed.packages())) {utils::install.packages("BiocManager")}
  if (!"flowCore" %in% rownames(utils::installed.packages())) {BiocManager::install("flowCore")}
  if (!"BiocGenerics" %in% rownames(utils::installed.packages())) {BiocManager::install("BiocGenerics")}

  if (!file.exists(file_path)) {
    stop(paste0(file_path, " not found."))
  }

  if (missing(kfactor_df)) {
    kfactor_df <- readRDS(system.file("extdata", "k_factors.rds", package = "fcexpr"))
  }

  ff_orig <- flowCore::read.FCS(file_path, truncate_max_range = F, emptyValue = F)

  if (!flowCore::keyword(ff_orig)[["$CYT"]] %in% names(kfactor_df)) {
    stop(paste0(flowCore::keyword(ff_orig)[["$CYT"]], " not found in names of kfactor_df. The fcs file had to be recorded at one of these machines: ", paste(names(kfactor_df), collapse = ", "), "."))
  }
  kfactor_df <- kfactor_df[[flowCore::keyword(ff_orig)[["$CYT"]]]]

  volts <- utils::stack(flowCore::keyword(ff_orig)[which(grepl("P[[:digit:]]{1,}V", names(flowCore::keyword(ff_orig))))])
  volts$ind <- gsub("V$", "", volts$ind)
  names(volts) <- c("volt", "ind")
  volts$volt <- as.numeric(volts$volt)
  pdata <- flowCore::pData(flowCore::parameters(ff_orig))
  volts <- merge(pdata, volts, by.x = "row.names", by.y = "ind")
  volts <- volts[which(!grepl("FSC|SSC", volts$name)), which(names(volts) %in% c("name", "volt"))]
  volts$channel <- as.character(gsub("-[[:alpha:]]{1}$", "", volts$name))
  volt_k_df <- merge(volts, kfactor_df, by = c("volt", "channel"))

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

  if (logicle_trans) {
    ff <- lgcl_trsfrm_ff(ff, channels = volt_k_df$name)
  }
  exprs <- flowCore::exprs(ff)
  exprs <- exprs[, which(colnames(exprs) %in% volt_k_df$name)]

  for (i in colnames(exprs)) {
    exprs[,i] <- exprs[,i]*volt_k_df[which(volt_k_df$name == i),"k"]
  }
  colnames(exprs) <- paste0(colnames(exprs), "_nPE")

  if (compensate && logicle_trans) {
    ext <- "_comp_lgcl"
  } else if (logicle_trans) {
    ext <- "_lgcl"
  } else if (compensate) {
    ext <- "_comp"
  }
  colnames(exprs) <- paste0(colnames(exprs), ext)

  new_p <- flowCore::parameters(ff_orig)[1,]
  new_kw <- flowCore::keyword(ff_orig)
  new_pars <- flowCore::parameters(ff_orig)
  for (z in 1:ncol(exprs)) {
    new_p_number <- as.integer(dim(ff_orig)[2]+z)
    rownames(new_p) <- c(paste0("$P", new_p_number))
    new_pars <- BiocGenerics::combine(new_pars, new_p)
    new_p_name <- colnames(exprs)[z]
    new_p_desc <- pdata[which(pdata$name == strsplit(colnames(exprs)[z], "_")[[1]][1]), "desc"]
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
    new_kw[paste0("$P",as.character(new_p_number),"DISPLAY")] <- "LIN"
    new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- min(exprs[,z])
    new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- max(exprs[,z])
  }
  ff_new <- methods::new("flowFrame", exprs = cbind(flowCore::exprs(ff_orig), exprs), parameters = new_pars, description = new_kw)

  if (!is.null(output_folder)) {
    dir.create(output_folder, recursive = T, showWarnings = F)
  } else {
    output_folder <- dirname(file_path)
  }

  flowCore::write.FCS(ff_new, file.path(output_folder, paste0(gsub("\\.fcs$", "", basename(file_path)), "_", new_file_suffix, ".fcs")))
  return(ff_new)
}

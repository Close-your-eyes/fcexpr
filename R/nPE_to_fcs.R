#' Convert fluorescence intensity (FI) to number of detected photo electrons (nPE)
#'
#' CAUTION: Very experimental function. Channels of selected flow cytometers have been analyzed with a QuantiFlash. By generation of very defined light pulses the
#' fluorescence intensity (FI) measured by a PMT can be correlated to an absolute number of detected photo electrons (nPE).
#' Currently this only works for fluorescence channels, not scatter channels. Also, strictly speaking, it is only valid for Height-channels (-H at the end)- For Area (-A) and Width (-W) it may not be correct.
#' Fluorescence signals acquired at different voltages and even by different flow cytometers become comparable on the nPE-scale. A good idea would be to test the conversion by analyzing the same sample with different
#' settings and/or at different machines.
#'
#' @param file_path character, path to a fcs file
#' @param kfactor_df data.frame, table with k-factors for every channel and voltage per machine, defaults to system.file("extdata", "k_factors.rds", package = "fcexpr")
#' @param output_folder character, optional, path to a folder where to save the newly generated fcs file. Default is dirname(file_path).
#' @param new_file_suffix character, the suffix to add to the the newly generated fcs file. Default is _nPE.
#' @param h_channels_only logical, have only the height channels converted to nPE; height channels are actually the only ones for which the conversion is correct.
#'
#' @return appended flowFrame which is also saved as fcs file to disk
#' @export
#'
#' @examples
#' \dontrun{
#' ff <- nPE_to_fcs(file_path = "myfolder/my_file.fcs")
#' }
nPE_to_fcs <- function(file_path,
                       h_channels_only = T,
                       kfactor_df = NULL,
                       output_folder = NULL,
                       new_file_suffix = "nPE") {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }
  if (!requireNamespace("BiocGenerics", quietly = T)){
    BiocManager::install("BiocGenerics")
  }
  if (!file.exists(file_path)) {
    stop(file_path, " not found.")
  }

  if (is.null(kfactor_df)) {
    kfactor_df <- readRDS(system.file("extdata", "k_factors.rds", package = "fcexpr"))
  } else {
    ## check the data frame
  }

  ff_orig <- flowCore::read.FCS(file_path, truncate_max_range = F, emptyValue = F)

  if (!flowCore::keyword(ff_orig)[["$CYT"]] %in% names(kfactor_df)) {
    # wayne and symphony are both called "Fortessa" - check for ambiguity!
    stop(flowCore::keyword(ff_orig)[["$CYT"]], " not found in names of kfactor_df. The fcs file had to be recorded at one of these machines: ", paste(names(kfactor_df), collapse = ", "), ".")
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
  volts$channel_type <- sapply(seq_along(volts$name), function (x) gsub(paste0(volts$channel[x], "-"), "", volts$name[x]))
  volt_k_df <- merge(volts, kfactor_df, by = c("volt", "channel"))
  if (all(volt_k_df$channel_type == "A")) {
    warning("Only Area-channels were detected in the fcs file. Conversion from FI to nPE is strictly valid only for Height-channels. Nevertheless, conversion will be done for Area channels as well. The more the width of events deviates the wronger the conversion becomes.")
    if (h_channels_only) {
      message("h_channels_only set to TRUE, so no channels will be converted.")
    }
  } else if (any(volt_k_df$channel_type != "H")) {
    warning("Conversion from FI to nPE is strictly valid only for Height-channels.")
  }
  if (h_channels_only) {
    volt_k_df <- volt_k_df[which(volt_k_df$channel_type == "H"),]
  }
  if (nrow(volt_k_df) == 0) {
    stop("No channels left to convert, consider setting h_channels_only to FALSE in order to have area and width channels converted. Converting these channels may be incorrect though.")
  }

'  if (compensate) {
    if (missing(compMat)) {
      compMat <- flowCore::keyword(ff_orig)[["SPILL"]]
      message("SPILL keyword used for compensation.")
    }
    ff <- flowCore::compensate(ff_orig, compMat)
  } else {
    ff <- ff_orig
    message("No compensation applied.")
  }
'

  exprs <- flowCore::exprs(ff_orig)
  exprs <- exprs[, which(colnames(exprs) %in% volt_k_df$name)]
  for (i in colnames(exprs)) {
    exprs[,i] <- exprs[,i]*volt_k_df[which(volt_k_df$name == i),"k"]
  }
  colnames(exprs) <- paste0(colnames(exprs), "_nPE")


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
    new_kw[paste0("$P",as.character(new_p_number),"DISPLAY")] <- "LOG"
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

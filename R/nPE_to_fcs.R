nPE_to_fcs <- function(file_path,
                       compensate = F,
                       compMat,
                       timeChannel,
                       logicle_trans = F,
                       kfactor_df) {

  if (!"BiocManager" %in% rownames(utils::installed.packages())) {utils::install.packages("BiocManager")}
  if (!"flowCore" %in% rownames(utils::installed.packages())) {BiocManager::install("flowCore")}
  if (!"BiocGenerics" %in% rownames(utils::installed.packages())) {BiocManager::install("BiocGenerics")}

  if (!file.exists(file_path)) {
    stop(paste0(file_path, " not found."))
  }

  if (missing(kfactor_df)) {
    kfactor_df <- readRDS(system.file("extdata", "k_factors.rds", package = "fcexpr"))
  }

  ff_orig <- flowCore::read.FCS(file_path, which.lines = which_lines, truncate_max_range = F, emptyValue = F)

  if (!flowCore::keyword(ff_orig)[["$CYT"]] %in% names(kfactor_df)) {
    stop(paste0(flowCore::keyword(ff_orig)[["$CYT"]], " not found in names of kfactor_df. The fcs file had to be recorded at one of these machines: ", paste(names(kfactor_df), collapse = ", "), "."))
  }

  volts <- stack(flowCore::keyword(ff_orig)[which(grepl("P[[:digit:]]{1,}V", names(flowCore::keyword(ff_orig))))])
  volts$ind <- gsub("V$", "", volts$ind)
  names(volts) <- c("volt", "ind")
  pdata <- flowCore::pData(flowCore::parameters(ff_orig))
  volts <- merge(pdata, volts, by.x = "row.names", by.y = "ind")
  volts <- volts[which(!grepl("FSC|SSC", volts$name)), which(names(volts) %in% c("name", "volt"))]

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


  if (logicle_trans) {
    ff <- lgcl_trsfrm_ff(ff, channels = channels)
  }
  exprs <- flowCore::exprs(ff)

  # filter scatter channels
  # get voltages



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



}

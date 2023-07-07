#' Write the compensation matrices generated in a flowjo workspace to the SPILL keywords of respective FCS files
#'
#' The compensation matrices generated or optimized in FlowJo only exist in the workspace. The underlying FCS files remain unchanged.
#' This is part of FlowJos philosophy, FCS remain original all the time. This is actually not a bad procedure but reduces flexibility in some cases.
#' This function saves the associated compensation matrices to the FCS files - the matrices are written into the SPILL keyword. That makes
#' them easily and unambiguously available for other calculations outside of FlowJo (e.g. in R) when compensated channels are required.
#'
#' @param ws path to flowjo workspace
#' @param ... additional arguments to fcexpr:::prep_spill()
#' @param groups which flowjo groups to get fcs files from
#' @param alt_FCS_file_folder if FCS files paths on disk are not the same as
#' in flowjos wsp file, then provide the correct path here
#'
#' @return no return but SPILL keyword updated in FCS files
#' @export
#'
#' @examples
#'\dontrun{
#' wsx_compMats_to_fcs(ws = "mypath/my.wsp")
#'}
wsx_compMats_to_fcs <- function(ws,
                                groups = NULL,
                                alt_FCS_file_folder = NULL,
                                ...) {

  if (!requireNamespace("CytoML", quietly = T)){
    BiocManager::install("CytoML")
  }
  if (!requireNamespace("flowWorkspace", quietly = T)){
    BiocManager::install("flowWorkspace")
  }
  ws <- check_ws(ws)

  ids <- wsx_get_groups(ws)
  if (is.null(groups)) {
    groups <- unique(ids[,"FlowJoGroup", drop=T])
  }
  ids <- ids[which(ids$FlowJoGroup %in% groups),"sampleID"]
  ss <- xml2::xml_find_all(xml2::xml_child(ws, "SampleList"), "Sample")
  ss <- ss[which(sapply(seq_along(ss), function(x) xml2::xml_attrs(xml2::xml_child(ss[[x]], "DataSet"))[["sampleID"]]) %in% ids)]

  compMats <- lapply(seq_along(ss), function(n) {
    sp <- xml2::xml_child(ss[[n]], "transforms:spilloverMatrix")
    # xml2::xml_children(sp)[2:length(xml2::xml_children(sp))]
    if (!is.na(sp)) {
      compMat <- t(do.call(cbind, lapply(xml2::xml_children(sp)[which(xml2::xml_name(xml2::xml_children(sp)) == "spillover")], function(x) {
        mat <- as.matrix(stats::setNames(as.numeric(xml2::xml_attr(xml2::xml_children(x), "value")), xml2::xml_attr(xml2::xml_children(x), "parameter")), nrow = 1, byrow = T)
        colnames(mat) <- xml2::xml_attr(x, "parameter")
        return(mat)
      })))
    } else {
      compMat <- NULL
    }
    return(compMat)
  })
  names(compMats) <- sapply(seq_along(ss), function(n) {
    gsub("^file:", "", xml2::xml_attr(xml2::xml_child(ss[[n]], "DataSet"), "uri"))
  })
  compMats <- compMats[which(!sapply(compMats, is.null))]


  if (!is.null(alt_FCS_file_folder)) {
    # replace path in compMats list with path to alt_FCS_file
    # match files via identities (the custom, own one)
    # get identity from keyword entries in wsp file below
    # make this a seperate functions somewhen

    fcs_identities <- .get_fcs_identities(kwl = wsx_get_keywords(ws = ws, return_type = "vector"))

    'k <- wsx_get_keywords(ws = ws, return_type = "vector")
    kk <- do.call(rbind, k)
    kk$FileName <- rep(names(k), sapply(k,nrow))
    kk <- dplyr::distinct(kk)
    ## taken from .check.FCS.files - made analogous
    filenames <- unique(kk$FileName)
    dd <- kk[which(kk$name == "$DATE"), "value"]
    tt <- kk[which(kk$name == "$BTIM"), "value"]
    et <- kk[which(kk$name == "$ETIM"), "value"]

    fil <- kk[which(kk$name == "$FIL"), "value"]
    tot <- kk[which(kk$name == "$TOT"), "value"]

    if (any(nchar(tt) - nchar(gsub(":", "", tt)) > 2)) {
      tt_fix_ind <- which(nchar(tt) - nchar(gsub(":", "", tt)) > 2)
      tt[tt_fix_ind] <- paste(rev(rev(strsplit(tt[tt_fix_ind], ":")[[1]])[-1]), collapse = ":")
    }
    datetime <- paste0(dd, "-", tt)
    sub <- ifelse(grepl("^2[[:digit:]]", tt) & grepl("^0[[:digit:]]", et), 86400, 0)
    datetime <- format(lubridate::parse_date_time(datetime, orders = c("%Y-%b-%d-%H:%M:%S", "%Y-%B-%d-%H:%M:%S", "%Y-%m-%d-%H:%M:%S", "%d-%b-%Y-%H:%M:%S",
                                                                       "%d-%m-%Y-%H:%M:%S", "%d-%B-%Y-%H:%M:%S", "%d-%b-%Y-%H:%M:%S")) - sub, "%Y.%m.%d-%H.%M.%S")
    if (any(is.na(datetime))) {
      warning("datetimes ", paste(paste0(dd, "-", tt)[which(is.na(datetime))], collapse = ", "), " could not be converted to a uniform format. Please, provide this to the package-maintainer.")
    }
    fcs.files <- stats::setNames(paste0(fil, "_-_", trimws(tot), "_-_", datetime), nm = filenames)
    if (length(unique(fcs.files)) != length(fcs.files)) {
      stop(paste0("Duplicate FCS files found. This is not allowed. Please, remove one of each duplicates. \n", paste(names(fcs.files[duplicated(fcs.files) |
                                                                                                                                       duplicated(fcs.files, fromLast = T)]), collapse = "\n")))
    }'

    alt_fcs_files <- .check.FCS.files(alt_FCS_file_folder)
    # replace paths with identity
    names(compMats) <- unname(fcs_identities[basename(names(compMats))])
    alt_fcs_files <- alt_fcs_files[which(alt_fcs_files %in% names(compMats))] # actually not needed
    compMats <- compMats[which(names(compMats) %in% alt_fcs_files)]
    alt_fcs_files_rev <- stats::setNames(names(alt_fcs_files), nm = alt_fcs_files)
    names(compMats) <- alt_fcs_files_rev[names(compMats)]
    message(length(names(compMats)), " FCS files in alt_FCS_file_folder matched with those from the workspace. CompMats will be written to those.")
  }

  for (i in seq_along(compMats)) {
    if (!file.exists(names(compMats)[i])) {
      warning("File ", names(compMats)[i], " not found. Did you change the names of FCS files and have not reconnected them to the workspace? Or did you copy the workspace from somewhere else? If so, open the workspace, reconnect the fcs files and save. Alternatively, provide alt_FCS_file_folder.")
    }
    ff <- flowCore::read.FCS(names(compMats)[i], truncate_max_range = F, emptyValue = F)
    sp <- flowCore::keyword(ff)[["SPILL"]]
    sp <- prep_spill(sp = sp, compMat = compMats[[i]], ...)
    flowCore::keyword(ff)[["SPILL"]] <- sp
    flowCore::write.FCS(ff, names(compMats)[i])
  }

}



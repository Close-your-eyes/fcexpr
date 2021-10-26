#' A modified version of flowCores read.FCS
#'
#' Please see ?flowCore::read.FCS for documentation. The additional feature of this modified version is
#' the attempt to read compensation-matrices (SPILL keyword) which have been written with a comma as decimal
#' separator instead of a dot. This error may occurs, when the language computer of BD Flow Cytometer is set
#' German. In this case the DIVA software writes an FCS file which does not obey the international standard
#' of using a dot as decimal separator. FlowJo by the way ignores this problem and just read the FCS file ignoring
#' the SPILL keyword. In this case information is lost but the file can be used for further analyses.
#' FlowCores read.FCS function though does not read these files.
#'
#' In essence the function txt2spillmatrix tries to to guess which comma is a decimal separator and which one is a
#' separator for different entries in the compensation matrix. One assumption is that the diagonal of the matrix has 1s only.
#'
#' @param filename see ?read.FCS
#' @param transformation see ?read.FCS
#' @param which.lines see ?read.FCS
#' @param alter.names see ?read.FCS
#' @param column.pattern see ?read.FCS
#' @param invert.pattern see ?read.FCS
#' @param decades see ?read.FCS
#' @param ncdf see ?read.FCS
#' @param min.limit see ?read.FCS
#' @param truncate_max_range see ?read.FCS
#' @param dataset see ?read.FCS
#' @param emptyValue see ?read.FCS
#' @param channel_alias see ?read.FCS
#' @param ... see ?read.FCS
#'
#' @return a flowFrame with corrected SPILL keyword
#' @export
#'
#' @examples
#' \dontrun{
#' # read the FCS file while the spill keyword is fixed (wrong comma replaced by dots)
#' ff <- read.FCS_spill.comma.sensitive(filename = "123.fcs")
#' # save the correct version as new file or overwrite the old one
#' # NOTE: the meta data reveal that the file was read and saved with flowCore
#' flowCore::write.FCS(x = ff, filename = "123_spill_corr.fcs")
#' }
read.FCS_comma <- function (filename, transformation = "linearize", which.lines = NULL,
                            alter.names = FALSE, column.pattern = NULL, invert.pattern = FALSE,
                            decades = 0, ncdf = FALSE, min.limit = NULL, truncate_max_range = TRUE,
                            dataset = NULL, emptyValue = TRUE, channel_alias = NULL,
                            ...) {

  if (!"BiocManager" %in% rownames(utils::installed.packages())) {utils::install.packages("BiocManager")}
  if (!"flowCore" %in% rownames(utils::installed.packages())) {BiocManager::install("flowCore")}
  if (!"flowWorkspace" %in% rownames(utils::installed.packages())) {BiocManager::install("flowWorkspace")}

  channel_alias <- flowCore:::check_channel_alias(channel_alias)
  if (ncdf)
    .Deprecated("'ncdf' argument is deprecated!Please use 'ncdfFlow' package for disk-based data structure.")
  if (!is.character(filename) || length(filename) != 1)
    stop("'filename' must be character scalar")
  if (!file.exists(filename))
    stop(paste("'", filename, "' is not a valid file", sep = ""))
  con <- file(filename, open = "rb")
  on.exit(close(con))
  fcsPnGtransform <- FALSE
  if (is.logical(transformation) && transformation || !is.null(transformation) && transformation == "linearize") {
    transformation <- TRUE
    scale <- FALSE
  }
  else if (!is.null(transformation) && transformation == "scale") {
    transformation <- TRUE
    scale <- TRUE
  }
  else if (!is.null(transformation) && transformation == "linearize-with-PnG-scaling") {
    transformation <- TRUE
    scale <- FALSE
    fcsPnGtransform <- TRUE
  }
  else if (is.null(transformation) || is.logical(transformation) && !transformation) {
    transformation <- FALSE
    scale <- FALSE
  }
  offsets <- flowCore:::findOffsets(con, emptyValue = emptyValue, dataset = dataset, ...)
  txt <- flowCore:::readFCStext(con, offsets, emptyValue = emptyValue, ...)
  if (fcsPnGtransform)
    txt[["flowCore_fcsPnGtransform"]] <- "linearize-with-PnG-scaling"
  if ("transformation" %in% names(txt) && txt[["transformation"]] %in% c("applied", "custom"))
    transformation <- FALSE
  mat <- flowCore:::readFCSdata(con, offsets, txt, transformation, which.lines,
                                scale, alter.names, decades, min.limit, truncate_max_range,
                                channel_alias)

  matRanges <- attr(mat, "ranges")
  id <- paste("$P", 1:ncol(mat), sep = "")
  zeroVals <- as.numeric(sapply(strsplit(txt[paste(id, "E", sep = "")], ","), function(x) x[2]))

  #absMin <- colMins(mat, , na.rm = TRUE)
  absMin <- apply(mat, 2, min, na.rm = T)

  realMin <- pmin(zeroVals, pmax(-111, absMin, na.rm = TRUE), na.rm = TRUE)
  keep_idx <- seq_along(colnames(mat))
  remove_idx <- NULL
  if (!is.null(column.pattern)) {
    n <- colnames(mat)
    keep_idx <- grep(column.pattern, n, invert = invert.pattern)
    remove_idx <- setdiff(seq_along(colnames(mat)), keep_idx)
    cols <- names(attr(mat, "dimnames")[[2]])
    mat <- mat[, keep_idx, drop = FALSE]
    matRanges <- matRanges[keep_idx]
    names(attr(mat, "dimnames")[[2]]) <- cols[keep_idx]
    attr(mat, "ranges") <- matRanges
    absMin <- absMin[keep_idx]
    realMin <- realMin[keep_idx]
    zeroVals <- zeroVals[keep_idx]
    id <- id[keep_idx]
  }
  if ("transformation" %in% names(txt) && txt[["transformation"]] ==
      "custom") {
    for (i in seq_along(colnames(mat))) {
      realMin[i] <- as.numeric(txt[[sprintf("flowCore_$P%sRmin", keep_idx[i])]])
    }
  }
  params <- flowCore:::makeFCSparameters(colnames(mat), txt, transformation, scale, decades, realMin, id = keep_idx)
  fix_pnr_idx <- which(is.na(params@data[, "maxRange"]))
  if (length(fix_pnr_idx) > 0) {
    fix_pnr_vals <- matRanges[fix_pnr_idx]
    params@data[fix_pnr_idx, "maxRange"] <- fix_pnr_vals
    params@data[fix_pnr_idx, "range"] <- fix_pnr_vals + 1
  }
  if (is.null(which.lines)) {
    total_number_of_events <- as.integer(flowCore:::readFCSgetPar(txt, "$TOT"))
    if (total_number_of_events != nrow(mat))
      stop("file", filename, "seems to be corrupted. \n The actual number of cells in data section (", nrow(mat), ") is not consistent with keyword '$TOT' (", total_number_of_events, ")")
  }
  txt[["FILENAME"]] <- filename
  if (transformation == TRUE) {
    txt[["transformation"]] <- "applied"
    for (p in seq_along(flowWorkspace::pData(params)$name)) {
      txt[[sprintf("$P%sE", p)]] <- sprintf("0,%g", 0)
      txt[[sprintf("flowCore_$P%sRmax", keep_idx[p])]] <- matRanges[p] + 1
      txt[[sprintf("flowCore_$P%sRmin", keep_idx[p])]] <- realMin[p]
    }
    txt[["$DATATYPE"]] <- "F"
  }
  if (offsets["FCSversion"] <= 2) {
    description <- strsplit(txt, split = "\n")
    names(description) <- names(txt)
  }
  else {
    description <- strsplit(txt, split = NA)
  }
  if (length(fix_pnr_idx) > 0)
    description[paste0("$P", fix_pnr_idx, "R")] <- fix_pnr_vals + 1
  if (!is.null(remove_idx) && length(remove_idx) > 0) {
    remove_regex <- paste0("\\$P", remove_idx, "[A-Z]+")
    remove_keys <- lapply(remove_regex, function(rx) grep(rx, names(description), value = TRUE))
    description <- description[!names(description) %in% do.call(c, remove_keys)]
  }

  for (sn in flowCore:::.spillover_pattern) {
    sp <- description[[sn]]
    if (!is.null(sp)) {
      #sp <- flowCore:::txt2spillmatrix(sp, cpp = F)
      sp <- txt2spillmatrix_comma(sp, cpp = F)
      if (is.matrix(sp)) {
        cnames <- colnames(sp)
        cnames <- flowCore:::update_channel_by_alias(cnames, channel_alias, silent = TRUE)
        if (alter.names)
          cnames <- make.names(cnames)
        colnames(sp) <- cnames
        description[[sn]] <- sp
      }
    }
  }
  tmp <- methods::new("flowFrame", exprs = mat, description = description,
                      parameters = params)
  flowCore::identifier(tmp) <- basename(flowCore::identifier(tmp))
  return(tmp)
}


txt2spillmatrix_comma <- function (txt, cpp = TRUE) {
  if (cpp) {
    flowCore:::string_to_spill(txt)
  }
  else {
    splt <- strsplit(txt, ",")[[1]]
    nrCols <- as.numeric(splt[1])

    if (!is.na(nrCols) && nrCols > 0) {
      if ((length(splt)-1) %% nrCols != 0) {
        ## assume there are ones the diagonal
        rows <- strsplit(paste(strsplit(sp, ",")[[1]][-1], collapse = ","), ",1,|,1$")[[1]]
        new.rows <- lapply(seq_along(rows)[-1], function(i) {
          split.row <- strsplit(rows[i], ",")[[1]]
          new.row <- c()
          j<-1
          while(j <= length(split.row)) {
            if (j < length(split.row) && split.row[j] == 0 && split.row[j+1] != 0) {
              new.row <- c(new.row, paste0(split.row[j], ".", split.row[j+1]))
              j<-j+2
            } else {
              new.row <- c(new.row, split.row[j])
              j<-j+1
            }
          }
          return(c("1", new.row))
        })
        txt <- paste(nrCols, rows[1], paste(do.call(c, new.rows), collapse = ","), "1", sep = ",")
        splt <- strsplit(txt, ",")[[1]]
        nrCols <- as.numeric(splt[1])
      }
      cnames <- splt[2:(nrCols + 1)]
      vals <- as.numeric(splt[(nrCols + 2):length(splt)])
      matrix(vals, ncol = nrCols, byrow = TRUE, dimnames = list(NULL, cnames))
    }
    else txt
  }
}

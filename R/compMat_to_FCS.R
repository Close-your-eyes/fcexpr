#' Title
#'
#' @param fcs_file_path
#' @param compMat_file_path
#' @param max_match_dist
#' @param skip_check
#'
#' @return
#' @export
#'
#' @examples
compMat_to_fcs <- function(fcs_file_path, compMat_file_path, max_match_dist = 1, skip_check = F) {

  if (!file.exists(compMat_file_path)) {
    stop("compMat not found.")
  }
  if (!file.exists(fcs_file_path)) {
    stop("fcs_file not found.")
  }
  if (class(compMat_file_path) != "character" || rev(strsplit(basename(compMat_file_path), "\\.")[[1]])[1] != "csv") {
    stop("compMat has to be character (path to a csv file).")
  }

  compMat <- utils::read.csv(compMat_file_path, header = T, row.names = 1, check.names = F)
  if (!identical(rownames(compMat), colnames(compMat))) {
    stop("colnames and rownames of compMat have to be equal.")
  }
  ff <- flowCore::read.FCS(fcs_file_path, truncate_max_range = F, emptyValue = F)
  sp <- flowCore::keyword(ff)[["SPILL"]]
  rownames(sp) <- colnames(sp)

  if (!all(colnames(compMat) %in% colnames(sp))) {
    print("Not all colnames of compMat found in those of the SPILLOVER keyword matrix from the FCS file.")
    # match channel names
    if (any(apply(adist(colnames(compMat), colnames(sp)), 1, min) > max_match_dist)) {
      stop("Too big string distances between channel names of compMat and FCS file. Please, check the column names or make sure you provide the correct compensation matrix.")
    }
    match_ind <- apply(adist(colnames(compMat), colnames(sp)), 1, which.min)
    if (any(duplicated(match_ind))) {
      stop("Channel names from compMat not uniquely matched to channel names from FCS file.")
    }
    print("Matched channel names:")
    print(data.frame(compMat = colnames(compMat), FCS = colnames(sp)[match_ind]))
    colnames(compMat) <- colnames(sp)[match_ind]
    rownames(compMat) <- colnames(sp)[match_ind]
    if (!skip_check) {
      if (interactive()) {
        choice <- utils::menu(c("Yes", "No"), title = "Channel names matched correctly - Continue?")
        if (choice == 2) {
          return(NULL)
        }
      }
    }
  }


  for (rr in rownames(compMat)) {
    for (cc in colnames(compMat)) {
      sp[which(rownames(sp) == rr), which(colnames(sp) == cc)] <- compMat[which(rownames(compMat) == rr), which(colnames(compMat) == cc)]
    }
  }

  rownames(sp) <- NULL
  flowCore::keyword(ff)[["SPILL"]] <- sp
  flowCore::write.FCS(ff, fcs_file_path)
  print(fcs_file_path)
}


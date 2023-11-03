#' Make an existing compensation matrix in a FCS neutral
#'
#' This will overwrite an existing compensation matrix in a FCS file so that
#' there are ones on the diagonal and else zeros. This is handy in order to
#' not having to type this manually in flow jo or so.
#'
#' @param fcs_file_path character, file path to the fcs file
#' @param compMat_keyword character, name of the keyword which holds the matrix
#' @param channel_names enter channels for comp matrix; only applied when matrix under compMat_keyword is missing of size 0
#'
#' @return no return, but FCS changed on disk
#' @export
#'
#' @examples
#'\dontrun{
#'compMap_neutral_to_fcs(fcs_file_path = "my/path/to/FCSfile.fcs)
#' }
compMat_neutral_to_fcs <- function(fcs_file_path,
                                   compMat_keyword = "SPILL",
                                   channel_names = NULL) {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }
  if (!file.exists(fcs_file_path)) {
    stop("fcs_file not found.")
  }


  ff <- flowCore::read.FCS(fcs_file_path, truncate_max_range = F, emptyValue = F)

  new_key <- F
  if (!compMat_keyword %in% names(flowCore::keyword(ff))) {
    candidates <- grep("Spill", names(flowCore::keyword(ff)), value = T, ignore.case = T)
    if (length(candidates > 0)) {
      message("compMat_keyword not found in keywords of FCS file. You may try one of these: ", paste(candidates, collapse = ","))
    } else {
      message("compMat_keyword not found in keywords of FCS file and also no other Keyword with 'Spill' in it. Will assign a new keyword for comp matrix: 'SPILL' and will create matrix from channel names in FCS parameters.")
      new_key <- T
      compMat_keyword <- "SPILL"
      flowCore::keyword(ff)[[compMat_keyword]] <- matrix(nrow = 0, ncol = 0)
    }
  }

  sp <- flowCore::keyword(ff)[[compMat_keyword]]
  if (!methods::is(sp, "matrix") || nrow(sp) == 0) {
    if (!methods::is(sp, "matrix")) {
      message(compMat_keyword, " in FCS file is not a matrix. Will create matrix from channel names in FCS parameters.")
    } else if (nrow(sp) == 0 && !new_key) {
      message("CompMat in FCS file has zero rows. Will create matrix from channel names in FCS parameters.")
    }
    if (is.null(channel_names)) {
      channel_names <- flowCore::parameters(ff)@data$name
      channel_names <- channel_names[which(!grepl("HDR|Time|FSC|SSC|BSC", channel_names, ignore.case = T))]
      message("Channels included in matrix are: ", paste(channel_names, collapse = ", "), ". If wrong or some are missing, please provide manually via the channel_names argument.")
    } else {
      temp <- channel_names[!which(channel_names %in% flowCore::parameters(ff)@data$name)]
      if (length(temp) > 0) {
        message("These channel_names were not found: ", paste(temp, collapse = ","), ". They are not included in the newly generated comp matrix.")
      }
      channel_names <- channel_names[which(channel_names %in% flowCore::parameters(ff)@data$name)]
    }
    sp_neutral <- diag(1, nrow = length(channel_names), ncol = length(channel_names))
    colnames(sp_neutral) <- unname(channel_names)
  } else {
    sp_neutral <- diag(1, nrow = nrow(sp), ncol = ncol(sp))
    colnames(sp_neutral) <- colnames(sp)
  }

  flowCore::keyword(ff)[[compMat_keyword]] <- sp_neutral
  flowCore::write.FCS(ff, fcs_file_path)
  message(fcs_file_path)
}

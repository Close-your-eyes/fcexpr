#' Make an existing compensation matrix in a FCS neutral
#'
#' This will overwrite an existing compensation matrix in a FCS file so that
#' there are ones on the diagonal and else zeros. This is handy in order to
#' not having to type this manually in flow jo or so.
#'
#' @param fcs_file_path character, file path to the fcs file
#' @param compMat_keyword character, name of the keyword which holds the matrix
#'
#' @return no return, but FCS changed on disk
#' @export
#'
#' @examples
#'\dontrun{
#'compMap_neutral_to_fcs(fcs_file_path = "my/path/to/FCSfile.fcs)
#' }
compMat_neutral_to_fcs <- function(fcs_file_path, compMat_keyword = "SPILL") {

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

  if (!compMat_keyword %in% names(flowCore::keyword(ff))) {
    stop("compMat_keyword not found in keywords of FCS file.")
  }

  sp <- flowCore::keyword(ff)[[compMat_keyword]]
  if (!methods::is(sp, "matrix")) {
    stop(compMat_keyword, " in FCS file is not a matrix.")
  }

  sp_neutral <- diag(1 ,nrow = nrow(sp), ncol = ncol(sp))
  colnames(sp_neutral) <- colnames(sp)

  flowCore::keyword(ff)[[compMat_keyword]] <- sp_neutral
  flowCore::write.FCS(ff, fcs_file_path)
  message(fcs_file_path)
}

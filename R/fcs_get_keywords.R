#' Get keywords from all FCS files in a folder or subfolder
#'
#' The returned data frame can be used to filter for relevant keywords which may then be joined to the sampledescription.
#' Note that if you have defined your keywords in flowjo you have to export the FCS files in order to hard-code the keyword into the FCS file. Otherwise it only exists in the flowjo workspace. In order to pull out keywords from there another function exists.
#'
#' @param FCS.file.folder path to root folder which contains FCS files; if missing file.path(getwd(), 'FCS_files') is assumed
#'
#' @return a data frame of keywords as columns and FCS files as rows
#' @export
#'
#' @examples
#' \dontrun{
#' key_df <- fcs_get_keywords()
#' key_df <- key_df[,c("FileName", "$OP")]
#' # sd is data frame of sampledescription
#' sd <- dplyr::left_join(sd, key_df, by = "FileName")
#' # save sd
#' }
fcs_get_keywords <- function(FCS.file.folder = file.path(getwd(), "FCS_files")) {

  if (!file.exists(FCS.file.folder)) {
    stop(paste0(FCS.file.folder, " not found."))
  }

  keys <- lapply(list.files(FCS.file.folder, "\\.fcs$", ignore.case = T, recursive = T, full.names = T), function(x) {
    k <- flowCore::keyword(flowCore::read.FCS(x, which.lines = 1, emptyValue = T, truncate_max_range = F))
    kdf <- suppressWarnings(utils::stack(k[which(!names(k) == "SPILL")]))
    names(kdf) <- c("value", "keyword")
    kdf[,"FileName"] <- basename(x)
    kdf[,"FilePath"] <- dirname(x)
    return(kdf)
  })
  keys <- tidyr::pivot_wider(do.call(rbind, keys), names_from = keyword, values_from = value)

  return(as.data.frame(keys))
}







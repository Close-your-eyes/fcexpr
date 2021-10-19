#' Get keywords from all FCS files in a folder or subfolder
#'
#' The returned data frame can be used to filter for relevant keywords which may then be joined to the sampledescription.
#' Note that if you have defined your keywords in flowjo you have to export the FCS files in order to hard-code the keyword into the FCS file. Otherwise it only exists in the flowjo workspace. In order to pull out keywords from there another function exists.
#'
#' @param file_path path to the fcs file
#'
#' @return a data frame of keywords as columns and FCS files as rows
#' @export
#'
#' @examples
#' \dontrun{
#' # get FCS files in a folder
#' files <- list.files(FCS.files.folder, recursive = T, full.names = T)
#' # loop through files and extract keywords; only SPILL is
#' keys <- lapply(files, function(x) {
#' fcs_get_keywords(x)
#' })
#' # convert the lists of data frames to one wide data frame
#' keys <- tidyr::pivot_wider(do.call(rbind, keys), names_from = keyword, values_from = value)
#' }
fcs_get_keywords <- function(file_path) {

  if (!file.exists(file_path)) {
    stop(paste0(file_path, " not found."))
  }

  k <- flowCore::keyword(flowCore::read.FCS(file_path, which.lines = 1, emptyValue = T, truncate_max_range = F))
  if ("SPILL" %in% names(k)) {
    k[["SPILL"]] <- paste(as.character(k[["SPILL"]]), collapse = ",")
  }
  kdf <- suppressWarnings(utils::stack(k))
  names(kdf) <- c("value", "keyword")
  kdf[,"FileName"] <- basename(file_path)
  kdf[,"FilePath"] <- dirname(file_path)

  return(kdf)
}







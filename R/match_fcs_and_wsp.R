#' Find out in which FlowJo workspaces FCS files have been analyzed.
#'
#' @param FCS.file.folder one or more absolute paths to folder containing FCS files;
#' will be scanned for .fcs files recursively (subfolder as well)
#' @param FCS.file.paths absolute paths to to FCS files
#' @param exclude.folders folders to exclude; if any of these strings appears
#' in the absolute path to a fcs file this file will be excluded
#' @param wsp one or more absolute paths to flowjo workspaces on disk
#'
#' @return a data.frame of matched (joined) information about fcs files on disk
#' and fcs files referred to in flowjo workspaced
#' @export
#'
#' @examples
#'\dontrun{
#' match_fcs_and_wsp(FCS.file.folder = file.path(wd, "FCS_files"), ws = "mypath/my.wsp")
#'}
match_fcs_and_wsp <- function(FCS.file.folder = NULL,
                              FCS.file.paths = NULL,
                              exclude.folders = NULL,
                              wsp) {

  if (missing(wsp)) {
    stop("wsp missing. Please provide a vector of paths to wsp files.")
  }

  if (!is.null(FCS.file.folder)) {
    fcs.file.paths1 <- list.files(path = FCS.file.folder, pattern = "\\.fcs", full.names = T, recursive = T, ignore.case = T)
    if (length(fcs.file.paths1) == 0) {
      message("No FCS files found in FCS.file.folder or its subfolders.")
    }
  }

  if (!is.null(FCS.file.paths)) {
    FCS.file.paths <- unique(FCS.file.paths)
    FCS.file.paths <- FCS.file.paths[which(file.exists(FCS.file.paths))]
    if (length(FCS.file.paths) == 0) {
      message("FCS files in FCS.file.paths not found. Did you provide full paths? If not, do so, please.")
    }
  }


  FCS.file.paths <- c(FCS.file.paths, fcs.file.paths1)
  if (!is.null(exclude.folders)) {
    FCS.file.paths <- FCS.file.paths[which(!grepl(paste0(tolower(exclude.folders), collapse = "|"), tolower(FCS.file.paths)))]
  }
  if (length(FCS.file.paths) == 0) {
    stop("No FCS files found or remaining after exclusion of 'exclude.folders' in FCS.file.folder and FCS.file.paths.")
  }

  wsp <- wsp[which(file.exists(wsp))]
  wsp <- wsp[which(grepl("\\.wsp$", x = wsp, ignore.case = T))]
  if (length(wsp) == 0) {
    stop("No wsps found in wsp. Did you provide full paths? If not, do so, please.")
  }


  fcs_file_idents <-  utils::stack(.get_fcs_identities(kwl = flowCore::read.FCSheader(FCS.file.paths, emptyValue = F)))
  fcs_file_idents <- dplyr::rename(fcs_file_idents, "identity" = values, "disk_file_path" = ind)
  fcs_file_idents$disk_file_path <- as.character(fcs_file_idents$disk_file_path)
  fcs_file_idents$disk_file_name <- basename(fcs_file_idents$disk_file_path)

  fcs_file_idents_wsp <- purrr::map_df(stats::setNames(wsp, wsp), function(x) utils::stack(.get_fcs_identities(kwl = wsx_get_keywords(ws = x, return_type = "vector"))), .id = "wsp_path")
  fcs_file_idents_wsp <- dplyr::rename(fcs_file_idents_wsp, "identity" = values, "flowjo_file_path" = ind)
  fcs_file_idents_wsp$flowjo_file_path <- as.character(fcs_file_idents_wsp$flowjo_file_path)
  fcs_file_idents_wsp$wsp_file <- basename(fcs_file_idents_wsp$wsp_path)
  fcs_file_idents_wsp$flowjo_file_name <- basename(fcs_file_idents_wsp$flowjo_file_path)

  matched_df <- dplyr::distinct(dplyr::full_join(fcs_file_idents_wsp, fcs_file_idents, by = "identity"))
  # if wsp_path is NA, then the fcs file does not appear in any wsp
  # if disk_file_path is NA, then there is a fcs file analyzed in a wsp which does not exist on disk (FCS.file.folder and FCS.file.paths)
  ## add groups??

  return(matched_df)
}

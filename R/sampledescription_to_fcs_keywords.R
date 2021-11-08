#' Write columns from a sampledescription as keyword into FCS files
#'
#' ATTENTION: This function will alter the meta data of FCS files
#' such that one can notice it was once opened and saved with flowCore.
#'
#' @param sampledescription data frame
#' @param columns vector of column names from sampledescription to write as keywords
#' @param FCS.file.folder folder of FCS files
#'
#' @return no return but keywords written to FCS files
#' @export
#'
#' @examples
#' \dontrun{
#' wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
#' setwd(wd)
#' sd <- openxlsx::read.xlsx(sampledescription.xlsx)
#' # if only a subset of files should be considered, select respective rows
#' sampledescription_to_fcs_keywords(sampledescription = sd, columns = c("Patient", "ExpPart"),FCS.file.folder = "FCS_files)
#' }
sampledescription_to_fcs_keywords <- function(sampledescription,
                                              columns,
                                              FCS.file.folder) {

  if (!"BiocManager" %in% rownames(utils::installed.packages())) {install.packages("BiocManager")}
  if (!"flowCore" %in% rownames(utils::installed.packages())) {BiocManager::install("flowCore")}

  fcs_files <- .check.FCS.files(FCS.file.folder)
  fcs_files <- stats::setNames(names(fcs_files), fcs_files)
  sd <- as.data.frame(sampledescription, stringsAsFactors = F)
  if (!"identity" %in% names(sd)) {
    stop("identity column not found in sd.")
  }
  fcs_files <- fcs_files[sd$identity]
  if (length(fcs_files) == 0) {
    stop("No matching FCS files found.")
  }

  if (length(columns[which(columns %in% names(sd))]) == 0) {
    stop("No matching column names found in sd.")
  }
  print(paste0(length(columns[which(columns %in% names(sd))]), " of ", length(columns), " columns found in sd: " ))
  print(columns[which(columns %in% names(sd))])
  columns <- columns[which(columns %in% names(sd))]

  for (x in seq_along(fcs_files)) {
    if (any(!is.na(sd[which(sd$identity == names(fcs_files[x])),columns]))) {
      f <- flowCore::read.FCS(fcs_files[x], truncate_max_range = F, emptyValue = F)
      for (k in columns) {
        if (!is.na(is.na(sd[which(sd$identity == names(fcs_files[x])),k]))) {
          flowCore::keyword(f)[[k]] <- sd[which(sd$identity == names(fcs_files[x])),k]
        }
      }
      flowCore::write.FCS(f, fcs_files[x])
    }
  }


}




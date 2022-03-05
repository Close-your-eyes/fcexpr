#' Write columns from sampledescription as keyword into FCS files
#'
#' Given a sampledescription with annotation of fcs files such meta data may be
#' written as keyword into respective fcs files.
#'
#' @param sampledescription data frame
#' @param columns vector of column names from sampledescription to write as keywords
#' @param FCS.file.folder path to folder of FCS files
#'
#' @return no return but keywords written to FCS files
#' @export
#'
#' @examples
#' \dontrun{
#' wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
#' sd <- openxlsx::read.xlsx(file.path(wd, sampledescription.xlsx))
#' # if only a subset of files should be considered, select respective rows
#' sampledescription_to_fcs_keywords(sampledescription = sd, columns = c("Patient", "ExpPart"), FCS.file.folder = file.path(wd, "FCS_files"))
#' }
sampledescription_to_fcs_keywords <- function(sampledescription,
                                              columns,
                                              FCS.file.folder) {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }

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
  message(length(columns[which(columns %in% names(sd))]), " of ", length(columns), " columns found in sd: " )
  message(columns[which(columns %in% names(sd))])
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



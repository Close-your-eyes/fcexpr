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
#' sampledescription_to_fcs_keywords(sampledescription = sd,
#' columns = c("Patient", "ExpPart"),
#' FCS.file.folder = "FCS_files)
#' }
sampledescription_to_fcs_keywords <- function(sampledescription,
                                              columns,
                                              FCS.file.folder) {

  if (!"flowCore" %in% rownames(installed.packages())) {BiocManager::install("flowCore")}

  fcs_files <- .check.FCS.files(FCS.file.folder)
  fcs_files <- stats::setNames(names(fcs_files), fcs_files)
  sampledescription <- as.data.frame(sampledescription)
  if (!"identity" %in% names(sampledescription)) {
    stop("identity column not found in sampledescription.")
  }
  fcs_files <- fcs_files[sampledescription$identity]
  if (length(fcs_files) == 0) {
    stop("No matching FCS files found.")
  }

  if (length(columns[which(columns %in% names(sampledescription))]) == 0) {
    stop("No matching column names found in sampledescription.")
  }
  print(paste0(length(columns[which(columns %in% names(sampledescription))]), " of ", length(columns), " columns found in sampledescription: " ))
  print(columns[which(columns %in% names(sampledescription))])
  columns <- columns[which(columns %in% names(sampledescription))]

  lapply(fcs_files, function(x) {
    if (any(!is.na(sampledescription[which(sampledescription$identity == names(x)),columns]))) {
      f <- flowCore::read.FCS(x, truncate_max_range = F, emptyValue = F)
      for (k in columns) {
        if (!is.na(is.na(sampledescription[which(sampledescription$identity == names(x)),k]))) {
          flowCore::keyword(f)[[k]] <- sampledescription[which(sampledescription$identity == names(x)),k]
        }
      }
      flowCore::write.FCS(f, x)
    }
  })

}






#' Get compensation matrices from flowjo workspace
#'
#' From xml keywords the compensation matrices are obtained. Keyword is SPILL, usually.
#' The matrix first comes as one string which is converted to a matrix with
#' fcexpr:::compmat_key_to_mat. Matrices created in flowjo only exist in the xml (wsp) file.
#' To write them to fcs files, use fcexpr::compMat_to_FCS.
#'
#' @param ws path to flowjo workspace
#' @param comp_keyword SPILL or SPILLOVER
#'
#' @return data frame with FileNames and prepared matrices
#' @export
#'
#' @examples
#' \dontrun{
#' compmat_df <- fcexpr::wsx_get_compmat(ws = "path/to/flowjo/workspace.wsp")
#' }
wsx_get_compmat <- function(ws, comp_keyword = "SPILL") {
  compmat_df <- tibble::as_tibble(dplyr::bind_rows(fcexpr::wsx_get_keywords(
    ws = ws,
    return = "data.frame",
    keywords = comp_keyword
  ), .id = "FileName"))
  if (nrow(compmat_df) == 0) {
    stop("keyword not found.")
  }
  compmat_df[["mat"]] <- lapply(compmat_df$value, compmat_key_to_mat)
  return(compmat_df)
}

compmat_key_to_mat <- function(character) {
  if (length(character) > 1) {
    stop("one character only.")
  }
  comp <- strsplit(character, ",")[[1]]
  nc <- as.numeric(comp[1])
  comp <- comp[-1]
  cnames <- comp[1:nc]
  cnum <- as.numeric(comp[(nc+1):length(comp)])
  nr <- length(cnum)/nc
  comp <- matrix(cnum, nrow = nr, ncol = nc, byrow = T, dimnames = list(NULL, cnames))
  return(comp)
}

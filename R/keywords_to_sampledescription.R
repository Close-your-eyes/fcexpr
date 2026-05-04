#' Add keywords from FCS files to sampledescription
#'
#' FCS files are identified by identity and searched for in
#' FCS.file.folder in dirname(sd_path).
#'
#' @param sd_path path to sampledescription.xlsx
#' @param keywords keywords to read from FCS files and write to sd
#' @param FCS.file.folder folder name with FCS files in dirname(sd_path)
#'
#' @return nothing, sampledescription file is updated on disk
#' @export
#'
#' @examples
#'\dontrun{
#' keywords_to_sampledescription(sd_path = "path/sampledescription.xlsx",
#'                               keywords = c("$TOT", "$OP"))
#'}
keywords_to_sampledescription <- function(sd_path,
                                          keywords,
                                          FCS.file.folder = "FCS_files") {

  if (!requireNamespace("brathering", quietly = T)) {
    pak::pak("Close-your-eyes/brathering")
  }

  xlsx <- openxlsx::read.xlsx(sd_path)
  fcs_files <- list.files(file.path(dirname(sd_path), FCS.file.folder), pattern = ".fcs", full.names = T, ignore.case = T)
  keys <- flowCore::read.FCSheader(fcs_files, keyword = keywords)
  keys <- purrr::map(keys, ~.x[which(!is.na(.x))])
  keys <- purrr::discard(keys, purrr::negate(length))
  if (!length(keys)) {
    message("no keyworrds found.")
    return(NULL)
  }
  names(keys) <- get_fcs_identities(flowCore::read.FCSheader(fcs_files))
  keys_df <-
    purrr::map_dfr(keys, stack, .id = "identity") |>
    tidyr::pivot_wider(names_from = ind, values_from = values)
  xlsx <- brathering::coalesce_join(xlsx, keys_df, by = "identity")
  openxlsx::write.xlsx(xlsx, file = sd_path)
  message(sd_path)
}

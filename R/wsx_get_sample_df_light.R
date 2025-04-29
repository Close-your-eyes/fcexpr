#' Get information on samples in a flowjo workspace
#'
#' Info on group belonging, some keywords which form the identity of FCS files
#' and location on disk as stored in the ws. Th summary is a collapsed version
#' which can be used to check whether there are sample in one FlowJoGroup that
#' have the same $FIL keyword(=name). Because CytoML::flowjo_to_gatingset may
#' complain that one of these files does not match, if either one is to be
#' imported. In comparison to wsx_get_sample_df, local FCS file existence (as on
#' disk on the current machine is not checked).
#'
#' @param ws path to a flowjo workspace
#'
#' @return list of data frames
#' @export
#'
#' @examples
wsx_get_sample_df_light <- function(ws) {

  paths <- wsx_get_fcs_paths(ws, split = F, filter_AllSamples = T)
  kwlist <- wsx_get_keywords(ws, return = c("data.frame", "vector"), keywords = c("$DATE", "$FIL", "$TOT", "$ETIM", "$BTIM"))
  fcs_ident <-
    stack(fcexpr:::.get_fcs_identities(kwlist[["vec"]], allow_duplicates = T)) |>
    dplyr::rename("identity" = values, "FlowJoFileName" = ind) |>
    dplyr::mutate(FlowJoFileName = as.character(FlowJoFileName))
  keys <-
    kwlist[["df2"]] |>
    dplyr::left_join(fcs_ident, by = "FlowJoFileName")

  sampledf <-
    paths |>
    dplyr::left_join(keys, by = c("FileName" = "FlowJoFileName")) |>
    dplyr::mutate(sampleID = as.numeric(sampleID)) |>
    dplyr::mutate(FileFolder = basename(dirname(FilePath)), .after = FilePath) |>
    dplyr::arrange(sampleID) |>
    dplyr::rename("FlowJoFilePath" = FilePath, "FlowJoFileName" = FileName, "FlowJoFileFolder" = "FileFolder")

  fil_summary <-
    sampledf |>
    dplyr::group_by(FlowJoGroup, `$FIL`) |>
    tidyr::nest(.key = "sample_info") |>
    dplyr::rowwise() |>
    dplyr::mutate(n_samples = nrow(sample_info)) |>
    dplyr::mutate(n_TOT_levels = length(unique(sample_info$`$TOT`))) |>
    dplyr::mutate(n_FileFolders = length(unique(sample_info$FlowJoFileFolder))) |>
    dplyr::ungroup()

  return(list(sampledf = sampledf, summary = fil_summary))
}

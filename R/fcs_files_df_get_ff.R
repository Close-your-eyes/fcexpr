#' Convenience function to read many FCS files as flowframes into memory
#'
#' Read FCS files as flowframes into memory and have them compensated with compensation
#' matrices from the respective flowjo workspace. This function combines flowCore::read.FCS,
#' fcexpr::wsx_get_compmat and fcexpr::ff_apply_compensation. The flowframes my then be
#' subsetted wit index matrices (ind_mat) as returned from fcexpr::wsp_get_ff or
#' fcexpr::wsp_get_indices.
#'
#' @param df dafa frame with information to FCS files; as returned from fcexpr::wsp_get_gs or fcexpr::wsp_get_ff
#' @param compensate apply compensation with matrices from flowjo workspaces, obtained with fcexpr::wsx_get_compmat
#' @param comp_keyword the keyword under which the compensation is saved, usually SPILL or SPILLOVER;
#' passed to fcexpr::wsx_get_compmat
#' @param match_channels if channel names in compmat do not exactly match those in flowframes
#' (due to conversion of special characters), try to find the best matching channel name
#' and change it in the compmat; passed to ff_apply_compensation
#'
#' @return a list of flowframes
#' @export
#'
#' @examples
#' \dontrun{
#' # read subsetted flowframes (by population) via flowjo workspace
#' ffdata <- fcexpr::wsp_get_ff(wsp = "mypath/to/workspace.wsp", population = c("CD4", "CD8"), groups = "flowjogroup", FCS.file.folder = "local/path/to/fcsfiles")
#' # do the equivalent separately obtained comp matrix, indices (which rows belong to which population) and fcs files on disk
#' ffs <- fcexpr::fcs_files_df_get_ff(df = ffdata[["FCS_files"]])
#' # subset with ind mat
#' ff1 <- flowCore::Subset(ffs[[1]], ffdata[["ind_mats"]][[1]][,"/Time/Lymphocytes/Single Cells/CD3+/CD8"])
#' }
fcs_files_df_get_ff <- function(df,
                                compensate = T,
                                comp_keyword = "SPILL",
                                match_channels = T) {
  message("reading flowframes into memory:")
  ffs <- purrr::map(stats::setNames(df$FilePathUse, df$FlowJoFileName), function(x) {
    message(basename(x))
    flowCore::read.FCS(x, truncate_max_range = F)
  })

  if (compensate) {
    compmats <- purrr::map_dfr(unique(df$wsp), wsx_get_compmat, comp_keyword = comp_keyword)
    compmats <- compmats |> dplyr::filter(FileName %in% names(ffs))
    compmats <- compmats[order(match(compmats$FileName, names(ffs))), "mat", drop = T]
    ffs <- purrr::map2(ffs, compmats, ~ff_apply_compensation(.x, .y), match_channels = match_channels)
  }
  return(ffs)
}

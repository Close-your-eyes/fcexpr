#' Read data from a flowjo wsp file
#'
#' @param ws path to flowjo workspace or a parsed xml-document (xml2::read_xml(ws))
#' @param return_stats logical whether to return statistics
#' @param groups vector of flowjo group names to consider
#' @param invert_groups logical whether to exclude the selected groups
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply
#' or parallel::mclapply are suggested; only for legacy method
#' @param strip_data remove extensive internal data of wsp from count data frame
#' @param skip_legacy do not run the legacy method?
#' @param ... additional argument to the lapply function;
#' mainly mc.cores when parallel::mclapply is chosen
#'
#' @return list of data frames of counts and stats
#' @export
#'
#' @examples
#' \dontrun{
#' # find workspaces
#' ws <- list.files(path = wd, pattern = '\\.wsp$', recursive = T, full.names = T)
#' # read counts
#' counts <- wsx_get_popstats(ws = ws[[1]])[["counts"]]
#' }
wsx_get_popstats <- function(ws,
                             groups = NULL,
                             invert_groups = F,
                             return_stats = F,
                             lapply_fun = lapply,
                             strip_data = T,
                             skip_legacy = F,
                             ...) {

  dots <- list(...)

  tryCatch(expr = {
    out <- wsx_get_popstats2(ws = ws,
                             groups = groups,
                             invert_groups = invert_groups,
                             return_stats = return_stats,
                             strip_data = strip_data)
  },
  error = function(err) {
    message("error in wsx_get_popstats2:")
    print(err)
    message("try legacy method.")
    skip_legacy <- F
  })

  if (skip_legacy) {
    return(out)
  }

  suppressMessages(expr = {
    out_legacy <- wsx_get_popstats_legacy(ws = ws,
                                          groups = groups,
                                          invert_groups = invert_groups,
                                          return_stats = return_stats,
                                          strip_data = strip_data,
                                          lapply_fun = lapply_fun,
                                          ...)
  })

  cols <- c("FileName", "PopulationFullPath", "Population", "Count", "ParentCount", "identity")
  res <- waldo::compare(x <- dplyr::select(out_legacy[["counts"]], dplyr::all_of(cols)) |> dplyr::arrange(FileName,identity,PopulationFullPath,Count),
                        y <- dplyr::select(out[["counts"]], dplyr::all_of(cols)) |> dplyr::arrange(FileName,identity,PopulationFullPath,Count))

  # out1 <- out[["counts"]] |> dplyr::filter(grepl("0041", FileName))
  # out2 <- out_legacy[["counts"]] |> dplyr::filter(grepl("0041", FileName))

  # zz <- dplyr::left_join(dplyr::filter(x, grepl("0022", FileName)),
  #                        dplyr::filter(y, grepl("0022", FileName)), by = c("FileName", "identity", "Count"))

  if (!length(res)) {
    return(out)
  } else {
    message("Found differences between new and legacy method. Returning legacy.
            You may report that workspace to vonskopnik@pm.me.")
    return(out_legacy)
  }
}

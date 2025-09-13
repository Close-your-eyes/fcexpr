#' Filter observations from flowframes
#'
#' Remove extreme values that may disturb dimension reduction or plotting
#'
#' @param ff flow frame or list of such
#' @param channels channels to use for filtering, goes into ff_get_channels
#' @param channels_select select channels from filtering? (remove other)
#' @param limits absolute limits for filtering
#' @param limits_include > or >=
#' @param quantiles quantiles to use for filtering (relative by channel and ff)
#'
#' @returns reduced ff
#' @export
#'
#' @examples
#' \dontrun{
#' fflist <- fcexpr::wsp_get_ff(w1,
#'                              FCS.file.folder = file.path(wd, "FCS_files"),
#'                              samples = grep("full", basename(f1[[1]]), value = T),
#'                              population = "Time")
#' ffs <- purrr::map(fflist[[1]], ~.x[[1]][[1]])
#' ffs <- ff_filter_obs(ffs, quantiles = c(0.001,0.999))
#' }
ff_filter_obs <- function(ff,
                          channels = NULL,
                          channels_select = F,
                          limits = c(-Inf, Inf),
                          limits_include = c(F, F),
                          quantiles = c(0, 1)) {
  if (!is.list(ff)) {
    ff <- list(ff)
  }

  compareMin <- if (limits_include[which.min(limits)]) `>=` else `>`
  compareMax <- if (limits_include[which.max(limits)]) `<=` else `<`

  summary <- purrr::map_dfr(ff, ~ data.frame(n_before = nrow(.x@exprs)), .id = "ff")

  ff <- purrr::map(ff, function(f) {
    chann <- ff_get_channels(ff = f, channels = channels)
    for (i in chann) {
      f@exprs <- f@exprs[which(compareMin(f@exprs[, i], min(limits)) & compareMax(f@exprs[, i], max(limits))), ]
    }

    f@exprs <- brathering::quantile_filter2(
      x = f@exprs,
      columns = chann,
      columns_select = channels_select,
      quantiles = quantiles
    )

    return(f)
  })


  ff <- ff_update(ff)


  summary <- summary |>
    dplyr::left_join(purrr::map_dfr(ff, ~ data.frame(n_after = nrow(.x@exprs)), .id = "ff"), by = "ff") |>
    dplyr::mutate(diff_abs = n_before - n_after) |>
    dplyr::mutate(diff_rel = round(diff_abs / n_before, 3)) |>
    tibble::as_tibble()

  print(summary, n = Inf)
  attr(ff, "summary") <- summary
  return(ff)
}

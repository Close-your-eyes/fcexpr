#' Narrow channels of stained antibody
#'
#' @param popstats_df data frame from wsx_get_popstats
#'
#' @returns data frame
#' @export
#'
#' @examples
#' df <- wsx_get_popstats(ws_path, strip_data = F)[[1]]
#' df2 <- guess_marker_channels(df)
guess_marker_channels <- function(popstats_df) {
  ps2 <- popstats_df |>
    dplyr::select(Population, xChannel, yChannel) |>
    tidyr::pivot_longer(-Population) |>
    tidyr::drop_na() |>
    dplyr::filter(!grepl("SSC|FSC|HDR-T|Time", value, ignore.case = T)) |>
    dplyr::distinct() |>
    dplyr::mutate(marker = strsplit(Population, "[^A-Za-z0-9]+")) |>
    tidyr::unnest(marker) |>
    dplyr::select(-Population) |>
    dplyr::distinct() |>
    dplyr::add_count(value)
  ps21 <- ps2 |>
    dplyr::filter(n == 1)
  ps22 <- ps2 |>
    dplyr::anti_join(ps21, by = "marker")
  ps3 <- dplyr::bind_rows(ps21, ps22) |>
    dplyr::mutate(ind = dplyr::row_number()) |>
    dplyr::select(-name)
  return(ps3)
}

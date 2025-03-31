#' Channel names / description lookup vector
#'
#' All gating hierarchies (gh) in gs are checked for equal channel descriptions.
#' For channel names with unique desc across all gh, desc are returned,
#' otherwise the channel name.
#'
#' @param gs a gatingset
#'
#' @return a named vector
#' @export
#'
#' @examples
gs_get_descname_lookup <- function(gs) {

  pars <- purrr::map_dfr(stats::setNames(flowWorkspace::sampleNames(gs), flowWorkspace::sampleNames(gs)), function(x) {
    cf <- flowWorkspace::gh_pop_get_data(gs[[x]])
    pars <- tibble::as_tibble(flowCore::parameters(cf)@data[,c(1,2)])
    pars <- dplyr::mutate(pars, dplyr::across(dplyr::everything(), as.character))
    return(pars)
  }, .id = "identity")
  uniquedesc_names <-
    pars |>
    dplyr::filter(!grepl("fsc|ssc|time", name, ignore.case = T)) |>
    dplyr::group_by(name) |>
    dplyr::summarise(ndesc = dplyr::n_distinct(desc)) |>
    dplyr::filter(ndesc == 1) |>
    dplyr::pull(name)
  parsselect <-
    pars |>
    dplyr::select(name, desc) |>
    dplyr::filter(name %in% uniquedesc_names) |>
    dplyr::distinct()

  psel1 <- parsselect |> dplyr::filter(!is.na(desc))
  psel2 <- parsselect |> dplyr::filter(is.na(desc))
  psel3 <- grep("fsc|ssc|time", unique(pars$name), ignore.case = T, value = T)
  descname_lookup <- c(stats::setNames(psel3, psel3),
                       stats::setNames(psel1$desc, psel1$name),
                       stats::setNames(psel2$name, psel2$name))

  return(descname_lookup)
}



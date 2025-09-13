#' Transform a flow frame
#'
#' Intended to switch between inverse.transform = T or F as returned by
#' flowWorkspace::gh_pop_get_data.
#'
#' @param ff flow frame or a list of flow frames
#' @param trafolist matching transformation list
#'
#' @return flow frame
#' @export
#'
#' @examples
#' \dontrun{
#' fflist <- fcexpr::wsp_get_ff(ws)
#' ff <- fflist[["flowframes"]][[1]][[1]][["transformed"]]
#' fflist[["flowframes"]][[1]][[1]][["untransformed"]] <- ff_transform(ff = ff, trafollist = attr(ff, "trafolistinv"))
#' ff_un <- fflist[["flowframes"]][[1]][[1]][["untransformed"]]
#' fflist[["flowframes"]][[1]][[1]][["transformed2"]] <- ff_transform(ff = ff_un, trafollist = attr(ff, "trafolist"))
#' }
ff_transform <- function(ff, trafolist) {

  if (is.list(ff)) {
    if (!all(purrr::map_lgl(ff, ~trafolist %in% names(attributes(.x))))) {
      stop("trafolist not found in all attributes of ff list.")
    }

    ff <- mapply(FUN = ff_transform,
                 ff = ff,
                 trafolist = lapply(ff, attr, which = trafolist))
    return(ff)
  } else {

    # from internal flowWorkspace::gh_pop_get_data
    trafolist[[1]]@transforms <- trafolist[[1]]@transforms[colnames(ff@exprs)]
    cf <- flowWorkspace:::flowFrame_to_cytoframe(ff)
    cs <- flowWorkspace::cytoset()
    flowWorkspace::cs_add_cytoframe(cs, names(trafolist), cf)
    cs2 <- flowWorkspace::gs_cyto_data(flowWorkspace:::gs_clone(cs))
    ff <- flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::transform(cs2, trafolist)[[1]])
    return(ff)
  }
}

#' Try to merge gatingsets by removing non-common populations/gates
#'
#' Uses flowWorkspace::merge_list_to_gs
#'
#' @param gs_list list of gatingsets
#'
#' @returns gatingset
#' @export
#'
#' @examples
gslist_try_merge <- function(gs_list) {
  message("Trying to merge ", length(gs_list), " gs by removing non-common populations.\n")

  gs_pop <- lapply(gs_list, function(x) flowWorkspace::gh_get_pop_paths(x[[1]]))
  common_pop <- Reduce(intersect, gs_pop)
  rm_pop <- purrr::map(gs_pop, setdiff, common_pop)

  message("common pops:\n", paste(common_pop, collapse = "\n"))
  message("\n")
  message("Non-common pops to be removed: \n", paste(names(rm_pop), rm_pop, collapse = "\n\n", sep = "\n"))
  message("\n")


  purrr::map2(gs_list, rm_pop, function(x,y) {
    # go reverse, because rm of a parent gate will rm its children
    for (i in rev(y)) {
      try(expr = {
        flowWorkspace::gs_pop_remove(x, i)
      })
    }
  })

  gs_list <- list(flowWorkspace::merge_list_to_gs(gs_list))
  return(gs_list)
}

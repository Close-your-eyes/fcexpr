#' Title
#'
#' @param file_paths
#'
#' @return
#' @export
#'
#' @examples
fcs_get_keywords2 <- function(file_paths,
                              keywords = NULL) {
  kwl <- tryCatch(
    expr = {
      # unknown error in c++ with some fcs files once
      # setting emptyvalue to F avoided errors
      # or specifically defining which keywords needed
      flowCore::read.FCSheader(files = file_paths,
                               keyword = keywords,
                               cpp = T,
                               emptyValue = F)
    },
    error = function(e) {

      message("Error in reading FCS headers from a FCS files at once:")
      print(e)
      message("Trying to read headers one by one in order to skip faulty FCS files.")
      kwl <- purrr::map(stats::setNames(file_paths, file_paths), function(x) {
        tryCatch(
          expr = {
            flowCore::read.FCSheader(x)[[1]]
          },
          error = function(e) {
            #print(e)
            return(NULL)
          }
        )
      }, .progress = T)
      kwl[which(purrr::map_lgl(kwl, ~ !is.null(.x)))]
      #return(kwl)
    }
  )
  return(kwl)
}

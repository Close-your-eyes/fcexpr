#' Obtain keywords from FCS files
#'
#' Simple wrapper around flowCore::read.FCSheader with some error handling.
#'
#' @param file_paths paths to fcs files
#' @param keywords selection of keywords to extract
#' @param return return format
#'
#' @return data frame of list of named vectors
#' @export
#'
#' @examples
fcs_get_keywords <- function(file_paths,
                             keywords = NULL,
                             return = c("data.frame", "vector")) {

  return <- rlang::arg_match(return)

  kwl <- tryCatch(
    expr = {
      # unknown error in c++ with some fcs files once
      # setting emptyvalue to F avoided errors
      # or specifically defining which keywords needed
      flowCore::read.FCSheader(files = file_paths,
                               keyword = keywords,
                               cpp = T, # the cpp function (default) is much quicker
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
  kwl <- purrr::map(kwl, trimws)

  if (return == "data.frame") {
    kwl <-
      purrr::map_dfr(kwl, stack, .id = "FilePath") |>
      dplyr::mutate(FileName = basename(FilePath)) |>
      dplyr::rename("value" = values, "name" = ind) |>
      tidyr::pivot_wider()
  }

  # if ("SPILL" %in% names(k)) {
  #   k[["SPILL"]] <- paste(as.character(k[["SPILL"]]), collapse = ",")
  # }
  return(kwl)
}

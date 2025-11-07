#' Read FCS file with flowCore and non_modifying defaults
#'
#' This simplest wrapper to imagine just to pass default arguments which
#' avoid any data manipulation
#'
#' @param file_path path to to fcs file
#' @param truncate_max_range see flowCore::read.FCS
#' @param transformation see flowCore::read.FCS
#' @param ... args to flowCore::read.FCS
#'
#' @returns
#' @export
#'
#' @examples
read_fcs <- function(file_path,
                     truncate_max_range = F,
                     transformation = F,
                     ...) {


  flowCore::read.FCS(filename = file_path,
                     truncate_max_range = truncate_max_range,
                     transformation = transformation,
                     ...)

}

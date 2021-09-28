#' Title
#'
#' @param path character
#' @param name character
#' @param date_prefix logical
#'
#' @return No return value. Instead a folder is written to disk.
#' @export
#'
#' @examples
#' \dontrun{
#' new_exp(path = '/Users/CMS/Documents/experiments', name = 'CD3_titration', date_prefix = T)
#' }
new_exp <- function(path = NULL, name = NULL, date_prefix = T) {

  if (base::is.null(path)) {
    stop("Please provide a directory (path) to create the folder in.")
  }

  if (date_prefix) {
    if (base::is.null(name)) {
      name <- base::paste0(base::gsub("-", "", base::Sys.Date()), "_experiment")
    } else {
      name <- base::paste0(base::gsub("-", "", base::Sys.Date()), "_", name)
    }
  } else {
    if (base::is.null(name)) {
      name <- "experiment"
    }
  }

  utils::untar(base::system.file("extdata", "template_folder.tgz", package = "fcexpr"), exdir = path)
  base::file.rename(base::file.path(path, "template_folder"), base::file.path(path, name))

  base::print(base::paste0(base::file.path(path, name), " created."))
}

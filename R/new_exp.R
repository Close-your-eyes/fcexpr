#' Create a template folder for flow cytometry experiments
#'
#' @param path path to the parent directory where the experiment folder is to be initiated
#' @param name name of the experiment folder
#' @param date_prefix logical if the folder name should be prefixed with the current date
#'
#' @return No return value. Instead a folder with template files is written to disk.
#' @export
#'
#' @examples
#' \dontrun{
#' new_exp(path = '/Users/CMS/Documents/experiments', name = 'CD3_titration')
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

  if (dir.exists(base::file.path(path, name))) {
    stop(paste0(base::file.path(path, name), " already exists."))
  }

  utils::untar(base::system.file("extdata", "template_folder.tgz", package = "fcexpr"), exdir = path)
  base::file.rename(base::file.path(path, "template_folder"), base::file.path(path, name))

  base::print(base::paste0(base::file.path(path, name), " created."))
}

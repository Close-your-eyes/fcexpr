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
new_exp <- function(path = getwd(), name = "fcexp", date_prefix = T) {

  if (is.null(path) || path == "") {
    stop("Please provide a directory (path) to create the folder in.")
  }
  if (is.null(name)) {
    stop("Please provide a name.")
  }

  path <- path.expand(path)

  if (date_prefix) {
    if (is.null(name)) {
      suffix <- "_experiment"
    } else {
      suffix <- paste0("_", name)
    }
    name <- paste0(gsub("-", "", format(Sys.Date(), "%y-%m-%d")), suffix)
  }

  if (dir.exists(file.path(path, name))) {
    stop(file.path(path, name), " already exists.")
  }

  utils::untar(system.file("extdata", "template_folder.tgz", package = "fcexpr"), exdir = path)
  file.rename(file.path(path, "template_folder"), file.path(path, name))

  if (getOS() == "Windows") {
    files <- list.files(file.path(path, name), all.files = T, recursive = T, full.names = T)
    file.remove(files[which(grepl("^\\.", basename(files)))])
  }

  message(file.path(path, name), " created.")
  invisible(file.path(path, name))
}


getOS <- function() {
  machine <- switch(Sys.info()[["sysname"]],
                    Windows= "Windows",
                    Linux  = "Linux",
                    Darwin = "Mac")
  return(machine)
}


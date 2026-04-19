#' Generate Unique FCS Identities from Keyword Lists
#'
#' Constructs unique identifiers for FCS files based on keyword metadata
#' (e.g., acquisition date, time, file name, and total events).
#' The function standardizes date-time formats and concatenates relevant
#' fields into a single identity string per file.
#'
#' @param kwl A **named list** of keyword lists, typically obtained from
#'   \code{flowCore::read.FCSheader()} or workspace-derived keyword extraction.
#'   Each element must contain at least the following keys:
#'   \code{"$DATE"}, \code{"$BTIM"}, \code{"$ETIM"}, \code{"$FIL"}, \code{"$TOT"}.
#'
#' @param allow_duplicates Logical. Whether duplicate FCS identities are allowed.
#' @param folder_path alternative to kwl: provide folder with fcs files
#' @param fcs_file_paths alternative to kwl: provide fcs files paths on disk
#' @param exclude_folders exclude sub-folders in folder_path, can be char vector
#' @param recursive scan folder_path recursively for fcs files?
#'
#' @returns A named character vector of FCS identities. Each element corresponds
#'   to an input file and is formatted as:
#'   \preformatted{
#'   <FIL>_-_<TOT>_-_<YYYY.MM.DD-HH.MM.SS>
#'   }
#'
#' @details
#' The function:
#' \itemize{
#'   \item Validates that \code{kwl} is a named list with unique names
#'   \item Extracts required FCS keywords
#'   \item Normalizes time strings (handles malformed time formats)
#'   \item Parses and standardizes date-time using \pkg{lubridate}
#'   \item Adjusts for overnight acquisitions (end time past midnight)
#'   \item Concatenates metadata into a unique identity string
#' }
#'
#' If date-time parsing fails for any entry, a warning is issued.
#'
#' @seealso \code{\link[flowCore]{read.FCSheader}}
#'
#' @examples
#' \dontrun{
#' library(flowCore)
#'
#' # Example: read FCS headers
#' fcs_files <- list.files(path = "path/to/fcs", pattern = ".fcs$", full.names = TRUE)
#' kwl <- flowCore::read.FCSheader(fcs_files, emptyValue = FALSE)
#'
#' ids <- get_fcs_identities(kwl = kwl)
#' head(ids)
#'
#' # Example from workspace-derived keywords
#' # kwl <- wsx_get_keywords(ws = ws, return = "vector")
#' # ids <- get_fcs_identities(kwl = kwl)
#' }
#'
#' @export
get_fcs_identities <- function(kwl = NULL,
                               folder_path = NULL,
                               fcs_file_paths = NULL,
                               exclude_folders = NULL,
                               recursive = T,
                               allow_duplicates = T) {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }

  # list_fcs_files()
  # read_fcs_headers()
  # compute_fcs_identities()
  # validate_fcs_identities()

  if (!is.null(folder_path)) {
    fcs_file_paths <- list_fcs_files(folder_path = folder_path[1],
                                     exclude_folders = exclude_folders,
                                     recursive = recursive)

  }

  if (!is.null(fcs_file_paths)) {
    if (any(!file.exists(fcs_file_paths))) {
      stop("some fcs_file_paths not found.")
    }
    kwl = flowCore::read.FCSheader(fcs_file_paths, emptyValue = F)
  }

  fcs_identities <- get_idents_from_kwl(kwl = kwl)

  if (!allow_duplicates) {
    check_fcs_ident_uniqueness(x = fcs_identities)
  }

  # if (!allow_duplicates && length(unique(fcs_identities)) != length(fcs_identities)) {
  #   stop(paste0("Duplicate FCS files found. This is not allowed. Please, remove one of each duplicates. \n", paste(names(fcs_identities[duplicated(fcs_identities) |
  #                                                                                                                                    duplicated(fcs_identities, fromLast = T)]), collapse = "\n")))
  # } else if (allow_duplicates && length(unique(fcs_identities)) != length(fcs_identities)) {
  #   message(paste0("Duplicate FCS files found. Be careful. \n", paste(names(fcs_identities[duplicated(fcs_identities) |
  #                                                                                                                                         duplicated(fcs_identities, fromLast = T)]), collapse = "\n")))
  #
  # }
  return(fcs_identities)
}

get_idents_from_kwl <- function(kwl) {
  if (is.null(kwl)) {
    stop("provide kwl, folder_path or fcs_file_paths.")
  }
  if (!methods::is(kwl, "list")) {
    stop("keyword list not a list. Could be made a list but then names are missing. Try to fix.")
  }
  if (is.null(names(kwl))) {
    stop("keyword list needs names")
  }
  if (any(duplicated(names(kwl)))) {
    stop("names of keyword list are not unique")
  }

  dd <- sapply(kwl, "[", "$DATE")
  tt <- sapply(kwl, "[", "$BTIM")
  et <- sapply(kwl, "[", "$ETIM")
  fil <- sapply(kwl, "[", "$FIL")
  tot <- sapply(kwl, "[", "$TOT")

  if (any(nchar(tt) - nchar(gsub(":", "", tt)) > 2)) {
    tt_fix_ind <- which(nchar(tt) - nchar(gsub(":", "", tt)) > 2)
    tt[tt_fix_ind] <- paste(rev(rev(strsplit(tt[tt_fix_ind], ":")[[1]])[-1]), collapse = ":")
  }
  datetime <- paste0(dd, "-", tt)
  ## fix datetime when analysis started before 00:00 and ended after 00:00
  sub <- ifelse(grepl("^2[[:digit:]]", tt) & grepl("^0[[:digit:]]", et), 86400, 0)
  datetime <- format(lubridate::parse_date_time(datetime, orders = c("%Y-%b-%d-%H:%M:%S", "%Y-%B-%d-%H:%M:%S", "%Y-%m-%d-%H:%M:%S", "%d-%b-%Y-%H:%M:%S",
                                                                     "%d-%m-%Y-%H:%M:%S", "%d-%B-%Y-%H:%M:%S", "%d-%b-%Y-%H:%M:%S")) - sub, "%Y.%m.%d-%H.%M.%S")
  if (any(is.na(datetime))) {
    warning("datetimes ", paste(paste0(dd, "-", tt)[which(is.na(datetime))], collapse = ", "), " could not be converted to a uniform format. Please, provide this to the package-maintainer.")
  }

  fcs_identities <- stats::setNames(paste0(fil, "_-_", trimws(tot), "_-_", datetime), nm = names(kwl))
  return(fcs_identities)
}

check_fcs_ident_uniqueness <- function(x) {

  df <-
    utils::stack(x) |>
    dplyr::mutate(ind = basename(as.character(ind))) |>
    dplyr::rename("identity" = values, "FileName" = ind) |>
    dplyr::add_count(identity, name = "same_identity") |>
    dplyr::add_count(identity, FileName, name = "same_identity_and_filename") |>
    tibble::as_tibble()

  if (any(df$same_identity > 1)) {
    message("Duplicate identities found. FCS files were duplicated on disk.")
    print(dplyr::filter(df, same_identity > 1), n = Inf)
  }

  if (any(df$same_identity_and_filename > 1)) {
    print(dplyr::filter(df, same_identity_and_filename > 1), n = Inf)
    stop("Duplicate identities with equal FileName are not allowed. Rename files to make them unique.")
  }
}

list_fcs_files <- function(folder_path,
                           exclude_folders,
                           recursive) {
  if (!dir.exists(folder_path[1])) {
    stop("folder_path not found.")
  }
  fcs_file_paths <- list.files(
    path = folder_path[1],
    pattern = "\\.fcs", #w/o $ is fine
    full.names = T,
    recursive = recursive,
    ignore.case = T
  )
  if (length(fcs_file_paths) == 0) {
    stop("No FCS files found.")
  }

  if (!is.null(exclude_folders)) {
    exclude.pattern <- paste0(tolower(exclude_folders), collapse = "|")
    fcs_file_paths <- fcs_file_paths[!grepl(exclude.pattern, tolower(fcs_file_paths))]
  }

  if (length(fcs_file_paths) == 0) {
    stop("No FCS files left after filtering for exclusion folders.")
  }
  return(fcs_file_paths)
}

#' Obtain PMT-voltages from a FCS file
#'
#' If there is a large number of fcs files and one wants to make sure that all files
#' have been acquired with the same settings this functions may help by retrieving the
#' voltage of PMTs.
#'
#' @param file_path path to the fcs file
#'
#' @return a data.frame with different columns depending on the machine the fcs was acquired with;
#' usually the 'V'-column indicates PMT voltages
#' @export
#'
#' @examples
#' \dontrun{
#' # When the script is saved to R_scripts in the experiment folder,
#' get the absolute path to the folder
#' wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
#' # Find workspaces
#' ws <- list.files(wd, pattern = '\\.wsp$', recursive = T, full.names = T)
#' # pull out paths to fcs files, optionally filter for relevant groups or so
#' paths <- wsx_get_fcs_paths(ws[1])
#' # iterate (loop) through the paths to get the voltage of each fcs file
#' volts <- do.call(rbind, lapply(paths, function(x) {
#' v <- voltage_from_fcs(x)
#' v$file <- basename(x)
#' return(v)
#' }))
#' }
fcs_get_voltages <- function(file_path) {
    if (!requireNamespace("BiocManager", quietly = T)){
        utils::install.packages("BiocManager")
    }
    if (!requireNamespace("flowCore", quietly = T)){
        BiocManager::install("flowCore")
    }

    if (!any(file.exists(file_path))) {
        warning("Not all files found.")
        file_path <- file_path[which(file.exists(file_path))]
        if (length(file_path) == 0) {
            stop("None of files found.")
        }
    }

    ff <- flowCore::read.FCSheader(file_path)
    names(ff) <- basename(file_path)
    out <- do.call(rbind, lapply(names(ff), function(x) {
        y <- ff[[x]]
        y <- y[which(names(y) != "SPILL")]
        f <- suppressWarnings(utils::stack(y))

        f[, "ind"] <- as.character(f[, "ind"])
        f <- f[which(grepl("\\$P[[:digit:]]", f[, "ind"]) & !grepl("flowCore", f[, "ind"])), ]
        f[, "ind"] <- gsub("\\$", "", f[, "ind"])
        f[, "type"] <- sapply(sapply(strsplit(f[, "ind"], ""), rev), "[", 1)
        f[, "ind"] <- gsub("[[:alpha:]]$", "", f[, "ind"])

        f <- tidyr::pivot_wider(f, names_from = "type", values_from = "values")
        f <- f[which(f[, "N"] != "Time"), ]
        f[,"FileName"] <- x
        f <- f[order(as.numeric(gsub("[^0-9.]", "", f$ind))),]
        return(f)
    }))

    return(as.data.frame(out))
}

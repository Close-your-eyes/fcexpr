#' Get voltages of PMTs from a fcs file
#'
#' If there is a large number of fcs files and one wants to make sure that all files
#' have been acquired with the same settings this functions helps to pull out the
#' voltage of PMTs.
#'
#' @param file_path character of path to the fcs file
#'
#' @return a data.frame with different columns depending on the machine the fcs was acquired with; usually the 'V'-column indicates PMT voltages
#' @export
#'
#' @examples
#' \dontrun{
#' # When the script is saved to R_scripts in the experiment folder,
#' get the absolute path to the folder
#' wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
#' # Find workspaces
#' ws <- list.files(wd, pattern = '\\.wsp$', recursive = T, full.names = T)
#' # Get groups from the first ws
#' gr <- lapply(ws[1], function(x) {
#' unique(as.character(CytoML::fj_ws_get_sample_groups(CytoML::open_flowjo_xml(x))$groupName))
#' })[[1]]
#' # pull out paths to fcs files of the first group
#' paths <- fcs_paths_from_ws(ws[1], gr[1])
#' # iterate (loop) through the paths to get the voltage of each fcs file
#' volts <- do.call(rbind, lapply(paths, function(x) {
#' v <- voltage_from_fcs(x)
#' v$file <- basename(x)
#' return(v)
#' }))
#' }
fcs_get_voltages <- function(file_path) {

    if (!file.exists(file_path)) {
        stop(paste0(file_path, " not found."))
    }

    ff <- flowCore::read.FCS(file_path, which.lines = 1, emptyValue = F, truncate_max_range = F)
    f <- suppressWarnings(utils::stack(flowCore::keyword(ff)[names(flowCore::keyword(ff)) != "SPILL"]))
    f[, "ind"] <- as.character(f[, "ind"])
    f <- f[which(grepl("\\$P[[:digit:]]", f[, "ind"]) & !grepl("flowCore", f[, "ind"])), ]
    f[, "ind"] <- gsub("\\$", "", f[, "ind"])
    f[, "type"] <- sapply(sapply(strsplit(f[, "ind"], ""), rev), "[", 1)
    f[, "ind"] <- gsub("[[:alpha:]]$", "", f[, "ind"])

    f <- tidyr::pivot_wider(f, names_from = "type", values_from = "values")
    f <- f[which(f[, "N"] != "Time"), ]

    return(as.data.frame(f))
}

#' Get the paths to fcs files within a group of a workspace
#'
#' @param ws character of path to flowjo workspaces
#' @param gr character of group name to consider
#' @param FCS.file.folder path to root folder which contains FCS files; if missing file.path(getwd(), 'FCS_files') is assumed
#'
#' @return character vector of absolute paths to fcs files
#' @export
#'
#' @examples
#' \dontrun{
#' # When the script is saved to R_scripts in the experiment folder,
#' # get the absolute path to the folder
#' wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
#' # find workspaces
#' ws <- list.files(wd, pattern = '\\.wsp$', recursive = T, full.names = T)
#' # get groups from the first ws
#' gr <- lapply(ws[1], function(x) {
#' unique(as.character(CytoML::fj_ws_get_sample_groups(CytoML::open_flowjo_xml(x))$groupName))
#' })[[1]]
#' # pull out paths to fcs files of the first group
#' paths <- fcs_paths_from_ws(ws[1], gr[1])
#' }
ws_get_fcs_paths <- function(ws, gr, FCS.file.folder) {

    if (missing(FCS.file.folder)) {
        FCS.file.folder <- file.path(getwd(), "FCS_files")
    }
    if (!dir.exists(FCS.file.folder)) {
        stop(paste0(FCS.file.folder, " not found."))
    }

    gs <- CytoML::flowjo_to_gatingset(CytoML::open_flowjo_xml(ws), name = gr, path = FCS.file.folder, execute = T, which.lines = 1)

    paths <- sapply(1:length(gs), function(z) {
        flowCore::keyword(flowWorkspace::gh_pop_get_data(gs[[z]]))[["FILENAME"]]
    })
    if (all(grepl("^\\.", paths))) {
        paths <- file.path(getwd(), gsub("^\\.\\/", "", paths))
    }

    return(paths)
}

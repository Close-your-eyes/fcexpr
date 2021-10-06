#' Retrieve MFIs from populations gated in one or more flowjo workspaces
#'
#' If inverse.transform is set to TRUE values as shown by flowjo will be used
#' for calculation. This may omit to use 'geo_mean' for mean.fun as negative
#' values are not allowed. Inverse.transform = F leads to data transformation by
#' the logicle transform.
#'
#' @param ws character vector of paths to flowjo workspaces
#' @param gr character vector of flowjo groups to import or a list of those vectors, one for each ws
#' @param population character vector of which populations to calculate values for; if omitted calculation is done for all population which increases computational time
#' @param FCS.file.folder path to root folder which contains FCS files; if missing file.path(getwd(), 'FCS_files') is assumed
#' @param inverse.transform logical indicating if fluorescence values are to be inversely transformed
#' @param mean.fun character name of the function name to use to calculate the average FI
#' @param variance.fun  character name of the function name to use to calculate the variacnce of FI
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' # When the script is saved to R_scripts in the experiment folder,
#' # get the absolute path to the folder
#' wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
#' # find workspaces
#' ws <- list.files(wd, pattern = '\\.wsp$', recursive = T, full.names = T)
#' # get groups
#' gr <- lapply(ws, function(x) {
#' unique(as.character(CytoML::fj_ws_get_sample_groups(CytoML::open_flowjo_xml(x))$groupName))
#' })[[1]]
#' # read the compensated values from flowjo and calculate MFIs
#' ws_get_MFIs(ws = ws, gr = gr)
#' }
ws_get_MFIs <- function(ws,
                        gr,
                        population = NULL,
                        FCS.file.folder = file.path(getwd(), "FCS_files"),
                        mean.fun = "median",
                        variance.fun = "sd",
                        inverse.transform = T) {

    if (missing(ws) || class(ws) != "character") {
        stop("Please provide a vector of paths to flowjo workspaces.")
    }
    if (missing(gr)) {
        stop("Please provide a vector of list of vectors of groups to import.")
    }
    if (class(gr) == "list" && (length(gr) != 1 & length(gr) != length(ws))) {
        stop("list of gr has to have length 1 or the same length as ws.")
    }
    if (class(gr) == "character") {
        gr <- rep(list(gr), length(ws))
    } else if (class(gr) == "list" && length(gr) == 1) {
        gr <- rep(gr, length(ws))
    }
    if (class(population) == "list" && (length(population) != 1 & length(population) != length(ws))) {
        stop("list of population has to have the same length as ws.")
    }
    if (class(population) == "character") {
        population <- rep(list(population), length(ws))
    } else if (class(population) == "list" && length(population) == 1) {
        population <- rep(population, length(ws))
    }
    if (!dir.exists(FCS.file.folder)) {
        stop(paste0(FCS.file.folder, " not found."))
    }

    mean.fun.fun <- match.fun(mean.fun)
    variance.fun.fun <- match.fun(variance.fun)

    MFI.table <- do.call(rbind, lapply(seq_along(ws), function(x) {
        wsp <- CytoML::open_flowjo_xml(ws[x])
        gr.wsp <- base::intersect(unique(CytoML::fj_ws_get_sample_groups(wsp)[, "groupName"]), gr[[x]])
        out <- do.call(rbind, lapply(gr.wsp, function(y) {
            gs <- CytoML::flowjo_to_gatingset(ws = wsp, name = y, path = FCS.file.folder, emptyValue = F, truncate_max_range = F, additional.keys = c())
            if (is.null(population)) {
                population.wsp <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")
            } else {
                population.wsp <- base::intersect(flowWorkspace::gs_get_pop_paths(gs, path = "auto"), population[[x]])
            }
            out <- do.call(rbind, lapply(population.wsp, function(p) {
                out <- do.call(rbind, lapply(seq_along(gs), function(g) {
                  dat <- flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[g]], p, inverse.transform = inverse.transform))
                  t <- cbind(utils::stack(apply(dat, 2, mean.fun.fun)), utils::stack(apply(dat, 2, variance.fun.fun))[, 1])
                  names(t) <- c(mean.fun, "channel", variance.fun)
                  t[, "channel.desc"] <- flowCore::parameters(flowWorkspace::gh_pop_get_data(gs[[g]], p))[["desc"]]
                  t[, "FileName"] <- basename(flowCore::keyword(flowWorkspace::gh_pop_get_data(gs[[g]]))[["FILENAME"]])
                  return(t)
                }))
                out[, "Population"] <- p
                return(out)
            }))
            out[, "group"] <- y
            return(out)
        }))
        out[, "ws"] <- basename(ws[x])
        return(out)
    }))
    MFI.table <- MFI.table[, c(5, 6, 2, 4, 1, 3, 7, 8)]
    return(MFI.table)
}

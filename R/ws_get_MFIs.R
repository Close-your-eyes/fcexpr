#' Retrieve MFIs from populations gated in one or more flowjo workspaces
#'
#' If inverse.transform is set to TRUE values as shown by flowjo will be used
#' for calculation. This may omit to use 'geo_mean' for mean.fun as negative
#' values are not allowed. Inverse.transform = F leads to data transformation by
#' the logicle transform.
#'
#' @param ws character vector of paths to flowjo workspaces
#' @param gr character vector of flowjo groups to import or a list of those vectors, each for one ws
#' @param FCS.file.folder path to root folder which contains FCS files; if missing file.path(getwd(), "FCS_files") is assumed
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
#' ws <- list.files(wd, pattern = "\\.wsp$", recursive = T, full.names = T)
#' # get groups
#' gr <- lapply(ws, function(x) {
#' unique(as.character(CytoML::fj_ws_get_sample_groups(CytoML::open_flowjo_xml(x))$groupName))
#' })
#' # read the compensated values from flowjo and calculate MFIs
#' ws_get_MFIs(ws = ws, gr = gr)
#' }
ws_get_MFIs <- function(ws,
                        gr,
                        FCS.file.folder,
                        mean.fun = "median",
                        variance.fun = "sd",
                        inverse.transform = T) {

  if (missing(ws) || class(ws) != "character") {
    stop("Please provide a vector of paths to flowjo workspaces.")
  }
  if (missing(gr)) {
    stop("Please provide a vector of list of vectors of groups to import.")
  }
  if (class(gr) == "list" && length(gr) != length(ws)) {
    stop("list of gr has to have the same length as ws. Alternatively pass a vector gr to use for all workspace.")
  }
  if (class(gr) == "character") {
    gr <- rep(list(gr), length(ws))
  }
  if (missing(FCS.file.folder)) {
    FCS.file.folder <- file.path(getwd(), "FCS_files")
  }
  if (!dir.exists(FCS.file.folder)) {
    stop(paste0(FCS.file.folder, " not found."))
  }

  mean.fun.fun <- match.fun(mean.fun)
  variance.fun.fun <- match.fun(variance.fun)

  MFI.table <- do.call(rbind, lapply(seq_along(ws), function(x) {
    wsp <- CytoML::open_flowjo_xml(ws[x])
    gr.wsp <- base::intersect(unique(CytoML::fj_ws_get_sample_groups(wsp)[,"groupName"]), gr[[x]])
    out <- do.call(rbind, lapply(gr.wsp, function(y) {
      gs <- CytoML::flowjo_to_gatingset(ws = wsp, name = y, path = FCS.file.folder, emptyValue = F, truncate_max_range = F)
      out <- do.call(rbind, lapply(flowWorkspace::gs_get_pop_paths(gs), function(p) {
        out <- do.call(rbind, lapply(1:length(gs), function(g) {
          t <- cbind(utils::stack(apply(flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[g]], p, inverse.transform = inverse.transform)), 2, mean.fun.fun)),
                     utils::stack(apply(flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[g]], p, inverse.transform = inverse.transform)), 2, variance.fun.fun))[,1])
          names(t) <- c(mean.fun, "channel", variance.fun)
          t[,"channel.desc"] <- flowCore::parameters(flowWorkspace::gh_pop_get_data(gs[[g]], p))[["desc"]]
          t[,"FileName"] <- basename(flowCore::keyword(flowWorkspace::gh_pop_get_data(gs[[g]]))[["FILENAME"]])
          return(t)
        }))
        out[,"PopulationFullPath"] <- p
        return(out)
      }))
      out[,"group"] <- y
      return(out)
    }))
    out[,"ws"] <- basename(ws[x])
    return(out)
  }))
  MFI.table <- MFI.table[,c(5,6,2,4,1,3,7,8)]
  return(MFI.table)
}

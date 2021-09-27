#' Title
#'
#' @param ws character
#' @param gr character
#' @param FCS.file.folder character
#' @param groupwise logical
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' ws_get_popstats(ws = "fj.wsp", gr = "full.stain", FCS.file.folder = "FCS.files", groupwise = T)
#' }
ws_get_popstats <- function(ws,
                            gr,
                            FCS.file.folder,
                            groupwise = T) {

  if (missing(ws) || class(ws) != "character") {stop("Please provide a vector of paths to flowjo workspaces.")}
  if (missing(gr)) {stop("Please provide a vector of list of vectors of groups to import.")}
  if (class(gr) == "list" && length(gr) != length(ws)) {stop("list of gr has to have the same length as ws. Alternatively pass a vector gr to use for all workspace.")}
  if (class(gr) != "list") {gr <- rep(list(gr), length(ws))}
  if (missing(FCS.file.folder)) {FCS.file.folder <- file.path(getwd(), "FCS_files")}
  if (!dir.exists(FCS.file.folder)) {stop(paste0(FCS.file.folder, " not found."))}

  popStats.raw <- do.call(rbind, lapply(seq_along(ws), function(x) {
    wsp <- CytoML::open_flowjo_xml(ws[x], sample_names_from = "sampleNode")
    gr.wsp <- base::intersect(unique(CytoML::fj_ws_get_sample_groups(wsp)[,"groupName"]), gr[[x]])
    ps <- if (!groupwise) {
      do.call(rbind, lapply(gr.wsp, function(y) {
        group_id <- unique(CytoML::fj_ws_get_sample_groups(wsp)[which(CytoML::fj_ws_get_sample_groups(wsp)$groupName == y), "groupID"]) + 1
        fcs.files <- CytoML::fj_ws_get_samples(wsp, group_id = group_id)[,"name"]
        ps <- data.frame()
        for (i in 1:length(fcs.files)) {
          ps <- as.data.frame(rbind(ps, flowWorkspace::gs_pop_get_count_fast(CytoML::flowjo_to_gatingset(ws = wsp, name = y, subset = i, path = FCS.file.folder, execute = F, emptyValue = F), path = "full", xml = T)))
        }
        ps[,"Group"] <- y
        return(ps)
      }))
    } else {
      do.call(rbind, lapply(gr.wsp, function(y) {
        gs <- CytoML::flowjo_to_gatingset(wsp, name = y, path = FCS.file.folder, execute = F, emptyValue = F)
        ps <- as.data.frame(flowWorkspace::gs_pop_get_count_fast(gs, path = "full", xml = T))
        ps[,"Population"] <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")[-1]
        ps[,"Group"] <- y
        return(ps)
      }))
    }
    names(ps)[which(names(ps) == "name")] <- "FileName"
    ps[,"PopulationFullPath"] <- gsub("^root", "", paste(ps[,"Parent"], sapply(sapply(base::strsplit(ps[,"Population"], "/"), rev), "[", 1), sep = "/"), "^root")
    ps[,"FractionOfParent"] <- ps[,"Count"]/ps[,"ParentCount"] * 100
    ps[,"wsp"] <- basename(ws[x])
    ps <- ps[,c(1,7,2,3,4,5,8,9,6)]
    return(ps)
  }))

  return(popStats.raw)
}

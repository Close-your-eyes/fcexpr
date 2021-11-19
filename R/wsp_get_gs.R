#' Title
#'
#' @param wsps
#' @param groups
#' @param FCS.file.folder
#'
#' @return
#' @export
#'
#' @examples
wsp_get_gs <- function(wsps,
                       groups,
                       FCS.file.folder) {

  if (!"CytoML" %in% rownames(installed.packages())) {BiocManager::install("CytoML")}
  if (!"flowWorkspace" %in% rownames(installed.packages())) {BiocManager::install("flowWorkspace")}


  if (missing(wsps) || class(wsps) != "character") {stop("Please provide a vector of paths to flowjo workspaces.")}
  if (missing(groups)) {stop("Please provide a vector of list of groups to import.")}
  if (class(groups) == "list" && length(groups) != length(wsps)) {stop("list of groups has to have the same length as wsps. Alternatively pass a vector groups to use for all workspace.")}
  if (class(groups) != "list") {groups <- rep(list(groups), length(wsps))}
  if (!dir.exists(FCS.file.folder)) {stop(paste0(FCS.file.folder, " not found."))}


  wsp.groups <- unlist(lapply(wsps, function(x) {
    wsp <- CytoML::open_flowjo_xml(x)
    sample.groups <- base::intersect(unique(CytoML::fj_ws_get_sample_groups(wsp)$groupName), groups)
    return(setNames(sample.groups, rep(x,length(sample.groups))))
  }))

  gs.list <- pblapply(seq_along(wsp.groups), function(x) {
    gs <- CytoML::flowjo_to_gatingset(CytoML::open_flowjo_xml(names(wsp.groups)[x]), name = wsp.groups[x], path = FCS.file.folder, truncate_max_range = F)
    flowWorkspace::sampleNames(gs) <- sapply(seq_along(gs), function(z) {basename(keyword(flowWorkspace::gh_pop_get_data(gs[[z]]))[["FILENAME"]])})
    return(gs)
  })

  return(setNames(gs.list, paste0(basename(names(wsp.groups)), "_-_", wsp.groups)))
}

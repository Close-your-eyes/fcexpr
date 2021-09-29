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
#' ws_get_popstats(ws = 'fj.wsp', gr = 'full.stain', FCS.file.folder = 'FCS.files', groupwise = T)
#' }
ws_get_popstats <- function(ws, gr, FCS.file.folder, groupwise = T) {

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
    } else {
        stop("gr has to be a character vector.")
    }
    if (missing(FCS.file.folder)) {
        FCS.file.folder <- file.path(getwd(), "FCS_files")
    }
    if (!dir.exists(FCS.file.folder)) {
        stop(paste0(FCS.file.folder, " not found."))
    }

    boic_subversion <- as.numeric(strsplit(as.character(BiocManager::version()), "\\.")[[1]][2])

    if (!groupwise && boic_subversion < 12) {
        stop("Non-groupwise import not possible with the old Bioconductor version.")
    }

    ps <- if (boic_subversion < 12) {

        do.call(rbind, lapply(seq_along(ws), function(x) {
            wsp <- CytoML::open_flowjo_xml(ws[x])
            gr.wsp <- base::intersect(unique(CytoML::fj_ws_get_sample_groups(wsp)[,"groupName"]), gr[[x]])
            do.call(rbind, lapply(gr.wsp, function(y) {
                gs <- CytoML::flowjo_to_gatingset(wsp, name = y, path = FCS.file.folder, execute = T, emptyValue = F, which.lines = 1)
                flowWorkspace::sampleNames(gs) <- sapply(1:length(gs), function(z) {basename(flowWorkspace::keyword(flowWorkspace::gh_pop_get_data(gs[[z]]))[["FILENAME"]])})
                ps <- as.data.frame(flowWorkspace::gs_pop_get_count_fast(gs, path = "full", xml = T))
                ps[,"Population"] <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")[-1]
                ps[,"Group"] <- y
                ps[, "wsp"] <- basename(ws[x])
                return(ps)
            }))
        }))

    } else {

        do.call(rbind, lapply(seq_along(ws), function(x) {
            wsp <- CytoML::open_flowjo_xml(ws[x], sample_names_from = "sampleNode")
            gr.wsp <- base::intersect(unique(CytoML::fj_ws_get_sample_groups(wsp)[, "groupName"]), gr[[x]])
            if (!groupwise) {
                do.call(rbind, lapply(gr.wsp, function(y) {
                    group_id <- unique(CytoML::fj_ws_get_sample_groups(wsp)[which(CytoML::fj_ws_get_sample_groups(wsp)$groupName == y), "groupID"]) +
                        1
                    fcs.files <- CytoML::fj_ws_get_samples(wsp, group_id = group_id)[, "name"]
                    ps <- data.frame()
                    for (i in 1:length(fcs.files)) {
                        ps <- as.data.frame(rbind(ps, flowWorkspace::gs_pop_get_count_fast(CytoML::flowjo_to_gatingset(ws = wsp, name = y, subset = i,
                                                                                                                       path = FCS.file.folder, execute = F, emptyValue = F), path = "full", xml = T)))
                    }
                    ps[, "Group"] <- y
                    ps[, "wsp"] <- basename(ws[x])
                    return(ps)
                }))
            } else {
                do.call(rbind, lapply(gr.wsp, function(y) {
                    gs <- CytoML::flowjo_to_gatingset(wsp, name = y, path = FCS.file.folder, execute = F, emptyValue = F)
                    ps <- as.data.frame(flowWorkspace::gs_pop_get_count_fast(gs, path = "full", xml = T))
                    ps[, "Population"] <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")[-1]
                    ps[, "Group"] <- y
                    ps[, "wsp"] <- basename(ws[x])
                    return(ps)
                }))
            }
        }))

    }

    names(ps)[which(names(ps) == "name")] <- "FileName"
    ps[, "PopulationFullPath"] <- gsub("^root", "", paste(ps[, "Parent"], sapply(sapply(base::strsplit(ps[, "Population"], "/"), rev), "[",
                                                                                 1), sep = "/"), "^root")
    ps[, "FractionOfParent"] <- ps[, "Count"]/ps[, "ParentCount"] * 100
    ps <- ps[,c(1,8,2,3,4,5,9,6,7)]

    return(ps)
}

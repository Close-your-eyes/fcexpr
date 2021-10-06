#' Convenient function to import population counts from one or more flowjo workspaces
#'
#' At least one gate per group in flowjo is required to import.
#' When groupwise is set to TRUE the gating tree for all samples within a group
#' has to be equal (same gating tree structure, though gates can be at different
#' positions for different samples). If this cannot be managed, groupwise has to be set to FALSE which
#' increases computational time to import a bit as samples are imported one by one.
#'
#'
#' @param ws character vector of paths to flowjo workspaces
#' @param gr character vector of flowjo groups to import or a list of those vectors, one for each
#' @param FCS.file.folder path to root folder which contains FCS files; if missing file.path(getwd(), 'FCS_files') is assumed
#' @param groupwise logical indicating if samples are to be imported groupwise or one by one
#'
#' @return returns a data.frame with population counts
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
#' # import the population counts:
#' ws_get_popstats(ws = ws, gr = gr)
#' }
ws_get_popstats <- function(ws, gr, FCS.file.folder, groupwise = T) {

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
    }
    if (class(gr) == "list" && length(gr) == 1) {
        gr <- rep(gr, length(ws))
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
            gr.wsp <- base::intersect(unique(CytoML::fj_ws_get_sample_groups(wsp)[, "groupName"]), gr[[x]])
            do.call(rbind, lapply(gr.wsp, function(y) {
                gs <- CytoML::flowjo_to_gatingset(wsp, name = y, path = FCS.file.folder, execute = F, emptyValue = F, which.lines = 1, additional.keys = c())
                flowWorkspace::sampleNames(gs) <- sapply(1:length(gs), function(z) {
                  basename(flowWorkspace::keyword(flowWorkspace::gh_pop_get_data(gs[[z]]))[["FILENAME"]])
                })
                ps <- as.data.frame(flowWorkspace::gs_pop_get_count_fast(gs, path = "full", xml = T))
                ps[, "Population"] <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")[-1]
                ps[, "group"] <- y
                ps[, "ws"] <- basename(ws[x])
                return(ps)
            }))
        }))

    } else {

        do.call(rbind, lapply(seq_along(ws), function(x) {
            wsp <- CytoML::open_flowjo_xml(ws[x], sample_names_from = "sampleNode")
            gr.wsp <- base::intersect(unique(CytoML::fj_ws_get_sample_groups(wsp)[, "groupName"]), gr[[x]])
            if (!groupwise) {
                do.call(rbind, lapply(gr.wsp, function(y) {
                  group_id <- unique(CytoML::fj_ws_get_sample_groups(wsp)[which(CytoML::fj_ws_get_sample_groups(wsp)$groupName == y), "groupID"]) + 1
                  fcs.files <- CytoML::fj_ws_get_samples(wsp, group_id = group_id)[, "name"]
                  ps <- do.call(rbind, lapply(fcs.files, function(z) {
                    gs <- flowjo_to_gatingset_CMS(ws = wsp, name = y, subset = z, path = FCS.file.folder, execute = F, emptyValue = F, which.lines = 1, additional.keys = c())
                    ps <- as.data.frame(flowWorkspace::gs_pop_get_count_fast(gs, path = "full", xml = T))
                    ps[, "Population"] <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")[-1]
                    return(ps)
                  }))
                  ps[, "group"] <- y
                  ps[, "ws"] <- basename(ws[x])
                  return(ps)
                }))
            } else {
                do.call(rbind, lapply(gr.wsp, function(y) {
                  gs <- CytoML::flowjo_to_gatingset(ws = wsp, name = y, path = FCS.file.folder, execute = F, emptyValue = F, which.lines = 1, additional.keys = c())
                  ps <- as.data.frame(flowWorkspace::gs_pop_get_count_fast(gs, path = "full", xml = T))
                  ps[, "Population"] <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")[-1]
                  ps[, "group"] <- y
                  ps[, "ws"] <- basename(ws[x])
                  return(ps)
                }))
            }
        }))
    }

    names(ps)[which(names(ps) == "name")] <- "FileName"
    ps[, "PopulationFullPath"] <- gsub("^root", "", paste(ps[, "Parent"], sapply(sapply(base::strsplit(ps[, "Population"], "/"), rev), "[", 1), sep = "/"), "^root")
    ps[, "FractionOfParent"] <- ps[, "Count"]/ps[, "ParentCount"] * 100
    ps <- ps[, c(1, 8, 2, 3, 4, 5, 9, 6, 7)]

    return(ps)
}

flowjo_to_gatingset_CMS <- function(ws, name = NULL, subset = list(), execute = TRUE, path = "", cytoset = NULL, backend_dir = tempdir(), backend = flowWorkspace::get_default_backend(),
    includeGates = TRUE, additional.keys = "$TOT", additional.sampleID = FALSE, keywords = character(), keywords.source = "XML", keyword.ignore.case = FALSE,
    extend_val = 0, extend_to = -4000, channel.ignore.case = FALSE, leaf.bool = TRUE, include_empty_tree = FALSE, skip_faulty_gate = FALSE, compensation = NULL,
    transform = TRUE, fcs_file_extension = ".fcs", greedy_match = FALSE, mc.cores = 1, ...) {
    if (is.null(cytoset))
        cytoset <- cytoset()
    backend <- match.arg(backend, c("h5", "tile"))
    g <- CytoML::fj_ws_get_sample_groups(ws)
    groups <- g[!duplicated(g$groupName), ]
    groups <- groups[order(groups$groupID), "groupName"]
    if (is.null(name)) {
        groupInd <- utils::menu(groups, graphics = FALSE, "Choose which group of samples to import:")
    } else if (is.numeric(name)) {
        if (length(groups) < name)
            stop("Invalid sample group index.")
        groupInd <- name
    } else if (is.character(name)) {
        if (is.na(match(name, groups)))
            stop("Invalid sample group name.")
        groupInd <- match(name, groups)
    }

    if (methods::is(subset, "character")) {
        subset <- list(name = subset)
    }
    if (!methods::is(subset, "list")) {
        stop("invalid 'subset' argument!")
    }
    if (is.null(additional.keys))
        additional.keys <- character(0)
    if (is.null(path))
        path <- ""
    args <- list(...)
    if (!is.null(args[["isNcdf"]])) {
        warning("'isNcdf' argument is deprecated!Data is always stored in h5 format by default!")
        args[["isNcdf"]] <- NULL
    }
    if (is.null(compensation)) {
        compensation <- list()
    } else {
        if (is.list(compensation) && !is.data.frame(compensation)) {
            compensation <- sapply(compensation, flowWorkspace:::check_comp, simplify = FALSE)
        } else compensation <- flowWorkspace:::check_comp(compensation)
    }
    args <- list(ws = ws@doc, group_id = groupInd - 1, subset = subset, execute = execute, path = suppressWarnings(normalizePath(path)), cytoset = cytoset@pointer,
        backend_dir = suppressWarnings(normalizePath(backend_dir)), backend = backend, includeGates = includeGates, additional_keys = additional.keys, additional_sampleID = additional.sampleID,
        keywords = keywords, is_pheno_data_from_FCS = keywords.source == "FCS", keyword_ignore_case = keyword.ignore.case, extend_val = extend_val, extend_to = extend_to,
        channel_ignore_case = channel.ignore.case, leaf_bool = leaf.bool, include_empty_tree = include_empty_tree, skip_faulty_gate = skip_faulty_gate, comps = compensation,
        transform = transform, fcs_file_extension = fcs_file_extension, greedy_match = greedy_match, fcs_parse_arg = args, num_threads = mc.cores)
    p <- do.call(CytoML:::parse_workspace, args)
    gs <- methods::new("GatingSet", pointer = p)
    gslist <- suppressMessages(flowWorkspace::gs_split_by_tree(gs))
    if (length(gslist) > 1) {
        msg <- "GatingSet contains different gating tree structures and must be cleaned before using it!\n "
        if (grepl("all samples", groups[groupInd], ignore.case = TRUE)) {
            msg <- c(msg, "It seems that you selected the 'All Samples' group,", " which is a generic group and typically contains samples with different gating schemes attached.",
                "Please choose a different sample group and try again.")
        }
        warning(msg)
    }
    gs
}

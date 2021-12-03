#' Write the compensation matrices generated in a FlowJo workspace to the SPILL keywords of respective FCS files
#'
#' The compensation matrices generated or optimized in FlowJo only exist in the workspace. The underlying FCS files remain unchanged.
#' This is part of FlowJos philosophy, FCS remain original all the time. This is actually not a bad procedure but reduces flexibility in some cases.
#' This function saves the associated compensation matrices to the FCS files - the matrices are written into the SPILL keyword. That makes
#' them easily and unambiguously available for other calculations outside of FlowJo (e.g. in R) when compensated channels are required.
#'
#' @param ws path to flowjo workspace
#' @param ... additional arguments to fcexpr:::prep_spill(): max_match_dist, skip_check, verbose
#'
#' @return no return but SPILL keyword updated in FCS files
#' @export
#'
#' @examples
wsx_compMats_to_fcs <- function(ws, ...) {

<<<<<<< HEAD
  if (!"BiocManager" %in% rownames(utils::installed.packages())) {utils::install.packages("BiocManager")}
  if (!"flowCore" %in% rownames(utils::installed.packages())) {BiocManager::install("flowCore")}

  ws <- fcexpr:::check_ws(ws)
=======
  if (!requireNamespace("CytoML", quietly = T)){
    BiocManager::install("CytoML")
  }
  if (!requireNamespace("flowWorkspace", quietly = T)){
    BiocManager::install("flowWorkspace")
  }
  ws <- check_ws(ws)
>>>>>>> dev

  ss <- xml2::xml_find_all(xml2::xml_child(ws, "SampleList"), "Sample")
  compMats <- lapply(seq_along(ss), function(n) {
    sp <- xml2::xml_child(ss[[n]], "transforms:spilloverMatrix")
    if (!is.na(sp)) {
      compMat <- t(do.call(cbind, lapply(xml2::xml_children(sp)[2:length(xml2::xml_children(sp))], function(x) {
        mat <- as.matrix(stats::setNames(as.numeric(xml2::xml_attr(xml2::xml_children(x), "value")), xml2::xml_attr(xml2::xml_children(x), "parameter")), nrow = 1, byrow = T)
        colnames(mat) <- xml2::xml_attr(x, "parameter")
        return(mat)
      })))
    } else {
      compMat <- NULL
    }
    return(compMat)
  })
  names(compMats) <- sapply(seq_along(ss), function(n) {
    gsub("^file:", "", xml2::xml_attr(xml2::xml_child(ss[[n]], "DataSet"), "uri"))
  })
  compMats <- compMats[which(!sapply(compMats, is.null))]

  for (i in seq_along(compMats)) {
    if (!file.exists(names(compMats)[i])) {
      print(paste0("File ", names(compMats)[i]), " not found. Did you change the names of FCS files and have not reconnected them to the workspace? Or did you copy the workspace from somewhere else? If so, open the workspace, reconnect the fcs files and save.")
    }
    ff <- flowCore::read.FCS(names(compMats)[i], truncate_max_range = F, emptyValue = F)
    sp <- flowCore::keyword(ff)[["SPILL"]]
    sp <- prep_spill(sp = sp, compMat = compMats[[i]], ...)
    flowCore::keyword(ff)[["SPILL"]] <- sp
    flowCore::write.FCS(ff, names(compMats)[i])
  }

}


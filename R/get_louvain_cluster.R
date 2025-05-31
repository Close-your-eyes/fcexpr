#' Calculate louvain clusters based on expression matrix
#'
#' @param exprs numeric matrix with channels as columns
#' @param FindNeighbors_args arguments to Seurat::FindNeighbors
#' @param FindClusters_args arguments to Seurat::FindClusters
#' @param mc.cores cores to use to calculate multiple resolutions in parallel
#'
#' @return matrix with cluster identities
#' @export
#'
#' @examples
get_louvain_cluster <- function(exprs,
                                FindNeighbors_args = list(),
                                FindClusters_args = list(resolution = c(0.1, 0.2, 0.3, 0.4)),
                                mc.cores = 1) {
  mc.cores <- min(mc.cores, parallel::detectCores()-1)
  rownames(exprs) <- 1:nrow(exprs)
  snn <- suppressMessages(suppressWarnings(do.call(Seurat::FindNeighbors, args = c(list(object = exprs), FindNeighbors_args))))

  if (any(grepl("resolution", names(FindClusters_args), ignore.case = T)) && length(FindClusters_args[["resolution"]]) > 1) {
    clust_idents <- do.call(cbind, parallel::mclapply(FindClusters_args[["resolution"]], function(x) {
      apply(do.call(Seurat::FindClusters, args = c(list(object = snn$snn, resolution = x, verbose = F),
                                                   FindClusters_args[which(names(FindClusters_args) != "resolution")])), 2, as.numeric)
    }, mc.cores = mc.cores))
  } else {
    clust_idents <- apply(do.call(Seurat::FindClusters, args = c(list(object = snn$snn), FindClusters_args)), 2, as.numeric)
  }
  return(clust_idents)
}

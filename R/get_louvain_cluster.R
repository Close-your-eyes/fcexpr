#' Calculate louvain clusters based on expression matrix
#'
#' @param exprs numeric matrix with channels as columns
#' @param FindNeighbors_args arguments to Seurat::FindNeighbors
#' @param FindClusters_args arguments to Seurat::FindClusters
#' @param mc.cores cores to use to calculate multiple resolutions in parallel
#' @param return_as
#' @param pad_len
#' @param unique_only
#' @param one_res_per_cluster_num_only
#' @param min_cluster
#'
#' @return matrix with cluster identities
#' @export
#'
#' @examples
#' \dontrun{
#' # try options(future.globals.maxSize = 8000 * 1024^2)
#' }
get_louvain_cluster <- function(exprs,
                                FindNeighbors_args = list(),
                                FindClusters_args = list(resolution = c(0.1), verbose = T),
                                mc.cores = 1,
                                return_as = c("character", "numeric"),
                                pad_len = "overall_max",
                                unique_only = T,
                                one_res_per_cluster_num_only = F,
                                min_cluster = 1) {

  return_as <- rlang::arg_match(return_as)

  mc.cores <- min(mc.cores, parallel::detectCores()-1)
  if (is.null(rownames(exprs))) {
    rownames(exprs) <- 1:nrow(exprs)
  }

  snn <- suppressMessages(suppressWarnings(do.call(Seurat::FindNeighbors, args = c(list(object = exprs), FindNeighbors_args))))

  if (any(grepl("resolution", names(FindClusters_args), ignore.case = T)) && length(FindClusters_args[["resolution"]]) > 1) {
    clust_idents <- do.call(cbind, parallel::mclapply(FindClusters_args[["resolution"]], function(x) {
      apply(do.call(Seurat::FindClusters, args = c(list(object = snn$snn, resolution = x, verbose = F),
                                                   FindClusters_args[which(names(FindClusters_args) != "resolution")])), 2, as.numeric)
    }, mc.cores = mc.cores))
  } else {
    clust_idents <- apply(do.call(Seurat::FindClusters, args = c(list(object = snn$snn), FindClusters_args)), 2, as.numeric)
  }

  if (return_as == "character") {
    if (pad_len == "overall_max") {
      pad_len <- nchar(max(clust_idents))
    }
    clust_idents <- apply(clust_idents, 2, as.character)
    clust_idents <- apply(clust_idents, 2, function(x) brathering::pad_num_in_str(x = x, len = pad_len))
  }

  if (unique_only) {
    hashes <- apply(clust_idents, 2, digest::digest)
    cl_unique <- hashes[!duplicated(hashes)]
    clust_idents <- clust_idents[,names(cl_unique), drop = F]
  }

  if (one_res_per_cluster_num_only) {
    nn <- apply(clust_idents, 2, function(x) length(unique(x)))
    clust_idents <- clust_idents[,names(nn[!duplicated(nn)])]
  }

  clust_idents <- clust_idents[,which(apply(clust_idents, 2, function(x) length(unique(x))) >= min_cluster)]
  if (ncol(clust_idents) == 0) {
    message("no clusterings left. decrease min_cluster?")
  } else {
    rownames(clust_idents) <- rownames(exprs)
  }

  return(clust_idents)
}



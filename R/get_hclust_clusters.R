#' Calculate clusters with hclust
#'
#' stats::dist takes quite some time, reduction of columns with pca may help
#'
#' @param exprs numeric matrix with channels as columns
#' @param distance distance metric for stats::dist
#' @param method method for stats::hclust
#' @param k number of clusters to generate by stats::cutree
#'
#' @return matrix with cluster annotations as columns
#' @export
#'
#' @examples
get_hclust_clusters <- function(exprs,
                                distance = c("euclidean", "correlation", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                                method = c("complete", "ward.D", "ward.D2", "ward", "single", "average", "mcquitty", "median", "centroid"),
                                k = c(5,10,20),
                                mc.cores = 1) {


  mc.cores <- min(mc.cores, parallel::detectCores()-1)
  method <- rlang::arg_match(method)
  distance <- rlang::arg_match(distance)

  if (distance == "correlation") {
    d = as.dist(1 - cor(t(exprs)))
  } else{
    d <- do.call(stats::dist, args = c(list(x = exprs, method = distance)))
  }

  h <- do.call(stats::hclust, args = c(list(d = d, method = method)))
  ks <- do.call(cbind, parallel::mclapply(k, function(x) stats::cutree(tree = h, k = x), mc.cores = mc.cores))


  # make sure that cluster 1 is the largest and so on
  ks <- cluster_ordering(ks = ks)
  colnames(ks) <- paste0("cutree_", k)

  return(ks)
}

cluster_ordering <- function(ks) {
  if (methods::is(ks, "matrix")) {
    ks <- apply(ks, 2, function (x) {
      new_order <- stats::setNames(names(table(x)), nm = names(sort(table(x), decreasing = T)))
      return(as.numeric(new_order[as.character(x)]))
    })
  } else if (methods::is(ks, "integer")) {
    new_order <- stats::setNames(names(table(ks)), nm = names(sort(table(ks), decreasing = T)))
    ks <- as.numeric(new_order[as.character(ks)])
  } else {
    new_order <- stats::setNames(names(table(ks)), nm = names(sort(table(ks), decreasing = T)))
    ks <- unname(new_order[as.character(ks)])
  }
  return(ks)
}

# nothing was necessarily faster than dist
# mat <- expr[1:10000,]
# system.time(d1 <- dist(mat))
# system.time(d2 <- as.dist(sqrt(outer(rowSums(mat^2), rowSums(mat^2), `+`) - 2 * tcrossprod(mat))) # slower)
# system.time(d3 <- as.dist(Rfast::Dist(mat)))
# system.time(d4 <- as.dist(proxy::dist(mat)))
# system.time(d4 <- parallelDist::parDist(mat, threads = 4))

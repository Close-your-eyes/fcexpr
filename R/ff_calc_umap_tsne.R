#' Calculate UMAP or tSNE dimension reduction based on flow frame
#'
#' @param ff flow frame
#' @param channels which channels to use, see ff_get_channels for selection
#' options
#' @param exprs expression matrix, either provide ff or exprs
#' @param fun_args arguments to fun
#' @param seed random seed
#' @param names names of new dimensions to add as channels, when auto
#' adjusted to UMAP or tSNE with trailing dimension numbers
#' @param return return UMAP object or flow frame with umap channels added
#' @param ... arguments to ff_get_channels
#' @param fun uwot::umap or Rtsne::Rtsne or Seurat::RunTSNE or Seurat::RunUMAP
#' or brathering::fft_tsne
#'
#' @return flow frame or umap object
#' @export
#'
#' @examples
#'\dontrun{
#' mat <- hdos::random_space(n_spheres = 40, n_tori = 0, n_cuboids = 0,
#'                           n_samples = 2000, max_dim = 20)
#'
#' umap <- ff_calc_umap_tsne(exprs = as.matrix(mat[["coord"]]),
#'                           fun_args = list(n_neighbors = 20, metric = "cosine", verbose = T, scale = T),
#'                           fun = uwot::umap)
#' brathering::plot2(umap)
#'
#' # fast tsne: fft or fit
#' tsne1 <- ff_calc_umap_tsne(exprs = scale(as.matrix(mat[["coord"]])),
#'                            fun_args = list(),
#'                            fun = brathering::fft_tsne)
#' brathering::plot2(tsne1)
#'
#' # slow tsne: barnes hut
#' tsne2 <- ff_calc_umap_tsne(exprs = as.matrix(mat[["coord"]])[1:4000,],
#'                            fun_args = list(theta = 0.8, verbose = T, pca_scale = T),
#'                            fun = Rtsne::Rtsne)
#' brathering::plot2(tsne2)
#'
#' # seurat wrapper, barnes hut
#' tsne3 <- ff_calc_umap_tsne(exprs = as.matrix(mat[["coord"]])[1:4000,],
#'                            fun_args = list(theta = 0.8, verbose = T, pca_scale = T, tsne.method = "Rtsne"),
#'                            fun = Seurat::RunTSNE)
#' brathering::plot2(tsne3)
#'
#' # seurat wrapper, fit
#' tsne4 <- ff_calc_umap_tsne(exprs = as.matrix(mat[["coord"]]),
#'                            fun_args = list(theta = 0.01, tsne.method = "FIt-SNE"),
#'                            fun = Seurat::RunTSNE)
#' brathering::plot2(tsne4)
#'
#' ## try 3D
#' # tsne: not supported
#' tsne3d <- ff_calc_umap_tsne(exprs = scale(as.matrix(mat[["coord"]])),
#'                             fun_args = list(dims = 3),
#'                             fun = brathering::fft_tsne)
#' # umap
#' umap3d <- ff_calc_umap_tsne(exprs = as.matrix(mat[["coord"]]),
#'                             fun_args = list(n_neighbors = 20, metric = "cosine", verbose = T, scale = T, n_components = 3),
#'                             fun = uwot::umap)
#'
#' brathering::plot3d(umap3d)
#' brathering::plot3d(umap3d, backend = "rgl")
#' }
ff_calc_umap_tsne <- function(ff,
                              channels = NULL,
                              exprs = NULL,
                              fun_args = list(n_neighbors = 20, metric = "cosine", verbose = T, scale = T),
                              seed = 42,
                              names = "..auto..",
                              return = c("coords", "ff"),
                              fun = uwot::umap,
                              ...) {

  return <- rlang::arg_match(return)
  if (return == "ff" && missing(ff)) {
    return <- "coords"
  }

  if (is.null(exprs)) {
    channels <- ff_get_channels(ff, channels = channels, ...)
    stopifnot("umap/tSNE: less than 3 channels provided." = length(channels) > 2)
    exprs <- flowCore::exprs(ff)[,channels]
    message("Channels for umap/tSNE: ", paste(channels, collapse = ", "))
  } else {
    stopifnot("exprs must be matrix." = is.matrix(exprs))
  }



  set.seed(seed)
  if (identical(fun, uwot::umap) || identical(fun, brathering::fft_tsne)) {
    coords <- Gmisc::fastDoCall(fun, args = c(list(X = exprs), fun_args))
    if (identical(fun, brathering::fft_tsne)) {
      what <- "tSNE"
    } else {
      what <- "UMAP"
    }
  } else if (identical(fun, Rtsne::Rtsne)) {
    coords <- Gmisc::fastDoCall(fun, args = c(list(X = exprs), fun_args))[["Y"]]
    what <- "tSNE"
  } else if (identical(fun, Seurat::RunTSNE)) {
    if (is.null(rownames(exprs))) {
      rownames(exprs) <- as.character(1:nrow(exprs))
    }
    coords <- Gmisc::fastDoCall(fun, args = c(list(object = exprs, assay = "RNA"), fun_args))@cell.embeddings
    what <- "tSNE"
  } else if (identical(fun, Seurat::RunUMAP)) {
    if (is.null(rownames(exprs))) {
      rownames(exprs) <- as.character(1:nrow(exprs))
    }
    coords <- Gmisc::fastDoCall(fun, args = c(list(object = exprs, assay = "RNA"), fun_args))@cell.embeddings
    what <- "UMAP"
  } else {
    stop("fun not supported.")
  }

  ## prep and check names to avoid error below
  if (names == "..auto..") {
    names <- paste(what, 1:ncol(coords), sep = "_")
  } else {
    if (length(names) < ncol(coords)) {
      names <- paste("Dim", 1:ncol(coords), sep = "_")
    }
  }
  names <- make.unique(names)[1:ncol(coords)]
  colnames(coords) <- names

  if (return == "coords") {
    return(coords)
  }
  if (return == "ff") {
    return(ff_add_columns(ff, coords))
  }
}

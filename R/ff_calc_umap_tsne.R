#' Calculate UMAP or tSNE dimension reduction based on flow frame
#'
#' @param ff A flow frame. Required when exprs is NULL. If omitted, coordinates
#' are returned even when return is "ff".
#' @param channels Channels to use when extracting expression values from ff;
#' see ff_get_channels for selection options. Ignored when exprs is supplied.
#' @param exprs An expression matrix with observations in rows and channels in
#' columns. If NULL, values are extracted from ff.
#' @param fun_args A list of arguments passed to fun. The default is NULL,
#' which uses list(n_neighbors = 20, metric = "cosine", verbose = TRUE,
#' scale = TRUE) for uwot::umap and an empty list for all other backends,
#' leaving their defaults in effect. Supply list() to use uwot's own defaults.
#' An explicit list replaces the wrapper defaults rather than merging with them.
#' @param seed Random seed, default 42. Also passed as rand_seed to
#' brathering::fft_tsne or seed.use to the Seurat backends, unless that argument
#' is already supplied in fun_args. Explicit backend seed arguments take precedence.
#' @param names Names of the output dimensions. The default, "..auto..", uses
#' UMAP or tSNE followed by an underscore and the dimension number. If too few
#' names are supplied, Dim names with trailing dimension numbers are used.
#' Duplicate names are made unique and excess names are discarded.
#' @param return Either "coords" (the default) for a coordinate matrix or "ff"
#' for the input flow frame with the dimensions added as channels.
#' @param ... Additional arguments passed to ff_get_channels when exprs is NULL.
#' @param fun Dimension reduction function: uwot::umap (the default),
#' Rtsne::Rtsne, Seurat::RunTSNE, Seurat::RunUMAP, or brathering::fft_tsne.
#'
#' @return A coordinate matrix with one row per observation, or a flow frame
#' with the dimensions added as channels. For uwot list outputs, only the
#' embedding is returned; model and other auxiliary results are discarded.
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
                              fun_args = NULL,
                              seed = 42,
                              names = "..auto..",
                              return = c("coords", "ff"),
                              fun = uwot::umap,
                              ...) {

  fcexpr:::.ensure_packages("Gmisc")

  if (missing(fun)) {
    fcexpr:::.ensure_packages("uwot")
  }
  if (!is.function(fun) || !isNamespace(environment(fun))) {
    stop("fun not supported.")
  }
  backend <- getNamespaceName(environment(fun))
  supported <- switch(backend,
                      uwot = "umap",
                      brathering = "fft_tsne",
                      Rtsne = "Rtsne",
                      Seurat = c("RunTSNE", "RunUMAP"),
                      character())
  matches <- vapply(supported, function(name) {
    identical(fun, get(name, envir = environment(fun), inherits = FALSE))
  }, logical(1))
  if (!any(matches)) {
    stop("fun not supported.")
  }
  method <- supported[matches][1]
  fcexpr:::.ensure_packages(backend)

  if (is.null(fun_args)) {
    fun_args <- if (backend == "uwot") {
      list(n_neighbors = 20, metric = "cosine", verbose = TRUE, scale = TRUE)
    } else {
      list()
    }
  }

  return <- rlang::arg_match(return)
  if (return == "ff" && missing(ff)) {
    return <- "coords"
  }

  if (is.null(exprs) || return == "ff") {
    fcexpr:::.ensure_packages("flowCore")
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
  if (backend == "brathering" && !"rand_seed" %in% names(fun_args)) {
    fun_args$rand_seed <- seed
  }
  if (backend == "Seurat" && !"seed.use" %in% names(fun_args)) {
    fun_args$seed.use <- seed
  }

  if (backend == "Seurat") {
    if (is.null(rownames(exprs))) {
      rownames(exprs) <- as.character(seq_len(nrow(exprs)))
    }
    coords <- Gmisc::fastDoCall(fun, c(list(object = exprs, assay = "RNA"), fun_args))@cell.embeddings
  } else {
    coords <- Gmisc::fastDoCall(fun, c(list(X = exprs), fun_args))
    if (backend == "Rtsne") {
      coords <- coords[["Y"]]
    } else if (backend == "uwot" && is.list(coords)) {
      coords <- coords[["embedding"]]
    }
  }
  stopifnot("Backend must return a coordinate matrix." = is.matrix(coords))
  what <- if (method %in% c("umap", "RunUMAP")) "UMAP" else "tSNE"

  ## prep and check names to avoid error below
  if (names[1] == "..auto..") {
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

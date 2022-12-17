#' Gates from gatingset for plotting with ggcyto
#'
#' @param gs gs a gatingset, e.g. made with fcexpr::wsp_get_gs
#' @param n_bins number of bins in total, will be used equally in x and y direction, bin size is adjusted to ranges in x and y direction
#' @param quantile_lim_filter min,max quantiles of signal intensities to set axis limits to; use quantiles to exclude extreme values
#' @param min_max_vals min,max required signal intensitiy for fluorescence channels of one event in order to condider it for axis limit calculation (to filter extreme values);
#' in logicle transformation
#' @param min_max_vals_scatter min,max required signal intensitiy for scatter channels of one event in order to condider it for axis limit calculation (to filter extreme values);
#' in inverse transformation which is equal to logicle transformation (for scatter channels)
#' @param x_statpos_name x-position for gate name labels
#' @param y_statpos_name y-position for gate name labels
#' @param x_statpos_pct x-position for gate percent labels
#' @param y_statpos_pct y-position for gate percent labels
#' @param statsize_name size of name label
#' @param statsize_pct size of percent label
#'
#' @return a data frame to loop over and produce plots with ggcyto
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' gates <- gs_get_gates(gs)
#' }
gs_get_gates <- function(gs,
                         n_bins = 50^2,
                         quantile_lim_filter = c(0.0001, 0.9999),
                         min_max_vals = c(0, 300),
                         min_max_vals_scatter = c(0, 250000),
                         x_statpos_name = 0.9,
                         y_statpos_name = 0.9,
                         x_statpos_pct = 0.9,
                         y_statpos_pct = 0.1,
                         statsize_name = 4,
                         statsize_pct = 4) {

  if (!requireNamespace("flowWorkspace", quietly = T)){
    utils::install.packages("flowWorkspace")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }

  if (!is.null(min_max_vals_scatter)) {
    if (!is.numeric(min_max_vals_scatter) || length(min_max_vals_scatter) != 2) {
      stop("min_max_vals_scatter has to be a numeric vector of length 2.")
    }
  }

  if (!is.null(min_max_vals)) {
    if (!is.numeric(min_max_vals) || length(min_max_vals) != 2) {
      stop("min_max_vals has to be a numeric vector of length 2.")
    }
  }

  if (is.null(quantile_lim_filter)) {
    quantile_lim_filter <- c(0,1)
  }

  gates <-
    data.frame(PopulationFullPath = gsub("^/", "", flowWorkspace::gs_get_pop_paths(gs)),
               Population = flowWorkspace::gs_get_pop_paths(gs, path = "auto"),
               GateLevel = nchar(flowWorkspace::gs_get_pop_paths(gs)) - nchar(gsub("/", "", flowWorkspace::gs_get_pop_paths(gs)))) %>%
    dplyr::filter(GateLevel > 0) %>%
    dplyr::mutate(Parent = gsub("^/", "", dirname(PopulationFullPath))) %>%
    dplyr::mutate(Parent = ifelse(Parent == ".", "root", Parent)) %>%
    dplyr::mutate(x_statpos_name = x_statpos_name,
                  y_statpos_name = y_statpos_name,
                  statsize_name = statsize_name,
                  x_statpos_pct = x_statpos_pct,
                  y_statpos_pct = y_statpos_pct,
                  statsize_pct = statsize_pct)

  gates$dims <- sapply(gates$PopulationFullPath, function(x) {
    y <- unname(flowCore::parameters({flowWorkspace::gs_pop_get_gate(gs[[1]], x)[[1]]}))
    return(y)
    # stupid handling of Not-gate. other booleans may require similar specific treatment
'    if (length(y) == 0) {
      y <- list(unname(flowCore::parameters({flowWorkspace::gs_pop_get_gate(gs[[1]], gsub("^!", "", zz@deparse))[[1]]})))
    } else {
      return(list(y))
    }'
  }, simplify = F)
  ## filter boolean gates - test further ... # boolean are not easy to handle (e.g. their children)
  gates <- gates[which(lengths(gates$dims) > 0),]

  gates$x <- unname(sapply(gates$dims, function(x) {unlist(x)[1]}))
  gates$y <- unname(sapply(gates$dims, function(x) {unlist(x)[2]}))
  gates$x_lab <- unname(sapply(gates$dims, function(x) {unlist(x)[1]}))
  gates$y_lab <- unname(sapply(gates$dims, function(x) {unlist(x)[2]}))
  gates$marginalFilter <- ifelse(grepl("fsc|ssc", gates$x, ignore.case = T) & grepl("fsc|ssc", gates$y, ignore.case = T), T, F)

  lims <- lapply(split(gates, 1:nrow(gates)), function(y) {
    parent <- dirname(y$PopulationFullPath)
    if (parent == ".") {
      parent <- "root"
    }

    out <- flowWorkspace::cytoset_to_flowSet(flowWorkspace::gs_pop_get_data(gs, y = parent, truncate_max_range = F))
    out_names <- names(out@frames)
    min_max_vals <- sort(min_max_vals)
    min_max_vals_scatter <- sort(min_max_vals_scatter)
    tempfun <- function(x, z) {
      if (grepl("fsc|ssc", x, ignore.case = T)) {
        flowCore::exprs(out[[z]])[,x] > min_max_vals_scatter[1] & flowCore::exprs(out[[z]])[,x] < min_max_vals_scatter[2]
      } else {
        flowCore::exprs(out[[z]])[,x] > min_max_vals[1] & flowCore::exprs(out[[z]])[,x] < min_max_vals[2]
      }
    }

    ## currently focus is on 2D-gates only
    quants <- do.call(rbind, lapply(out_names, function(z) {
      # rel are rows for which all values above or below min_max_vals; not 100 % correct as outliers in one column are also removed for all columns
      temp <- sapply(c(y$x, y$y), tempfun, z = z)
      if (all(lengths(temp) > 0)) {
        temp <- as.matrix(temp)
        if (ncol(temp) == 1) {
          temp <- t(temp)
        }
        rel <- apply(temp, 1, all)
        if (length(rel) > 0) {
          apply(flowCore::exprs(out[[z]])[rel,c(y$x, y$y),drop=F], 2, stats::quantile, quantile_lim_filter)
        } else {
          NULL
        }
      } else {
        NULL
      }
    }))
    # get min and max from all flowFrames
    return(c(apply(quants, 2, min), apply(quants, 2, max)))
  })

  # order is known
  gates$x_lowlim <- sapply(lims, "[", 1)
  gates$x_uplim <- sapply(lims, "[", 3)
  gates$y_lowlim <- sapply(lims, "[", 2)
  gates$y_uplim <- sapply(lims, "[", 4)

  mat <- cbind((gates$x_uplim - gates$x_lowlim)/sqrt(n_bins), (gates$y_uplim - gates$y_lowlim)/sqrt(n_bins))
  gates$binwidths <- split(t(mat), rep(1:nrow(mat), each = ncol(mat)))

  gates$facet_strip <- c(T, rep(F, nrow(gates)-1))

  return(gates)
}


# https://stackoverflow.com/questions/24299171/function-to-split-a-matrix-into-sub-matrices-in-r
matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
}


#' Gates from gatingset for plotting with ggcyto
#'
#' @param gs gs a gatingset, e.g. made with fcexpr::wsp_get_gs
#' @param n_bins number of bins in total, will be used equally in x and y direction, bin size is adjusted to ranges in x and y direction
#' @param quantile_lim_filter min,max quantiles of signal intensities to set axis limits to; use quantiles to exclude extreme values
#' @param min_max_vals min,max required signal intensity for fluorescence channels of one event in order to consider it for axis limit calculation (to filter extreme values);
#' in logicle transformation
#' @param min_max_vals_scatter min,max required signal intensity for scatter channels of one event in order to consider it for axis limit calculation (to filter extreme values);
#' in inverse transformation which is equal to logicle transformation (for scatter channels)
#' @param x_statpos_name x-position for gate name labels
#' @param y_statpos_name y-position for gate name labels
#' @param x_statpos_pct x-position for gate percent labels
#' @param y_statpos_pct y-position for gate percent labels
#' @param statsize_name size of name label
#' @param statsize_pct size of percent label
#' @param reduce_n_bins_below
#' @param reduce_factor
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
                         n_bins = 200^2,
                         quantile_lim_filter = c(0.0001, 0.9999),
                         min_max_vals = c(0, 300),
                         min_max_vals_scatter = c(0, 250000),
                         x_statpos_name = 0.9,
                         y_statpos_name = 0.9,
                         x_statpos_pct = 0.9,
                         y_statpos_pct = 0.1,
                         statsize_name = 4,
                         statsize_pct = 4,
                         reduce_n_bins_below = 0,
                         reduce_factor = 2) {

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
  min_max_vals <- sort(min_max_vals)
  min_max_vals_scatter <- sort(min_max_vals_scatter)

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

  # for 1D gates: dimension becomes x-dim, always. does not matter. or does it?
  gates$dims <- lapply(gates$PopulationFullPath, function(x) {
    y <- unname(flowCore::parameters(flowWorkspace::gs_pop_get_gate(gs[[1]], x)[[1]]))
    return(y)
    # stupid handling of Not-gate. other booleans may require similar specific treatment
    '    if (length(y) == 0) {
      y <- list(unname(flowCore::parameters({flowWorkspace::gs_pop_get_gate(gs[[1]], gsub("^!", "", zz@deparse))[[1]]})))
    } else {
      return(list(y))
    }'
  })

  ## filter boolean gates - test further ... # boolean are not easy to handle (e.g. their children)
  gates <- gates[which(lengths(gates$dims) > 0),]

  gates$x <- sapply(gates$dims, "[", 1)
  gates$y <- sapply(gates$dims, "[", 2)
  gates$x_lab <- gates$x
  gates$y_lab <- gates$y
  gates$marginalFilter <- ifelse(grepl("fsc|ssc", gates$x, ignore.case = T) & grepl("fsc|ssc", gates$y, ignore.case = T), T, F)

  lims_count <- purrr::map(split(gates, 1:nrow(gates)), function(gate) {

    fs <- flowWorkspace::cytoset_to_flowSet(flowWorkspace::gs_pop_get_data(gs, y = gate$Parent, truncate_max_range = F))

    ## currently focus is on 2D-gates only and 1D
    quants_count <- purrr::map(names(fs@frames),
                               ~ get_quantiles_and_count(gate = gate,
                                                         fs = fs,
                                                         fs_name = .x,
                                                         min_max_vals_scatter = min_max_vals_scatter,
                                                         min_max_vals = min_max_vals,
                                                         quantile_lim_filter = quantile_lim_filter))

    quants <- do.call(rbind, purrr::discard(sapply(quants_count, "[[", "quants", simplify = F), is.null))
    counts <- purrr::discard(sapply(quants_count, "[[", "count", simplify = F), is.null)

    # get min and max from all flowFrames
    if (is.null(quants)) {
      return(NULL)
    }
    quants <- c(apply(quants, 2, min), apply(quants, 2, max))
    if (is.na(gate$x)) {
      quants <- c(NA, quants[1], NA, quants[2])
    }
    if (is.na(gate$y)) {
      quants <- c(quants[1], NA, quants[2], NA)
    }
    return(list(quants = quants, counts = counts))
  })

  lims <- sapply(lims_count, "[[", "quants", simplify = F)
  lims[which(sapply(lims, is.null))] <- NA

  counts <- sapply(lims_count, "[[", "counts", simplify = F)
  counts <- sapply(counts, unlist)

  # order is known
  gates$x_lowlim <- sapply(lims, "[", 1)
  gates$x_uplim <- sapply(lims, "[", 3)
  gates$y_lowlim <- sapply(lims, "[", 2)
  gates$y_uplim <- sapply(lims, "[", 4)

  n_bins <- rep(n_bins, nrow(gates))
  below <- sapply(counts, function(x) any(x<reduce_n_bins_below))
  n_bins[which(below)] <- n_bins[which(below)] / reduce_factor

  mat <- cbind((gates$x_uplim - gates$x_lowlim)/sqrt(n_bins), (gates$y_uplim - gates$y_lowlim)/sqrt(n_bins))
  gates$binwidths <- split(t(mat), rep(1:nrow(mat), each = ncol(mat)))

  gates$facet_strip <- c(T, rep(F, nrow(gates)-1))

  tryCatch(
    expr = {
      descnames <- gs_get_descname_lookup(gs)
      gates$x_lab <- descnames[gates$x]
      gates$y_lab <- descnames[gates$y]
    }, error = function(e) {
      message("gs_get_descname_lookup failed.")
      print(e)
    }
  )
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

get_inds <- function(channel, fs, fs_name, min_max = c(0, 300)) {
  if (is.na(channel)) {return(NULL)}
  if (nrow(flowCore::exprs(fs[[fs_name]])) == 0) {return(NULL)}
  inds_in_range <- dplyr::between(flowCore::exprs(fs[[fs_name]])[,channel], min_max[1], min_max[2])
  return(inds_in_range)
}
get_quantiles_and_count <- function(gate,
                                    fs,
                                    fs_name,
                                    min_max_vals_scatter,
                                    min_max_vals,
                                    quantile_lim_filter) {

  # inds are rows for which all values above or below min_max_vals; not 100 % correct as outliers in one column are also removed for all columns
  xy_channel <- stats::setNames(c(gate$x, gate$y), c("x", "y"))
  xy_channel <- xy_channel[which(!is.na(xy_channel))]
  #xy_channel <- stats::setNames(c(gate$x, gate$y), c(gate$x, gate$y))
  inds_xy <- purrr::map(xy_channel,
                        ~get_inds(channel = .x,
                                  fs = fs,
                                  fs_name = fs_name,
                                  min_max = if (grepl("fsc|ssc", .x, ignore.case = T)) min_max_vals_scatter else min_max_vals))
  if (length(inds_xy) == 2) {
    # 2D gate: no NA channel
    inds <- inds_xy[[1]] & inds_xy[[2]] # event within limits in both dimension?
  } else {
    inds <- inds_xy[[1]]
  }

  if (any(inds)) {
    quants <- apply(flowCore::exprs(fs[[fs_name]])[inds, xy_channel, drop=F],
                    MARGIN = 2,
                    FUN = stats::quantile,
                    probs = quantile_lim_filter)
    count <- sum(inds)
  } else {
    return(NULL)
  }



  return(list(quants = quants, count = count))

}

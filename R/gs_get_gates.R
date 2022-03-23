#' Gates from gatingset for plotting with ggcyto
#'
#' @param gs gatingset
#' @param n_bins number of bins in total, will be used equally in x and y direction, bin size is adjusted to range in x and y direction
#' @param quantile_lim_filter quantiles of signals to set axis limits to
#' @param min_max_vals minimum and/or maximum required signal of one event in order to condider it for axis limit calculation (to filter extreme values)
#' @param scatter_lim manual limits for scatter channel, set to NULL to get the actual limits (min and max)
#' @param x_statpos x-position of percent labels for gates
#' @param y_statpos y-position of percent labels for gates
#' @param stat_size size of percent labels for gates
#'
#' @return a data frame to loop over and produce plots with ggcyto
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' gates <- gs_get_gates(gs)

#' out <- purrr::flatten(lapply(unique(gates$gate.level), function (z) {
#'   g <- gates[which(gates[,"gate.level"] == z),]
#'   p <- lapply(split(g, paste(g$gate.level, g$subset, g$x, g$y, sep = "__")), function(gg) {
#'     my.filter <- if (gg[1,"marginalFilter"]) {ggcyto::marginalFilter} else {NULL}
#'     p <- ggcyto::ggcyto(gs, subset = gg[1,"subset"], filter = my.filter, aes(!!sym(gg[1,"x"]), !!sym(gg[1,"y"])), max_nrow_to_plot = 5e4) +
#'       geom_hex(binwidth = gg[1,"binwidths"][[1]]) +
#'       xlab(gg[1,"x_lab"]) +
#'       ylab(gg[1,"y_lab"]) +
#'       ggcyto::ggcyto_par_set(limits = list(x = c(gg[1,"x_lowlim"], gg[1,"x_uplim"]), y = c(gg[1,"y_lowlim"], gg[1,"y_uplim"]))) +
#'      scale_fill_gradientn(colours = scexpr::col_pal("spectral"), trans = "pseudo_log") +
#'       theme(legend.position = "none", axis.text = element_blank()) +
#'       facet_wrap(vars(FileName), nrow = 1, labeller = label_wrap_gen(multi_line=FALSE)) # order inline: facet_grid(cols = vars(PBMC.donor.short, factor(IL15.pre.stim.conc.ug.ml, levels=c("0", "0.12", "0.37", "1.11", "3.33", "10"))), rows = vars(RPTECs, RPTEC.IFNg.pre.stim))
#'
#'     if (all(!gg$facet.strip)) {
#'       p <- p + theme(strip.background = element_blank(), strip.text = element_blank())
#'     }
#'
#'     ## loop through rows of gg for gates
#'     for (i in 1:nrow(gg)) {
#'       p <-
#'         p +
#'         ggcyto::geom_gate(gg[i,"gate.path.full"], colour = "black") +
#'         ggcyto::geom_stats(gg[i,"gate.path.full"], type = "gate_name", size = gg[i,"stat_size"], color = "black", adjust = c(gg[i,"x_statpos"], gg[i,"y_statpos"]), fill = alpha(c("white"),0.5))
#'     }
#'     return(p)
#'   })
#'   return(p)
#' }))
#' out.grid <- cowplot::plot_grid(plotlist = lapply(out, function(x) ggcyto::as.ggplot(x)), ncol = 1, align = "hv", axis = "tblr")
#' ggsave(out.grid, filename = paste0("plot.png"), device = "png", path = im_path, dpi = "retina", width = 7, height = 24)
#' }
#'
gs_get_gates <- function(gs,
                         n_bins = 50^2,
                         quantile_lim_filter = c(0.0001, 0.9999),
                         min_max_vals = c(0, 300),
                         scatter_lim = c(0, 250000),
                         x_statpos = 0.8,
                         y_statpos = 0.2,
                         stat_size = 4) {

  if (!requireNamespace("flowWorkspace", quietly = T)){
    utils::install.packages("flowWorkspace")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }

  if (!is.null(scatter_lim)) {
    if (!is.numeric(scatter_lim) || length(scatter_lim) != 2) {
      stop("scatter_lim has to be a numeric vector of length 2.")
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
    data.frame(gate.path.full = flowWorkspace::gs_get_pop_paths(gs),
               gate.path.auto = flowWorkspace::gs_get_pop_paths(gs, path = "auto"),
               gate.level = nchar(flowWorkspace::gs_get_pop_paths(gs)) - nchar(gsub("/", "", flowWorkspace::gs_get_pop_paths(gs)))) %>%
    dplyr::filter(gate.level > 0) %>%
    dplyr::mutate(subset = gsub("^/", "", dirname(gate.path.full))) %>%
    dplyr::mutate(subset = ifelse(subset == "", "root", subset)) %>%
    dplyr::mutate(x_statpos = x_statpos, y_statpos = y_statpos, stat_size = stat_size)

  gates$dims <- sapply(gates$gate.path.full, function(x) {
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
    parent <- dirname(y$gate.path.full)
    if (parent == "/") {
      parent <- "root"
    }

    out <- flowWorkspace::cytoset_to_flowSet(flowWorkspace::gs_pop_get_data(gs, y = parent, truncate_max_range = F))
    out_names <- names(out@frames)
    min_max_vals <- sort(min_max_vals)
    tempfun <- function(x) x > min_max_vals[1] & x < min_max_vals[2]

    ## currently focus is on 2D-gates only
    quants <- do.call(rbind, lapply(out_names, function(z) {
      # rel are rows for which all values above or below min_max_vals; not 100 % correct as outliers in one column are also removed for all columns
      temp <- apply(flowCore::exprs(out[[z]])[,c(y$x, y$y),drop=F], 2, tempfun)
      if (length(temp) > 0) {
        temp <- as.matrix(temp)
        if (ncol(temp) == 1) {
          temp <- t(temp)
        }
        rel <- apply(temp, 1, all)
        if (length(rel) > 0) {
          apply(flowCore::exprs(out[[z]])[rel,c(y$x, y$y),drop=F], 2, quantile, quantile_lim_filter)
        } else {
          NULL
        }
      } else {
        NULL
      }
    }))

    # get min and max from all flowFrames
    min_max_quants <- c(apply(quants, 2, min), apply(quants, 2, max))
    # order is known
    min_max_quants[which(grepl("fsc", names(min_max_quants), ignore.case = T))] <- scatter_lim
    min_max_quants[which(grepl("ssc", names(min_max_quants), ignore.case = T))] <- scatter_lim
    return(min_max_quants)
  })

  # order is known
  gates$x_lowlim <- sapply(lims, "[", 1)
  gates$x_uplim <- sapply(lims, "[", 3)
  gates$y_lowlim <- sapply(lims, "[", 2)
  gates$y_uplim <- sapply(lims, "[", 4)
  mat <- cbind((gates$x_uplim - gates$x_lowlim)/sqrt(n_bins), (gates$y_uplim - gates$y_lowlim)/sqrt(n_bins))
  gates$binwidths <- split(t(mat), rep(1:nrow(mat), each = ncol(mat)))

  if (!is.null(scatter_lim)) {
    scatter_lim <- sort(scatter_lim)
    gates[which(grepl("fsc|ssc", gates$x, ignore.case = T)),"x_lowlim"] <- scatter_lim[1]
    gates[which(grepl("fsc|ssc", gates$x, ignore.case = T)),"x_uplim"] <- scatter_lim[2]
    gates[which(grepl("fsc|ssc", gates$y, ignore.case = T)),"y_lowlim"] <- scatter_lim[1]
    gates[which(grepl("fsc|ssc", gates$y, ignore.case = T)),"y_uplim"] <- scatter_lim[2]
  }

  gates$facet.strip <- c(T, rep(F, nrow(gates)-1))

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


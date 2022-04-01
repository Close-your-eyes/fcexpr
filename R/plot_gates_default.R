#' Plot a default layout of gated populations of a set of FCS files (e.g. from one experiment)
#'
#'
#'
#' @param gs a gatingset, e.g. made with fcexpr::wsp_get_gs
#' @param gates_df a data frame with informaion of how to plot gates, made with fcexpr::gs_get_gates
#' @param facetting which facetting (ggplot2) to apply, facet_wrap or facet_grid with respective arguments, check flowCore::pData(gs) for valid columns
#' @param plot_gates logical whether to plot gates
#' @param plot_gate_names logical whether to plot gate names
#' @param plot_gate_pct logical whether to plot gate percentages (fraction of parent)
#' @param ... arguments passed to ggplot2::theme; set a global default theme with ggplot2:theme_set() and ggplot2::theme_update()
#'
#' @return
#' @export
#'
#' @examples
plot_gates <- function(gs,
                       gates_df,
                       facetting = NULL, # complete ggplot2::facet_wrap or ggplot2::facet_grid
                       plot_gates = T,
                       plot_gate_names = T,
                       plot_gate_pct = T,
                       ...) {
  # order inline: facet_grid(cols = vars(PBMC.donor.short, factor(IL15.pre.stim.conc.ug.ml, levels=c("0", "0.12", "0.37", "1.11", "3.33", "10"))), rows = vars(RPTECs, RPTEC.IFNg.pre.stim))

  dots <- list(...)

  out <- purrr::flatten(lapply(unique(gates_df$GateLevel), function (z) {
    g <- gates_df[which(gates_df[,"GateLevel"] == z),]

    p <- lapply(split(g, paste(g$GateLevel, g$Parent, g$x, g$y, sep = "__")), function(gg) {

      my.filter <- if (gg[1,"marginalFilter"]) {ggcyto::marginalFilter} else {NULL}

      p <- ggcyto::ggcyto(gs, subset = gg[1,"Parent"], filter = my.filter, aes(!!sym(gg[1,"x"]), !!sym(gg[1,"y"])), max_nrow_to_plot = 5e4) +
        ggplot2::geom_hex(binwidth = gg[1,"binwidths"][[1]]) +
        ggplot2::xlab(gg[1,"x_lab"]) +
        ggplot2::ylab(gg[1,"y_lab"]) +
        ggcyto::ggcyto_par_set(limits = list(x = c(gg[1,"x_lowlim"], gg[1,"x_uplim"]), y = c(gg[1,"y_lowlim"], gg[1,"y_uplim"]))) +
        ggplot2::scale_fill_gradientn(colours = rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"), interpolate = "linear")(100)), trans = "pseudo_log")


      if (length(dots) > 0) {
        p <- p + do.call(ggplot2::theme, args = dots[which(names(dots) %in% names(formals(ggplot2::theme)))])
      }

      if (!is.null(facetting)) {
        p <- p + facetting
      }

      if (all(!gg$facet_strip)) {
        p <- p + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_blank())
      }

      if (plot_gates) {
        for (i in 1:nrow(gg)) {
          p <- p + ggcyto::geom_gate(gg[i,"PopulationFullPath"], colour = "black")
        }
      }
      if (plot_gate_names) {
        for (i in 1:nrow(gg)) {
          p <- p + ggcyto::geom_stats(gg[i,"PopulationFullPath"], type = "gate_name", size = gg[i,"statsize_name"], color = "black", adjust = c(gg[i,"x_statpos_name"], gg[i,"y_statpos_name"]), fill = alpha(c("white"),0.5))
        }
      }
      if (plot_gate_pct) {
        for (i in 1:nrow(gg)) {
          p <- p + ggcyto::geom_stats(gg[i,"PopulationFullPath"], type = "percent", size = gg[i,"statsize_pct"], color = "black", adjust = c(gg[i,"x_statpos_pct"], gg[i,"y_statpos_pct"]), fill = alpha(c("white"),0.5))
        }
      }
      return(p)
    })
    return(p)
  }))


  return(lapply(out, function(x) ggcyto::as.ggplot(x)))
}

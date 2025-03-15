#' Plot a default layout of gated populations of a set of FCS files (e.g. from one experiment)
#'
#' Getting an overview of the results of a flow cytometry experiment can be a generic task if gating and analysis strategy is fixed.
#' Here, ggcyto is used to create an arrangement of pseudocolor plots. A number of default settings has been selected which may garuantee
#' good looking plots for most cases. Upon problems, the function may be improved to handle edge cases. Currently, customization is limited to a few rather technical
#' aspects of the plots. This may be subject to expansion but generally the range of customization is too large to fit it into a generic function with a finite number of arguments.
#' Hence, for very specific and detailed requriements manual plotting and fiddling is unavoidable. For that, the code of this function may serve as
#' a template.
#'
#' @param gs a gatingset, e.g. made with fcexpr::wsp_get_gs
#' @param gates_df a data frame with informaion of how to plot gates, made with fcexpr::gs_get_gates
#' @param facetting which facetting (ggplot2) to apply, facet_wrap or facet_grid with respective arguments, check flowCore::pData(gs) for valid columns;
#' facet_null will completely remove facets; by default facetting is done across fcs each individual fcs file
#' @param plot_gates logical whether to plot gates
#' @param plot_gate_names logical whether to plot gate names
#' @param plot_gate_pct logical whether to plot gate percentages (fraction of parent)
#' @param inverse_trans logical wheter to inverse transform axes numbers; if TRUE this will make axes look like in flowjo
#' @param geom how to plot data; recommendation is hex; hex = geom_hex taking the binwidths column of gates_df into account, pointdensity = ggpointdensity::geom_pointdensity (see ...
#' and optionally provide max_nrow_to_plot (passed to ggcyto::ggcyto) to limit the number of dots per plot, defaults to 2000; be careful may increase time for plotting a lot),
#' scattermore = scattermore::geom_scattermore which is a high performance dot plot version for quickly plotting millions of points (only black-white currently)
#' @param gate_stats_color font color of gate statistics
#' @param gate_color line color of gates
#' @param plot_contours logical whether to plot contour lines on top with ggplot2::geom_density_2d
#' @param pct_digits how many digits after comma to print; passed to 'digits' of ggcyto::geom_stats
#' @param col_pal color palette to use for color gradient generation
#' @param col_pal_trans argument passed as 'trans' in scale_fill_gradientn
#'
#' @return a list of ggplot2 objects, one for every gating level; each list index contains respective plots for every fcs file
#' @export
#'
#' @examples
#' \dontrun{
#' ## read gatingset
#' gs <- fcexpr::wsp_get_gs(wsp = ws, groups = "Group1")
#'
#' ## write meta data to pData of gs; sd is sampledescription
#' p.df <-
#' flowCore::pData(gs) %>%
#' tibble::rownames_to_column("FileName") %>%
#' dplyr::left_join(sd) %>%
#' tibble::column_to_rownames("FileName")
#' p.df$FileName <- rownames(p.df)
#' flowCore::pData(gs) <- p.df
#'
#' ## get the gates_df, optionally select relevant gates, and modify
#' gates <- fcexpr::gs_get_gates(gs, n_bins = 100^2)
#' gates <- gates[which(gates$Population %in% c("CD8+", "CD8-")),]
#' gates$facet_strip <- T
#'
#' ## selected gates; to order facetted plot the inline factor level
#' ## ordering is required as flowCore::pData(gs) cannot contain factors
#' ## axis.text = element_blank() is part of ... and will omit axis numbers (passed to ggplot2::theme)
#' plotlist <-
#' fcexpr::plot_gates(gs = gs,
#' gates_df = gates,
#' facetting = facet_grid(cols = vars(factor(diluton_factor, levels = c(unique(p.df$diluton_factor)))), rows = vars(CD8_biotin_batch)),
#' axis.text = element_blank())
#'
#' ## paste plots together with patchwork and save
#' ## patchwork is superior to cowplot as is will completely ignore ommitted facet_strips
#' ggsave(patchwork::wrap_plots(plotlist, ncol = 1), filename = paste0("facsplots.png"), device = "png", path = im_path, dpi = "retina", width = 18, height = 7)
#' }
plot_gates <- function(gs,
                       gates_df,
                       facetting = NULL,
                       plot_gates = T,
                       plot_gate_names = T,
                       plot_gate_pct = T,
                       inverse_trans = F,
                       geom = c("hex", "pointdensity", "scattermore"),
                       gate_color = "black",
                       gate_stats_color = "black",
                       pct_digits = 1,
                       plot_contours = F,
                       contour_args = list(fill = "white",
                                           geom = "polygon",
                                           color = "black",
                                           contour_var = "ndensity",
                                           breaks = seq(0.05,0.95,0.1),
                                           alpha = 0.8,
                                           linewidth = 0.2),
                       col_pal = RColorBrewer::brewer.pal(11, "Spectral"),
                       col_pal_trans = "pseudo_log",
                       theme = ggplot2::theme_bw(),
                       theme_args = list(panel.grid = element_blank(),
                                         strip.background = element_rect(fill = "grey95", color = "white"),
                                         axis.text = element_blank(),
                                         axis.ticks = element_blank(),
                                         legend.position = "none"),
                       ...) {

  if (!requireNamespace("grDevices", quietly = T)) {
    utils::install.packages("grDevices")
  }
  if (!requireNamespace("RColorBrewer", quietly = T)) {
    utils::install.packages("RColorBrewer")
  }
  if (!requireNamespace("BiocManager", quietly = T)) {
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("ggcyto", quietly = T)) {
    BiocManager::install("ggcyto")
  }
  if (!requireNamespace("purrr", quietly = T)) {
    utils::install.packages("purrr")
  }

  dots <- list(...)

  geom <- match.arg(geom, c("hex", "pointdensity", "scattermore"))

  if (plot_contours) {
    message("Caution: Contour lines are affected across multiple facets.")
  }

  out <- purrr::flatten(lapply(unique(gates_df$GateLevel), function (z) {
    g <- gates_df[which(gates_df[,"GateLevel"] == z),]

    p <- lapply(split(g, paste(g$GateLevel, g$Parent, g$x, g$y, sep = "__")), function(gg) {

      my.filter <- if (gg[1,"marginalFilter"]) {ggcyto::marginalFilter} else {NULL}

      if ("max_nrow_to_plot" %in% names(dots)) {
        max_nrow_to_plot <- dots[["max_nrow_to_plot"]]
      } else {
        max_nrow_to_plot <- switch(geom,
                                   "hex" = 5e4,
                                   "pointdensity" = 2000,
                                   "scattermore" = 2e6)

      }

      p <- ggcyto::ggcyto(gs, subset = gg[1,"Parent"], filter = my.filter, ggplot2::aes(!!rlang::sym(gg[1,"x"]), !!rlang::sym(gg[1,"y"])), max_nrow_to_plot = max_nrow_to_plot)

      if (geom == "hex") {
        p <-
          p +
          ggplot2::geom_hex(binwidth = gg[1,"binwidths"][[1]]) +
          ggplot2::scale_fill_gradientn(colours = rev(grDevices::colorRampPalette(col_pal, interpolate = "linear")(100)), trans = col_pal_trans)
      }

      if (geom == "pointdensity") {
        if (!requireNamespace("ggpointdensity", quietly = T)) {
          utils::install.packages("ggpointdensity")
        }
        temp_dots <- dots[which(grepl("^pointdensity__", names(dots), ignore.case = T))]
        names(temp_dots) <- gsub("^pointdensity__", "", names(temp_dots), ignore.case = T)
        if (!"adjust" %in% names(temp_dots)) {
          temp_dots <- c(temp_dots, list(adjust = 5))
        }
        if (!"size" %in% names(temp_dots)) {
          temp_dots <- c(temp_dots, list(size = 0.3))
        }
        p <-
          p +
          do.call(ggpointdensity::geom_pointdensity, args = temp_dots) +
          ggplot2::scale_color_gradientn(colours = rev(grDevices::colorRampPalette(col_pal, interpolate = "linear")(100)), trans = col_pal_trans)
      }

      if (geom == "scattermore") {
        if (!requireNamespace("devtools", quietly = T)) {
          utils::install.packages("devtools")
        }
        if (!requireNamespace("scattermore", quietly = T)) {
          devtools::install_github("exaexa/scattermore")
        }

        temp_dots <- dots[which(grepl("^scattermore__", names(dots), ignore.case = T))]
        names(temp_dots) <- gsub("^scattermore__", "", names(temp_dots), ignore.case = T)
        p <- p + do.call(scattermore::geom_scattermore, args = temp_dots)
      }


      if (plot_contours) {
        p <- p + do.call(ggplot2::stat_density_2d, args = contour_args)
      }

      p <-
        p +
        ggplot2::xlab(gg[1,"x_lab"]) +
        ggplot2::ylab(gg[1,"y_lab"]) +
        ggcyto::ggcyto_par_set(limits = list(x = c(gg[1,"x_lowlim"], gg[1,"x_uplim"]), y = c(gg[1,"y_lowlim"], gg[1,"y_uplim"])))

      if (inverse_trans) {
        p <-
          p +
          ggcyto::axis_x_inverse_trans() +
          ggcyto::axis_y_inverse_trans()
      }

      p <- p + theme
      p <- p + do.call(ggplot2::theme, args = theme_args)

      if (!is.null(facetting)) {
        p <- p + facetting
      }

      if (all(!gg$facet_strip)) {
        p <- p + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_blank())
      }

      if (plot_gates) {
        for (i in 1:nrow(gg)) {
          p <- p + ggcyto::geom_gate(gg[i,"PopulationFullPath"], colour = gate_color)
        }
      }
      if (plot_gate_names) {
        for (i in 1:nrow(gg)) {
          p <- p + ggcyto::geom_stats(gg[i,"PopulationFullPath"], type = "gate_name", size = gg[i,"statsize_name"], colour = gate_stats_color, adjust = c(gg[i,"x_statpos_name"], gg[i,"y_statpos_name"]), fill = scales::alpha(c("white"),0.5))
        }
      }
      if (plot_gate_pct) {
        for (i in 1:nrow(gg)) {
          p <- p + ggcyto::geom_stats(gg[i,"PopulationFullPath"], digits = pct_digits, type = "percent", size = gg[i,"statsize_pct"], colour = gate_stats_color, adjust = c(gg[i,"x_statpos_pct"], gg[i,"y_statpos_pct"]), fill = scales::alpha(c("white"),0.5))
        }
      }
      return(p)
    })
    return(p)
  }))


  return(lapply(out, function(x) ggcyto::as.ggplot(x)))
}


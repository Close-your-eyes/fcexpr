#' Plot a default layout of gated populations of a set of FCS files
#'
#' Getting an overview of the results of a flow cytometry experiment can be a
#' generic task if gating and analysis strategy is fixed.
#' Here, ggcyto is used to create an arrangement of pseudocolor plots. A number
#' of default settings has been selected which may guarantee
#' good looking plots for most cases. Upon problems, the function may be
#' improved to handle edge cases. Currently, customization is limited to a few
#' rather technical aspects of the plots. This may be subject to expansion but
#' generally the range of customization is too large to fit it into a generic
#' function with a finite number of arguments. Hence, for very specific and
#' detailed requirements manual plotting and fiddling is unavoidable. For that,
#' the code of this function may serve as a template.
#'
#' @param gs a gatingset, e.g. made with fcexpr::wsp_get_gs
#' @param gates_df a data frame with information of how to plot gates,
#' made with fcexpr::gs_get_gates
#' @param facetting which facetting to apply, facet_wrap or facet_grid
#' with respective arguments, check flowCore::pData(gs) for valid columns;
#' facet_null will completely remove facets; by default facetting is done across
#' each fcs file
#' @param plot_gates logical whether to plot gates
#' @param plot_gate_names logical whether to plot gate names
#' @param plot_gate_pct logical whether to plot gate percentages
#' (fraction of parent)
#' @param inverse_trans logical whether to inverse transform axes numbers
#' if TRUE this will make axes look like in flowjo
#' @param geom how to plot data; recommendation is hex; hex = geom_hex taking
#' the binwidths column of gates_df into account,
#' pointdensity = ggpointdensity::geom_pointdensity
#' scattermore = scattermore::geom_scattermore (only black-white currently)
#' @param gate_stats_color font color of gate statistics
#' @param gate_color line color of gates
#' @param plot_contours logical whether to plot contour lines on top with
#' ggplot2::geom_density_2d
#' @param pct_digits how many digits after comma to print;
#' passed to 'digits' of ggcyto::geom_stats
#' @param col_pal color palette to use for color gradient generation
#' @param col_pal_trans argument passed as 'trans' in scale_fill_gradientn
#' @param plot_title do plot titles?
#' @param title title format
#' @param contour_args arguments to stat_density_2d
#' @param theme ggplot theme
#' @param max_nrow_to_plot to ggcyto::ggcyto; defaults are chosen with respect
#' to plotting time
#' @param geom_args arguments to geom
#' @param theme_args arguments to ggplot2::theme
#'
#' @return a list of ggplot2 objects, one for every gating level;
#' each list index contains respective plots for every fcs file
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
#' facetting = facet_grid(cols = vars(factor(dilution_factor, levels = c(unique(p.df$dilution_factor)))), rows = vars(CD8_biotin_batch)),
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
                       plot_title = T,
                       title = c("full_path", "short_path", "final_node"),
                       contour_args = list(fill = "white",
                                           geom = "polygon",
                                           color = "black",
                                           contour_var = "ndensity",
                                           breaks = seq(0.05,0.95,0.1),
                                           alpha = 0.8,
                                           linewidth = 0.2),
                       col_pal = colrr::col_pal("Spectral", direction = -1),
                       col_pal_trans = "pseudo_log",
                       theme = ggplot2::theme_bw(),
                       theme_args = list(panel.grid = ggplot2::element_blank(),
                                         strip.background = ggplot2::element_rect(fill = "grey95", color = "white"),
                                         axis.text = ggplot2::element_blank(),
                                         axis.ticks = ggplot2::element_blank(),
                                         axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 2, unit = "pt")),
                                         axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 2, unit = "pt")),
                                         legend.position = "none",
                                         strip.text.x = ggplot2::element_text(margin = ggplot2::margin(2,0,2,0, unit = "pt")),
                                         strip.text.y = ggplot2::element_text(margin = ggplot2::margin(0,2,0,2, unit = "pt")),
                                         plot.margin = grid::unit(c(1,1,1,1), "pt"),
                                         panel.spacing = grid::unit(2, "pt")),
                       max_nrow_to_plot = switch(geom,
                                                 "hex" = 5e4,
                                                 "pointdensity" = 2000,
                                                 "scattermore" = 2e6),
                       geom_args = list()) {

  if (!requireNamespace("ggcyto", quietly = T)) {
    BiocManager::install("ggcyto")
  }
  if (!requireNamespace("colrr", quietly = T)) {
    devtools::install_github("Close-your-eyes/colrr")
  }

  if (geom == "scattermore") {
    if (!requireNamespace("devtools", quietly = T)) {
      utils::install.packages("devtools")
    }
    if (!requireNamespace("scattermore", quietly = T)) {
      devtools::install_github("exaexa/scattermore")
    }
  }

  if (geom == "pointdensity") {
    if (!requireNamespace("ggpointdensity", quietly = T)) {
      utils::install.packages("ggpointdensity")
    }
  }

  geom <- rlang::arg_match(geom)
  title <- rlang::arg_match(title)

  if (plot_contours) {
    message("Caution: Contour lines are calculated across multiple facets.")
  }

  geom_fun <- switch(geom,
                     "hex" = ggplot2::geom_hex,
                     "pointdensity" = ggpointdensity::geom_pointdensity,
                     "scattermore" = scattermore::geom_scattermore)

  if (geom == "pointdensity") {
    if (!"adjust" %in% names(geom_args)) {
      geom_args <- c(geom_args, list(adjust = 5))
      message("pointdensity: adjust=5")
    }
    if (!"size" %in% names(geom_args)) {
      geom_args <- c(geom_args, list(size = 0.3))
      message("pointdensity: size=0.3")
    }
  }

  scale_fun <- switch(geom,
                      "hex" = ggplot2::scale_fill_gradientn,
                      "pointdensity" = ggplot2::scale_color_gradientn)


  out <- purrr::flatten(lapply(sort(unique(gates_df$GateLevel)), function (z) {
    g <- gates_df[which(gates_df[,"GateLevel"] == z),]

    p <- lapply(split(g, paste(g$GateLevel, g$Parent, g$x, g$y, sep = "__")), function(gg) {

      if (geom == "hex" && !"binwidth" %in% names(geom_args)) {
        geom_args <- c(geom_args, list(binwidth = gg[1,"binwidths"][[1]]))
      }

      p <- ggcyto::ggcyto(
        data = gs,
        subset = gg[1,"Parent"],
        filter = if (gg[1,"marginalFilter"]) {ggcyto::marginalFilter} else {NULL},
        mapping = ggplot2::aes(!!rlang::sym(gg[1,"x"]), !!rlang::sym(gg[1,"y"])),
        max_nrow_to_plot = max_nrow_to_plot) +
        do.call(what = geom_fun, args = geom_args) +
        theme +
        do.call(ggplot2::theme, args = theme_args)

      if (geom == "hex") {
        p$scales$scales <- list() # Remove any scales to avoid message of new fill scale
      }
      if (geom != "scattermore") {
        p <- p + do.call(scale_fun, args = list(colors = col_pal, trans = col_pal_trans))
      }

      if (plot_contours) {
        p <- p + do.call(ggplot2::stat_density_2d, args = contour_args)
      }
      # capture.output only to suppress text about coord system
      bin <- suppressMessages(capture.output(
        p <- p +
          ggplot2::xlab(gg[1,"x_lab"]) +
          ggplot2::ylab(gg[1,"y_lab"]) +
          ggcyto::ggcyto_par_set(limits = list(
            x = c(gg[1,"x_lowlim"], gg[1,"x_uplim"]),
            y = c(gg[1,"y_lowlim"], gg[1,"y_uplim"])
          ))
      ))

      if (inverse_trans) {
        p <- p +
          ggcyto::axis_x_inverse_trans() +
          ggcyto::axis_y_inverse_trans()
      }

      if (!is.null(facetting)) {
        p <- p + facetting
      }

      if (all(!gg$facet_strip)) {
        p <- p + ggplot2::theme(strip.background = ggplot2::element_blank(),
                                strip.text.x = ggplot2::element_blank(),
                                strip.text.y = ggplot2::element_blank())
      }

      if (plot_gates) {
        for (i in 1:nrow(gg)) {
          p <- p + ggcyto::geom_gate(
            gg[i,"PopulationFullPath"],
            colour = gate_color
          )
        }
      }
      if (plot_gate_names) {
        for (i in 1:nrow(gg)) {
          p <- p + ggcyto::geom_stats(
            gg[i,"PopulationFullPath"],
            type = "gate_name",
            size = gg[i,"statsize_name"],
            colour = gate_stats_color,
            adjust = c(gg[i,"x_statpos_name"], gg[i,"y_statpos_name"]),
            fill = scales::alpha(c("white"),0.5)
          )
        }
      }
      if (plot_gate_pct) {
        for (i in 1:nrow(gg)) {
          p <- p + ggcyto::geom_stats(
            gg[i,"PopulationFullPath"],
            digits = pct_digits,
            type = "percent",
            size = gg[i,"statsize_pct"],
            colour = gate_stats_color,
            adjust = c(gg[i,"x_statpos_pct"], gg[i,"y_statpos_pct"]),
            fill = scales::alpha(c("white"),0.5)
          )
        }
      }

      if (!plot_title) {
        p <- p + ggplot2::labs(title = NULL)
      } else {
        if (title == "short_path") {
          p <- p + ggplot2::labs(title = gg$Population)
        } else if (title == "final_node") {
          p <- p + ggplot2::labs(title = rev(strsplit(gg$Population, "/")[[1]])[1])
        }
      }
      attr(p, "Population") <- paste(gg$Population, collapse = "__")

      return(p)
    })
    return(p)
  }))

  popattr <- purrr::map_chr(out, attr, "Population")
  out <- purrr::map(out, ggcyto::as.ggplot)
  attr(out, "Population") <- popattr

  return(out)
}


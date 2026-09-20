#' Plot gated populations in a GatingSet
#'
#' Create ggcyto plots for populations described by `gates_df`. Gates sharing
#' a gating level, parent population, and x and y channels appear in one plot.
#' By default, ggcyto facets the samples in the GatingSet.
#'
#' @param gs A GatingSet, for example an element of the `gs_list` returned by
#'   [wsp_get_gs()].
#' @param gates_df A data frame from [gs_get_gates()]. Its columns control
#'   channels, limits, bin widths, gate display, labels, and facet strips.
#' @param facetting A ggplot2 facet specification, such as
#'   `ggplot2::facet_wrap()` or `ggplot2::facet_grid()`. Use columns from
#'   `flowCore::pData(gs)` for facet variables. `NULL` retains ggcyto's default
#'   sample facets; `ggplot2::facet_null()` removes facets.
#' @param plot_gates,plot_gate_names,plot_gate_pct Logical values that override
#'   the corresponding `plot_gate`, `plot_gate_name`, and `plot_gate_pct` columns
#'   of `gates_df` for every gate. The default `"gates_df"` uses each column's
#'   existing values. Percentages are relative to the parent population.
#' @param inverse_trans If `TRUE`, show inverse-transformed axis labels, as in
#'   FlowJo.
#' @param geom Event geometry: `"hex"` uses `ggplot2::geom_hex()` and the
#'   `binwidths` column of `gates_df`; `"pointdensity"` uses
#'   `ggpointdensity::geom_pointdensity()`; `"scattermore"` uses
#'   `scattermore::geom_scattermore()` without a color gradient.
#' @param gate_stats_color Color of gate name and percentage labels.
#' @param pct_digits Number of decimal places in gate percentages, passed to
#'   `ggcyto::geom_stats()`.
#' @param plot_contours If `TRUE`, add density contours. These are calculated
#'   across facets.
#' @param plot_title If `TRUE`, show a title based on the parent population.
#' @param title Title format: `"final_node"`, `"short_path"`, or `"full_path"`.
#' @param title_superscript If `TRUE`, render plus and minus signs in titles as
#'   superscripts.
#' @param contour_args Named list passed to `ggplot2::stat_density_2d()` when
#'   `plot_contours` is `TRUE`.
#' @param col_pal Colors for the hex fill or point-density color gradient.
#' @param col_pal_trans Transformation for that gradient. It does not apply to
#'   `"scattermore"`.
#' @param theme Base ggplot2 theme.
#' @param theme_args Named list passed to `ggplot2::theme()`.
#' @param theme_args_repl Named list that replaces matching entries in
#'   `theme_args` before `style_preset` is applied.
#' @param max_nrow_to_plot Maximum events passed to `ggcyto::ggcyto()` per plot.
#'   Defaults depend on `geom`. This limit does not apply when a gate uses
#'   `ggcyto::marginalFilter`.
#' @param geom_args Named list passed to the selected event geometry.
#' @param gate_args Named list passed to `ggcyto::geom_gate()`.
#' @param as_ggplot If `TRUE`, convert each result with `ggcyto::as.ggplot()`.
#' @param style_preset One of `"technical"` (show axis text and ticks),
#'   `"clean"` (hide axis text, ticks, and grid), or `"none"` (add no preset).
#'
#' @return A list of plots, one per combination of gating level, parent
#'   population, and x and y channels. Plots are ggcyto objects unless
#'   `as_ggplot = TRUE`. The list has a `"Population"` attribute with the
#'   population names represented by each plot.
#' @export
#'
#' @examples
#' \dontrun{
#' gs_data <- fcexpr::wsp_get_gs(wsp = "path/to/workspace.wsp")
#' gs <- gs_data$gs_list[[1]]
#' gates <- gs_data$gate_dfs[[1]]
#'
#' # Show gate outlines but omit gate statistics; set a smaller title.
#' plots <- fcexpr::plot_gates(
#'   gs = gs,
#'   gates_df = gates,
#'   plot_gate_names = FALSE,
#'   plot_gate_pct = FALSE,
#'   theme_args_repl = list(plot.title = ggplot2::element_text(size = 10))
#' )
#' print(plots[[1]])
#' }
plot_gates <- function(gs,
                       gates_df,
                       facetting = NULL,
                       plot_gates = "gates_df",
                       plot_gate_names = "gates_df",
                       plot_gate_pct = "gates_df",
                       inverse_trans = T,
                       geom = c("hex", "pointdensity", "scattermore"),
                       gate_stats_color = "black",
                       pct_digits = 1,
                       plot_contours = F,
                       plot_title = T,
                       title = c("final_node", "short_path", "full_path"),
                       title_superscript = F,
                       contour_args = list(fill = "white",
                                           geom = "polygon",
                                           color = "black",
                                           contour_var = "ndensity",
                                           breaks = seq(0.05,0.95,0.1),
                                           alpha = 0,
                                           linewidth = 0.2),
                       col_pal = colrr::col_pal("Spectral", direction = -1),
                       col_pal_trans = "pseudo_log",
                       theme = ggplot2::theme_bw(),
                       theme_args = list(strip.background = ggplot2::element_rect(color = NA),
                                         axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 2, unit = "pt")),
                                         axis.title.y = ggplot2::element_text(margin = ggplot2::margin(b = 2, unit = "pt")),
                                         strip.text.x = ggplot2::element_text(margin = ggplot2::margin(2,0,2,0, unit = "pt")),
                                         strip.text.y = ggplot2::element_text(margin = ggplot2::margin(0,2,0,2, unit = "pt")),
                                         plot.margin = ggplot2::margin(1,1,1,1, "pt"),
                                         plot.title = ggplot2::element_text(margin = ggplot2::margin(1,1,1,1, unit = "pt"), size = 12),
                                         panel.spacing = grid::unit(2, "pt"),
                                         legend.position = "none"),
                       theme_args_repl = list(),
                       max_nrow_to_plot = switch(geom,
                                                 "hex" = 5e4,
                                                 "pointdensity" = 2000,
                                                 "scattermore" = 2e6),
                       geom_args = list(),
                       gate_args = list(colour = "black",
                                        linewidth = 0.3),
                       as_ggplot = F,
                       style_preset = c("technical", "clean", "none")) {
  fcexpr:::.ensure_packages(c("colrr", "ggcyto", "ggplot2", "ggtext", "Gmisc", "scales"))

  geom <- rlang::arg_match(geom)
  title <- rlang::arg_match(title)
  style_preset <- rlang::arg_match(style_preset)

  for (i in names(theme_args_repl)) {
    theme_args[[i]] <- theme_args_repl[[i]]
  }


  if (geom == "scattermore") {
    fcexpr:::.ensure_package("scattermore")
  }

  if (geom == "pointdensity") {
    fcexpr:::.ensure_package("ggpointdensity")
  }

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

  if (is.logical(plot_gates)) {
    gates_df$plot_gate <- plot_gates
  }
  if (is.logical(plot_gate_names)) {
    gates_df$plot_gate_name <- plot_gate_names
  }
  if (is.logical(plot_gate_pct)) {
    gates_df$plot_gate_pct <- plot_gate_pct
  }

  conv <- stats::setNames(gates_df$Population, gates_df$PopulationFullPath)
  gates_df$Parent_short <- conv[gates_df$Parent]
  gates_df$Parent_short[which(is.na(gates_df$Parent_short))] <- ""

  if (any(gates_df$marginalFilter)) {
    message("max_nrow_to_plot does not apply to gates with ggcyto::marginalFilter enabled.")
  }

  if (style_preset == "technical") {
    repl <- list(panel.grid.minor = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(),
                 axis.text.y = ggplot2::element_text(),
                 axis.ticks.x = ggplot2::element_line(),
                 axis.ticks.y = ggplot2::element_line())
    for (i in names(repl)) {
      theme_args[[i]] <- repl[[i]]
    }
  } else if (style_preset == "clean") {
    repl <- list(panel.grid = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())
    for (i in names(repl)) {
      theme_args[[i]] <- repl[[i]]
    }
  }


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
      bin <- suppressMessages(utils::capture.output(
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
      } else {
        p <- p +
          ggplot2::scale_x_continuous(expand = ggplot2::expansion()) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion())
      }

      if (!is.null(facetting)) {
        p <- p + facetting
      }

      if (all(!gg$facet_strip)) {
        p <- p + ggplot2::theme(strip.background = ggplot2::element_blank(),
                                strip.text.x = ggplot2::element_blank(),
                                strip.text.y = ggplot2::element_blank())
      }

      for (i in 1:nrow(gg)) {

        if (gg[i,"plot_gate"]) {
          p <- p + Gmisc::fastDoCall(what = ggcyto::geom_gate,
                                     args = c(list(data = gg[i,"PopulationFullPath"]),
                                              gate_args))
        }
        if (gg[i,"plot_gate_name"]) {
          p <- p + ggcyto::geom_stats(
            gate = gg[i,"PopulationFullPath"],
            type = "gate_name",
            size = gg[i,"statsize_name"],
            colour = gate_stats_color,
            adjust = c(gg[i,"x_statpos_name"], gg[i,"y_statpos_name"]),
            fill = scales::alpha(c("white"),0.5)
          )
        }
        if (gg[i,"plot_gate_pct"]) {
          p <- p + ggcyto::geom_stats(
            gate = gg[i,"PopulationFullPath"],
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
        titlechr <- switch(title,
                           short_path = gg$Parent_short[1],# multiple boolean parents?
                           final_node = gsub("root", "", rev(strsplit(gg$Parent[1], "/")[[1]])[1]),
                           full_path = gsub("root", "", gg$Parent[1]))
        if (title_superscript && !trimws(titlechr) == "") {
          titlechr <- gsub("-", "<sup>-</sup>", titlechr, fixed = T)
          titlechr <- gsub("+", "<sup>+</sup>", titlechr, fixed = T)
          p <- p + ggplot2::theme(plot.title = ggtext::element_markdown(margin = ggplot2::margin(1,1,1,1, unit = "pt"), size = 12))
        }
        p <- p + ggplot2::labs(title = titlechr)
      }
      attr(p, "Population") <- paste(gg$Population, collapse = "__")


      p <- p +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion()) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion())
      return(p)
    })
    return(p)
  }))

  # browser()

  # dims <- attr(out[[1]][["data"]], "dims")
  # dims <- dims[axis != "order", ]
  # as.ggplot(out[[1]])


  popattr <- purrr::map_chr(out, attr, "Population")
  # spare that until fix comes https://github.com/RGLab/ggcyto/issues/108
  if (as_ggplot) {
    out <- purrr::map(out, ggcyto::as.ggplot)
  }
  attr(out, "Population") <- popattr

  return(out)
}

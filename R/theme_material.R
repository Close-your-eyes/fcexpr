#' ggplot2 theme with style of Rstudios color palette 'Material'
#'
#' @param base_size base font size, given in pts.
#' @param text_color color of various text elements
#' @param bg_color background color, not intended to be changed dramatically
#' @param bg_color2 background color 2, e.g. for facet strip
#' @param text_fun function for text elements, ggplot2::element_text or
#' ggtext::element_markdown
#' @param style appearance of axes and grid
#' @param white use a light color scale with white bg
#' @param legend_tight legend more compact and closer to plot panel
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
theme_material <- function(base_size = 14,
                           text_color = "white",
                           bg_color = "#273238",
                           bg_color2 = "#4F5C69",
                           text_fun = ggtext::element_markdown,
                           style = c("prism", "prismy", "prismx", "base"),
                           white = F,
                           legend_tight = F) {

  style <- rlang::arg_match(style)
  if (white) {
    text_color = "black"
    bg_color = "white"
    bg_color2 = "grey90"
    axes_color = "black"
    grid_color = "grey60"
  } else {
    axes_color <- "white"
    grid_color = bg_color2
  }

  thistheme <- ggplot2::theme_grey(base_size = base_size) +
    ggplot2::theme(
      # Text elements
      text = ggplot2::element_text(color = text_color),
      plot.title = text_fun(color = text_color, size = base_size * 1.2, face = "bold"),
      plot.subtitle = text_fun(color = text_color, size = base_size),
      plot.caption = text_fun(color = text_color),
      plot.tag = text_fun(color = text_color),
      axis.title = text_fun(color = text_color),
      axis.title.x = text_fun(color = text_color),
      axis.title.y = text_fun(color = text_color),
      axis.text = text_fun(color = text_color),
      axis.text.x = text_fun(color = text_color),
      axis.text.y = text_fun(color = text_color),
      legend.title = text_fun(color = text_color),
      legend.text = ggplot2::element_text(color = text_color), # this as element_markdown locks height of legend bar
      strip.text = text_fun(color = text_color),
      strip.text.x = text_fun(color = text_color),
      strip.text.y = text_fun(color = text_color),

      # Backgrounds
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      strip.background = ggplot2::element_rect(fill = bg_color2, color = "white"),
      legend.background = ggplot2::element_blank(),

      # axes and misc
      axis.ticks = ggplot2::element_line(color = axes_color),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid = ggplot2::element_line(linewidth = 0.3, color = grid_color),
      panel.border = ggplot2::element_rect(fill = NA, color = bg_color2, linewidth = 1.5),
      legend.key = ggplot2::element_blank(),

      plot.margin = ggplot2::margin(t=5, r=5, b=5, l=5, unit = "pt"),
      legend.margin = ggplot2::margin(t=5, r=5, b=5, l=5, unit = "pt"),
      legend.position = "right",
      legend.direction = "vertical"
    )

  if (legend_tight) {
    thistheme <-
      thistheme +
      ggplot2::theme(
        legend.key.size = ggplot2::unit(5, "pt"),
        legend.key.spacing.y = ggplot2::unit(2, "pt"),
        legend.key.spacing.x = ggplot2::unit(2, "pt"),
        # negative left: bring closer to panel
        legend.margin = ggplot2::margin(t=5, r=5, b=5, l=-10, unit = "pt")
      )
  }

  #guides(color = guide_legend(override.aes = list(size = 3)))

  if (style %in% c("prism", "prismy")) {
    thistheme <-
      thistheme +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(color = axes_color),
        axis.line.y = ggplot2::element_line(color = axes_color),
        panel.border = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = bg_color2, color = NA)
      )
    if (style == "prismy") {
      thistheme <-
        thistheme +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_line(color = grid_color, linewidth = 0.3)
        )
    }
    if (style == "prismx") {
      thistheme <-
        thistheme +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_line(color = grid_color, linewidth = 0.3)
        )
    }
  }
  return(thistheme)
}



#' Add default legend guide elements to ggplot
#'
#' Make default barheight = 10 and width = 1, legend title as
#' ggtext::element_markdown, point legends of size = 4. Valid for fill and color.
#'
#' @param colorbar plot color legend as bar
#' @param ... arguments to guide_colorbar or guide_legend
#'
#' @return
#' @export
#'
#' @examples
guides_default <- function(colorbar = F, ...) {
  # fill and color !!!
  ggplot2::guides(
    fill = if (colorbar) {
      ggplot2::guide_colorbar(barwidth = 1,
                              barheight = 10,
                              title.theme = ggtext::element_markdown(),
                              ...
      )
    } else {
      ggplot2::guide_legend(override.aes = list(size = 4),
                            title.theme = ggtext::element_markdown(),
                            ...
                            # label.theme = ggtext::element_markdown(), # blocks barheight, bug
                            # reverse = F,
                            # nrow = NULL,
                            # ncol = 1,
                            # label.position = "right",
                            # title.position = "top",
                            # order = 1,
                            # title.position = "top",
                            # title = "..auto..", # character or NULL or ..auto..
                            # title.hjust = 0.5
      )
    },
    color = if (colorbar) {
      ggplot2::guide_colorbar(barwidth = 1,
                              barheight = 10,
                              title.theme = ggtext::element_markdown(),
                              ...
      )
    } else {
      ggplot2::guide_legend(override.aes = list(size = 4),
                            title.theme = ggtext::element_markdown(),
                            ...
      )
    }
  )
}

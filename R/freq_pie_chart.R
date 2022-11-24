.freq_pie_chart <- function(x,
                            meta.col,
                            inset.text.size = 5,
                            inset.text.radius = 0.75,
                            legend.position = "right",
                            col_pal = fcexpr::col_pal("custom"),
                            order_pieces = T) {

  if (!requireNamespace("ggforce", quietly = T)) {
    utils::install.packages("ggforce")
  }
  if (!requireNamespace("farver", quietly = T)) {
    utils::install.packages("farver")
  }

  ## add option to repel label or to increase their position to furhter outside gradually

  # https://stackoverflow.com/questions/16184188/ggplot-facet-piechart-placing-text-in-the-middle-of-pie-chart-slices (ggforce)
  tab <- table(x[,meta.col], exclude = c())
  tab <- data.frame(frac = as.numeric(tab/sum(tab)), cluster = factor(names(tab), levels = names(tab)))
  if (order_pieces) {
    tab <- tab[order(tab$frac, decreasing = T), ]
  }
  tab$start_angle <- c(0,cumsum(tab$frac))[-(length(tab$frac) + 1)]*pi*2
  tab$end_angle <- c(cumsum(tab$frac))*pi*2
  tab$mid_angle <-  0.5*(tab$start_angle + tab$end_angle)

  if (length(col_pal) != length(unique(tab[,"cluster"]))) {
    if (is.null(names(col_pal))) {
      if (length(col_pal) < length(unique(tab[,"cluster"]))) {
        col_pal <- scales::hue_pal()(length(unique(tab[,"cluster"])))
        warning("Number of colors provided not sufficient for number of factor levels. Falling back to scales::hue_pal().")
      } else {
        col_pal <- col_pal[1:length(unique(tab[,"cluster"]))]
      }
    } else {
      if (length(col_pal) > length(unique(tab[,"cluster"])) && all(names(col_pal) %in% unique(tab[,"cluster"]))) {
        col_pal <- col_pal[unique(tab[,"cluster"])]
      } else {
        warning("Number of colors provided not matching the number of factor levels in meta.col. Falling back to scales::hue_pal().")
        col_pal <- scales::hue_pal()(length(unique(tab[,"cluster"])))
      }
    }
  } else {
    if (!is.null(names(col_pal)) && !all(names(col_pal) %in% unique(tab[,"cluster"]))) {
      warning("Not all names of col_pal found in factor levels of meta.col. Falling back to scales::hue_pal().")
      col_pal <- scales::hue_pal()(length(unique(tab[,"cluster"])))
    }
  }

  ggplot2::ggplot(tab, ggplot2::aes(x0 = 0, y0 = 0, r0 = 0.3, r = 1, start = start_angle, end = end_angle, fill = cluster)) +
    ggforce::geom_arc_bar(colour = "white") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = legend.position,
                   panel.grid = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
    ggplot2::geom_text(ggplot2::aes(color = farver::decode_colour(col_pal, to = "hcl")[,"l"] > 50,
                                    x = inset.text.radius*sin(mid_angle),
                                    y = inset.text.radius*cos(mid_angle),
                                    label = format(round(frac, 2), nsmall = 2)),
                       size = inset.text.size) +
    ggplot2::scale_fill_manual(values = col_pal) +
    ggplot2::scale_color_manual(guide = "none", values = c("white", "black")) +
    ggplot2::coord_fixed(ratio = 1)
}

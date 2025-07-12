#' Plot heatmap from data frame in long format
#'
#'
#'
#' @param df data frame in long format
#' @param groups column with groups
#' @param features column with features
#' @param values column with fill values for color scale
#' @param dotsizes column for dot sizes, if provided will cause plotting of
#' dots instead of tiles
#' @param fill color palette vector for fill of tiles or dots. when auto,
#' RColorBrewer::RdBu is used
#' @param color color of stroke (border) around tiles or dots; "NA" means no
#' stroke is plotted; NA has
#' to be put in quotation mark ("NA"), such that geom_point accepts it.
#' other choices may be black, white or any other color code; when "auto",
#' grey70 is used by default when !is.null(dotsizes) and a the number of
#' features is below 100.
#' @param scale how to scale values: not, zscore or from -1 to 1
#' @param omit_featurelabels do not plot axis text of features, e.g. when
#' their number is large
#' @param topn_features only plot the top n features (ordered by value) per
#' group
#' @param border_linewidth linewidth of borders (stroke) around tiles or dots
#' @param legendbreaks a single number, a vector of explicit breaks, or "auto"
#' for ggplot default or "minmidmax" for three breaks at minimum, middle and
#' maximum of value range
#' @param legendlabels labels for breaks, e.g. c("min", "mid", "max")
#' @param colorsteps NULL to have normal colorbar, auto for default colorsteps,
#' a single number or a vector of explicit steps; may not work with any number
#' when nice_colorsteps is TRUE
#' @param nice_colorsteps heuristic for pretty steps
#' @param flip_axes do flip them?
#' @param group_seplines plot lines that separate features belonging to
#' different groups
#' @param seplines_args arguments to geom_hline for seplines
#' @param legend_fill_args arguments to ggplot2::guide_colorsteps or
#' ggplot2::guide_colorbar, depend upon the color scale
#' @param legend_size_args arguments to ggplot2::guide_legend to modify the
#' size legend
#' @param ... arguments to heatmap_ordering
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' df <- readRDS(system.file("extdata", "heatmap_df.rds", package = "fcexpr"))
#' # everything default
#' fcexpr::heatmap_long_df(df = df,
#'                 groups = "cluster",
#'                 features = "channel",
#'                 values = "mean_cluster_scale")
#' # -log10(pvalues) as dot size
#' fcexpr::heatmap_long_df(df = df,
#'                 groups = "cluster",
#'                 features = "channel",
#'                 values = "mean_cluster_scale",
#'                 dotsizes = "pvalue2")
#' # continuuos color scale, no feature axis text, flipped axes
#' # only 4 features per group, lines to separate features
#' fcexpr::heatmap_long_df(df = df,
#'                 groups = "cluster",
#'                 features = "channel",
#'                 values = "mean_cluster_scale",
#'                 dotsizes = "pvalue2",
#'                 omit_featurelabels = T,
#'                 flip_axes = T,
#'                 topn_features = 4,
#'                 group_seplines = T,
#'                 colorsteps = NULL)
#' # scale in function, alter legend text
#' fcexpr::heatmap_long_df(df = df,
#'                 groups = "cluster",
#'                 features = "channel",
#'                 values = "mean_cluster",
#'                 scale = "zscore",
#'                 colorsteps = NULL,
#'                 legendlabels = c("min", "mid", "max"),
#'                 legendbreaks = "minmidmax")
heatmap_long_df <- function(df,
                            groups,
                            features,
                            values,
                            dotsizes = NULL,
                            fill = "auto",
                            color = "auto",
                            scale = c("none", "zscore", "1"),
                            omit_featurelabels = F,
                            topn_features = NULL,
                            border_linewidth = 0.2,
                            legendbreaks = "auto",
                            legendlabels = "auto",
                            colorsteps = "auto",
                            nice_colorsteps = T,
                            flip_axes = F,
                            group_seplines = F,
                            seplines_args = list(),
                            legend_fill_args = list(
                              barwidth = 1,
                              barheight = 8,
                              order = 1
                              #ticks.colour = "black",
                              #frame.colour = "black",
                              #frame.linewidth = 0.1
                            ),
                            legend_size_args = list(
                              order = 2,
                              ncol = NULL,
                              nrow = NULL,
                              override.aes = list(color = "black")
                            ),
                            ...) {

  if (!requireNamespace("brathering", quietly = T)) {
    devtools::install_github("Close-your-eyes/brathering")
  }

  stopifnot("df must be a data frame" = is.data.frame(df))
  scale <- rlang::arg_match(scale)

  # optional filter for top n features per group
  if (!is.null(topn_features)) {
    select <-
      df |>
      dplyr::group_by(!!rlang::sym(features)) |>
      dplyr::slice_max(order_by = !!rlang::sym(values), n = 1) |>
      dplyr::ungroup() |>
      dplyr::group_by(!!rlang::sym(groups)) |>
      dplyr::slice_max(order_by = !!rlang::sym(values), n = topn_features) |>
      dplyr::pull(!!rlang::sym(features))
    df <- df[which(df[[features]] %in% select),]
  }

  # optional scaling
  df <-
    df |>
    dplyr::group_by(!!rlang::sym(features)) |>
    dplyr::mutate(!!values := dplyr::case_when(
      scale == "zscore" ~ as.vector(scale(!!rlang::sym(values))),
      scale == "1" ~ brathering::scale2(!!rlang::sym(values), min = -1, max = 1),
      .default = !!rlang::sym(values)  # fallback (optional)
    )) |>
    dplyr::ungroup()

  # assign factors to features and groups
  df <- fcexpr::heatmap_ordering(
    df = df,
    features = features,
    groups = groups,
    values = values,
    ...
  )


  if (color[1] == "auto") { # catch if length(color) > 1
    if (!is.null(dotsizes) || nlevels(df[["features"]]) > 100) {
      color <- "NA"
    } else {
      color <- "grey70"
    }
  }


  if (length(fill) == 1 && fill == "auto") {
    fill <- colrr::col_pal(name = "RColorBrewer::RdBu", n = 11, direction = -1)
  } else if (length(fill) == 1) {
    fill <- colrr::col_pal(name = fill)
  }

  # start ggplot pipeline
  plot <- ggplot2::ggplot(df, ggplot2::aes(
    x = !!rlang::sym(groups),
    y = !!rlang::sym(features),
    fill = !!rlang::sym(values)
  ))
  if (flip_axes) {
    plot <- plot + ggplot2::coord_flip()
  }

  if (!is.null(dotsizes)) {
    plot <- plot + ggplot2::geom_point(
      ggplot2::aes(size = !!rlang::sym(dotsizes)),
      shape = 21,
      color = color,
      stroke = border_linewidth
    )
  } else {
    plot <- plot + ggplot2::geom_tile(
      color = color,
      linewidth = border_linewidth
    )
  }

  # check if values are z-scored
  dfmat <- brathering::df_long_to_mat(
    df,
    to_rows = groups,
    to_cols = features,
    values = values
  )
  values_zscored <- sum(apply(dfmat, 2, brathering::is_z_scored, verbose = F)) > 0.75*ncol(dfmat) # 0.75: arbitrary choice

  # decide for colorsteps or continuuous colorbar
  scale_fill <- colorscale_heuristic(colorscale_values = df[[values]],
                                     values_zscored = values_zscored,
                                     colorsteps = colorsteps,
                                     legendbreaks = legendbreaks,
                                     legendlabels = legendlabels,
                                     fill = fill,
                                     nice_colorsteps = nice_colorsteps)
  if (grepl("coloursteps", scale_fill[["guide"]])) {
    guide_fun <- ggplot2::guide_colorsteps
  } else {
    guide_fun <- ggplot2::guide_colorbar
  }

  plot <-
    plot +
    scale_fill +
    ggplot2::guides(fill = Gmisc::fastDoCall(guide_fun, args = legend_fill_args),
                    size = Gmisc::fastDoCall(ggplot2::guide_legend, args = legend_size_args))

  if (omit_featurelabels) {
    if (flip_axes) {
      plot <- plot + theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    } else {
      plot <- plot + theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
    }

  }

  if (group_seplines) {
    dfsort <- df |>
      dplyr::group_by(!!rlang::sym(features)) |>
      dplyr::slice_max(order_by = !!rlang::sym(values), n = 1) |>
      dplyr::arrange(!!rlang::sym(groups))
    hlines <- cumsum(rle(as.character(dfsort[[groups]]))[["lengths"]]) + 0.5
    hlines <- hlines[-length(hlines)]
    plot <- plot + Gmisc::fastDoCall(ggplot2::geom_hline, args = c(list(yintercept = hlines), seplines_args))
  }

  return(plot)
}




colorscale_heuristic <- function(colorscale_values,
                                 values_zscored,
                                 colorsteps,
                                 legendbreaks,
                                 legendlabels,
                                 fill,
                                 nice_colorsteps,
                                 type = c("fill", "color"),
                                 col_na = "grey50",
                                 qmin = 0,
                                 qmax = 1,
                                 scale.max = NULL,
                                 scale.min = NULL) {

  #qmin, qmax for featureplot2 from scexpr, for correct limits of colorsteps, colorsteps must be auto or vector
  # scale.min: provided from scexpr featureplot2 but exclude non expressers (=0)
  type <- rlang::arg_match(type)
  if (is.null(scale.max)) {
    scale.max <- as.numeric(format(brathering::floor2(max(colorscale_values), 0.1), nsmall = 1))
  }
  if (is.null(scale.min)) {
    scale.min <- as.numeric(format(brathering::ceiling2(min(colorscale_values), 0.1), nsmall = 1))
  }
  scale.mid <- ifelse(values_zscored, 0, as.numeric(format(round(scale.min + ((scale.max - scale.min) / 2), 1), nsmall = 1)))

  if (is.null(colorsteps)) {
    if (length(legendbreaks) == 1 && legendbreaks == "auto") {
      legendbreaks <- ggplot2::waiver()
    } else if (length(legendbreaks) == 1 && legendbreaks == "minmidmax") {
      legendbreaks <- c(scale.min, scale.mid, scale.max)
    } else if (length(legendbreaks) == 1) {
      legendbreaks <- seq(scale.min, scale.max, length.out = legendbreaks)
    } else {
      # legendbreaks is vector
    }
    if (length(legendlabels) == 1 && legendlabels == "auto") {
      legendlabels <- ggplot2::waiver()
    } else if (length(legendlabels) != length(legendbreaks)) {
      message("length(legendlabels) != length(legendbreaks), using ggplot2 default")
      legendlabels <- ggplot2::waiver()
    }

    if (type == "fill") {
      scalefun <- ggplot2::scale_fill_gradientn
    } else {
      scalefun <- ggplot2::scale_color_gradientn
    }
    scale_fill <-
      scalefun(values = scales::rescale(c(scale.min, scale.mid, scale.max)),
               colors = fill,
               breaks = legendbreaks,
               labels = legendlabels,
               na.value = col_na)
  } else {
    if (length(colorsteps) == 1 && colorsteps == "auto") {
      # colorsteps is auto
      if (values_zscored) {
        colorsteps <- seq(round(scale.min), round(scale.max), 1)
      } else {
        colorsteps <- seq(round(scale.min), round(scale.max), length.out = 6)
      }
    }
    if (type == "fill") {
      scalefun <- ggplot2::scale_fill_stepsn
    } else {
      scalefun <- ggplot2::scale_color_stepsn
    }
    if (length(colorsteps) == 1) {
      # colorsteps given as one number
      scale_fill <-
        scalefun(colors = fill,
                 values = scales::rescale(c(scale.min, scale.mid, scale.max)),
                 n.breaks = colorsteps,
                 limits = c(scale.min, scale.max),
                 show.limits = T,
                 nice.breaks = nice_colorsteps,
                 na.value = col_na)
    } else {
      # colorsteps given as a vector
      if (qmin > 0 || qmax < 1) {
        # only to alter limit labels
        colorsteps <- round(seq(scale.min, scale.max, length.out = 6),1)
        if (qmin > 0) {min.lab <- paste0(scale.min, " (q", round(qmin*100, 0), ")")} else {min.lab <- scale.min}
        if (qmax < 1) {max.lab <- paste0(scale.max, " (q", round(qmax*100, 0), ")")} else {max.lab <- scale.max}

       colorstepbreaks <- colorsteps
       if (dplyr::near(colorstepbreaks[1], scale.min)) {
         colorstepbreaks <- colorstepbreaks[-1]
       }
       if (dplyr::near(colorstepbreaks[length(colorstepbreaks)], scale.max)) {
         colorstepbreaks <- colorstepbreaks[-length(colorstepbreaks)]
       }
        scale_fill <-
          scalefun(colors = fill,
                   values = scales::rescale(c(scale.min, scale.mid, scale.max)),
                   breaks = c(scale.min, colorstepbreaks, scale.max), # manually add limits as breaks
                   limits = c(scale.min, scale.max),
                   labels = c(min.lab, colorstepbreaks, max.lab),
                   na.value = col_na)
      } else {
        scale_fill <-
          scalefun(colors = fill,
                   values = scales::rescale(c(scale.min, scale.mid, scale.max)),
                   breaks = colorsteps,
                   limits = c(scale.min, scale.max),
                   show.limits = T,
                   na.value = col_na)
      }
      #labels = function(x) round(x, legend.decimals))
    }
  }
  return(scale_fill)
}

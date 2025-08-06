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
#' grey70 is used by default when is.null(dotsizes) and a the number of
#' features is below 100.
#' @param scale how to scale values: not, zscore or from -1 to 1
#' @param features_topn only plot the top n features (ordered by value) per
#' group
#' @param color_linewidth linewidth of borders (stroke) around tiles or dots
#' @param legendbreaks a single number, a vector of explicit breaks, or "auto"
#' for ggplot default or "minmidmax" for three breaks at minimum, middle and
#' maximum of value range
#' @param legendlabels labels for breaks, e.g. c("min", "mid", "max")
#' @param colorsteps NULL to have normal colorbar, auto for default colorsteps,
#' a single number or a vector of explicit steps; may not work with any number
#' when colorsteps_nice is TRUE
#' @param colorsteps_nice heuristic for pretty steps
#' @param axes_flip do flip them?
#' @param group_seplines plot lines that separate features belonging to
#' different groups
#' @param seplines_args arguments to geom_hline for seplines
#' @param legend_fill_args arguments to ggplot2::guide_colorsteps or
#' ggplot2::guide_colorbar, depend upon the color scale
#' @param legend_size_args arguments to ggplot2::guide_legend to modify the
#' size legend; e.g. use override.aes = list(size = c(1,3,5)) to adjust dot size
#' in legend in contrast to dotsize_range, one number for each dot in legend
#' needed
#' @param ... arguments to heatmap_ordering like feature_order or group_order
#' @param featurelabels subset of feature labels to plot; NULL to plot all,
#' "" to plot none; can be a named vector with names being aliases to use for
#' plotting e.g. c("CD20" = "MS4A1", "CD3", "KLRG1") to only alter MS4A1;
#' auto to omit labels by default when more than 100 features are there
#' @param featurelabels_repel do repel feature axis labels?
#' @param featuresitalic shorthand to plot feature labels in italic
#' @param theme_args arguments to ggplot2::theme
#' @param repel_args fine tuning for featurelabel repelling
#' @param theme ggplot2 theme
#' @param dotsize_range range for dot size
#' @param topn_cols columns in df to use for ordering when selecting
#' features_topn with dplyr::slice_max. hint: p-values have to multiplied
#' by -1 to make this work properly
#' @param topn_ties break ties for features_topn? if TRUE, more then
#' features_topn may be plotted
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
#'                 featurelabels_omit = T,
#'                 axes_flip = T,
#'                 features_topn = 4,
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
                            dotsize_range = c(2,7),
                            fill = "auto",
                            color = "auto",
                            scale = c("none", "zscore", "1"),
                            features_topn = NULL,
                            topn_cols = values,
                            topn_ties = F,
                            featurelabels = "auto",
                            featurelabels_repel = F,
                            featuresitalic = F,
                            color_linewidth = 0.2,
                            legendbreaks = "auto",
                            legendlabels = "auto",
                            colorsteps = "auto",
                            colorsteps_nice = T,
                            axes_flip = F,
                            group_seplines = F,
                            seplines_args = list(),
                            theme = ggplot2::theme_classic(),
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
                                                  #size = c(2, 7))
                            ),
                            theme_args = list(
                              panel.grid = element_blank()
                            ),
                            repel_args = list(featurelabels_width = 0.2,
                                              featurelabels_nudhe_x = -1),
                            ...) {

  if (!requireNamespace("brathering", quietly = T)) {
    devtools::install_github("Close-your-eyes/brathering")
  }

  stopifnot("df must be a data frame" = is.data.frame(df))
  scale <- rlang::arg_match(scale)

  if (featuresitalic) {
    theme_args <- brathering::gg_inject_theme_element(theme_args = theme_args,
                                                      elem = if (axes_flip) "axis.text.x" else "axis.text.y",
                                                      elem_sub = "face",
                                                      value = "italic")
  }

  # optional filter for top n features per group
  if (!is.null(features_topn)) {
    select <-
      df |>
      # max group per feature
      dplyr::slice_max(order_by = !!rlang::sym(values), n = 1, by = !!rlang::sym(features)) |>
      # then best features per group
      # dplyr::slice_max(order_by = !!rlang::sym(values), n = features_topn, by = !!rlang::sym(groups)) |>
      # dplyr::slice_max(order_by = tibble::tibble(!!rlang::sym("auc"), !!rlang::sym("padj"), !!rlang::sym("logFC")), n = features_topn, by = !!rlang::sym(groups)) |>
      dplyr::slice_max(order_by = tibble::tibble(!!!rlang::syms(as.list(topn_cols))), n = features_topn, by = !!rlang::sym(groups), with_ties = topn_ties) |>
      dplyr::pull(!!rlang::sym(features))
    df <- df[which(df[[features]] %in% select),,drop = F]
  }

  # optional scaling
  df <- dplyr::mutate(df, !!values := dplyr::case_when(
    scale == "zscore" ~ as.vector(scale(!!rlang::sym(values))),
    scale == "1" ~ brathering::scale2(!!rlang::sym(values), min = -1, max = 1),
    .default = !!rlang::sym(values)  # fallback (optional)
  ), .by = !!rlang::sym(features))

  # assign factors to features and groups
  df <- heatmap_ordering(
    df = df,
    features = features,
    groups = groups,
    values = values,
    ...
  )


  if (color[1] == "auto") { # catch if length(color) > 1
    # dots: never with color
    if (!is.null(dotsizes)) {
      color <- "NA"
    } else {
      # tiles: it depends
      if (nlevels(df[[features]]) > 100) {
        color <- "NA"
      } else {
        color <- "grey70"
      }
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
  if (axes_flip) {
    plot <- plot + ggplot2::coord_flip()
  }

  if (!is.null(dotsizes)) {
    plot <- plot + ggplot2::geom_point(
      ggplot2::aes(size = !!rlang::sym(dotsizes)),
      shape = 21,
      color = color,
      stroke = color_linewidth
    ) + ggplot2::scale_size(range = dotsize_range)
  } else {
    plot <- plot + ggplot2::geom_tile(
      color = color,
      linewidth = color_linewidth
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
                                     colorsteps_nice = colorsteps_nice)
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

  if (group_seplines) {
    dfsort <- df |>
      dplyr::group_by(!!rlang::sym(features)) |>
      dplyr::slice_max(order_by = !!rlang::sym(values), n = 1) |>
      dplyr::arrange(!!rlang::sym(groups))
    hlines <- cumsum(rle(as.character(dfsort[[groups]]))[["lengths"]]) + 0.5
    hlines <- hlines[-length(hlines)]
    plot <- plot + Gmisc::fastDoCall(ggplot2::geom_hline, args = c(list(yintercept = hlines), seplines_args))
  }


  if (!is.null(featurelabels) && featurelabels[1] == "auto") {
    if (nrow(df) > 100) {
      featurelabels <- ""
      names(featurelabels) <- featurelabels
      theme_args[[if (determine_feature_axis(plot) == "x") "axis.ticks.x" else "axis.ticks.y"]] <- ggplot2::element_blank()
    } else {
      featurelabels <- stats::setNames(as.character(df[[features]]), as.character(df[[features]]))
    }
  } else {
    if (is.null(featurelabels)) {
      featurelabels <- stats::setNames(as.character(df[[features]]), as.character(df[[features]]))
    } else if (is.null(names(featurelabels))) {
      # also works for featurelabels = ""
      names(featurelabels) <- featurelabels
    } else if (any(names(featurelabels) == "") && length(featurelabels) > 1) {
      # replace missing names in case only some features got alt labels
      names(featurelabels)[which(names(featurelabels) == "")] <- featurelabels[which(names(featurelabels) == "")]
    }
  }

  plot <- plot + theme + Gmisc::fastDoCall(ggplot2::theme, args = theme_args)

  # breaks not present are ignored
  # axes_flip is incorporated
  if (featurelabels_repel) {
    plot <- repel_features(
      df = df,
      plot = plot,
      repel_args = repel_args,
      featurelabels = featurelabels,
      featuresitalic = featuresitalic)
  } else {
    plot <- plot + ggplot2::scale_y_discrete(breaks = featurelabels, labels = names(featurelabels))
  }


  return(plot)
}




colorscale_heuristic <- function(colorscale_values,
                                 values_zscored,
                                 colorsteps,
                                 legendbreaks,
                                 legendlabels,
                                 fill,
                                 colorsteps_nice,
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
                 nice.breaks = colorsteps_nice,
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


repel_features <- function(df, plot, repel_args, featurelabels, featuresitalic) {

  axis.df <- data.frame(
    y = 1:length(levels(df$feature)),
    feature = levels(df$feature))
  axis.df$label = stats::setNames(names(featurelabels), featurelabels)[axis.df$feature]
  axis <- ggplot2::ggplot(axis.df, ggplot2::aes(x = 0, y = y, label = feature)) +
    ggrepel::geom_text_repel(
      fontface = ifelse(featuresitalic, "italic", "plain"),
      data = axis.df[which(axis.df$feature %in% featurelabels),],
      ggplot2::aes(label = label),
      nudge_x = repel_args[["featurelabels_nudge_x"]],
      min.segment.length = 0,
      direction = "y"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-0.1, 0),
      expand = c(0, 0),
      breaks = NULL,
      labels = NULL,
      name = NULL
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, length(levels(df$feature)) + 0.5),
      expand = c(0, 0),
      breaks = NULL,
      labels = NULL,
      name = NULL
    ) +
    ggplot2::theme_void()
  plot <- plot + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  ) + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")) +
    ggplot2::scale_y_discrete(name = NULL)
  plot <- cowplot::plot_grid(
    axis,
    plot,
    align = "h",
    axis = "tb",
    nrow = 1,
    rel_widths = c(repel_args[["featurelabels_width"]],1)
  )
  return(plot)
}

determine_feature_axis <- function(plot) {
  if (length(unique(plot[["data"]][[rlang::quo_get_expr(plot[["mapping"]][["x"]])]])) > 100) {
    return("x")
  } else {
    return("y")
  }
}

#' Plot a heatmap from long-format data
#'
#' `heatmap_long_df()` plots one observation for each group-feature combination.
#' By default, observations are drawn as tiles whose fill represents `values`.
#' Supplying `dotsizes` switches to a dot heatmap in which fill represents
#' `values` and point size represents a second variable.
#'
#' @param df A data frame in long format. Each row should describe a
#'   group-feature combination. The columns named by `groups`, `features`, and
#'   `values` must be present; the value column must be numeric.
#' @param groups A single character string naming the column that defines the
#'   heatmap groups. If omitted, the function attempts to infer it from the
#'   first two character or factor columns in `df`.
#' @param features A single character string naming the feature column. If
#'   omitted, the function attempts to infer it from the first two character or
#'   factor columns in `df`, choosing the column with more unique values when
#'   both `groups` and `features` are omitted.
#' @param values A single character string naming the numeric column mapped to
#'   fill. If omitted, the first numeric column in `df` is used.
#' @param featuregroup Optional single character string naming a column that
#'   assigns features to higher-level feature groups.
#' @param featuregroup_style One or both of `"facet"` and `"color"`. Faceting
#'   separates feature groups into panels; coloring renders feature-axis labels
#'   in group-specific colors.
#' @param featuregroup_col_name Character string used as the feature-group
#'   color legend title.
#' @param featuregroup_col_pal Palette name passed to [colrr::col_pal()] for
#'   feature-group label colors.
#' @param dotsizes Optional single character string naming a numeric column to
#'   map to point size. When supplied, points are drawn instead of tiles.
#' @param dotsize_range Numeric vector of length two giving the minimum and
#'   maximum plotted point sizes.
#' @param fill A vector of colors or a palette name understood by
#'   [colrr::col_pal()]. The sentinel `"..auto.."` uses a reversed 11-color
#'   `"RdBu"` palette.
#' @param color Border color for tiles or points. Use `"NA"` for no border.
#'   With `"..auto.."`, dot borders are omitted; tile borders are `"grey70"`
#'   for at most 100 features and omitted otherwise.
#' @param scale Value transformation applied separately within each feature:
#'   `"none"` leaves values unchanged, `"zscore"` standardizes them, and
#'   `"1"` rescales them to `scale_range`.
#' @param scale_range Numeric vector of length two giving the output range used
#'   when `scale = "1"`.
#' @param features_topn Optional positive integer. If supplied, retain the top
#'   features selected per group before plotting.
#' @param topn_cols Character vector naming columns used, in order, to rank
#'   features for `features_topn` via [dplyr::slice_max()]. Transform columns
#'   first when smaller values should rank higher (for example, negate
#'   p-values).
#' @param topn_ties Logical; whether ties may cause more than `features_topn`
#'   features to be retained.
#' @param featurelabels Controls feature-axis labels. `NULL` labels all
#'   features, `""` labels none, and `"..auto.."` omits labels when there are
#'   more than 200 features. A character vector selects labels; a named vector
#'   uses its names as display labels and its values as feature names, for
#'   example `c("CD20" = "MS4A1", "CD3", "KLRG1")`.
#' @param featurelabels_repel Logical; whether to draw feature labels in a
#'   separate repelled-label panel using [ggrepel::geom_text_repel()].
#' @param featuresitalic Logical; whether to render feature labels in italics.
#' @param color_linewidth Numeric border width for tiles or points.
#' @param legendbreaks A numeric vector of fill-scale breaks, `"..auto.."` for
#'   automatic breaks, or `"minmidmax"` for breaks at the minimum, midpoint,
#'   and maximum of the value range.
#' @param legendlabels A character vector of labels corresponding to
#'   `legendbreaks`, or `"..auto.."` for automatic labels.
#' @param colorsteps Controls discretization of the fill guide. Use `NULL` for
#'   a continuous color bar, `"..auto.."` for automatic steps, a single number
#'   for the requested number of steps, or a numeric vector of explicit step
#'   boundaries.
#' @param colorsteps_nice Logical; whether to adjust color-step boundaries to
#'   visually convenient values. Some requested step counts may be adjusted.
#' @param color_trans_log Logical; whether to use a logarithmic transformation
#'   for the fill scale.
#' @param color_center_zero Logical; whether to center the fill scale at zero.
#' @param axes_flip Logical; whether to exchange the group and feature axes
#'   with [ggplot2::coord_flip()].
#' @param group_seplines Logical; whether to draw lines between runs of features
#'   assigned to different groups by their maximum value.
#' @param seplines_args Named list of additional arguments passed to
#'   [ggplot2::geom_hline()] when `group_seplines = TRUE`.
#' @param theme A complete or partial ggplot2 theme added to the plot.
#' @param legend_fill_args Named list of arguments passed to
#'   [ggplot2::guide_colorsteps()] or [ggplot2::guide_colorbar()], depending on
#'   the selected fill scale.
#' @param legend_size_args Named list of arguments passed to
#'   [ggplot2::guide_legend()] for the dot-size legend. For example,
#'   `override.aes = list(size = c(1, 3, 5))` customizes the sizes shown in the
#'   legend independently of `dotsize_range`.
#' @param theme_args Named list of arguments passed to [ggplot2::theme()] after
#'   `theme` is added.
#' @param repel_args Named list controlling the repelled feature-label panel.
#'   Supported entries include `featurelabels_width` and
#'   `featurelabels_nudge_x`.
#' @param heatmap_ordering_args Named list of arguments passed to
#'   [heatmap_ordering()], such as `feature_order` and `group_order`.
#' @param values_zscored Logical indicating whether `values` are already
#'   z-scored. If `NULL`, this is inferred from the plotted matrix with
#'   [brathering::is_z_scored()].
#' @param pvals Optional single character string naming a p-value column. When
#'   supplied, significant cells are annotated with symbols.
#' @param pval_features Optional vector of feature values eligible for p-value
#'   annotation. `NULL` makes all plotted features eligible.
#' @param pval_max Numeric significance threshold; only p-values less than or
#'   equal to this value are annotated.
#' @param pval_symnum_args Named list passed to [stats::symnum()] to convert
#'   p-values to annotation symbols.
#' @param pval_filter Rule used to choose cells for p-value annotation. `"top"`
#'   selects the highest-value cell for each feature; `"pos_fc"` selects rows
#'   with a positive value in the column named by `pval_logfc`.
#' @param pval_logfc Single character string naming the fold-change column used
#'   when `pval_filter = "pos_fc"`.
#' @param pval_text_args Named list of additional arguments passed to
#'   [ggplot2::geom_text()] for p-value symbols.
#' @param impute_missing_to Optional scalar used to replace missing `values`
#'   after completing all group-feature combinations. If `NULL`, rows with
#'   missing values are removed.
#' @param lower_tri Logical; whether to retain only the lower triangle of the
#'   group-by-feature value matrix. This is primarily useful when the groups
#'   and features describe the same entities.
#' @param col_na Color used for missing values. With `"..auto.."`, the plot
#'   background color is used.
#' @param ... Additional arguments intended for [heatmap_ordering()], including
#'   `feature_order` and `group_order`.
#'
#' @details
#' The function first optionally selects top features, completes the grid of
#' group-feature combinations, handles missing values, and scales values within
#' each feature. It then orders the axes with [heatmap_ordering()] and constructs
#' either a tile or dot heatmap. Optional layers add significance symbols,
#' feature-group facets or label colors, and separation lines.
#'
#' `groups`, `features`, `values`, `dotsizes`, `featuregroup`, `pvals`, and
#' `pval_logfc` use column names supplied as character strings rather than tidy
#' evaluation expressions.
#'
#' @return A ggplot2 plot object. When `featurelabels_repel = TRUE`, a cowplot
#'   object containing the label panel and heatmap is returned.
#'
#' @importFrom rlang :=
#'
#' @export
#'
#' @examples
#' df <- readRDS(system.file("extdata", "heatmap_df.rds", package = "fcexpr"))
#'
#' # Tile heatmap with default styling
#' heatmap_long_df(
#'   df = df,
#'   groups = "cluster",
#'   features = "channel",
#'   values = "mean_cluster_scale"
#' )
#'
#' # Dot heatmap: fill shows scaled mean and size shows -log10(p-value)
#' heatmap_long_df(
#'   df = df,
#'   groups = "cluster",
#'   features = "channel",
#'   values = "mean_cluster_scale",
#'   dotsizes = "pvalue2"
#' )
#'
#' # Show four top features per group, omit labels, flip axes, and use a
#' # continuous color bar
#' heatmap_long_df(
#'   df = df,
#'   groups = "cluster",
#'   features = "channel",
#'   values = "mean_cluster_scale",
#'   dotsizes = "pvalue2",
#'   featurelabels = "",
#'   axes_flip = TRUE,
#'   features_topn = 4,
#'   group_seplines = TRUE,
#'   colorsteps = NULL
#' )
#'
#' # Scale within features and label the fill legend at its range endpoints
#' # and midpoint
#' heatmap_long_df(
#'   df = df,
#'   groups = "cluster",
#'   features = "channel",
#'   values = "mean_cluster",
#'   scale = "zscore",
#'   colorsteps = NULL,
#'   legendbreaks = "minmidmax",
#'   legendlabels = c("min", "mid", "max")
#' )
heatmap_long_df <- function(df,
                            groups,
                            features,
                            values,
                            featuregroup = NULL,
                            featuregroup_style = c("facet", "color"),
                            featuregroup_col_name = "",
                            featuregroup_col_pal = "custom_light",
                            dotsizes = NULL,
                            dotsize_range = c(2,7),
                            fill = "..auto..",
                            color = "..auto..",
                            scale = c("none", "zscore", "1"),
                            scale_range = c(-1,1),
                            features_topn = NULL,
                            topn_cols = values,
                            topn_ties = F,
                            featurelabels = "..auto..",
                            featurelabels_repel = F,
                            featuresitalic = F,
                            color_linewidth = 0.2,
                            legendbreaks = "..auto..",
                            legendlabels = "..auto..",
                            colorsteps = "..auto..",
                            colorsteps_nice = T,
                            color_trans_log = F,
                            color_center_zero = T,
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
                              override.aes = list(color = "..auto..", fill = "..auto..")
                              #size = c(2, 7))
                            ),
                            theme_args = list(
                              panel.grid = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
                            ),
                            repel_args = list(featurelabels_width = 0.2,
                                              featurelabels_nudhe_x = -1),
                            heatmap_ordering_args = list(feature_order = "custom",
                                                         group_order = "hclust"),
                            values_zscored = NULL,
                            pvals = NULL,
                            pval_features = NULL,
                            pval_max = 0.01,
                            pval_symnum_args = list(cutpoints = c(0, 0.0001, Inf),
                                                    symbols = c("*", "ns")),
                            pval_filter = c("top", "pos_fc"),
                            pval_logfc = "logFC",
                            pval_text_args = list(size = 5, vjust = 0.75),
                            impute_missing_to = NULL,
                            lower_tri = F,
                            col_na = "..auto..",
                            ...) {
  fcexpr:::.ensure_packages(c("brathering", "colrr", "ggplot2", "ggtext", "Gmisc", "scales"))



  stopifnot("df must be a data frame" = is.data.frame(df))

  dots <- list(...)

  if ("feature_order" %in% names(dots)) {
    heatmap_ordering_args[[feature_order]] <- dots[[feature_order]]
  }
  if ("group_order" %in% names(dots)) {
    heatmap_ordering_args[[group_order]] <- dots[[group_order]]
  }

  if (missing(values)) {
    # first numeric column becomes values
    values <- names(which(sapply(df, is.numeric)))[1]
    message("values: ", values)
  }

  if (!is.null(featuregroup) && !featuregroup %in% names(df)) {
    stop("featuregroup column not found in df.")
  }
  featuregroup_style <- rlang::arg_match(featuregroup_style, multiple = T)

  if (missing(features) || missing(groups)) {
    # first two char columns
    featgrou <- names(c(which(sapply(df, is.character)), which(sapply(df, is.factor))))[1:2]
    if (missing(features) && missing(groups)) {
      lens <- c(length(unique(df[[featgrou[1]]])), length(unique(df[[featgrou[2]]])))
      # features: more levels
      features <- featgrou[which.max(lens)[1]]
      groups <- setdiff(featgrou, features)
      message("features: ", features)
      message("groups: ", groups)
    } else if (missing(features)) {
      features <- setdiff(featgrou, groups)
      message("features: ", features)
    } else if (missing(groups)) {
      groups <- setdiff(featgrou, features)
      message("features: ", groups)
    }
  }

  scale <- rlang::arg_match(scale)
  pval_filter <- rlang::arg_match(pval_filter)

  if (featuresitalic) {
    theme_args <- brathering::gg_inject_theme_element(theme_args = theme_args,
                                                      elem = if (axes_flip) "axis.text.x" else "axis.text.y",
                                                      elem_sub = "face",
                                                      value = "italic")
  }

  if (!is.null(pvals) && !pvals %in% names(df)) {
    message("pvals column not found.")
    pvals <- NULL
  }
  if (!is.null(pvals) && pval_filter == "pos_fc" && !pval_logfc %in% names(df)) {
    message("pval_logfc column not found.")
    pvals <- NULL
  }



  if (!is.null(pvals)) {
    if (is.null(pval_features)) {
      pval_features <- unique(df[[features]])
    } else {
      pval_features <- intersect(pval_features, df[[features]])
      if (!length(pval_features)) {
        pvals <- NULL
        message("none of pval_features found.")
      }
    }
  }

  # optional filter for top n features per group
  if (!is.null(features_topn)) {
    select <- df |>
      # max group per feature
      dplyr::slice_max(order_by = !!rlang::sym(values), n = 1, by = !!rlang::sym(features)) |>
      # then best features per group
      # dplyr::slice_max(order_by = !!rlang::sym(values), n = features_topn, by = !!rlang::sym(groups)) |>
      # dplyr::slice_max(order_by = tibble::tibble(!!rlang::sym("auc"), !!rlang::sym("padj"), !!rlang::sym("logFC")), n = features_topn, by = !!rlang::sym(groups)) |>
      dplyr::slice_max(order_by = tibble::tibble(!!!rlang::syms(as.list(topn_cols))), n = features_topn, by = !!rlang::sym(groups), with_ties = topn_ties) |>
      dplyr::pull(!!rlang::sym(features))
    df <- df[which(df[[features]] %in% select),,drop = F]
  }


  df <- tidyr::complete(df, !!rlang::sym(groups), !!rlang::sym(features))
  if (anyNA(df[[values]])) {
    if (!is.null(impute_missing_to)) {
      message("missing values imputed to ", impute_missing_to)
      df[[values]][which(is.na(df[[values]]))] <- impute_missing_to
    } else {
      message("missing values found!")
      df <- dplyr::filter(df, !is.na(!!rlang::sym(values)))
    }
  }

  # optional scaling
  if (scale != "none") {
    df <- df |>
      #tidyr::complete(!!rlang::sym(groups), !!rlang::sym(features)) |>
      #dplyr::mutate(!!values := ifelse(is.na(!!rlang::sym(values)), 0, !!rlang::sym(values)))
      dplyr::mutate(!!values := dplyr::case_when(
        scale == "zscore" ~ as.vector(scale(!!rlang::sym(values))),
        scale == "1" ~ scales::rescale(!!rlang::sym(values), to = scale_range),
        .default = !!rlang::sym(values)  # fallback (optional)
      ), .by = !!rlang::sym(features))
  }


  # assign factors to features and groups
  df <- Gmisc::fastDoCall(fcexpr::heatmap_ordering,
                          args = c(list(df = df,
                                        features = features,
                                        groups = groups,
                                        values = values),
                                   heatmap_ordering_args))

  if (color[1] == "..auto..") { # catch if length(color) > 1
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


  if (length(fill) == 1 && fill == "..auto..") {
    fill <- colrr::col_pal(name = "RdBu", n = 11, direction = -1)
  } else if (length(fill) == 1) {
    fill <- colrr::col_pal(name = fill)
  }

  df <- dplyr::arrange(df, !!rlang::sym(features))

  dfmat <- brathering::df_long_to_mat(
    df,
    to_rows = groups,
    to_cols = features,
    values = values)

  if (lower_tri) {
    # if (ncol(corr_mat) != nrow(corr_mat)) {
    #   message("Correlation matrix is not quadratic. Returning the lower triangle may not yield intended results.")
    # }
    df_lowtri <- dfmat
    df_lowtri[which(!lower.tri(df_lowtri))] <- NA
    df_lowtri <- brathering::mat_to_df_long(df_lowtri,
                                            rownames_to = groups,
                                            colnames_to = features,
                                            values_to = values) |>
      dplyr::left_join(df[,-which(names(df) == values)], by = c(groups, features))
    df_lowtri[[groups]] <- factor(df_lowtri[[groups]], levels = levels(df[[groups]]))
    df_lowtri[[features]] <- factor(df_lowtri[[features]], levels = levels(df[[features]]))
    df <- df_lowtri
  }

  # start ggplot pipeline
  plot <- ggplot2::ggplot(df, ggplot2::aes(
    x = !!rlang::sym(groups),
    y = !!rlang::sym(features),
    fill = !!rlang::sym(values)))
  if (axes_flip) {
    plot <- plot + ggplot2::coord_flip()
  }

  if (!is.null(dotsizes)) {
    plot <- plot + ggplot2::geom_point(
      ggplot2::aes(size = !!rlang::sym(dotsizes)),
      shape = 21,
      color = color,
      stroke = color_linewidth) +
      ggplot2::scale_size(range = dotsize_range)
  } else {
    plot <- plot + ggplot2::geom_tile(
      color = color,
      linewidth = color_linewidth)
  }

  # check if values are z-scored
  if (is.null(values_zscored)) {
    values_zscored <- sum(apply(dfmat, 2, brathering::is_z_scored, verbose = F, tol = 0.05)) > 0.9*ncol(dfmat) # 0.9: arbitrary choice
    if (values_zscored) {
      message("values interpreted as z-scored.")
    }
  }


  if (!is.null(pvals)) {
    if (pval_filter == "max") {
      df_pval <- dplyr::slice_max(df, order_by = !!rlang::sym(values), n = 1, by = c(!!rlang::sym(features)))
    } else if (pval_filter == "pos_fc") {
      df_pval <- dplyr::filter(df, !!rlang::sym(pval_logfc) > 0)
    }
    df_pval <- df_pval |>
      dplyr::filter(!!rlang::sym(pvals) <= pval_max) |>
      dplyr::filter(!!rlang::sym(features) %in% pval_features)
    df_pval$pval_sym <- do.call(stats::symnum, c(list(x = df_pval[[pvals]]), pval_symnum_args))

    plot <- plot + Gmisc::fastDoCall(ggplot2::geom_text, args = c(list(data = df_pval,
                                                                       mapping = ggplot2::aes(label = pval_sym)),
                                                                  pval_text_args))
  }

  col_na <- col_na[1]
  if (col_na == "..auto..") {
    col_na <- brathering::gg_get_theme_element(plot, element = "plot.background")@fill
  }

  # decide for colorsteps or continuous colorbar
  scale_fill <- colrr::get_scale_fill_fun(values = df[[values]],
                                          zscored = values_zscored,
                                          steps = colorsteps,
                                          legendbreaks = legendbreaks,
                                          legendlabels = legendlabels,
                                          palette = fill,
                                          steps_nice = colorsteps_nice,
                                          trans_log = color_trans_log,
                                          col_na = col_na,
                                          center_zero = color_center_zero)

  if (grepl("coloursteps", scale_fill[["guide"]])) {
    guide_fun <- ggplot2::guide_colorsteps
  } else {
    guide_fun <- ggplot2::guide_colorbar
  }


  plot <- plot + theme + Gmisc::fastDoCall(ggplot2::theme, args = theme_args)
  ## auto legend size fill and color ?
  if ("override.aes" %in% names(legend_size_args)) {
    if ("color" %in% names(legend_size_args[["override.aes"]])) {
      if (legend_size_args[["override.aes"]][["color"]] == "..auto..") {
        legend_size_args[["override.aes"]][["color"]] <- brathering::bw_txt(brathering::gg_get_theme_element(plot, "plot.background")@fill)
      }
    }
    if ("fill" %in% names(legend_size_args[["override.aes"]])) {
      if (legend_size_args[["override.aes"]][["fill"]] == "..auto..") {
        legend_size_args[["override.aes"]][["fill"]] <- brathering::bw_txt(brathering::gg_get_theme_element(plot, "plot.background")@fill)
      }
    }
  }

  plot <- plot +
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


  if (!is.null(featurelabels) && featurelabels[1] == "..auto..") {
    if (length(unique(df[[features]])) > 200) {
      message("featurelabels omitted as n>200. Set to NULL to plot all.")
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

  # plot <- plot + theme + Gmisc::fastDoCall(ggplot2::theme, args = theme_args)

  if (!is.null(featuregroup)) {
    if ("facet" %in% featuregroup_style) {
      if (axes_flip) {
        facet_args <- list(scales = "free_x",
                           space = "free_x",
                           nrow = 1,
                           ncol = NULL,
                           strip.position = "top")
      } else {
        facet_args <- list(scales = "free_y",
                           space = "free_y",
                           ncol = 1,
                           nrow = NULL,
                           strip.position = "right")
      }

      plot <- plot +
        Gmisc::fastDoCall(ggplot2::facet_wrap,
                          args = c(list(facets = ggplot2::vars(!!rlang::sym(featuregroup))),
                                   facet_args))
    }

    if ("color" %in% featuregroup_style) {
      marker_df <- dplyr::distinct(df, !!rlang::sym(features), !!rlang::sym(featuregroup))
      color_conv <- colrr::col_pal(featuregroup_col_pal, n = marker_df[[featuregroup]], return = "char")
      marker_df$color <- color_conv[marker_df[[featuregroup]]]
      colman <- stats::setNames(marker_df$color, marker_df[[featuregroup]])

      plot <- plot +
        ggplot2::geom_point(
          data = data.frame(
            x = NA,
            y = NA,
            yaxis = marker_df[[featuregroup]]
          ) |> dplyr::mutate(!!values := 0),
          ggplot2::aes(x = x, y = y, color = yaxis, fill = !!rlang::sym(values)),
          inherit.aes = F
        ) +
        ggplot2::scale_color_manual(
          name = featuregroup_col_name,
          values = colman[!duplicated(colman)]
        ) +
        ggplot2::scale_x_discrete(na.translate = F)
      if (axes_flip) {
        plot <- plot + ggplot2::theme(axis.text.x = ggtext::element_markdown())
      } else {
        plot <- plot + ggplot2::theme(axis.text.y = ggtext::element_markdown())
      }
    }
  }


  # breaks not present are ignored
  # axes_flip is incorporated
  if (featurelabels_repel) {
    if (!is.null(featuregroup)) {
      warning("grouped y axis and feature repel combination is not tested or established.")
    }
    plot <- repel_features(
      df = df,
      plot = plot,
      repel_args = repel_args,
      featurelabels = featurelabels,
      featuresitalic = featuresitalic)
  } else {
    if (!is.null(featuregroup) && "color" %in% featuregroup_style) {
      plot <- plot +
        ggplot2::scale_y_discrete(na.translate = F,
                                  labels = function(y) color_labels(y, stats::setNames(marker_df[["color"]],
                                                                                       marker_df[[features]])),
                                  breaks = featurelabels)
    } else {
      plot <- plot +
        ggplot2::scale_y_discrete(breaks = featurelabels, labels = names(featurelabels))
    }
  }

  return(plot)
}


repel_features <- function(df, plot, repel_args, featurelabels, featuresitalic) {
  fcexpr:::.ensure_packages(c("cowplot", "ggplot2", "ggrepel"))

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


round_auto_any <- function(x,
                           start_at = 100,               # start magnitude-based rounding at |x| >= this
                           method = c("nearest","up","down")) {
  method <- match.arg(method)
  ax <- abs(x)

  # base = 10^(floor(log10(|x|))) when |x| >= start_at; otherwise 1
  base <- ifelse(ax >= start_at, 10^floor(log10(ax)), 1)

  # helpers to do up/down as "away/toward zero"
  scale <- x / base
  up_scaled   <- ifelse(scale >= 0, ceiling(scale), floor(scale))   # away from zero
  down_scaled <- ifelse(scale >= 0, floor(scale), ceiling(scale))   # toward zero

  res <- switch(method,
                nearest = round(scale) * base,
                up      = up_scaled * base,
                down    = down_scaled * base
  )
  # keep NAs and zeros as-is
  res[is.na(x)] <- NA
  res
}

color_labels <- function(labels, col_map) {
  sapply(labels, function(x) {
    if (x %in% names(col_map)) {
      paste0("<span style='color:", col_map[x], ";'>", x, "</span>")
    } else {
      x
    }
  })
}

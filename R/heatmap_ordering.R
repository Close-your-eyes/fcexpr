#' Order columns and rows of a long data frame for an appealing heatmap
#'
#' Heatmaps only become meaningful when rows and columns are ordered properly.
#' If features or groups column of df are factors, then this remains unchanged
#' if order is 'none'. That way either order can be fixed as desired.
#' If both methods are set to 'hclust' this corresponds
#' to the procedure in pheatmap::pheatmap.'Custom' is a another heurisitc to
#' order features, preferentially. The default settings provide a good
#' basis, you may pre-set an order for columns by defining factors and leave
#' row_order 'custom'.
#'
#' @param df data frame in long format
#' @param features column name with features
#' @param groups column name with groups
#' @param values column name that becomes values in heatmap
#' @param feature_order how to order features
#' @param group_order how to order groups
#'
#' @return data frame with updated factor levels for rows and cols
#' @export
#'
#' @examples
heatmap_ordering <- function(df,
                             features = "channel",
                             groups = "cluster",
                             values = "mean_cluster_scale",
                             feature_order = c("custom", "hclust", "none"),
                             group_order = c("hclust", "none", "custom")) {
  #browser()

  cols <- groups
  rows <- features
  row_order_method <- rlang::arg_match(feature_order)
  col_order_method <- rlang::arg_match(group_order)

  if (group_order == "none" && feature_order == "none") {
    return(df)
  }

  mat <- df_long_to_mat(df = df,
                        to_rows = rows,
                        to_cols = cols,
                        values = values)

  clust_rows <- cluster_mat(mat = mat)
  clust_cols <- cluster_mat(mat = t(mat))


  # the default
  # if rows or cols are factors, keep that order
  # otherwise unique
  if (is.factor(df[[rows]])) {
    row_order <- levels(df[[rows]])
  } else {
    row_order <- unique(df[[rows]])
  }
  if (is.factor(df[[cols]])) {
    col_order <- levels(df[[cols]])
  } else {
    col_order <- unique(df[[cols]])
  }
  # pass orders already here to have them in df for potential order_custom_for_heatmap
  df[[rows]] <- factor(df[[rows]], levels = row_order)
  df[[cols]] <- factor(df[[cols]], levels = col_order)

  for (i in 1:2) {
    # i not used, just run this 2 times
    # in case column order has changed, re-ordering of rows must react to this
    if (row_order_method == "hclust") {
      row_order <- rownames(mat)[clust_rows$order]
    } else if (row_order_method == "custom") {
      row_order <- order_custom_for_heatmap(df = df,
                                            col_to_order = rows,
                                            order_metric = values,
                                            col_secondary = cols,
                                            secondary_order = col_order)
    }
    df[[rows]] <- factor(df[[rows]], levels = row_order)

    if (col_order_method == "hclust") {
      col_order <- colnames(mat)[clust_cols$order]
    } else if (col_order_method == "custom") {
      col_order <- order_custom_for_heatmap(df = df,
                                            col_to_order = cols,
                                            order_metric = values,
                                            col_secondary = rows,
                                            secondary_order = row_order)
    }
    df[[cols]] <- factor(df[[cols]], levels = col_order)
  }


  return(df)
}


cluster_mat = function(mat,
                       distance = c("euclidean", "correlation", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                       method = c("complete", "ward.D", "ward.D2", "ward", "single", "average", "mcquitty", "median", "centroid")) {
  method <- rlang::arg_match(method)
  distance <- rlang::arg_match(distance)

  if (distance == "correlation"){
    d = stats::as.dist(1 - cor(t(mat)))
  } else{
    d = stats::dist(mat, method = distance)
  }

  return(stats::hclust(d, method = method))
}

df_long_to_mat <- function(df, to_rows, to_cols, values) {

  if (any(dplyr::count(df, !!rlang::sym(to_rows), !!rlang::sym(to_cols))[["n"]] > 1)) {
    cat(paste0("dplyr::add_count(df, ", to_rows, ", ", to_cols, ")"))
    stop("pairs of ", to_rows, " and ", to_cols, " are not unique in your data frame. fix that. inspect n after running this command.")
  }


  mat <-
    df |>
    dplyr::select(dplyr::all_of(c(to_rows, to_cols, values))) |>
    tidyr::pivot_wider(names_from = !!rlang::sym(to_cols), values_from = !!rlang::sym(values)) |>
    tibble::column_to_rownames(to_rows) |>
    as.matrix()
  return(mat)
}

order_custom_for_heatmap <- function(df,
                                     col_to_order,
                                     order_metric,
                                     col_secondary,
                                     secondary_order = NULL) {

  # this is from scexpr::heatmap_pseudobulk

  if (is.null(secondary_order)) {
    if (is.factor(df[[col_secondary]])) {
      secondary_order <- levels(secondary_order)
    } else {
      secondary_order <- unique(df[[col_secondary]])
    }
  }

  df <-
    df %>%
    dplyr::group_by(!!rlang::sym(col_to_order)) %>%
    dplyr::slice_max(order_by = !!rlang::sym(order_metric), n = 1, with_ties = F) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(!!col_secondary := factor(!!rlang::sym(col_secondary), levels = unlist(secondary_order))) %>%
    dplyr::arrange(!!rlang::sym(col_secondary), !!rlang::sym(order_metric))

  return(df[[col_to_order]])
}

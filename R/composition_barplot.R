.composition_barplot <- function(x,
                                 x_cat,
                                 fill_cat,
                                 col_pal = fcexpr::col_pal("custom"),
                                 plot_labels = F) {


  if (!x_cat %in% names(x)) {
    stop("x_cat not found in x.")
  }
  if (!fill_cat %in% names(x)) {
    stop("fill_cat not found in x.")
  }

  table <-
    x %>%
    dplyr::count(!!rlang::sym(x_cat), !!rlang::sym(fill_cat)) %>%
    dplyr::left_join(dplyr::count(x, !!rlang::sym(x_cat), name = "total"), by = x_cat) %>%
    dplyr::mutate(rel = n/total) %>%
    tibble::as_tibble()

  plot <- ggplot2::ggplot(table, ggplot2::aes(x = !!rlang::sym(x_cat), y = rel, fill = !!rlang::sym(fill_cat))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = col_pal)

  if (plot_labels) {
    plot <- plot + ggplot2::geom_text(ggplot2::aes(label = round(rel,2)), position = ggplot2::position_stack(vjust = 0.5))
  }

  return(list(table = table, plot = plot))
}

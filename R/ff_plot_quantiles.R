#' Plot quantiles and their values channel-wise and flowframe-wise
#'
#' @param ff one flowframe or a list of such
#' @param channels channels to plot, goes into ff_get_channels
#' @param quantiles which quantils to calculate
#'
#' @returns list
#' @export
#'
#' @examples
ff_plot_quantiles <- function(ff,
                              channels = NULL,
                              quantiles = seq(0,1,0.01)) {

  if (!is.list(ff)) {
    ff <- list(ff)
  }



  qdf <- purrr::map(ff, function(x) {
    chann <- ff_get_channels(x, channels = channels)
    qu <- apply(
      flowCore::exprs(x)[,names(chann)],
      2,
      quantile,
      probs = quantiles,
      simplify = T
    )
    df <- brathering::mat_to_df_long(
      qu,
      rownames_to = "quantile",
      colnames_to = "channel"
    )
    df$quantile <- as.numeric(gsub("%", "", df$quantile))
    return(df)
  })

  plots <- purrr::map(qdf, function(df) {
    ggplot2::ggplot(df, ggplot2::aes(x=quantile, y=value)) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(ggplot2::vars(channel))
  })

  return(list(df = qdf, plot = plots))
}

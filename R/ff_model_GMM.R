#' Model expression data from flow frame with Gaussian mixture model (GMM)
#'
#' Functions from mclust package are used for modelling: mclust::mclustBIC and
#' mclust::mclustModel. GMMs_BIC is the result from mclust::mclustBIC showing
#' goodness of different model types. Model from the return list can be used by
#' fcexpr::ff_simulate to generate a flow frame with modelled exprs data.
#'
#' @param ff flow frame
#' @param channels channel selection passed to fcexpr::ff_get_channels
#' @param sample_n downsample evts to that number, set NULL for no sampling
#' @param leverage_sampling use leverage scores for sampling, see
#' brathering::sample_leverage
#' @param scale_for_sampling scale (z-score) before sampling (only relevant if
#' leverage_sampling)
#' @param scale_for_modelling scale (z-score) before modelling
#' @param mclustBIC_args arugments to mclust::mclustBIC
#' @param seed seed for reproducing random results
#' @param source_file name of ff to be added to model list, if $FIL then the
#' keyword from ff is used
#' @param ... arguments to fcexpr::ff_get_channels
#'
#' @return list of model data
#' @export
#'
#' @examples
ff_model_GMM <- function(ff,
                         channels = NULL,
                         sample_n = NULL,
                         leverage_sampling = T,
                         scale_for_sampling = T,
                         scale_for_modelling = T,
                         mclustBIC_args = list(model_names = NULL, G = NULL),
                         keywords = c("$BEGINANALYSIS", "$BEGINDATA", "$BEGINSTEXT",
                                      "$BYTEORD", "$DATATYPE", "$ENDANALYSIS", "$ENDDATA",
                                      "$ENDSTEXT", "$MODE", "$NEXTDATA", "$PAR", "$TOT",
                                      "$CYT", "$DATE", "$FIL", "$SRC", "$INST", "$OP",
                                      "$SYS", "$TIMESTEP", "$UNICODE", "$BTIM", "$ETIM",
                                      "SPILL", "SPILLOVER", "FCSversion"),
                         seed = 42,
                         source_file = "$FIL",
                         ...) {

  if (!requireNamespace("brathering", quietly = T)) {
    devtools::install_github("Close-your-eyes/brathering")
  }

  channels <- ff_get_channels(ff, channels = channels, rm_scatter = F, return = "data.frame", ...)

  if (source_file == "$FIL") {
    source_file <- flowCore::keyword(ff)[["$FIL"]]
  }

  # sample evts for modelling
  # use leverage: do not care about skewed population sizes in simulation
  inds <- do.call(brathering::sample_leverage,
                  args = c(list(x = scale(flowCore::exprs(ff)[,channels],
                                          scale = scale_for_sampling,
                                          center = scale_for_sampling),
                                leverage = leverage_sampling,
                                size = sample_n,
                                seed = seed)))
  # scale data before modelling, then rescale simulated results
  exprs <- scale(flowCore::exprs(ff)[,channels],
                 scale = scale_for_modelling,
                 center = scale_for_modelling)

  # model the data
  # fit different gaussian mix model to describe the data; BIC is used to find the best fit
  GMMs_BIC <- do.call(mclust::mclustBIC, args = c(list(data = exprs), mclustBIC_args))
  # select best model
  exprs_model <- mclust::mclustModel(data = exprs, BICvalues = GMMs_BIC)

  ## add info for rescaling after simulation
  # if attr does not exist: NULL
  exprs_model[["scale_center"]] <- attr(exprs, "scaled:center")
  exprs_model[["scale_scale"]] <- attr(exprs, "scaled:scale")

  ## add to exprs_model
  new_pars <- flowCore::parameters(ff)
  new_pars@data <- new_pars@data[which(new_pars@data$name %in% channels),]
  exprs_model[["new_pars"]] <- new_pars

  new_kw <- flowCore::keyword(ff)
  new_kw <- new_kw[keywords[which(keywords %in% names(new_kw))]]
  new_kw[c("$BTIM", "$ETIM", "$DATE", "$FIL", "$OP", "$TOT")]
  exprs_model[["new_kw"]] <- new_kw

  exprs_model[["source_file"]] <- source_file

  return(list(GMMs_BIC = GMMs_BIC, model = exprs_model))
}

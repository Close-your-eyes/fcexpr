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
#' @param mclustBIC_args arguments to mclust::mclustBIC
#' @param seed seed for reproducing random results
#' @param source_file name of ff to be added to model list, if $FIL then the
#' keyword from ff is used
#' @param ... arguments to fcexpr::ff_get_channels
#' @param keywords keywords to retain from ff
#'
#' @return list of model data
#' @export
#'
#' @examples
#' \dontrun{
#' library(fcexpr)
#' library(mclust)
#' ws <- "path/to/your/flowjo.wsp"
#'
#' # check samples, groups, populations
#' groups <- fcexpr::wsx_get_groups(ws)
#' paths <- fcexpr::wsx_get_fcs_paths(ws, filter_AllSamples = T)
#' samplenames <- basename(grep("AnKr", paths$Exp_part_1, value = T))
#' pops <- fcexpr::wsx_get_poppaths(ws, groups = "Exp_part_1")
#'
#' # read flow frames
#' fflist <- fcexpr::wsp_get_ff(wsp = ws,
#'                              FCS.file.folder = "path/to/your/FCS_files",
#'                              population = "Single Cells",
#'                              samples = samplenames[1])
#' ff <- purrr::list_flatten(fflist[["flowframes"]], name_spec = "{outer}")
#' ff <- purrr::list_flatten(ff, name_spec = "{outer}")
#'
#' # create model, define MmdelNames and/or G to speed modelling up
#' # mclustBIC_args <- list()
#' # mclustBIC_args <- list(modelNames = "VVV", G = 8:10)
#' modellist <- ff_model_GMM(ff = ff[[1]],
#'                           sample_n = 50000,
#'                           mclustBIC_args = mclustBIC_args,
#'                           source_file = names(ff)[1])
#' model <- modellist[[2]]
#' # save model
#' saveRDS(model, file.path("your/path", paste0(names(ff)[2], "_mclustmodel.rds")))
#'
#' flowframe <- ff_simulate(model = model)
#' flowCore::write.FCS(flowframe[[1]], "fcs/save/path.fcs")
#' }
ff_model_GMM <- function(ff,
                         channels = NULL,
                         sample_n = NULL,
                         leverage_sampling = T,
                         scale_for_sampling = T,
                         scale_for_modelling = T,
                         mclustBIC_args = list(modelNames = NULL, G = NULL),
                         keywords = c("$BEGINANALYSIS", "$BEGINDATA", "$BEGINSTEXT",
                                      "$BYTEORD", "$DATATYPE", "$ENDANALYSIS", "$ENDDATA",
                                      "$ENDSTEXT", "$MODE", "$NEXTDATA", "$PAR", "$TOT",
                                      "$CYT", "$DATE", "$FIL", "$SRC", "$INST", "$OP",
                                      "$SYS", "$TIMESTEP", "$UNICODE", "$BTIM", "$ETIM",
                                      "SPILL", "SPILLOVER", "FCSversion", "GUID"),
                         seed = 42,
                         source_file = "$FIL",
                         ...) {

  stopifnot("ff has to be a flowframe" = methods::is(ff, "flowFrame"))

  if (!requireNamespace("brathering", quietly = T)) {
    pak::pak("Close-your-eyes/brathering")
  }
  if (!requireNamespace("flowCore", quietly = T)) {
    BiocManager::install("flowCore")
  }

  channels <- ff_get_channels(ff, channels = channels, rm_scatter = F, return = "vector", ...)

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
  # downsample and optionally scale data before modelling, rescale simulated results later
  exprs <- scale(flowCore::exprs(ff)[inds[[1]], channels],
                 scale = scale_for_modelling,
                 center = scale_for_modelling)

  dtach <- !"mclust" %in% .packages()
  # model the data
  library(mclust) # necessary for mclust functions
  # fit different gaussian mix model to describe the data; BIC is used to find the best fit
  GMMs_BIC <- do.call(mclust::mclustBIC, args = c(list(data = exprs), mclustBIC_args))
  # select best model
  exprs_model <- mclust::mclustModel(data = exprs, BICvalues = GMMs_BIC)
  if (dtach) {
    detach("package:mclust", unload = T)
  }

  exprs_model[["original_data"]] <- exprs

  ## add info for rescaling after simulation
  # if attr does not exist: NULL
  exprs_model[["scale_center"]] <- attr(exprs, "scaled:center")
  exprs_model[["scale_scale"]] <- attr(exprs, "scaled:scale")

  ## add to exprs_model
  params <- flowCore::parameters(ff)
  params@data <- params@data[which(params@data$name %in% channels),]
  exprs_model[["params"]] <- params

  keywrd <- flowCore::keyword(ff)
  keywrd <- keywrd[keywords[which(keywords %in% names(keywrd))]]
  keywrd[c("$BTIM", "$ETIM", "$DATE", "$FIL", "$OP", "$TOT")] <- ""
  exprs_model[["keywrd"]] <- keywrd

  exprs_model[["source_file"]] <- source_file
  exprs_model[["time_channel"]] <- ifelse(length(flowCore:::findTimeChannel(ff)), flowCore:::findTimeChannel(ff), "Time")

  return(list(GMMs_BIC = GMMs_BIC, model = exprs_model))
}

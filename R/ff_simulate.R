#' Generate random flow frame from Gaussian mixture model (GMM)
#'
#' mclust::sim is used to generate an expression matrix based on a GMM
#' generated with fcexpr::ff_model_GMM. The matrix is decorated with
#' parameters and keywords to form a flow frame.
#'
#' @param model GMM and accessory information from fcexpr::ff_model_GMM
#' @param path path to folder with multiple models on disk as .rds files to
#' choose randomly from
#' @param n number of events to simulate
#' @param m number of flow frames to simulate, if path is given every flow frame
#' uses a randomly picked model
#' @param seed seed for mclust::sim
#'
#' @return list of flow frame(s)
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
#'
#' flowframe <- ff_simulate(model = modellist[[2]])
#' flowCore::write.FCS(flowframe[[1]], "fcs/save/path.fcs")
#' }
ff_simulate <- function(model,
                        path = NULL,
                        n = 50000,
                        m = 1,
                        seed = 42) {
  # model or path to folder with models to choose
  if (!is.null(path)) {
    modelfiles <- list.files(path, "\\.rds$", full.names = T)
    stopifnot("No rds files found in path." = length(modelfiles) != 0)
    models <- purrr::map(modelfiles, readRDS)
  }

  ff_list <- purrr::map(1:m, function(x) {
    if (!is.null(path)) {
      model <- sample(models, 1)[[1]]
    }
    seed <- seed+x

    dtach <- !"mclust" %in% .packages()
    library(mclust)
    y <- mclust::sim(modelName = model$modelName,
                     parameters = model$parameters,
                     seed = seed,
                     n = n)[,-1] # rm group column
    if (dtach) {
      detach("package:mclust", unload = T)
    }
    colnames(y) <- model[["params"]]@data[["name"]]

    # rescale
    if (!is.null(model[["scale_scale"]])) {
      y <- sweep(y, 2, model[["scale_scale"]], "*")
    }
    if (!is.null(model[["scale_center"]])) {
      y <- sweep(y, 2, model[["scale_center"]], "+")
    }

    # after optional sweeping:
    # add a time channel, random but equally spaced
    # BED[4] is the diff of BTIM and ETIM, so analysis duration
    # first evt random between 0.2 and 0.9 sec
    BED <- random_BTIM_ETIM_DATE(seed = seed)
    y <- cbind(y, matrix(seq(runif(1,0.2,0.9), as.numeric(BED[4]), length.out = n), ncol = 1,
                         dimnames = list(NULL, model[["time_channel"]])))

    kw_par <- get_kw_and_pars(exprs = y,
                              keywrd = model[["keywrd"]],
                              params = model[["params"]])

    # generate a few random keywords
    kw_par[["keywrd"]][["$OP"]] <- random_OP(seed = seed)
    kw_par[["keywrd"]][["$FIL"]] <- random_FIL(seed = seed)
    kw_par[["keywrd"]][["$BTIM"]] <- BED[1]
    kw_par[["keywrd"]][["$ETIM"]] <- BED[2]
    kw_par[["keywrd"]][["$DATE"]] <- BED[3]
    kw_par[["keywrd"]][["source_file"]] <- model[["source_file"]]

    y <- methods::new(
      "flowFrame",
      exprs = y,
      parameters = kw_par[["params"]],
      description = kw_par[["keywrd"]]
    )

    return(y)
  })

  names(ff_list) <- purrr::map_chr(ff_list, ~paste0(flowCore::keyword(.x)[["$FIL"]], "__", flowCore::keyword(.x)[["source_file"]]))

  return(ff_list)
}


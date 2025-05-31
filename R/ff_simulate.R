#' Title
#'
#' @param model
#' @param path
#' @param n
#' @param m
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
ff_simulate <- function(model,
                        path,
                        n = 50000,
                        m = 1,
                        seed = 42) {
  # model or path to folder with models to choose
  if (!missing(path)) {
    modelfiles <- list.files(path, "\\.rds$", full.names = T)
    models <- purrr::map(modelfiles, readRDS)
  }

  ff_list <- purrr::map(1:m, function(x) {
    if (!missing(path)) {
      model <- sample(models, 1)[[1]]
    }

    y <- mclust::sim(modelName = model$modelName,
                     parameters = model$parameters,
                     seed = seed,
                     n = n)[,-1] # rm group column
    colnames(y) <- model[["new_pars"]]@data$name

    # rescale
    y <- sweep(y, 2, model[["scale_scale"]], "*")
    y <- sweep(y, 2, model[["scale_center"]], "+")

    kw_par <- get_new_kw_and_pars(exprs = y,
                                  new_kw = model[["new_kw"]],
                                  new_pars = model[["new_pars"]])

    # generate a few random keywords
    kw_par[["new_kw"]][["$OP"]] <- random_OP()
    kw_par[["new_kw"]][["$FIL"]] <- random_FIL()
    BED <- random_BTIM_ETIM_DATE()
    kw_par[["new_kw"]][["$BTIM"]] <- BED[1]
    kw_par[["new_kw"]][["$ETIM"]] <- BED[2]
    kw_par[["new_kw"]][["$DATE"]] <- BED[3]

    y <- methods::new(
      "flowFrame",
      exprs = y,
      parameters = kw_par[["new_pars"]],
      description = kw_par[["new_kw"]]
    )

    return(y)
  })

  return(ff_list)
}


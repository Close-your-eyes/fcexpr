#' Add new columns to a flow frame
#'
#' @param ff flow frame
#' @param mat matrix with numeric columns
#' @param overwrite do overwrite existing columns in ff?
#'
#' @return flow frame
#' @export
#'
#' @examples
ff_add_columns <- function(ff,
                           mat = NULL,
                           overwrite = F) {

  stopifnot("mat must be a matrix" = is.matrix(mat),
            "row numbers of ff and mat do not match." = nrow(flowCore::exprs(ff)) == nrow(mat),
            "mat has colnames that exist in ff. fix this or set overwrite = T." = !(any(colnames(mat) %in% colnames(flowCore::exprs(ff))) && !overwrite))

  # if (any(colnames(mat) %in% colnames(flowCore::exprs(ff))) && !overwrite) {
  #   stop("mat has colnames that exist in ff. fix this or set overwrite = T.")
  # }

  if (is.null(mat)) {
    return(ff)
  }

  # removing existing colnames from exprs(ff) equals overwrite

  exprs <- cbind(
    flowCore::exprs(ff)[,which(!colnames(flowCore::exprs(ff)) %in% colnames(mat)), drop = F],
    mat
  )
  #exprs <- cbind(flowCore::exprs(ff), mat)

  kw_par <- get_kw_and_pars(exprs = exprs,
                            ff = ff)

  # save previous attr
  ff_default_attr <- c("exprs", "parameters", "description", "class")
  prev_attr <- attributes(ff)[which(!names(attributes(ff)) %in% ff_default_attr)]

  # new ff
  ff <- methods::new(
    "flowFrame",
    exprs = exprs,
    parameters = kw_par[["params"]],
    description = kw_par[["keywrd"]]
  )

  # add previous attributes
  if (length(prev_attr) > 0) {
    for (i in names(prev_attr)) {
      attr(ff, i) <- prev_attr[[i]]
    }
  }

  return(ff)
}

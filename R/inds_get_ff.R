#' Title
#'
#' @param ind_mat
#' @param population
#' @param alias_attr_name
#' @param path_attr_name
#' @param downsample
#' @param inverse_transform
#' @param lapply_fun
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
inds_get_ff <- function(ind_mat,
                        population,
                        alias_attr_name = "short_names",
                        path_attr_name = "FilePath",
                        downsample = 1,
                        inverse_transform = c(T,F),
                        lapply_fun = lapply,
                        ...) {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }
  lapply_fun <- match.fun(lapply_fun)


  check_in(wsp = "wsp", samples = NULL, groups = NULL, FCS.file.folder = NULL, inverse_transform = inverse_transform)
  ff.list <- lapply_fun(ind_mat,
                        get_ff2,
                        downsample = downsample,
                        population = population,
                        inverse_transform = inverse_transform)

## test naming when one ind_mat has errors and NULL is returned
  names(ff.list) <- unname(sapply(ind_mat, function(x) basename(attr(x, path_attr_name))))
  ffs <- lapply(seq_along(ff.list[[1]]), function(x) {
    sapply(ff.list, "[", x, simplify = T)
  })
  names(ffs) <- stats::setNames(c("inverse", "logicle"), c(T,F))[as.character(inverse_transform)]

  return(ffs)
}

get_ff2 <- function(x, downsample, population = population, inverse_transform) {

  if (!path_attr_name %in% names(attributes(x))) {
    print(paste0(path_attr_name, " not found in attributes."))
    return(NULL)
  }
  if (!file.exists(attr(x,path_attr_name))) {
    print(paste0(attr(x,path_attr_name), " not found."))
    return(NULL)
  }

  if (population %in% colnames(x)) {
    inds <- x[,which(colnames(x) == population)]
  } else if (alias_attr_name %in% names(attributes(x)) && all(names(attr(x,alias_attr_name)) == colnames(x)) && population %in% attr(x,alias_attr_name)) {
    inds <- x[,which(attr(x,alias_attr_name) == population)]
  } else {
    print(paste0("population not found for ", attr(x, path_attr_name)))
    return(NULL)
  }

  s <- if (downsample < 1) {
    sort(sample(which(inds), ceiling(length(which(inds))*downsample)))
  } else if (downsample > 1) {
    sort(sample(which(inds), min(length(which(inds)),downsample)))
  } else {
    which(inds)
  }
  inds[which(inds)[!which(inds) %in% s]] <- F

  ff <- flowCore::read.FCS(attr(x, path_attr_name), which.lines = which(inds), truncate_max_range = F, emptyValue = F)

  if (F %in% inverse_transform) {
    if (which(inverse_transform) == 1) {
      ff <- list(ff, fcexpr:::lgcl_trsfrm_ff(ff))
    } else {
      ff <- list(fcexpr:::lgcl_trsfrm_ff(ff), ff)
    }
  } else {
    ff <- fcexpr:::lgcl_trsfrm_ff(ff)
  }

  return(ff)
}


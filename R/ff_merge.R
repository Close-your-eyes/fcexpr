#' Combine/merge flow frames
#'
#' Default arguments assume that ff in ff_list were not inverse transformed
#' during import.
#'
#' @param ff_list list of flow frames
#' @param ident_col identity column to add
#' @param grouplist other sample info to add, see ff_add_groups
#' @param keywords_to_channels keywords to add as channels; $-signs are
#' removed
#' @param add_transformed_channels also add transformed channels to final ff;
#' all arguments below only become relevant when this is TRUE
#' @param trafoname transformation name to use , must be in attributes of ff
#' @param suffix suffix to add to channel names without trafo
#' @param suffix_trafo suffix to add to channel names that were transformed
#' @param exclude channels to exclude from transformation
#' @param exclude_from were to exclude those channels from, the transformed
#' or not-transformed channels
#' @param order_first which channels to appear first in ff
#'
#' @return flow frame
#' @export
#'
#' @examples
ff_merge <- function(ff_list,
                     ident_col = "ident",
                     grouplist = NULL,
                     keywords_to_channels = NULL,
                     add_transformed_channels = F,
                     trafoname = "trafolistinv",
                     suffix = "_trans",
                     suffix_trafo = "",
                     exclude = fcexpr:::non_fluo_channels,
                     exclude_from = c("raw", "trafo"),
                     order_first = c("trafo", "raw")) {

  if (add_transformed_channels) {
    # check trafoname in attr of ff.list
    if (!all(purrr::map_lgl(ff_list, ~trafoname %in% names(attributes(.x))))) {
      stop("trafolist not found in all attributes of ff_list.")
    }
    # check any suffix unequal ""
    if (suffix_trafo == "" && suffix == "") {
      stop("One suffix is needd, both cannot be empty.")
    }
    exclude_from <- rlang::arg_match(exclude_from)
    order_first <- rlang::arg_match(order_first)


    exprs <- do.call(rbind, lapply(ff_list, flowCore::exprs))
    if (exclude_from == "raw") {
      exprs <- exprs[,which(!colnames(exprs) %in% exclude)]
    }
    colnames(exprs) <- paste0(colnames(exprs), suffix)

    exprstrafo <- do.call(rbind, lapply(ff_transform(ff = ff_list, trafolist = trafoname), flowCore::exprs))
    if (exclude_from == "trafo") {
      exprstrafo <- exprstrafo[,which(!colnames(exprstrafo) %in% exclude)]
    }
    colnames(exprstrafo) <- paste0(colnames(exprstrafo), suffix_trafo)

    if (order_first == "trafo") {
      exprs <- cbind(exprstrafo, exprs)
    } else {
      exprs <- cbind(exprs, exprstrafo)
    }

  } else {
    exprs <- do.call(rbind, lapply(ff_list, flowCore::exprs))
  }

  if (!is.null(ident_col)) {
    exprs <- cbind(exprs, ident = rep(1:length(ff_list), sapply(ff_list, nrow)))
    colnames(exprs)[ncol(exprs)] <- ident_col
  }

  new_kw <- ff_get_common_kw(ff_list)
  new_desc <- ff_get_common_desc(ff_list)
  new_pars <- flowCore::parameters(ff_list[[1]])

  kw_par <- fcexpr:::get_kw_and_pars(exprs = exprs,
                                     new_kw = new_kw,
                                     new_desc = new_desc,
                                     new_pars = new_pars)


  # add desc to channels with suffix
  if (add_transformed_channels && any(!is.na(kw_par[["new_pars"]]@data$desc))) {
    dat <- kw_par[["new_pars"]]@data
    descs <- stats::setNames(dat$desc, dat$name)
    descs <- descs[which(!is.na(descs))]
    names_strip <- gsub(suffix_trafo, "", gsub(suffix, "", dat$name))
    dat$desc <- descs[names_strip]
    kw_par[["new_pars"]]@data <- dat
  }

  ## carry on attributes?

  ff <- methods::new(
    "flowFrame",
    exprs = exprs,
    parameters = kw_par[["new_pars"]],
    description = kw_par[["new_kw"]]
  )

  if (!is.null(ident_col) && !is.null(names(ff_list))) {
    ## add ident again to have it in grouplist attr by standard method
    ff <- ff_add_groups(ff, ident_col = ident_col, grouplist = list(ident = stats::setNames(seq_along(ff_list), names(ff_list))), overwrite = T)
  }
  if (!is.null(grouplist) && !is.null(ident_col)) {
    ff <- ff_add_groups(ff, ident_col = ident_col, grouplist = grouplist)
  }

  if (!is.null(keywords_to_channels)) {
    kw_list <- fcexpr:::ff_get_kw(ff_list = ff_list, keywords = keywords_to_channels)
    ff <- ff_add_groups(ff, ident_col = ident_col, grouplist = kw_list)
  }

  return(ff)
}

ff_get_kw <- function(ff_list,
                      keywords) {

  kw_df <-
    purrr::map_dfr(ff.list, ~stack(unlist(flowCore::keyword(.x))), .id = "FileName") |>
    dplyr::rename("value" = values, "name" = ind) |>
    tidyr::pivot_wider()

  # check for missing kw
  if (any(!keywords %in% names(kw_df))) {
    missingkw <- keywords[which(!keywords %in% names(kw_df))]
    message("Keywords not found: ", paste(missingkw, collapse = ", "))
    keywords <- keywords[which(!keywords %in% missingkw)]
    if (length(keywords) == 0) {
      message("No keywords left.")
      return(NULL)
    }
  }
  kw_df <- dplyr::select(kw_df, dplyr::all_of(keywords))

  # check for kw with NA (missing in some)
  NAkw <- apply(kw_df, 2, anyNA)
  if (any(NAkw)) {
    message("Keywords with NA: ", paste(names(NAkw[which(NAkw)]), collapse = ", "), ". Some FCS files may not have this keyword.")
    keywords <- keywords[which(!keywords %in% names(NAkw[which(NAkw)]))]
    if (length(keywords) == 0) {
      message("No keywords left.")
      return(NULL)
    }
  }
  kw_df <- dplyr::select(kw_df, dplyr::all_of(keywords))

  names(kw_df) <- gsub("\\$", "", names(kw_df))
  kw_list <- as.list(kw_df)
  numeric_cols <- suppressWarnings(apply(kw_df, 2, function(x) !anyNA(as.numeric(x))))
  kw_list[numeric_cols] <- lapply(kw_list[numeric_cols], as.numeric)
  return(kw_list)
}

lgcl_trsfrm_ff <- function(ff, m_max = 500, channels = NULL, ...) {
  # ... argument like .progress
  #future::plan(future::multisession, workers = 3)
  #future::plan(future::sequential)
  if (is.null(channels)) {
    channels <- colnames(flowCore::exprs(ff))
    channels <- channels[which(channels != flowCore:::findTimeChannel(ff))]
  }

  trfms <- furrr::future_map(channels, function(z) {
    m <- 4.5
    lgcl <- NULL
    while(is.null(lgcl) && m < m_max) {
      lgcl <- tryCatch(flowCore::estimateLogicle(ff, z, m = m), error = function(e) NULL)
      m <- m + 0.1
    }
    if (is.null(lgcl)) {
      warning("m reached ", m_max, ". No logicle trans found for channel ", z, ".")
    }
    return(lgcl)
  }, ...)
  trfms <- purrr::discard(trfms, is.null)
  trfms_list <- flowCore::transformList(purrr::map_chr(trfms, function(x) names(x@transforms)),
                                        purrr::map(trfms, function(x) x@transforms[[1]]@f))

  ff <- flowCore::transform(ff, trfms_list)
  #for (i in seq_along(trfms)) {ff <- flowCore::transform(ff, trfms[[i]])}

  return(ff)
}

wsp_xml_get_samples <- function(x) {

  if (is.character(x)) {
    x <- xml2::read_xml(x)
  }
  s <- as.data.frame(t(sapply(xml2::xml_children(xml2::xml_child(x, "SampleList")), function(y) {
    xml2::xml_attrs(xml2::xml_child(y, "DataSet"))
  })), stringsAsFactors = F)
  names(s) <- c("FilePath", "sampleID")
  s$FilePath <- gsub("file:", "", s$FilePath)
  s$FileName <- basename(s$FilePath)
  return(s)
}

shortest_unique_path <- function(p) {
  p_rev <- sapply(strsplit(p, "/"), rev)
  p_rev <- lapply(seq_along(p_rev), function(x) {
    i<-1
    while (any(sapply(p_rev[-x], function(y) {
      identical(p_rev[[x]][1:i], y[1:i])
    }))) {
      i<-i+1
    }
    return(p_rev[[x]][1:i])
  })
  p <- sapply(sapply(p_rev, rev), function(x) paste(x, collapse = "/"))
  return(p)
}

check_ws <- function(ws) {
  if (is.character(ws)) {
    if (!file.exists(ws)) {
      stop("ws not found.")
    }
    if (length(ws) > 1) {
      stop("Only provide one workspace (ws) at a time.")
    }
    if (!grepl("\\.", basename(ws))) {
      stop("Did you pass a directory as ws? Please pass the full path to the wsp-file.")
    }
    if (rev(strsplit(ws, "\\.")[[1]])[1] != "wsp") {
      stop("ws has to be a file path that ends with .wsp.")
    }
    ws <- xml2::read_xml(ws)
  }
  if (!any(class(ws) == "xml_document")) {
    stop("ws must be a xml-document or a character path to its location on disk")
  }
  return(ws)
}

get_ff <- function(x,
                   return_untransformed = T,
                   return_transformed = T,
                   population,
                   seed = 42,
                   channels = NULL,
                   leverage_score_for_sampling = F,
                   return_ind_mat_only = F) {

  if ("downsample" %in% names(attributes(x))) {
    downsample <- suppressWarnings(as.numeric(attr(x, "downsample"))) # when min --> NA
  } else {
    downsample <- NA
  }
  if (is.na(downsample)) {
    downsample <- 1 # set to neutral when downsample was "min", then sampling happens outside of get_ff
  }

  if (!return_untransformed && !return_transformed) {
    stop("At least one of return_untransformed or return_transformed has to be TRUE.")
  }

  if (downsample == 1 && leverage_score_for_sampling) {
    message("No downsampling with leverage_score_for_sampling = T is not meaningful. leverage_score_for_sampling set to F.")
    leverage_score_for_sampling <- F
  }

  if (downsample != 1 && leverage_score_for_sampling && (!requireNamespace("Seurat", quietly = T) || utils::packageDescription("Seurat")[["RemoteRef"]] != "feat/dictionary")) {
    if (!requireNamespace("remotes", quietly = T)) {
      utils::install.packages("remotes")
    }
    remotes::install_github("satijalab/seurat", "feat/dictionary")
  }

  if (!is.null(channels) && !leverage_score_for_sampling) {
    message("channels are only needed for leverage score aided sampling. since leverage_score_for_sampling = F channels are ignored.")
  }

  ind_mat <- flowWorkspace::gh_pop_get_indices_mat(gh = x, y = flowWorkspace::gh_get_pop_paths(x = x))
  attr(ind_mat, "short_names") <- stats::setNames(shortest_unique_path(colnames(ind_mat)), nm = colnames(ind_mat))
  attr(ind_mat, "ws") <- attr(x, "ws")
  attr(ind_mat, "FilePath") <- attr(x, "FilePath")

  if (return_ind_mat_only) {
    flowWorkspace::gs_cleanup_temp(x)
    return(ind_mat)
  }

  if (return_untransformed && !return_transformed) {
    inverse_transform <- stats::setNames(T, "untransformed")
  } else if (!return_untransformed && return_transformed) {
    inverse_transform <- stats::setNames(F, "transformed")
  } else if (return_untransformed && return_transformed) {
    inverse_transform <- stats::setNames(c(F,T), c("transformed", "untransformed"))
  }

  ff <- tryCatch({
    lapply(inverse_transform, function(y) flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gh_pop_get_data(x, inverse.transform = y)))
  }, error = function(err) {
    message(err)
    message("Do FCS files contain a valid compensation matrix?")
  })

  population <-
    population |>
    dplyr::filter(FileName == attr(x, "FlowJoFileName")) |>
    dplyr::pull(Population)

  # alter population here by an optional leading forward slash in order to not make changes to ind_mat construction which could
  # have effects elsewhere. Maybe find a better solution some when
  # wsx_get_poppaths(x) return population paths without leading fwd slash
  inds <- ind_mat[,ifelse(population %in% attr(ind_mat, "short_names"),
                          names(attr(ind_mat, "short_names"))[which(attr(ind_mat, "short_names") %in% population)],
                          ifelse(grepl("^/", population), population, paste0("/", population))),
                  drop=F]

  if (leverage_score_for_sampling) {
    message("Calculating leverage scores.")
    channels <- .get.channels(ff[[1]], channels = channels)
    lev_scores <- lapply(asplit(inds, 2), function(x) {
      Seurat::LeverageScore(object = t(flowCore::exprs(ff[[1]])[which(x),channels]), verbose = F, seed = seed)
    })
  } else {
    lev_scores <- lapply(asplit(inds, 2), function(x) rep(1, length(which(x))))
  }

  if (downsample != 1) {
    s <- mapply(sampling_fun, inds_col = asplit(inds, 2), lev_score = lev_scores)
  } else {
    s <- lapply(asplit(inds, 2), which)
  }

  inds <- mapply(inds = asplit(inds, 2), s = s, function(inds,s) {
    inds[which(inds)[!which(inds) %in% s]] <- F
    return(inds)
  }, SIMPLIFY = F)
  #browser()
  # pull multiple population from flowframe
  ff <- lapply(inds, function(x) {
    for (i in seq_along(ff)) {
      ff[[i]] <- flowCore::Subset(ff[[i]], as.logical(x)) # inds may become a 1d array, maybe if inds has 1 column only (when there is 1 population only)
    }
    return(ff)
  })
  names(ff) <- attr(ind_mat, "short_names")[names(ff)]


  # ind_mat <- ind_mat[which(inds),,drop=F]
  # attr(ind_mat, "short_names") <- stats::setNames(shortest_unique_path(colnames(ind_mat)), nm = colnames(ind_mat))
  # attr(ind_mat, "ws") <- x$wsp
  # attr(ind_mat, "FilePath") <- x$FilePath

  flowWorkspace::gs_cleanup_temp(x)
  return(list(ff = ff, ind_mat = ind_mat))
}

sampling_fun <- function(inds_col, lev_score) {
  set.seed(seed)
  s <- sort(sample(x = which(inds_col),
                   size = ifelse(downsample < 1,
                                 ceiling(length(which(inds_col))*downsample),
                                 min(c(length(which(inds_col)), downsample))),
                   prob = lev_score))
}

get_ff2 <- function(x,
                    population,
                    alias_attr_name = "short_names",
                    path_attr_name = "FilePath",
                    leverage_score_for_sampling = F,
                    channels = NULL,
                    seed = 42) {

  if (!path_attr_name %in% names(attributes(x))) {
    message(path_attr_name, " not found in attributes.")
    return(NULL)
  }
  if (!file.exists(attr(x,path_attr_name))) {
    message(attr(x,path_attr_name), " not found.")
    return(NULL)
  }

  if ("downsample" %in% names(attributes(x))) {
    downsample <- suppressWarnings(as.numeric(attr(x, "downsample"))) # when min --> NA
  } else {
    downsample <- NA
  }
  if (is.na(downsample)) {
    downsample <- 1 # set to neutral when downsample was "min", then sampling happens outside of get_ff2
  }


  if (downsample == 1 && leverage_score_for_sampling) {
    message("No downsampling with leverage_score_for_sampling = T is not meaningful. leverage_score_for_sampling set to F.")
    leverage_score_for_sampling <- F
  }

  if (leverage_score_for_sampling && (!requireNamespace("Seurat", quietly = T) || utils::packageDescription("Seurat")[["RemoteRef"]] != "feat/dictionary")) {
    if (!requireNamespace("remotes", quietly = T)) {
      utils::install.packages("remotes")
    }
    remotes::install_github("satijalab/seurat", "feat/dictionary")
  }

  if (!is.null(channels) && !leverage_score_for_sampling) {
    message("channels are only needed for leverage score aided sampling. leverage_score_for_sampling = F though, so channels are ignored.")
  }


  col_inds <- which(colnames(x) %in% population)
  if (alias_attr_name %in% names(attributes(x))) {
    col_inds <- c(col_inds, which(attr(x,alias_attr_name) %in% population))
  }
  col_inds <- unique(col_inds)

  if (length(col_inds) == 0) {
    message("no population found")
  } else if (length(col_inds) != length(population)) {
    message("not all population found.")
  }
  inds <- x[, col_inds]
  ff <- flowCore::read.FCS(attr(x, path_attr_name), truncate_max_range = F, emptyValue = F)

  if (leverage_score_for_sampling) {
    message("Calculating leverage scores.")
    channels <- .get.channels(ff[[1]], channels = channels)
    lev_scores <- lapply(asplit(inds, 2), function(x) {
      Seurat::LeverageScore(object = t(flowCore::exprs(ff[[1]])[which(x),channels]), verbose = F, seed = seed)
    })
  } else {
    lev_scores <- lapply(asplit(inds, 2), function(x) rep(1, length(which(x))))
  }

  if (downsample != 1) {
    s <- mapply(sampling_fun, inds_col = asplit(inds, 2), lev_score = lev_scores)
  } else {
    s <- lapply(asplit(inds, 2), which)
  }

  inds <- mapply(inds = asplit(inds, 2), s = s, function(inds,s) {
    inds[which(inds)[!which(inds) %in% s]] <- F
    return(inds)
  }, SIMPLIFY = F)

  # pull multiple population from flowframe
  ff <- lapply(inds, function(x) {
    # for (i in seq_along(ff)) {
    #   ff[[i]] <- flowCore::Subset(ff[[i]], as.logical(x)) # inds may become a 1d array, maybe if inds has 1 column only (when there is 1 population only)
    # }
    # inds may become a 1d array, maybe if inds has 1 column only (when there is 1 population only)
    flowCore::Subset(ff, as.logical(x))
  })

  if (alias_attr_name %in% names(attributes(x))) {
    names(ff) <- attr(x, alias_attr_name)[col_inds]
  }
  return(ff)
}

get_gs <- function(x,
                   remove_redundant_channels = F,
                   merge_to_gs = T) {

  # split(x, (seq(nrow(x))-1) %/% split_size
  message("tempdir: ", tempdir(), "\n")
  gs_list <- lapply(asplit(x,1), function(y) {
    message(y[["FlowJoFileName"]], ", ", format(as.numeric(y[["$TOT"]]), big.mark = ","), " evts, ", round(file.info(y[["FilePathUse"]])[["size"]]/1000/1000), " Mb")

    y[["equal_FileDirs"]] <- F # using filename as subset caused error
    gs <- CytoML::flowjo_to_gatingset(CytoML::open_flowjo_xml(y[["wsp"]]),
                                      name = unique(y[["FlowJoGroup"]]),
                                      path = dirname(y[["FilePathUse"]]),
                                      subset = if (y[["equal_FileDirs"]]) {basename(y[["FilePathUse"]])} else {`$FIL` %in% y[["$FIL"]] & `$TOT` %in% y[["$TOT"]] & `$ETIM` %in% y[["$ETIM"]] & `$BTIM` %in% y[["$BTIM"]]},
                                      truncate_max_range = F,
                                      keywords = c("$FIL", "$ETIM", "$BTIM", "$TOT"),
                                      additional.keys = c("$TOT", "$ETIM", "$BTIM"))
    if (length(gs) != 1) {
      message("number of samples in gs is not 1. this should not be.")
    }

    # validity check in flowWorkspace::merge_list_to_gs checks for eqaul names of gs. they have to be unique.
    # so using y[["FlowJoFileName"]] may not be sufficient to avoid ambiguities, use identity
    flowWorkspace::sampleNames(gs) <- y[["identity"]]
    attr(gs, "ws") <- y[["wsp"]]
    attr(gs, "FlowJoFileName") <- y[["FlowJoFileName"]]
    attr(gs, "FilePath") <- y[["FilePathUse"]]
    if ("downsample" %in% names(y)) {
      attr(gs, "downsample") <- y[["downsample"]]
    }

    #rownames(y) <- paste(y[["$FIL"]], y[["$TOT"]], y[["$ETIM"]], y[["$BTIM"]], sep = "_")
    #flowWorkspace::sampleNames(gs) <- y[flowWorkspace::sampleNames(gs),"FileName"]

    if (remove_redundant_channels) {
      gs <- suppressMessages(flowWorkspace::gs_remove_redundant_channels(gs))
    }

    # identity has become rownames of pData, name column may still contain duplicate fcs file names --> make them unique
    flowCore::pData(gs)$name <- rownames(flowCore::pData(gs))
    path <- list.files(flowWorkspace::gs_get_uri(gs), full.names = T)
    message("Written to temdir/", file.path(basename(dirname(path)), basename(path)), "\n")
    return(gs)
  })

  if (merge_to_gs) {
    return(flowWorkspace::merge_list_to_gs(gs_list))
  } else {
    return(stats::setNames(gs_list, basename(x$FilePathUse)))
  }

}

.get.channels <- function(ff,
                          timeChannel = NULL,
                          channels = NULL) {
  if (!is.null(timeChannel)) {
    timeChannel <- unlist(lapply(timeChannel, function(x) grep(paste0("^",x,"$"),
                                                               colnames(flowCore::exprs(ff)),
                                                               value = TRUE, ignore.case = TRUE)))
    if (all(is.na(timeChannel))) {
      stop("None of timeChannels not found in exprs of flowFrame.")
    }
  }

  if (is.null(channels)) {
    channels <- stats::setNames(flowCore::pData(flowCore::parameters(ff))$name, flowCore::pData(flowCore::parameters(ff))$desc)
    channels <- channels[which(!channels %in% timeChannel)]
  } else {
    channels <- trimws(channels)
    inds <- unique(c(which(flowCore::pData(flowCore::parameters(ff))$name %in% channels),
                     which(flowCore::pData(flowCore::parameters(ff))$desc %in% channels)))
    notfound <- channels[intersect(which(!channels %in% flowCore::pData(flowCore::parameters(ff))$name),
                                   which(!channels %in% flowCore::pData(flowCore::parameters(ff))$desc))]
    if (length(notfound) > 0) {
      warning(paste0(paste("These channels were not found in all flowFrames: ", notfound, collapse = ", "), "."))
    }

    channels_ff <- stats::setNames(flowCore::pData(flowCore::parameters(ff))$name[inds], nm = flowCore::pData(flowCore::parameters(ff))$desc[inds])
    channels_match_inds <- unique(c(which(channels %in% channels_ff),
                                    which(channels %in% names(channels_ff)),
                                    which(names(channels) %in% channels_ff),
                                    which(names(channels) %in% names(channels_ff))))
    channels <- channels_ff[channels_match_inds]
    na_inds <- which(is.na(names(channels)))
    names(channels)[na_inds] <- stats::setNames(names(channels_ff), nm = channels_ff)[channels[na_inds]]
    diff_inds <- which(!channels %in% channels_ff)
    if (length(diff_inds) > 0) {
      channels[diff_inds] <- channels_ff[names(channels[diff_inds])]
    }
    # order by ff, important!
    channels <- channels[order(match(channels, flowCore::pData(flowCore::parameters(ff))$name))]

  }
  if (length(channels) == 0) {
    stop("no channels matched to those in the flowFrame.")
  }
  return(channels)
}

.check.ff.list <- function(ff.list, channels = NULL, strict = T) {

  ## combine with .get.channels?
  ## check if untransformed and transformed ffs are equal
  if (length(ff.list) > 2) {
    stop("ff.list can not be larger than 2.")
  }

  if (length(ff.list) == 2) {
    if(any(unlist(purrr::map2(ff.list[[1]], ff.list[[2]], ~ length(unique(list(flowCore::pData(flowCore::parameters(.x))[,c("name", "desc")], flowCore::pData(flowCore::parameters(.y))[,c("name", "desc")]))) != 1)))) {
      stop("One or more paired flowframes (transformed and untransformed) do share the same pData.")
    }
  }

  if (strict) {
    #out <- purrr::map(.x = ff.list, .f = ~purrr::map_dfr(.x = .x, .f = ~flowCore::parameters(.x)$name)) ## change this somehow (transformed and untransformed are combined)
    out <- purrr::map(.x = ff.list[[1]], .f = ~flowCore::parameters(.x)$name)
    out <- purrr::pmap_lgl(out, ~length(unique(.x)) == 1)
    if (!all(out)) {
      warning("Channels of flowFrames do not have the same names. This cannot be handled. Will return data frame(s) of channel names.")
      return(purrr::map(.x = ff.list, .f = ~purrr::map(.x = .x, .f = ~flowCore::parameters(.x)$name)))
    }


    ## NA-columns are return without row names, hence set row names manually for binding to df
    #out <- purrr::map(.x = ff.list, .f = ~purrr::map_dfr(.x = .x, .f = ~stats::setNames(flowCore::pData(flowCore::parameters(.x))[,"desc"], flowCore::pData(flowCore::parameters(.x))[,"name"])))
    out <- purrr::map(.x = ff.list[[1]], .f = ~stats::setNames(flowCore::pData(flowCore::parameters(.x))[,"desc"], flowCore::pData(flowCore::parameters(.x))[,"name"]))
    out <- purrr::pmap_lgl(out, ~length(unique(.x)) == 1)
    if (!all(out)) {
      warning("Channel description are not equal across flowFrames.")
    }
    return(NULL)
  }

  if (!strict) {
    ## names
    #out <- purrr::map_dfr(.x = ff.list, .f = ~purrr::map_dfr(.x = ff.list[[1]], .f = ~flowCore::parameters(.x)$name))
    #out2 <- apply(out, 2, function(x) unique(x))

    out <- purrr::map(.x = ff.list[[1]], .f = ~flowCore::parameters(.x)$name)
    out2 <- purrr::pmap(out, ~unique(.x))
    if (any(purrr::map_int(out2, length) > 1)) {
      if (any(channels %in% unlist(out2[which(purrr::map_int(out2, length) > 1)]))) {
        warning("Channels of flowframes do not have the same names including one of selected channels.
        If this is intended try to select respective channels by equal channel descriptions.
                Modify flowframes accordingly before.
                Will now return data frame of channel names.")
        return(out)
      } else {
        warning("Channels of flowFrames do not have the same names. But non of selected channels is affected/included.")
      }
    }


    #descs
    out <- purrr::map(.x = ff.list[[1]], .f = ~stats::setNames(flowCore::pData(flowCore::parameters(.x))[,"desc"], flowCore::pData(flowCore::parameters(.x))[,"name"]))
    out <- purrr::map(out, function(x) x[which(!is.na(x))])
    channels_descs <- channels[which(channels %in% unique(unlist(out)))]

    if (length(unique(out)) != 1) {
      # check for uniqueness
      message("Channel descriptions are not equal across flowframes.")
      if (all(channels_descs %in% purrr::reduce(out, intersect))) {
        message("Selected channels are found in every flowframe though.")
        out2 <- purrr::map(.x = .x, .f = ~flowCore::pData(flowCore::parameters(.x))[,c("name", "desc")])
        out2 <- purrr::map(out2, tidyr::drop_na)
        if (length(unique(out2)) != 1) {
          message("Equal channel descriptions belong to different channels:")
          print(unique(out2))
        }
      } else {
        warning("At least one selected channel are affected. Please check and fix.
                Will return list of unique channel names and descriptions now.")
        return(unique(purrr::map(.x = .x, .f = ~flowCore::pData(flowCore::parameters(.x))[,c("name", "desc")])))
      }
    }
  }
  return(NULL)
}


min.max.normalization <- function (x, min.val = 0, max.val = 1) {
  if (is.matrix(x) || is.data.frame(x)) {
    if (is.data.frame(x)) {
      if (!all(apply(x, 2, is.numeric))) {
        stop("Please make sure that all columns of the data frame are numeric.")
      }
    }
    return(apply(x, 2, function (y) min.val + ((y- min(y)) * (max.val- min.val) / (max(y)-min(y)))))
  } else {
    return(min.val + ((x- min(x)) * (max.val- min.val) / (max(x)-min(x))))
  }
}


shift.to.positive <- function(x, rm.na = F) {
  if (rm.na) {
    x <- x[which(!is.na(x))]
  }
  if (min(x) <= 0) {
    return(x + abs(min(x - 1)))
  } else {
    return(x)
  }
  #one-liner: apply(x, 2, function(z){ if (min(z) < 0) {z + abs(min(z - 1))} else {z}})
}



lgcl_trsfrm_ff <- function(ff, channels = NULL) {

  if (is.null(channels)) {
    channels <- colnames(flowCore::exprs(ff))
    channels <- channels[which(channels != flowCore:::findTimeChannel(ff))]
  }

  trfms <- lapply(channels, function(z) {
    m <- 4.5
    lgcl <- NULL
    while(is.null(lgcl)) {
      lgcl <- tryCatch(flowCore::estimateLogicle(ff, z, m = m),
                       error = function(e) {
                         #message(m)
                         return(NULL)
                       }
      )
      m <- m + 0.1
    }
    return(lgcl)
  })

  for (i in seq_along(trfms)) {
    ff <- flowCore::transform(ff, trfms[[i]])
  }
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
    stop("x must be a xml-document or a character path to its location on disk")
  }
  return(ws)
}


get_smpl_df <- function(wsp,
                        groups,
                        invert_groups,
                        samples,
                        invert_samples,
                        FCS.file.folder) {

  smpl <- do.call(rbind, lapply(seq_along(wsp), function(x) {
    y <- wsx_get_fcs_paths(wsp[x], split = F)
    y$wsp <- wsp[x]
    y$FileName <- basename(y$FilePath)

    key <- sapply(wsx_get_keywords(wsp[x]), function(z) {
      z[which(z$name == "$FIL"),"value"]
    })
    y$FIL <- key[y$FileName]
    key <- sapply(wsx_get_keywords(wsp[x]), function(z) {
      z[which(z$name == "$TOT"),"value"]
    })
    y$TOT <- trimws(key[y$FileName])
    key <- sapply(wsx_get_keywords(wsp[x]), function(z) {
      z[which(z$name == "$BEGINDATA"),"value"]
    })
    y$BEGINDATA <- key[y$FileName]

    if (!is.null(groups)) {
      if (invert_groups) {
        y <- y[which(!y$group %in% groups[[x]]),]
      } else {
        y <- y[which(y$group %in% groups[[x]]),]
      }
    }
    if(nrow(y) == 0) {
      return(NULL)
    }

    if (!is.null(samples)) {
      if (invert_samples) {
        y <- y[which(!y$FileName %in% samples[[x]]),]
      } else {
        y <- y[which(y$FileName %in% samples[[x]]),]
      }
    }
    if (nrow(y) == 0) {
      return(NULL)
    }

    # remove All Samples group
    y <- do.call(rbind, lapply(unique(y$sampleID), function(zz) {
      if (length(y[which(y$sampleID== zz),"group"]) > 1) {
        y[intersect(which(y$sampleID == zz), which(y$group != "All Samples")), ]
      } else {
        y[which(y$sampleID == zz), ]
      }
    }))
    # remove other duplicates (multiple groups)
    y <- dplyr::arrange(y, group)
    y <- dplyr::distinct(y, FilePath, wsp, .keep_all = T)

    if (any(duplicated(y[,which(names(y) %in% c("FIL", "TOT", "BEGINDATA"))]))) {
      stop("Samples cannot be identified unambiguously.")
    }

    if (is.null(FCS.file.folder)) {
      y$FCS.file.folder <- NA
    } else {
      y$FCS.file.folder <- FCS.file.folder[x]
      y$FilePath <- sapply(y$FileName, function(z) {
        match_files <- list.files(path = FCS.file.folder[x], recursive = T, full.names = T, pattern = z)
        if (length(match_files) > 1) {
          message("Found multiple FCS files with equal names. Will select the one which matches best the keywords from flowjo workspace.")
          # match via keywords
          all_key_wsx <- wsx_get_keywords(wsp[x])[[z]]
          all_key_wsx <- all_key_wsx[which(!grepl("spill|^\\$P|^P[[:digit:]]{1,}", all_key_wsx$name, ignore.case = T)),]
          keysss <- lapply(match_files, function(match_file) {
            all_key_fcs <- utils::stack(flowCore::read.FCSheader(match_file)[[1]])
            names(all_key_fcs) <- names(all_key_wsx)[c(2,1)]
            all_key_fcs <- all_key_fcs[which(!grepl("spill|^\\$P|^P[[:digit:]]{1,}", all_key_fcs$name, ignore.case = T)),]
            all_key_fcs <- all_key_fcs[which(trimws(all_key_fcs$value) != ""),]
            return(all_key_fcs)
          })
          intersect_keys <- Reduce(intersect, c(list(all_key_wsx$name), unname(sapply(keysss, "[", "name"))))

          scores <- sapply(keysss, function(match_file_keys) {
            match_file_keys <- match_file_keys[which(match_file_keys$name %in% intersect_keys), ]
            all_key_wsx <- all_key_wsx[which(all_key_wsx$name %in% intersect_keys), ]
            all_key_wsx <- all_key_wsx[match(all_key_wsx$name, match_file_keys$name),]
            return(length(which(all_key_wsx$value == match_file_keys$value)))
          })
          # select best match
          return(match_files[which.max(scores)])
        } else {
          return(match_files)
        }
      })
    }
    return(y)
  }))

  return(smpl)
}

check_in <- function(wsp,
                     samples,
                     groups,
                     FCS.file.folder,
                     return_untransformed = NULL,
                     return_transformed = NULL) {

  if (missing(wsp) || class(wsp) != "character") {stop("Please provide a vector of paths to flowjo workspaces.")}
  if (!is.null(groups)) {
    if (class(groups) == "list" && length(groups) != length(wsp)) {stop("list of groups has to have the same length as wsp. Alternatively pass a vector groups to use for all workspace.")}
    if (class(groups) != "list") {groups <- rep(list(groups), length(wsp))}
  }
  if (!is.null(samples)) {
    if (class(samples) == "list" && length(samples) != length(wsp)) {stop("list of samples has to have the same length as wsp. Alternatively pass a vector samples to use for all workspace.")}
    if (class(samples) != "list") {samples <- rep(list(samples), length(wsp))}
  }
  if (!is.null(FCS.file.folder)) {
    if (any(!dir.exists(FCS.file.folder))) {stop(paste0(FCS.file.folder[which(!dir.exists(FCS.file.folder))], " not found."))}
    if (length(FCS.file.folder) == 1) {FCS.file.folder <- rep(FCS.file.folder, length(wsp))}
    if (length(FCS.file.folder) != length(wsp)) {stop("FCS.file.folder has to have the same length as wsp or 1.")}
  }

  if ((!is.null(return_untransformed) && !return_untransformed) && (!is.null(return_transformed) && !return_transformed)) {
    stop("At least one of return_transformed or return_untransformed has to be TRUE.")
  }

  #inverse_transform <- unique(inverse_transform)
  #if (!length(inverse_transform) %in% c(1,2)) {stop("inverse_transform must be of length 1 or 2, T or F or c(T,F) or c(F,T)")}

  return(list(groups = groups, samples = samples, FCS.file.folder = FCS.file.folder))
}

get_ff <- function(x,
                   return_untransformed = T,
                   return_transformed = T,
                   downsample = 1,
                   remove_redundant_channels = F,
                   population,
                   seed = 42,
                   channels = NULL,
                   leverage_score_for_sampling = F,
                   return_ind_mat_only = F) {

  # one file at a time avoids problems due to different gating trees, but this may leave unintentional different gating trees undetected
  # pass full path as attr and check consistency later?

  if (nrow(x) > 1) {
    stop("Only one fcs file at a time.")
  }
  if (!return_untransformed && !return_transformed) {
    stop("At least one of return_untransformed or return_transformed has to be TRUE.")
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

  if (is.na(x$FCS.file.folder)) {
    path <- dirname(x$FilePath)
    if (!file.exists(path)) {
      stop(paste0(path, " not found. Was the workspace saved on another computer? If so, reconnect FCS files in flowjo or provide the FCS.file.folder(s) on the current computer."))
    }
  } else {
    path <- x$FCS.file.folder
  }

  gs <- CytoML::flowjo_to_gatingset(ws = CytoML::open_flowjo_xml(x$wsp),
                                    name = x$group,
                                    path = path,
                                    subset = `$FIL` == x$FIL & `$TOT` == x$TOT & `$BEGINDATA` == x$BEGINDATA, # not && !
                                    truncate_max_range = F,
                                    transform = T,
                                    keywords = c("$FIL", "$TOT", "$BEGINDATA"),
                                    additional.keys = c("$TOT", "$BEGINDATA"))

  ind_mat <- flowWorkspace::gh_pop_get_indices_mat(gs[[1]], y = flowWorkspace::gh_get_pop_paths(gs[[1]]))
  attr(ind_mat, "short_names") <- stats::setNames(shortest_unique_path(colnames(ind_mat)), nm = colnames(ind_mat))
  attr(ind_mat, "ws") <- x$wsp
  attr(ind_mat, "FilePath") <- x$FilePath

  if (return_ind_mat_only) {
    return(ind_mat)
  }

  if (remove_redundant_channels) {
    gs <- suppressMessages(flowWorkspace::gs_remove_redundant_channels(gs))
  }

  if (return_untransformed && !return_transformed) {
    inverse_transform <- stats::setNames(T, "untransformed")
  } else if (!return_untransformed && return_transformed) {
    inverse_transform <- stats::setNames(F, "transformed")
  } else if (return_untransformed && return_transformed) {
    inverse_transform <- stats::setNames(c(F,T), c("transformed", "untransformed"))
  }

  ex <- lapply(inverse_transform, function(y) flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gh_pop_get_data(gs[[1]], inverse.transform = y)))

  inds <- ind_mat[,ifelse(population %in% attr(ind_mat, "short_names"),
                          names(which(attr(ind_mat, "short_names") == population)),
                          population),drop=T]

  ## overwrite downsample argument if provided as attr in x
  if ("downsample" %in% names(x)) {
    downsample <- x$downsample
  }

  if (leverage_score_for_sampling) {
    message("Calculating leverage scores.")
    channels <- .get.channels(ex[[1]], channels = channels)
    lev_scores <- Seurat::LeverageScore(object = t(flowCore::exprs(ex[[1]])[which(inds),channels]), verbose = F, seed = seed)
  } else {
    lev_scores <- rep(1, nrow(flowCore::exprs(ex[[1]])))
  }

  if (downsample != 1) {
    set.seed(seed)
    s <- sort(sample(x = which(inds),
                     size = ifelse(downsample < 1, ceiling(length(which(inds))*downsample), min(c(length(which(inds)), downsample))),
                     prob = lev_scores))
  } else {
    s <- which(inds)
  }
  inds[which(inds)[!which(inds) %in% s]] <- F

  for (i in seq_along(ex)) {
    ex[[i]] <- subset(ex[[i]], inds)
  }

  ind_mat <- ind_mat[which(inds),,drop=F]
  attr(ind_mat, "short_names") <- stats::setNames(shortest_unique_path(colnames(ind_mat)), nm = colnames(ind_mat))
  attr(ind_mat, "ws") <- x$wsp
  attr(ind_mat, "FilePath") <- x$FilePath

  flowWorkspace::gs_cleanup_temp(gs)
  return(list(expr = ex, ind_mat = ind_mat))
}


get_ff2 <- function(x,
                    downsample = 1,
                    population,
                    alias_attr_name,
                    path_attr_name,
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

  if (length(population) > 1) {
    stop("Only provide one population.")
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

  if (population %in% colnames(x)) {
    inds <- x[,which(colnames(x) == population)]
  } else if (alias_attr_name %in% names(attributes(x)) && all(names(attr(x,alias_attr_name)) == colnames(x)) && population %in% attr(x,alias_attr_name)) {
    inds <- x[,which(attr(x,alias_attr_name) == population)]
  } else {
    message("population not found for ", attr(x, path_attr_name))
    return(NULL)
  }

  ## overwrite downsample argument if provided as attr in x
  if ("downsample" %in% names(attributes(x))) {
    downsample <- attr(x, "downsample")
  }

  ff <- flowCore::read.FCS(attr(x, path_attr_name), truncate_max_range = F, emptyValue = F)

  if (leverage_score_for_sampling) {
    channels <- .get.channels(ff, channels = channels)
    lev_scores <- Seurat::LeverageScore(object = t(flowCore::exprs(ff)[which(inds),channels]), verbose = F, seed = seed)
  } else {
    lev_scores <- rep(1, nrow(flowCore::exprs(ff)))
  }

  if (downsample != 1) {
    set.seed(seed)
    s <- sample(x = which(inds),
                size = ifelse(downsample < 1, ceiling(length(which(inds))*downsample), min(c(length(which(inds)),downsample))),
                prob = lev_scores)
  } else {
    s <- which(inds)
  }
  inds[which(inds)[!which(inds) %in% s]] <- F

  return(stats::setNames(list(subset(ff, inds)), "untransformed"))
}

get_gs <- function(x,
                   remove_redundant_channels) {

  gs <- CytoML::flowjo_to_gatingset(CytoML::open_flowjo_xml(unique(x$wsp)),
                                    name = unique(x$group),
                                    path = unique(x$FCS.file.folder),
                                    subset = `$FIL` %in% x$FIL & `$BEGINDATA` %in% x$BEGINDATA & `$TOT` %in% x$TOT,
                                    truncate_max_range = F,
                                    keywords = c("$FIL", "$BEGINDATA", "$TOT"),
                                    additional.keys = c("$TOT", "$BEGINDATA"))

  rownames(x) <- paste(x$FIL, x$TOT, x$BEGINDATA, sep = "_")
  flowWorkspace::sampleNames(gs) <- x[flowWorkspace::sampleNames(gs),"FileName"]

  if (remove_redundant_channels) {
    gs <- suppressMessages(flowWorkspace::gs_remove_redundant_channels(gs))
  }

  return(gs)
}

.get.channels <- function(ff,
                          timeChannel = NULL,
                          channels = NULL) {
  if (!is.null(timeChannel)) {
    if (!timeChannel %in% colnames(flowCore::exprs(ff))) {
      stop("timeChannel not found in exprs of flowFrame.")
    }
  } else {
    timeChannel <- flowCore:::findTimeChannel(ff)
    message("time channel detected: ", timeChannel)
  }

  if (is.null(channels)) {
    channels <- stats::setNames(flowCore::pData(flowCore::parameters(ff))$name, flowCore::pData(flowCore::parameters(ff))$desc)
    channels <- channels[which(channels != timeChannel)]
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



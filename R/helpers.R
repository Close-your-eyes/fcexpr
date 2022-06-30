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
      stop("Only one ws at a time.")
    }
    if (!grepl("\\.", basename(ws))) {
      stop("Did you pass a directory as ws? Please pass the path the wsp-file.")
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
                     return_logicle_transformed = NULL) {

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

  if ((!is.null(return_untransformed) && !return_untransformed) && (!is.null(return_logicle_transformed) && !return_logicle_transformed)) {
    stop("At least one of return_logicle_transformed or return_untransformed has to be TRUE.")
  }

  #inverse_transform <- unique(inverse_transform)
  #if (!length(inverse_transform) %in% c(1,2)) {stop("inverse_transform must be of length 1 or 2, T or F or c(T,F) or c(F,T)")}

  return(list(groups = groups, samples = samples, FCS.file.folder = FCS.file.folder))
}

get_inds <- function(x) {

  if (nrow(x) > 1) {
    stop("Only one fcs file at a time.")
  }

  if (is.na(x$FCS.file.folder)) {
    #path <- x[,which(names(x) %in% c("sampleID", "FilePath")),drop=F]
    #names(path)[which(names(path) == "FilePath")] <- "file"
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
                                    keywords = c("$FIL", "$TOT", "$BEGINDATA"),
                                    additional.keys = c("$TOT", "$BEGINDATA"))

  ind_mat <- flowWorkspace::gh_pop_get_indices_mat(gs[[1]], y = flowWorkspace::gh_get_pop_paths(gs[[1]]))
  attr(ind_mat, "short_names") <- stats::setNames(shortest_unique_path(colnames(ind_mat)), nm = colnames(ind_mat))
  attr(ind_mat, "ws") <- x$wsp
  attr(ind_mat, "FilePath") <- x$FilePath

  flowWorkspace::gs_cleanup_temp(gs)
  return(ind_mat)
}


get_ff <- function(x,
                   return_untransformed = T,
                   return_logicle_transformed = T,
                   downsample,
                   remove_redundant_channels,
                   population) {

  # one file at a time avoids problems due to different gating trees, but this may leave unintentional different gating trees undetected
  # pass full path as attr and check consistency later?

  if (nrow(x) > 1) {
    stop("Only one fcs file at a time.")
  }

  if (is.na(x$FCS.file.folder)) {
    #path <- x[,which(names(x) %in% c("sampleID", "FilePath")),drop=F]
    #names(path)[which(names(path) == "FilePath")] <- "file"
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
                                    keywords = c("$FIL", "$TOT", "$BEGINDATA"),
                                    additional.keys = c("$TOT", "$BEGINDATA"))

  if (remove_redundant_channels) {
    gs <- suppressMessages(flowWorkspace::gs_remove_redundant_channels(gs))
  }

  if (return_untransformed && !return_logicle_transformed) {
    inverse_transform <- T
  } else if (!return_untransformed && return_logicle_transformed) {
    inverse_transform <- F
  } else if (return_untransformed && return_logicle_transformed) {
    inverse_transform <- c(F,T)
  }
  ex <- lapply(inverse_transform, function (y) flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gh_pop_get_data(gs[[1]], inverse.transform = y)))

  inds <- flowWorkspace::gh_pop_get_indices(gs[[1]], y = population)

  ## overwrite downsample argument if provided as attr in x
  if ("downsample" %in% names(x)) {
    downsample <- x$downsample
  }

  s <- if (downsample < 1) {
    sort(sample(which(inds), ceiling(length(which(inds))*downsample)))
  } else if (downsample > 1) {
    sort(sample(which(inds), min(length(which(inds)), downsample)))
  } else {
    which(inds)
  }
  inds[which(inds)[!which(inds) %in% s]] <- F

  for (i in seq_along(ex)) {
    ex[[i]] <- subset(ex[[i]], inds)
  }

  flowWorkspace::gs_cleanup_temp(gs)
  return(list(ex, inds))
}


get_ff2 <- function(x,
                    downsample,
                    population = population,
                    return_untransformed = return_untransformed,
                    return_logicle_transformed = return_logicle_transformed,
                    alias_attr_name,
                    path_attr_name) {

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

  s <- if (downsample < 1) {
    sort(sample(which(inds), ceiling(length(which(inds))*downsample)))
  } else if (downsample > 1) {
    sort(sample(which(inds), min(length(which(inds)),downsample)))
  } else {
    which(inds)
  }
  inds[which(inds)[!which(inds) %in% s]] <- F

  # which.lines with which(inds) argument is much slower!
  ff <- subset(flowCore::read.FCS(attr(x, path_attr_name), truncate_max_range = F, emptyValue = F), inds)

  if (return_logicle_transformed) {
    ff <- list(ff, lgcl_trsfrm_ff(ff))
  }

  if (!return_untransformed) {
    ## test!
    ff <- list(ff[[2]])
  }

  return(ff)
}

get_gs <- function(x,
                   remove_redundant_channels) {

  gs <- CytoML::flowjo_to_gatingset(CytoML::open_flowjo_xml(unique(x$wsp)),
                                    name = unique(x$group),
                                    path = unique(x$FCS.file.folder),
                                    subset = `$FIL` == x$FIL & `$BEGINDATA` == x$BEGINDATA & `$TOT` == x$TOT,
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

.check.ff.list <- function(ff.list) {

  sapply(ff.list, function (ff) {
    if(!all(apply(sapply(ff, function(x) {flowCore::parameters(x)$name}), 1, function(x) length(unique(x))) == 1)) {
      print(sapply(ff, function(x) {flowCore::parameters(x)$name}))
      stop("Channels of flowFrames do not have the same names. This cannot be handled.")
    }
  })

  sapply(ff.list, function (ff) {
    if(!all(apply(sapply(ff, function(x) {flowCore::parameters(x)$desc}), 1, function(x) length(unique(x))) == 1)) {
      print(sapply(ff, function(x) {flowCore::parameters(x)$desc}))
      warning("Channel description are not equal across flowFrames.")
    }
  })
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



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


get_smpl_df <- function(wsp,
                        groups,
                        invert_groups,
                        samples,
                        invert_samples,
                        FCS.file.folder
) {
  # FlowJoGroup
  # one wsp at a time
  smpl <- do.call(rbind, lapply(seq_along(wsp), function(x) {
    y <- wsx_get_fcs_paths(wsp[x], split = F)
    names(y)[3] <- "FlowJoFilePath"
    names(y)[4] <- "FlowJoFileName"
    y$FlowJoFilePath <- URLdecode(y$FlowJoFilePath) # nice !!!
    y$wsp <- wsp[x]
    y$FlowJoFileName <- basename(y$FlowJoFilePath) # redo this after urldecode

    kwl <- wsx_get_keywords(wsp[x], return = "vector", keywords = c("$DATE", "$BTIM", "$ETIM", "$TOT", "$FIL"))
    kwldf <- dplyr::bind_rows(wsx_get_keywords(wsp[x], return = "data.frame", keywords = c("$DATE", "$BTIM", "$ETIM", "$TOT", "$FIL")),
                              .id = "FlowJoFileName")
    kwldf <- tidyr::pivot_wider(kwldf, names_from = name, values_from = value)
    kwldf$`$TOT` <- trimws(kwldf$`$TOT`)
    fcs_ident <- stack(fcexpr:::.get_fcs_identities(kwl, allow_duplicates = T))
    names(fcs_ident) <- c("identity", "FlowJoFileName")
    fcs_ident$FlowJoFileName <- as.character(fcs_ident$FlowJoFileName)
    y <- dplyr::left_join(y, fcs_ident, by = "FlowJoFileName")
    y <- dplyr::left_join(y, kwldf, by = "FlowJoFileName")

    if (!is.null(groups)) {
      if (invert_groups) {
        y <- y[which(!y$FlowJoGroup %in% groups[[x]]),]
      } else {
        y <- y[which(y$FlowJoGroup %in% groups[[x]]),]
      }
    }
    if (nrow(y) == 0) {
      message("No FCS files left after filtering for groups. Available groups:\n", paste("'", unique(y$FlowJoGroup), "'", collapse = ", ", sep = ""))
      return(NULL)
    }

    if (!is.null(samples)) {
      if (invert_samples) {
        y <- y[which(!y$FlowJoFileName %in% samples[[x]]),]
      } else {
        y <- y[which(y$FlowJoFileName %in% samples[[x]]),]
      }
    }
    if (nrow(y) == 0) {
      message("No FCS files left after filtering for samples. Available samples:\n", paste("'", unique(y$FlowJoFileName), "'", collapse = ", ", sep = ""))
      return(NULL)
    }

    # remove All Samples group
    if ("All Samples" %in% y$FlowJoGroup) {
      y <- do.call(rbind, lapply(unique(y$sampleID), function(zz) {
        if (length(y[which(y$sampleID== zz),"FlowJoGroup"]) > 1) {
          y[intersect(which(y$sampleID == zz), which(y$FlowJoGroup != "All Samples")),,drop=F]
        } else {
          y[which(y$sampleID == zz),,drop=F]
        }
      }))
    }
    if (nrow(y) == 0) {
      message("No FCS files left after removing for 'All Samples' group. That's a bug.")
      return(NULL)
    }

    # remove duplicate fcs files, e.g. same file in multiple groups
    y <- dplyr::distinct(y, identity, .keep_all = T)

    # check if files exist under flowjos path
    y$FlowJoFilePathExists <- file.exists(y$FlowJoFilePath)

    # if (anyDuplicated(y$identity)) {
    #   message("Duplicate fcs identites from FlowJo in one wsp. This should not be. Filter somehow?")
    #   print(tibble::as_tibble(y[,c("FileName", "identity")]) |> group_by(FileName, identity) |> dplyr::mutate(group = dplyr::cur_group_id()), n = Inf)
    #   return(NULL)
    # }

    if (!is.null(FCS.file.folder)) {
      local_fcs_files <- list.files(path = FCS.file.folder[x], recursive = T, full.names = T, pattern = "\\.fcs$", ignore.case = T)
      local_fcs_files <- fcexpr:::.get_fcs_identities(kwl = flowCore::read.FCSheader(local_fcs_files))
      local_fcs_files_df <- stack(local_fcs_files)
      names(local_fcs_files_df) <- c("identity", "LocalFilePath")
      local_fcs_files_df$LocalFilePath <- as.character(local_fcs_files_df$LocalFilePath)
      local_fcs_files_df$LocalFileName <- basename(local_fcs_files_df$LocalFilePath)
      # only consider fcs files that are present in flowjo
      local_fcs_files_df <-
        local_fcs_files_df |>
        dplyr::filter(identity %in% y$identity)
      ##test
      #local_fcs_files_df <- rbind(local_fcs_files_df,local_fcs_files_df)
      if (anyDuplicated(local_fcs_files_df$identity)) {
        dups <- table(local_fcs_files_df$identity)[which(table(local_fcs_files_df$identity) > 1)]
        message("Some local FCS files have same identities. They have been copied by the user.")# Will match with files in FlowJo based on identity and FileName.")
        print(tibble::as_tibble(local_fcs_files_df) |> dplyr::filter(identity %in% names(dups)) |> dplyr::arrange(identity), n = Inf)

        # check for duplicate FileName within groups of same identity
        dup_FileName <-
          local_fcs_files_df |>
          dplyr::group_by(identity, LocalFileName) |>
          dplyr::filter(dplyr::n() > 1) |>
          dplyr::ungroup() |>
          dplyr::arrange(identity)
        if (nrow(dup_FileName) < 1) {
          message("FCS files with equal identity also have same FileNames. Will pick one file arbitrarily for reading in those cases.")
          message("Other FCS files are matched to FlowJo by identity and FileName.")
          local_fcs_files_df <-
            local_fcs_files_df |>
            dplyr::distinct(identity, LocalFileName, .keep_all = T)
        } else {
          message("FCS files are matched to FlowJo by identity and FileName.")
        }
        # when y and local_fcs_files_df are joined now, additional rows should not arise
        ## test
        #local_fcs_files_df <- rbind(local_fcs_files_df, local_fcs_files_df[4,] |> dplyr::mutate(LocalFileName = "othername")) # filtered below by joining
      }
      nrowbefore <- nrow(y)
      y <- dplyr::left_join(y, local_fcs_files_df, by = c("identity" = "identity", "FlowJoFileName" = "LocalFileName"))
      if (nrowbefore != nrow(y)) {
        message("Joining flowjo df and local df of fcs files caused new rows. This is a bug. Check.")
      }
      y$equal_FilePaths <- y$FlowJoFilePath == y$LocalFilePath
      y$equal_FileDirs <- dirname(y$FlowJoFilePath) == dirname(y$LocalFilePath)

      # equal paths means: flowjo dir was passed as fcs.file.folder and hence, it is the very same fcs file - ignore those below
      y[which(y$equal_FilePaths), "FlowJoFilePathExists"] <- F

      if (any(inds <- !y$FlowJoFilePathExists & is.na(y$LocalFilePath))) {
        message(sum(inds), " of ", nrow(y), " FCS files from FlowJo could not be matched to any local FCS file and cannot be found by the FilePath in Flowjo.")
        message(unique(basename(y$wsp)))
        message(paste("'", y[which(inds), "FlowJoFileName"], "'", collapse = ", ", sep = ""))
        y <- y[-which(inds),,drop=F]
        if (nrow(y) == 0) {
          message("No FCS files left after checking for matching local FCS files.")
          return(NULL)
        }
      }
      if (any(inds <- y$FlowJoFilePathExists & !is.na(y$LocalFilePath))) {
        message(sum(inds), " of ", nrow(y), " FCS files from FlowJo have matching FilePath in FlowJo and locally (in FCS.file.folder).")
        message(unique(basename(y$wsp)))
        message(paste("'", y[which(inds), "FlowJoFileName"], "'", collapse = ", ", sep = ""))
        y[inds,"FlowJoFilePathExists"] <- F #prefer the path from fcs.file.folder
        # ysub <- y[inds,]
        # if (any(ysub$equal_FileDirs)) {
        #   inds2 <- which(ysub$equal_FileDirs)
        #   message("At least some of these files are in the same directory. In this case the FCS file that appears first by name will be read.
        #           To force which file to read, modify one of their filename on disk and, for instance, add a prefix underscore to make sure this file comes first.")
        #   message(paste("'", ysub[inds2, "FlowJoFileName"], "'", collapse = ", ", sep = ""))
        # } else {
        #   message("Will use the path from FCS.files.folder in this case.")
        # }
        message("Will use the path from FCS.files.folder in this case.")
      }
      if (any(inds <- y$FlowJoFilePathExists & is.na(y$LocalFilePath))) {
        message(sum(inds), " of ", nrow(y), " FCS files from FlowJo have matching FilePath in FlowJo but not locally (in FCS.file.folder). Will use the FlowJo path here.")
        message(unique(basename(y$wsp)))
        message(paste("'", y[which(inds), "FlowJoFileName"], "'", collapse = ", ", sep = ""))
      }
      y$FilePathUse <- ifelse(y$FlowJoFilePathExists, y$FlowJoFilePath, y$LocalFilePath)
    } else {
      if (any(!y$FlowJoFilePathExists)) {
        message(sum(!y$FlowJoFilePathExists), " of ", nrow(y), " FCS files were not found under the path in FlowJo wsp. Without another path hint on the current machine (using the FCS.file.folder argument) these files cannot be read.")
        y <- y[which(y$FlowJoFilePathExists),,drop = F]
        if (nrow(y) == 0) {
          message("No FCS files left.")
          return(NULL)
        }
      }
      y$FilePathUse <- y$FlowJoFilePath
    }

    return(y)
  }))

  if (any(table(smpl$identity) > 1)) {
    message("Same FCS files found in multiple workspaces. This may cause problems. Provide samples argument as a list of vectors, one for each wsp and filter the duplicate files.")
    message(paste(names(table(smpl$identity))[which(table(smpl$identity) > 1)], collapse = ", "))
  }

  return(as.data.frame(smpl))
}

check_in <- function(wsp,
                     samples = NULL,
                     groups = NULL,
                     FCS.file.folder = NULL,
                     return_untransformed = NULL,
                     return_transformed = NULL) {

  # check invert samples and groups herein

  if (missing(wsp) || !methods::is(wsp, "character")) {stop("Please provide a vector of paths to flowjo workspaces.")}

  wsp_data <- lapply(wsp, wsx_get_fcs_paths, split = F, filter_AllSamples = T)

  if (!is.null(samples)) {
    if (length(samples) == 0) {
      stop("No samples provided.")
    }
    if (methods::is(samples, "list") && length(samples) != length(wsp)) {stop("list of samples has to have the same length as wsp. Alternatively pass a vector samples to use for all workspace.")}
    if (!methods::is(samples, "list")) {samples <- rep(list(samples), length(wsp))}
  }

  if (!is.null(groups)) {
    if (methods::is(groups, "list") && length(groups) != length(wsp)) {stop("list of groups has to have the same length as wsp. Alternatively pass a vector groups to use for all workspace.")}
    if (!methods::is(groups, "list")) {groups <- rep(list(groups), length(wsp))}
  } else if (is.null(samples)) {
    groups <- rep(list("All Samples"), length(wsp))
  } else {
    groups <- purrr::map2(wsp_data, samples, function(x,y) {
      dplyr::filter(x, FileName  %in% y) |> dplyr::distinct(FlowJoGroup) |> dplyr::pull(FlowJoGroup)
    })
  }

  if (is.null(samples)) {
    samples <- purrr::map2(wsp_data, groups, function(x,y) {
      dplyr::filter(x, FlowJoGroup  %in% y) |> dplyr::distinct(FileName) |> dplyr::pull(FileName)
    })
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
                   population,
                   seed = 42,
                   channels = NULL,
                   leverage_score_for_sampling = F,
                   return_ind_mat_only = F) {

  downsample <- suppressWarnings(as.numeric(attr(x, "downsample")))
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

  # pull multiple population from flowframe
  ff <- lapply(inds, function(x) {
    for (i in seq_along(ff)) {
      ff[[i]] <- flowCore::Subset(ff[[i]], x)
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
    lev_scores <- NULL
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

  return(stats::setNames(list(flowCore::Subset(ff, inds)), "untransformed"))
}

get_gs <- function(x,
                   remove_redundant_channels,
                   merge_to_gs = T) {

  # split(x, (seq(nrow(x))-1) %/% split_size
  message("tempdir: ", tempdir())
  gs_list <- lapply(asplit(x,1), function(y) {
    message(y[["FlowJoFileName"]], ", ", format(as.numeric(y[["$TOT"]]), big.mark = ","), " evts, ", round(file.info(y[["FilePathUse"]])[["size"]]/1000/1000), " Mb")

    gs <- CytoML::flowjo_to_gatingset(CytoML::open_flowjo_xml(y[["wsp"]]),
                                      name = unique(y[["FlowJoGroup"]]),
                                      path = dirname(y[["FilePathUse"]]),
                                      subset = if (y[["equal_FileDirs"]]) {basename(y[["FilePathUse"]])} else {`$FIL` %in% y[["$FIL"]] & `$TOT` %in% y[["$TOT"]] & `$ETIM` %in% y[["$ETIM"]] & `$BTIM` %in% y[["$BTIM"]]},
                                      truncate_max_range = F,
                                      keywords = c("$FIL", "$ETIM", "$BTIM", "$TOT"),
                                      additional.keys = c("$TOT", "$ETIM", "$BTIM"))

    flowWorkspace::sampleNames(gs) <- y[["FlowJoFileName"]]
    attr(gs, "ws") <- y[["wsp"]]
    attr(gs, "FlowJoFileName") <- y[["FlowJoFileName"]]
    attr(gs, "FilePath") <- y[["FilePathUse"]]
    attr(gs, "downsample") <- y[["downsample"]]
    #rownames(y) <- paste(y[["$FIL"]], y[["$TOT"]], y[["$ETIM"]], y[["$BTIM"]], sep = "_")
    #flowWorkspace::sampleNames(gs) <- y[flowWorkspace::sampleNames(gs),"FileName"]

    if (remove_redundant_channels) {
      gs <- suppressMessages(flowWorkspace::gs_remove_redundant_channels(gs))
    }
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



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
  p <- unique(p)
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
  names(p_rev) <- p
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

get_ff <- function(gs,
                   return_untransformed = F,
                   return_transformed = T,
                   population,
                   seed = 42,
                   channels = NULL,
                   leverage_score_for_sampling = F,
                   return_ind_mat = F) {

  if ("downsample" %in% names(attributes(gs))) {
    downsample <- suppressWarnings(as.numeric(attr(gs, "downsample"))) # when min --> NA
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

  if (downsample != 1 && leverage_score_for_sampling && !requireNamespace("Seurat", quietly = T)) {
    utils::install.packages("Seurat")
  }

  if (!is.null(channels) && !leverage_score_for_sampling) {
    message("channels are only needed for leverage score aided sampling. since leverage_score_for_sampling = F channels are ignored.")
  }

  ## little unhandy to avoid memory is populated by ind_mat if not needed
  population2 <-
    population |>
    dplyr::filter(FileName == attr(gs, "FlowJoFileName")) |>
    dplyr::pull(Population)
  poppaths <- flowWorkspace::gh_get_pop_paths(x = gs)
  shortpaths <- stats::setNames(shortest_unique_path(poppaths), nm = poppaths)
  # alter population here by an optional leading forward slash in order to not make changes to ind_mat construction which could
  # have effects elsewhere. Maybe find a better solution some when
  # wsx_get_poppaths(x) return population paths without leading fwd slash
  population3 <- ifelse(population2 %in% shortpaths,
                        names(shortpaths)[which(shortpaths %in% population2)],
                        ifelse(grepl("^/", population2), population2, paste0("/", population2)))
  if (return_ind_mat) {
    poppaths_return <- poppaths
    shortpaths_return <- shortpaths
  } else {
    poppaths_return <- population3
    shortpaths_return <- shortpaths[population3]
  }
  ind_mat <- flowWorkspace::gh_pop_get_indices_mat(gh = gs, y = poppaths_return)
  attr(ind_mat, "short_names") <- shortpaths_return
  attr(ind_mat, "ws") <- attr(gs, "ws")
  attr(ind_mat, "FilePath") <- attr(gs, "FilePath")

  # inds: actually needed populations
  inds <- ind_mat[,population3,drop=F]


  if (return_untransformed && !return_transformed) {
    inverse_transform <- stats::setNames(T, "untransformed")
  } else if (!return_untransformed && return_transformed) {
    inverse_transform <- stats::setNames(F, "transformed")
  } else if (return_untransformed && return_transformed) {
    inverse_transform <- stats::setNames(c(F,T), c("transformed", "untransformed"))
  }

  ff <- tryCatch({
    lapply(inverse_transform, function(y) flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gh_pop_get_data(gs, inverse.transform = y)))
  }, error = function(err) {
    message(err)
    stop("Do FCS files contain a valid compensation matrix?")
  })

  attr(ff[[1]], "trafolist") <- flowWorkspace:::gs_get_transformlists(gs, inverse = F)
  attr(ff[[1]], "trafolistinv") <- flowWorkspace:::gs_get_transformlists(gs, inverse = T)

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
    s <- mapply(
      sampling_fun,
      inds_col = split_mat(inds, colnames(inds), byrow = F),
      lev_score = lev_scores,
      seed = seed,
      downsample = downsample,
      SIMPLIFY = F
    )
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
      ff[[i]] <- flowCore::Subset(ff[[i]], as.logical(x)) # inds may become a 1d array, maybe if inds has 1 column only (when there is 1 population only)
    }
    return(ff)
  })
  names(ff) <- attr(ind_mat, "short_names")[names(ff)]


  # ind_mat <- ind_mat[which(inds),,drop=F]
  # attr(ind_mat, "short_names") <- stats::setNames(shortest_unique_path(colnames(ind_mat)), nm = colnames(ind_mat))
  # attr(ind_mat, "ws") <- x$wsp
  # attr(ind_mat, "FilePath") <- x$FilePath

  flowWorkspace::gs_cleanup_temp(gs)
  if (!return_ind_mat) {
    ind_mat <- NULL
  }
  return(list(ff = ff, ind_mat = ind_mat))
}

sampling_fun <- function(inds_col, lev_score, seed = 42, downsample = 1) {

  if (length(unique(lev_score)) == 1) {
    prob <- NULL
  } else {
    # sampling with prob becomes very slow for large x
    # dplyr::slice_sample is not much faster
    message("Sampling with probabilities is slow for large x.")
    lev_score <- lev_score / sum(lev_score)
    prob <- lev_score
  }

  set.seed(seed)
  s <- sort(sample(x = which(inds_col),
                   size = ifelse(downsample < 1,
                                 ceiling(length(which(inds_col))*downsample),
                                 min(c(length(which(inds_col)), downsample))),
                   prob = prob))
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

  if (leverage_score_for_sampling && (!requireNamespace("Seurat", quietly = T))) { # || utils::packageDescription("Seurat")[["RemoteRef"]] != "feat/dictionary")) {
    if (!requireNamespace("remotes", quietly = T)) {
      utils::install.packages("remotes")
    }
    install.packages("Seurat")
    #remotes::install_github("satijalab/seurat", "feat/dictionary")
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
    s <- mapply(
      sampling_fun,
      inds_col = split_mat(inds, colnames(inds), byrow = F),
      lev_score = lev_scores,
      seed = seed,
      downsample = downsample
    )
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
                   dir = tempdir(),
                   merge_to_gs = T) {

  # split(x, (seq(nrow(x))-1) %/% split_size
  message("tempdir: ", dir, "\n")
  gs_list <- lapply(asplit(x,1), function(y) {
    message(y[["FlowJoFileName"]], ", ", format(as.numeric(y[["$TOT"]]), big.mark = ","), " evts, ", round(file.info(y[["FilePathUse"]])[["size"]]/1000/1000), " Mb")

    y[["equal_FileDirs"]] <- F # using filename as subset caused error
    gs <- CytoML::flowjo_to_gatingset(CytoML::open_flowjo_xml(y[["wsp"]]),
                                      name = unique(y[["FlowJoGroup"]]),
                                      path = dirname(y[["FilePathUse"]]),
                                      subset = if (y[["equal_FileDirs"]]) {basename(y[["FilePathUse"]])} else {`$FIL` %in% y[["$FIL"]] & `$TOT` %in% y[["$TOT"]] & `$ETIM` %in% y[["$ETIM"]] & `$BTIM` %in% y[["$BTIM"]]},
                                      truncate_max_range = F,
                                      keywords = c("$FIL", "$ETIM", "$BTIM", "$TOT"),
                                      additional.keys = c("$TOT", "$ETIM", "$BTIM"),
                                      backend_dir = dir)
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


get_kw_and_pars <- function(exprs,
                            ff = NULL,
                            keywrd = list(),
                            params = NULL,
                            insert_neutral_spill = T) {

  if (!is.null(ff)) {
    # provide the flowframe which was basis for creation of modified/extended exprs
    keywrd <- flowCore::keyword(ff)
    params <- flowCore::parameters(ff)
  }

  if (is.null(params)) {
    # creation of empty params:
    # params@data <- data.frame()
    # saveRDS(params, "/Users/vonskopnik/Documents/R_packages/fcexpr/inst/extdata/ff_params_empty.rds")
    params <- readRDS(system.file("extdata", "ff_params_empty.rds", package = "fcexpr"))
    desc <- NA
  } else {
    channelmarker <- stats::setNames(as.character(params@data$desc), as.character(params@data$name))
    desc <- unname(channelmarker[colnames(exprs)])
  }

  max <- as.integer(ceiling(apply(exprs, 2, max)))
  # newdata could be generated from scratch with exprs only except for desc from params@data
  params@data <- data.frame(name = colnames(exprs),
                            desc = desc,
                            minRange = as.integer(floor(apply(exprs, 2, min))),
                            maxRange = max,
                            range = max,
                            row.names = paste0("$P", seq(1,ncol(exprs))))
  #   | Keyword | Meaning                                                              |
  #   | ------- | -------------------------------------------------------------------- |
  #   | `$PnN`  | **Name** of the parameter (e.g., `FSC-A`, `CD4`, `FL1-H`).           |
  #   | `$PnS`  | **Short description** or label of the parameter (optional).          |
  #   | `$PnB`  | **Bit width** of the parameter data (e.g., `16`, `32`).              |
  #   | `$PnR`  | **Range** of values for the parameter (e.g., `1024` means 0–1023).   |
  #   | `$PnE`  | **Amplification type** and exponent (e.g., `0,0` or `4,1`).          |
  #   | `$PnG`  | **Gain** applied during data acquisition (optional).                 |
  #   | `$PnV`  | **Detector voltage** applied to the channel (optional, if recorded). |
  #

  for (z in rownames(params@data)) {
    keywrd[[paste0(z, "N")]] <- params@data[z,"name"]
    keywrd[[paste0(z, "S")]] <- params@data[z,"desc"]
    keywrd[[paste0(z, "R")]] <- as.character(params@data[z,"range"])
    keywrd[[paste0(z, "E")]] <- "0,0"
    if (is.null(keywrd[[paste0(z, "G")]])) {
      keywrd[[paste0(z, "G")]] <- "1"
    }
    if (is.null(keywrd[[paste0(z, "B")]])) {
      keywrd[[paste0(z, "B")]] <- "32"
    }
    if (is.null(keywrd[[paste0(z, "V")]])) {
      keywrd[[paste0(z, "V")]] <- "1"
    }
  }

  keywrd[["$PAR"]] <- as.character(ncol(exprs))
  keywrd[["$TOT"]] <- as.character(nrow(exprs))

  # correctly added by flowCore::write.FCS()
  keywrd[["$BEGINDATA"]] <- ""
  keywrd[["$ENDDATA"]] <- ""

  # create or fix spill mat
  # when a spill mat is there but channels were compensated already (Comp- as leading txt in channel name)
  # then creating a matched spill mat does not make sense as repeated application of spill by flowjo
  # would be a second round of compensation
  try(
    expr = {
      fluo_channels <- get_fluo_channels(colnames(exprs))
      # check cytof - somewhen
      if (is.null(keywrd[["SPILL"]]) || insert_neutral_spill) {
        # create new spill kw
        spill_mat <- diag(nrow = length(fluo_channels))
        colnames(spill_mat) <- fluo_channels
        keywrd[["SPILL"]] <- spill_mat
      } else if (!is.null(keywrd[["SPILL"]])) {
        matches <- match_channels_and_spill(spillcols = colnames(keywrd[["SPILL"]]), channels = fluo_channels, strict = F)
        if (ncol(keywrd[["SPILL"]]) != length(fluo_channels)) {
          # fix spill: remove cols and rows if respective channels were removed
          spillcol_to_rm <- which(matches$lv_dist > 1)
          keywrd[["SPILL"]] <- keywrd[["SPILL"]][-spillcol_to_rm,-spillcol_to_rm]
        }
        # fix colnames of spill keyword
        colnames(keywrd[["SPILL"]]) <- matches[which(matches$lv_dist <= 1),"channels"]
        # make spill neutral when compensated channels are found
        if (any(grepl("^Comp-", channels))) {
          spill_mat <- diag(nrow = length(keywrd[["SPILL"]]))
          colnames(spill_mat) <- colnames(keywrd[["SPILL"]])
          keywrd[["SPILL"]] <- spill_mat
        }
      }
    }, silent = T
  )

  return(list(keywrd = keywrd, params = params))
}



#fix_spill_kw <- function()

get_fluo_channels <- function(channels, ff = NULL) {

  if (is.null(ff)) {
    # use channels argument
    channels <- channels[which(!grepl("FSC|SSC|Time|HDR-T", channels, ignore.case = T))]
  } else {
    # use ff
    channels <- flowCore::pData(flowCore::parameters(ff))[["name"]]
    if (is.null(flowCore::keyword(ff)[["SPILL"]])) {
      message("inferring fluo channel names from flowCore::pData(flowCore::parameters(ff)).")
      channels <- channels[which(!grepl("FSC|SSC|Time|HDR-T", channels, ignore.case = T))]
    } else {
      spillcols <- colnames(flowCore::keyword(ff2)[["SPILL"]])
      matches <- match_spill_and_channels(spillcols = spillcols, channels = channels)
      if (is.null(matches)) {
        message("inferring fluo channel names from SPILL keyword.")
        channels <- matches[["channels"]]
      } else {
        message("inferring fluo channel names from flowCore::pData(flowCore::parameters(ff)).")
        channels <- channels[which(!grepl("FSC|SSC|Time|HDR-T", channels, ignore.case = T))]
      }

    }
  }
  return(channels)
}

match_spill_and_channels <- function(spillcols, channels) {
  # optionally removing leading comp was neccessary for adist to work properly
  # these two modification to spillcols and channels should allow a near-perfect match
  channels <- get_fluo_channels(channels)
  dists <- adist(gsub("/", "_", spillcols), gsub("^Comp-", "", channels))
  min_inds <- apply(dists, 1, which.min)
  min_dist <- apply(dists, 1,  min)
  if (length(unique(min_inds)) != length(min_inds)) {
    message("spillcols and channels not match unambiguously.")
    return(NULL)
  }
  return(data.frame(spillcols = spillcols[min_inds], channels = channels, lv_dist = min_dist))
}


match_channels_and_spill <- function(spillcols, channels, strict = T) {
  # optionally removing leading comp was neccessary for adist to work properly
  # these two modification to spillcols and channels should allow a near-perfect match
  channels <- get_fluo_channels(channels)
  dists <- adist(gsub("/", "_", spillcols), gsub("^Comp-", "", channels))
  min_inds <- apply(dists, 1, which.min)
  min_dist <- apply(dists, 1,  min)
  if (length(unique(min_inds)) != length(min_inds) && strict) {
    message("spillcols and channels not match unambiguously.")
    return(NULL)
  }
  return(data.frame(spillcols = spillcols, channels = channels[min_inds], lv_dist = min_dist))
}

get_new_kw_and_pars <- function(exprs,
                                new_kw,
                                new_desc = NULL,
                                new_pars) {


  if (is.null(new_desc)) {
    new_desc <- new_pars@data$desc
  }

  n <- nrow(new_pars)+1
  if (n > ncol(exprs)) {
    # happens when no new column was added to exprs
    # but one may have been overwritten
    n <- ncol(exprs)
  }
  new_p <- new_pars[1,]
  new_kw[["$PAR"]] <- as.character(ncol(exprs))

  ## new channels
  browser()
  if (n < ncol(exprs)) { # or <=
    # < needed by ff_simulate
    for (z in n:ncol(exprs)) {

      rownames(new_p) <- c(paste0("$P", z))
      new_pars <- BiocGenerics::combine(new_pars, new_p)
      new_p_name <- colnames(exprs)[z]
      new_p_desc <- new_desc[z]

      flowCore::pData(new_pars)$name[z] <- new_p_name
      flowCore::pData(new_pars)$desc[z] <- new_p_desc
      flowCore::pData(new_pars)$minRange[z] <- as.integer(round(min(exprs[, z])))
      flowCore::pData(new_pars)$maxRange[z] <- as.integer(round(max(exprs[, z])))
      flowCore::pData(new_pars)$range[z] <- as.integer(round(max(exprs[, z])))

      #### WRITE KEYWORD WITH 2 BRACKETS!! OTHERWISE FLOWJO DOES NOT READ IT #### ????????????
      new_kw[[paste0("$P", z, "N")]] <- new_p_name
      new_kw[[paste0("$P", z, "S")]] <- new_p_desc
      new_kw[[paste0("$P", z, "E")]] <- "0,0"
      new_kw[[paste0("$P", z, "G")]] <- "1"
      new_kw[[paste0("$P", z, "B")]] <- new_kw[["$P1B"]]
      new_kw[[paste0("$P", z, "R")]] <- as.integer(round(max(exprs[, z])))
      new_kw[[paste0("$P", z, "V")]] <- "1"
      new_kw[[paste0("$P", z, "DISPLAY")]] <- "LIN"
      new_kw[[paste0("flowCore_$P", z, "Rmin")]] <- as.integer(round(min(exprs[, z])))
      new_kw[[paste0("flowCore_$P", z, "Rmax")]] <- as.integer(round(max(exprs[, z])))
    }
  }

  ## previous channels
  for (z in 1:n) {
    flowCore::pData(new_pars)$minRange[z] <- as.integer(round(min(exprs[, z])))
    flowCore::pData(new_pars)$maxRange[z] <- as.integer(round(max(exprs[, z])))
    flowCore::pData(new_pars)$range[z] <- as.integer(round(max(exprs[, z])))
    if (!grepl("FSC|SSC", new_kw[paste0("$P", as.character(z), "N")])) {
      new_kw[[paste0("$P", as.character(z), "R")]] <- as.integer(round(max(exprs[, z])))
    }
    new_kw[[paste0("flowCore_$P", as.character(z), "Rmin")]] <- as.integer(round(min(exprs[, z])))
    new_kw[[paste0("flowCore_$P", as.character(z), "Rmax")]] <- as.integer(round(max(exprs[, z])))
  }

  new_kw[["$TOT"]] <- nrow(exprs)
  return(list(new_kw = new_kw, new_pars = new_pars))
}


.get.channels <- function(ff,
                          timeChannel = NULL,
                          channels = NULL) {
  if (!is.null(timeChannel)) {
    timeChannel <- trimws(timeChannel)
    timeChannel <- unlist(lapply(timeChannel, function(x) grep(paste0("^",x,"$"),
                                                               colnames(flowCore::exprs(ff)),
                                                               value = TRUE, ignore.case = TRUE)))
    if (all(is.na(timeChannel))) {
      message("None of timeChannels not found in exprs of flowFrame. Trying flowCore:::findTimeChannel(ff).")
      timeChannel <- flowCore:::findTimeChannel(ff)
    }
  } else {
    timeChannel <- flowCore:::findTimeChannel(ff)
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
      message(paste("These channels were not found in all flowFrames: ", notfound, collapse = ", "), ".")
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


.cluster_ordering <- function(ks) {
  if (methods::is(ks, "matrix")) {
    ks <- apply(ks, 2, function (x) {
      new_order <- stats::setNames(names(table(x)), nm = names(sort(table(x), decreasing = T)))
      return(as.numeric(new_order[as.character(x)]))
    })
  } else if (methods::is(ks, "integer")) {
    new_order <- stats::setNames(names(table(ks)), nm = names(sort(table(ks), decreasing = T)))
    ks <- as.numeric(new_order[as.character(ks)])
  } else {
    new_order <- stats::setNames(names(table(ks)), nm = names(sort(table(ks), decreasing = T)))
    ks <- unname(new_order[as.character(ks)])
  }

  return(ks)
}

.calc.pairwise.cluster.marker <- function(dat, cluster, levels = NULL, mc.cores = 1) {
  mc.cores <- min(mc.cores, parallel::detectCores() - 1)

  dat_split <- split_mat(x = dat, f = cluster)

  if (is.null(levels)) {
    levels <- sort(unique(cluster))
  }

  #all_pairs <- utils::combn(levels, 2, simplify = T)
  ## redundant calculation is easier for subsequent marker analysis; or one has to think about the logic doing it the non-redundant way
  all_pairs <- expand.grid(levels, levels)
  all_pairs <- all_pairs[which(all_pairs$Var1 != all_pairs$Var2),]
  all_pairs <- t(all_pairs)

  out <- dplyr::bind_rows(parallel::mcmapply(x = as.character(all_pairs[1,]),
                                             y = as.character(all_pairs[2,]),
                                             function(x,y) {

                                               out <-
                                                 presto::wilcoxauc(X = cbind(t(dat_split[[x]]),t(dat_split[[y]])), y = c(rep("y", length(which(as.character(cluster) == x))),
                                                                                                                         rep("z", length(which(as.character(cluster) == y))))) %>%
                                                 dplyr::filter(group == "y") %>%
                                                 dplyr::select(feature, pval) %>%
                                                 dplyr::rename("pvalue" = pval, "channel" = feature)

                                               out[,"mean_1"] <- round(matrixStats::colMeans2(dat_split[[x]]), 2)
                                               out[,"mean_2"] <- round(matrixStats::colMeans2(dat_split[[y]]), 2)
                                               out[,"mean_diff"] <- round(out[,"mean_1"] - out[,"mean_2"], 2)
                                               out[,"diptest_pvalue_1"] <- suppressWarnings(apply(dat_split[[x]], 2, function(z) diptest::dip.test(x = if(length(z) > 71999) {sample(z,71999)} else {z})[["p.value"]]))
                                               out[,"diptest_pvalue_2"] <- suppressWarnings(apply(dat_split[[y]], 2, function(z) diptest::dip.test(x = if(length(z) > 71999) {sample(z,71999)} else {z})[["p.value"]]))
                                               out[,"cluster_1"] <- as.character(x) #sapply(strsplit(out$cluster12, "_____"), "[", 1, simplify = T)
                                               out[,"cluster_2"] <- as.character(y) #sapply(strsplit(out$cluster12, "_____"), "[", 2, simplify = T)
                                               out[,"diff_sign"] <- ifelse(out[,"mean_diff"] == 0, "+/-", ifelse(out[,"mean_diff"] > 0, "+", "-"))

                                               #cluster_sizes <- utils::stack(table(cluster)) %>% dplyr::mutate(ind = as.character(ind))
                                               out <-
                                                 out %>%
                                                 #dplyr::left_join(cluster_sizes, by = c("cluster_1" = "ind")) %>%
                                                 #dplyr::rename("cluster_1_events" = "values") %>%
                                                 #dplyr::left_join(cluster_sizes, by = c("cluster_2" = "ind")) %>%
                                                 #dplyr::rename("cluster_2_events" = "values") %>%
                                                 #, cluster_1_events, cluster_2_events
                                                 dplyr::select(channel, cluster_1, cluster_2, pvalue, mean_1, mean_2, mean_diff, diff_sign, diptest_pvalue_1, diptest_pvalue_2) %>%
                                                 dplyr::arrange(pvalue)

                                               'tryCatch({
      out <- suppressWarnings(matrixTests::col_wilcoxon_twosample(dat_split[[x]],
                                                                  dat_split[[y]])) %>%
        dplyr::select(pvalue) %>%
        tibble::rownames_to_column("channel")
    }, error=function(err) {
      message("Ran matrixTests::col_wilcoxon_twosample with error in level : ", x, " vs ", y, ": ")
      message(err)
      message("Trying presto::wilcoxauc.")
      out <-
        presto::wilcoxauc(cbind(t(dat_split[[x]]),t(dat_split[[y]])), c(rep("y", length(which(as.character(cluster) == x))),
                                                                        rep("z", length(which(as.character(cluster) == y))))) %>%
        dplyr::filter(group == "y") %>%
        dplyr::select(feature, pval) %>%
        dplyr::rename("pvalue" = pval, "channel" = feature)
    }, warning = function(war) {
      message("Ran matrixTests::col_wilcoxon_twosample with warning in level : ", x, " vs ", y, ": ")
      message(war)
      message("Trying presto::wilcoxauc.")
      out <-
        presto::wilcoxauc(cbind(t(dat_split[[x]]),t(dat_split[[y]])), c(rep("y", length(which(as.character(cluster) == x))),
                                                                        rep("z", length(which(as.character(cluster) == y))))) %>%
        dplyr::filter(group == "y") %>%
        dplyr::select(feature, pval) %>%
        dplyr::rename("pvalue" = pval, "channel" = feature)
    })'


                                               return(out)
                                             }, mc.cores = mc.cores, SIMPLIFY = F))

  out$cluster_1 <- factor(out$cluster_1, levels = levels)
  out$cluster_2 <- factor(out$cluster_2, levels = levels)
  return(out)
}

.calc.global.cluster.marker <- function(dat, cluster, levels = NULL, mc.cores = 1) {

  #method = c("presto", "matrixTests")
  #method <- match.arg(method, c("presto", "matrixTests"))

  mc.cores <- min(mc.cores, parallel::detectCores() - 1)

  # cluster is ident for each row in dat
  if (nrow(dat) != length(cluster)) {
    stop("nrow(dat) != length(cluster)")
  }
  if (is.null(levels)) {
    levels <- sort(unique(cluster))
  }

  #levels_out <<- levels
  #dat_out <<- dat
  #cluster_out <<- cluster

  ## try matrixStats first and on error run presto which requires transposation though
  ## matrixStats caused an error once (Integer Overflow with large matrices (approx. 1e6 cells as initial input))
  out <- dplyr::bind_rows(parallel::mclapply(levels, function(x) {
    ' tryCatch({
      out <- matrixTests::col_wilcoxon_twosample(dat[which(cluster == x),,drop = F],
                                                 dat[which(cluster != x),,drop = F])  |> |>
        dplyr::select(pvalue) |>
        tibble::rownames_to_column("channel")
    }, error=function(err) {
      message("Ran matrixTests::col_wilcoxon_twosample with error in level : ", x, ": ")
      message(err)
      message("Trying presto::wilcoxauc.")
      out <-
        presto::wilcoxauc(cbind(t(dat[which(cluster == x),,drop = F]),
                                t(dat[which(cluster != x),,drop = F])), c(rep("y", length(which(cluster == x))),
                                                                          rep("z", length(which(cluster != x))))) |>
        dplyr::filter(group == "y") |>
        dplyr::select(feature, pval) |>
        dplyr::rename("pvalue" = pval, "channel" = feature)
    }, warning = function(war) {
      message("Ran matrixTests::col_wilcoxon_twosample with warning in level : ", x, ": ")
      message(war)
      message("Trying presto::wilcoxauc.")
      out <-
        presto::wilcoxauc(cbind(t(dat[which(cluster == x),,drop = F]),
                                t(dat[which(cluster != x),,drop = F])), c(rep("y", length(which(cluster == x))),
                                                                          rep("z", length(which(cluster != x))))) |>
        dplyr::filter(group == "y") |>
        dplyr::select(feature, pval) |>
        dplyr::rename("pvalue" = pval, "channel" = feature)
    })'


    out <-
      presto::wilcoxauc(X = cbind(t(dat[which(cluster == x),,drop = F]),
                                  t(dat[which(cluster != x),,drop = F])), y = c(rep("y", length(which(cluster == x))),
                                                                                rep("z", length(which(cluster != x))))) |>
      dplyr::filter(group == "y") |>
      dplyr::select(feature, pval) |>
      dplyr::rename("pvalue" = pval, "channel" = feature)

    out[,"mean_cluster"] <- round(matrixStats::colMeans2(dat[which(cluster == x),,drop = F]), 2)
    out[,"mean_notcluster"] <- round(matrixStats::colMeans2(dat[which(cluster != x),,drop = F]), 2)
    out[,"mean_diff"] <- round(out[,"mean_cluster"] - out[,"mean_notcluster"], 2)
    out[,"diptest_pvalue_cluster"] <- suppressWarnings(apply(dat[which(cluster == x),,drop = F], 2, function(z) diptest::dip.test(x = if(length(z) > 71999) {sample(z,71999)} else {z})[["p.value"]]))
    out[,"diptest_pvalue_notcluster"] <- suppressWarnings(apply(dat[which(cluster != x),,drop = F], 2, function(z) diptest::dip.test(x = if(length(z) > 71999) {sample(z,71999)} else {z})[["p.value"]]))
    out[,"cluster"] <- as.character(x)
    out[,"diff_sign"] <- ifelse(out[,"mean_diff"] == 0, "+/-", ifelse(out[,"mean_diff"] > 0, "+", "-"))

    out <-
      out |>
      #dplyr::left_join(utils::stack(table(cluster)) |> dplyr::mutate(ind = as.character(ind)), by = c("cluster" = "ind")) |>
      #dplyr::rename("cluster_events" = "values") |>
      #cluster_eventsdf <-
      dplyr::select(channel, cluster, pvalue, mean_cluster, mean_notcluster, mean_diff, diff_sign, diptest_pvalue_cluster, diptest_pvalue_notcluster) |>
      dplyr::arrange(pvalue)
    return(out)
  }, mc.cores = mc.cores))

  out <-
    out |>
    dplyr::group_by(channel) |>
    dplyr::mutate(mean_cluster_scale = as.vector(scale(mean_cluster))) |>
    dplyr::ungroup() |>
    dplyr::mutate(cluster = factor(cluster, levels = levels))

  out <- heatmap_ordering(df = out,
                          feature_order = "custom",
                          group_order = "hclust")

  return(out)
}



random_FIL <- function(prefix = "Specimen_001_", seed = 1) {
  set.seed(seed)
  num1 <- sample(1:9,1)
  num2 <- sample(1:999,1)
  num2 <- sprintf("%03d", num2)
  FIL <- paste0(prefix, num1, "_", num2)
  return(FIL)
}

random_BTIM_ETIM_DATE <- function(seed = 1) {

  set.seed(seed)
  # random begin and end time
  base_time <- as.POSIXct(Sys.time())
  random_seconds <- runif(1, min = 0, max = 86400) # 24h
  BTIM <- base_time + random_seconds
  random_seconds <- as.integer(runif(1, min = 60, max = 600))
  ETIM <- BTIM + random_seconds

  # random date
  start_date <- as.Date("2015-01-01")
  end_date <- as.Date(Sys.Date())
  random_date <- as.Date(runif(1, min = as.numeric(start_date), max = as.numeric(end_date)), origin = "1970-01-01")
  DATE <- format(random_date, "%d-%b-%Y")
  DATE <- toupper(random_date)


  return(c(format(BTIM, "%H:%M:%S"), format(ETIM, "%H:%M:%S"), DATE, random_seconds))
}

random_OP <- function(seed = 1) {
  set.seed(seed)
  return(sample(fcexpr:::random_operators, 1))
}




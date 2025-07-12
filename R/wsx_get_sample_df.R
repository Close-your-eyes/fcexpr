#' Get information on samples as in ws and local FCS files that match those in ws
#'
#' This function is used within wsp_get_ff, wsp_get_gs and wsp_get_indices. Usage
#' by the user may be followed by fcexpr::fcs_files_df_get_ff to import the FCS files
#' as flow frames easily. With ind_mats from fcexpr::wsp_get_indices the flowframes
#' can be subsetted by gated populations. See wsp_get_ff, wsp_get_gs and wsp_get_indices
#' for explanation of parameters (function arguments).
#'
#' @param ws vector of paths to flowjo workspace
#' @param groups flowjo groups to consider, one vector for all ws or a list with one vector for each ws
#' @param invert_groups consider all other flowjo groups found in ws but groups
#' @param samples samples (FCS filenames) to consider, one vector for all ws or a list with one vector for each ws
#' @param invert_samples consider all other samples found in ws but samples
#' @param FCS.file.folder one path to the top local FCS file folder or a vector
#'
#' @return data frame with sample data
#' @export
#'
#' @examples
wsx_get_sample_df <- function(ws,
                              groups = NULL,
                              invert_groups = F,
                              samples = NULL,
                              invert_samples = F,
                              FCS.file.folder) {

  checked_in <- check_wsp_for_samples_and_groups(
    wsp = ws,
    groups = groups,
    samples = samples,
    FCS.file.folder = FCS.file.folder,
    invert_groups = invert_groups,
    invert_samples = invert_samples
  )

  smpl <- get_smpl_df(wsp = ws,
                      groups = checked_in[["groups"]],
                      invert_groups = F,
                      samples = checked_in[["samples"]],
                      invert_samples = F,
                      FCS.file.folder = checked_in[["FCS.file.folder"]])

  return(smpl)
}

check_wsp_for_samples_and_groups <- function(wsp,
                                             samples = NULL,
                                             groups = NULL,
                                             FCS.file.folder = NULL,
                                             invert_groups = F,
                                             invert_samples = F) {

  # check invert samples and groups herein

  if (missing(wsp) || !methods::is(wsp, "character")) {stop("Please provide a vector of paths to flowjo workspaces.")}

  wsp_data <- purrr::map(wsp, wsx_get_fcs_paths, split = F, filter_AllSamples = T)

  if (!is.null(samples)) {
    if (length(samples) == 0) {
      stop("No samples provided.")
    }

    if (methods::is(samples, "list") && length(samples) != length(wsp) && length(samples) != 1) {
      stop("list of samples has to have the same length as wsp or 1. Alternatively pass a vector of samples to use for all workspaces.")
    } else if ((methods::is(samples, "list") && length(samples) == 1) || !methods::is(samples, "list")) {
      if (!methods::is(samples, "list")) {
        samples <- rep(list(samples), length(wsp)) # samples is vector
      } else {
        samples <- rep(samples, length(wsp))
      }
      if (invert_samples) {
        samples <- purrr::map2(wsp_data, samples, ~ .x$FileName[which(!.x$FileName %in% .y)])
      } else {
        not_found_samples <- unlist(purrr::map2(wsp_data, samples, ~ .y[which(!.y %in% .x$FileName)]))
        if (length(not_found_samples) > 0) {
          message("samples not found: ", paste(not_found_samples, collapse = ", "))
        }
        samples <- purrr::map2(wsp_data, samples, ~ .x$FileName[which(.x$FileName %in% .y)])
        # empty elements are checked for below
      }
    }
  }

  if (!is.null(groups)) {
    if (methods::is(groups, "list") && length(groups) != length(wsp) && length(groups) != 1) {
      stop("list of groups has to have the same length as wsp or 1. Alternatively pass a vector of groups to use for all workspaces.")
    } else if ((methods::is(groups, "list") && length(groups) == 1) || !methods::is(groups, "list")) {
      if (!methods::is(groups, "list")) {
        groups <- rep(list(groups), length(wsp)) # groups is vector
      } else {
        groups <- rep(groups, length(wsp))
      }
      if (invert_groups) {
        groups <- purrr::map2(wsp_data, groups, ~ .x$FlowJoGroup[which(!.x$FlowJoGroup %in% .y)])
      } else {
        not_found_groups <- unlist(purrr::map2(wsp_data, groups, ~ .y[which(!.y %in% .x$FlowJoGroup)]))
        if (length(not_found_groups) > 0) {
          message("groups not found: ", paste(not_found_groups, collapse = ", "))
        }
        groups <- purrr::map2(wsp_data, groups, ~ .x$FlowJoGroup[which(.x$FlowJoGroup %in% .y)])
        # empty elements are checked for below
      }
    }
  } else if (is.null(samples)) {
    groups <- purrr::map(wsp_data, ~unique(.x$FlowJoGroup))
  } else {
    # pick groups based on samples
    groups <- purrr::map2(wsp_data, samples, function(x,y) {
      dplyr::filter(x, FileName %in% y) |> dplyr::distinct(FlowJoGroup) |> dplyr::pull(FlowJoGroup)
    })
  }

  # group now has values in any case, so samples can be populated
  if (is.null(samples)) {
    samples <- purrr::map2(wsp_data, groups, function(x,y) {
      dplyr::filter(x, FlowJoGroup %in% y) |> dplyr::distinct(FileName) |> dplyr::pull(FileName)
    })
  }

  if (!is.null(FCS.file.folder)) {
    if (any(!dir.exists(FCS.file.folder))) {stop(paste(FCS.file.folder[which(!dir.exists(FCS.file.folder))], collapse = ", "), " not found.")}
    if (length(FCS.file.folder) == 1) {FCS.file.folder <- rep(FCS.file.folder, length(wsp))}
    if (length(FCS.file.folder) != length(wsp)) {stop("FCS.file.folder has to have the same length as wsp or 1.")}
  }

  # check zero length vector in list
  # remove respective groups and wsp
  if (any(lengths(samples) == 0)) {
    inds <- lengths(samples) == 0
    message("samples not found or left from wsp(s): ", paste(wsp[inds], collapse = ","))
    wsp <- wsp[-inds]
    groups <- groups[-inds]
    FCS.file.folder <- FCS.file.folder[-inds]
  }
  if (any(lengths(groups) == 0)) {
    inds <- lengths(groups) == 0
    message("groups not found or left from wsp(s): ", paste(wsp[inds], collapse = ","))
    wsp <- wsp[-inds]
    samples <- samples[-inds]
    FCS.file.folder <- FCS.file.folder[-inds]
  }

  if (length(wsp) == 0) {
    stop("No workspace or samples left or found for import.")
  }

  return(list(groups = groups, samples = samples, FCS.file.folder = FCS.file.folder))
}

get_smpl_df <- function(wsp,
                        groups,
                        invert_groups,
                        samples,
                        invert_samples,
                        FCS.file.folder) {

  fcskeywords <- c("$DATE", "$BTIM", "$ETIM", "$TOT", "$FIL")

  # FlowJoGroup
  # one wsp at a time
  smpl <- do.call(rbind, lapply(seq_along(wsp), function(x) {
    y <- wsx_get_fcs_paths(wsp[x], split = F)
    names(y)[3:4] <- c("FlowJoFilePath", "FlowJoFileName")
    y$FlowJoFilePath <- URLdecode(y$FlowJoFilePath) # nice !!!
    y$wsp <- wsp[x]
    y$FlowJoFileName <- basename(y$FlowJoFilePath) # redo this after urldecode

    # get keywords from wsp and generate fcs identities
    kwlist <- wsx_get_keywords(wsp[x], return = c("data.frame", "vector"), keywords = fcskeywords)
    fcs_ident <-
      stack(fcexpr:::get_fcs_identities(kwlist[["vec"]], allow_duplicates = T)) |>
      dplyr::rename("identity" = values, "FlowJoFileName" = ind) |>
      dplyr::mutate(FlowJoFileName = as.character(FlowJoFileName))
    y <-
      y |>
      dplyr::left_join(fcs_ident, by = "FlowJoFileName") |>
      dplyr::left_join(kwlist[["df2"]], by = "FlowJoFileName")

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

    if (!is.null(FCS.file.folder)) {
      local_fcs_files <- list.files(path = FCS.file.folder[x], recursive = T, full.names = T, pattern = "\\.fcs$", ignore.case = T)
      kwl <- fcs_get_keywords(file_paths = local_fcs_files,
                              keywords = fcskeywords,
                              return = "vector")

      local_fcs_files_df <-
        stack(fcexpr:::get_fcs_identities(kwl = kwl, allow_duplicates = T)) |>
        dplyr::rename("identity" = values, "LocalFilePath" = ind) |>
        dplyr::filter(identity %in% y$identity) |>  # only consider fcs files that are present in flowjo
        dplyr::mutate(LocalFilePath = as.character(LocalFilePath)) |>
        dplyr::mutate(LocalFileName = basename(LocalFilePath))

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
      y <- dplyr::left_join(y, local_fcs_files_df, by = c("identity" = "identity")) # , "FlowJoFileName" = "LocalFileName"
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

  smpl <-
    smpl |>
    dplyr::mutate(sampleID = as.numeric(sampleID)) |>
    dplyr::group_by(wsp) |>
    dplyr::arrange(sampleID, .by_group = T) |>
    dplyr::ungroup()

  return(as.data.frame(smpl))
}

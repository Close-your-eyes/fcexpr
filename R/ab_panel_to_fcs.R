#' Write antibody related information into a fcs file
#'
#' This writes channel descriptions and optional keywords to fcs files.
#'
#' @param sampledescription data frame of sampledescription
#' @param FileNames vector of which fcs files to consider
#' @param FCS.file.folder path to the folder containing the fcs files specified in FileNames
#' @param channel_conjugate_match_file path to an xlsx file holding information about fluorochromes matched to channels (better leave as it is for now)
#' @param AbCalcFile_col column name in sampledescription indicating the file containing the antibody panel calculation
#' @param AbCalcSheet_col column name in sampledescription indicating the respective sheet name in AbCalcFile
#' @param conjugate_to_desc alter/update the channel description of fcs files with stained molecule and optionally the fluorochrome
#' @param other_keywords column names in AbCalcSheet of which keywords to write to fcs files
#' @param AbCalcFile_folder path to the folder containing the AbCalcFile; if NULL then AbCalcFile_col must contain the full, absolute path of AbCalcFile
#' @param clear_previous clear all previous entries (channel descriptions and keywords) in fcs files?
#' @param ignore_duplicate_ag whether or not to ignore duplicate antigen entries in antibody calculation sheet
#' @param machine name of the Flow Cytometer the FCS were acquired on. Only necessary when channels cannot be matched unambiguously.
#' @param manual_df manually provide a data frame with AbCalcFile_col, AbCalcSheet_col and FilePath columns for every FCS file
#' to write channel information to; sampledescription, FileNames, FCS.file.folder are not required when !is.null(manual_df)
#'
#' @return no return, but updated fcs files on disk
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#'fcexpr::ab_panel_to_fcs(sampledescription = sd,
#'FileNames = sd$FileName,
#'FCS.file.folder = file.path(wd, "FCS_files"),
#'other_keywords = c("Isotype", "Clone", "totalDF", "Vendor", "Cat", "Lot"))
#' }
ab_panel_to_fcs <- function(sampledescription,
                            FileNames,
                            FCS.file.folder,
                            channel_conjugate_match_file = system.file("extdata", "channel_conjugate_matches.xlsx", package = "fcexpr"),
                            AbCalcFile_col = "AbCalcFile",
                            AbCalcSheet_col = "AbCalcSheet",
                            AbCalcFile_folder = file.path(dirname(FCS.file.folder), "Protocols"),
                            conjugate_to_desc = T,
                            other_keywords = c(),
                            clear_previous = F,
                            ignore_duplicate_ag = F,
                            machine = NULL,
                            manual_df = NULL) {

  # how to handle non-fluorochrome conjugates?
  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("CytoML", quietly = T)){
    BiocManager::install("CytoML")
  }
  if (!requireNamespace("flowWorkspace", quietly = T)){
    BiocManager::install("flowWorkspace")
  }

  ccm <- .check.and.get.ccm(ccm = channel_conjugate_match_file)

  if (is.null(manual_df)) {
    if (missing(FCS.file.folder)) {
      stop("FCS.file.folder missing. Please provide.")
    }
    if (missing(FileNames)) {
      stop("FileNames missing. Please provide.")
    }

    ## check sd for columns
    sd <- as.data.frame(sampledescription, stringsAsFactors = F)
    if (!"FileName" %in% names(sd)) {
      stop("Column 'FileName' has to exist in sampledescription.")
    }
    if (!AbCalcFile_col %in% names(sd)) {
      stop(AbCalcFile_col, " not found in sampledescription columns.")
    }
    if (!AbCalcSheet_col %in% names(sd)) {
      stop(AbCalcSheet_col, " not found in sampledescription columns.")
    }

    # select cols
    sd <- sd[which(sd[,"FileName"] %in% FileNames), c("FileName", AbCalcFile_col, AbCalcSheet_col)]
    if (nrow(sd) == 0) {
      stop("Non of FileNames found in sampledescription.")
    }

    # check for empty cells
    if (length(unique(c(which(is.na(sd[,AbCalcFile_col])), which(trimws(sd[,AbCalcFile_col]) == "")),
                      c(which(is.na(sd[,AbCalcSheet_col])), which(trimws(sd[,AbCalcSheet_col]) == "")))) > 0) {
      message(paste0("AbCalcFile_col and/or AbCalcSheet_col empty or missing for: ", paste(sd[unique(c(which(is.na(sd[,AbCalcFile_col])), which(trimws(sd[,AbCalcFile_col]) == "")),
                                                                                                     c(which(is.na(sd[,AbCalcSheet_col])), which(trimws(sd[,AbCalcSheet_col]) == ""))), "FileName"], collapse = ", ")))

      sd <- sd[unique(c(which(!is.na(sd[,AbCalcFile_col])), which(trimws(sd[,AbCalcFile_col]) != "")),
                      c(which(!is.na(sd[,AbCalcSheet_col])), which(trimws(sd[,AbCalcSheet_col]) != ""))),]
    }

    sd[,"FilePath"] <- sapply(sd[,"FileName"], function(x) list.files(FCS.file.folder, pattern = stringr::str_replace_all(x, c("\\+" = "\\\\+", "\\." = "\\\\.",
                                                                                                                               "\\|" = "\\\\|", "\\(" = "\\\\(",
                                                                                                                               "\\)" = "\\\\)", "\\[" = "\\\\[",
                                                                                                                               "\\{" = "\\\\{", "\\*" = "\\\\*",
                                                                                                                               "\\?" = "\\\\?")), recursive = T, full.names = T))
  } else {
    if (!AbCalcFile_col %in% names(manual_df) || !AbCalcSheet_col %in% names(manual_df) || !"FilePath" %in% names(manual_df)) {
      stop("AbCalcFile_col, AbCalcSheet_col, FilePaths have to be in names of manual_df")
    }
    sd <- manual_df
  }

  sd[,"config"] <- dplyr::coalesce(sapply(flowCore::read.FCSheader(sd[,"FilePath"], emptyValue = F), "[", "CYTOMETER CONFIG NAME"),
                                   sapply(flowCore::read.FCSheader(sd[,"FilePath"], emptyValue = F), "[", "$CYT"))

  if (any(is.na(sd[,"config"]))) {
    stop("Config keyword could not be retrieved from all FCS files. Another FCS keyword is needed - fix the function.")
  }
  #glob2rx(x)

  # sd grouped by sheet and cytometers
  sd <-
    sd %>%
    dplyr::group_by(!!rlang::sym(AbCalcFile_col), !!rlang::sym(AbCalcSheet_col), config) %>%
    #dplyr::summarise(FileNames = list(FileName), .groups = "drop") %>%
    dplyr::summarise(FilePaths = list(FilePath), .groups = "drop") %>%
    as.data.frame()

  # loop though ab info files
  out <- lapply(split(sd, 1:nrow(sd)), function(x) {
    if (!is.null(AbCalcFile_folder)) {
      file <- file.path(AbCalcFile_folder, x[,AbCalcFile_col])
    } else {
      file <- x[,AbCalcFile_col]
    }

    if (!file.exists(file)) {
      warning("File '", file, "' not found.")
    } else {
      if (!x[,AbCalcSheet_col] %in% openxlsx::getSheetNames(file)) {
        warning("Sheet '", x[,AbCalcSheet_col], "' not found in ", file, ".")
      } else {
        sh <- openxlsx::read.xlsx(file, sheet = x[,AbCalcSheet_col])
        # check columns
        if (any(!c("Conjugate", "Antigen") %in% names(sh))) {
          warning("Conjugate and/or Antigen columns not found in Sheet '", x[,AbCalcSheet_col], "' in '", x[,AbCalcFile_col], "'.")
        } else {
          # check other_keywords
          if (any(!other_keywords %in% names(sh))) {
            message(paste(other_keywords[which(!other_keywords %in% names(sh))], collapse = ", "), " columns not found in AbCalcSheet. Those will not be written to FCS files.")
            other_keywords <- intersect(other_keywords, names(sh))
          }

          # do this after checking for LiveDeadMarker - this will allow to have no antibody entered w/o error
          #sh <- sh[which(!is.na(sh[,"Antigen"])),unique(c("Antigen", "Conjugate", other_keywords, "LiveDeadMarker"))]
          if (is.na(sh[1,"LiveDeadMarker"])) {
            message("No LiveDeadMarker entered.")
          } else {
            ## allow more than one LiveDeadMarker, separated by comma
            rrr <- nrow(sh)+1
            for (z in trimws(strsplit(sh[1,"LiveDeadMarker"], ",")[[1]])) {
              sh[rrr,"Conjugate"] <- z
              sh[rrr,"Antigen"] <- "LiveDead"
              rrr <- rrr + 1
            }
            sh[1,"LiveDeadMarker"] <- NA
          }

          # check for "channel" in other_keywords
          sh <- sh[which(!is.na(sh[,"Antigen"])),unique(c("Antigen", "Conjugate", other_keywords, "LiveDeadMarker"))]

          sh <- sh[,colSums(is.na(sh))<nrow(sh)] # also removes LiveDeadMarker column
          other_keywords <- intersect(other_keywords, names(sh))
          # exclude LiveDead from duplicate check as it may be in 2 channels
          if (length(unique(sh[which(sh[,"Antigen"] != "LiveDead"),"Antigen"])) != length(sh[which(sh[,"Antigen"] != "LiveDead"),"Antigen"]) && !ignore_duplicate_ag) {
            warning("Duplicate Antigen found in Ab.calc.sheet (",  x[,AbCalcSheet_col], "). ",  "Please check. To ignore this warning and write information in FCS files anyway, pass 'ignore_duplicate_ag = TRUE'.")
          } else {
            ## read FCS here, then call the ccm function
            # for loop as sh is modified by the first FCS file
            for (j in x[,"FilePaths"][[1]]) {
              'fcs.path <- list.files(FCS.file.folder, pattern = stringr::str_replace_all(j, c("\\+" = "\\\\+", "\\." = "\\\\.",
                                                                                              "\\|" = "\\\\|", "\\(" = "\\\\(",
                                                                                              "\\)" = "\\\\)", "\\[" = "\\\\[",
                                                                                              "\\{" = "\\\\{", "\\*" = "\\\\*",
                                                                                              "\\?" = "\\\\?")), recursive = T, full.names = T)'

              ff <- flowCore::read.FCS(j, truncate_max_range = F, emptyValue = F)
              ff <- .check_comp_mat(ff)
              channels <- flowCore::pData(flowCore::parameters(ff))[,"name"]
              channels.inv <- stats::setNames(names(channels),channels)

              # special treatment on first index
              if (j == x[,"FilePaths"][[1]][1]) {
                matches <- conjugate_to_channel(conjugates = sh[,"Conjugate"], channels = channels, channel_conjugate_match_file = ccm, machine = machine)
                if (any(duplicated(matches))) {
                  message("At least on channel used by more than one marker.")
                }
                sh[,"channel"] <- matches[match(sh[,"Conjugate"], names(matches))]
                if (any(is.na(matches))) {
                  warning("Conjugates not found in ccm: ", paste(sh[which(is.na(sh[,"channel"])),"Conjugate"], collapse = ","))
                  sh <- sh[which(!is.na(sh[,"channel"])),]
                }
                sh[,"Antigen.Conjugate"] <- paste(sh[,"Antigen"], sh[,"Conjugate"], sep = "-")
                # collapse multiple information for same channel
                collapse.fun <- function(x) paste0(x, collapse = ", ")
                sh <-
                  sh %>%
                  dplyr::group_by(channel) %>%
                  dplyr::summarise(dplyr::across(c(Antigen, Conjugate, Antigen.Conjugate, dplyr::all_of(other_keywords)), collapse.fun), .groups = "drop")
                sh <- as.data.frame(sh)

                # same order
                sh <- sh[order(match(sh[,"channel"], flowCore::pData(flowCore::parameters(ff))[,"name"])),]
                sh[,"channel.name"] <- stringr::str_extract(channels.inv[sh[,"channel"]], "P[:digit:]{1,}")
                print(sh)
              }

              if (clear_previous && any(!is.na(flowCore::pData(flowCore::parameters(ff))[,"desc"]))) {
                flowCore::pData(flowCore::parameters(ff))[,"desc"] <- NA
              }
              if (conjugate_to_desc) {
                flowCore::pData(flowCore::parameters(ff))[,"desc"][which(flowCore::pData(flowCore::parameters(ff))[,"name"] %in% sh[,"channel"])] <- sh[,"Antigen.Conjugate"]
              } else {
                flowCore::pData(flowCore::parameters(ff))[,"desc"][which(flowCore::pData(flowCore::parameters(ff))[,"name"] %in% sh[,"channel"])] <- sh[,"Antigen"]
              }

              ## to do: clear previous keywords

              # other meta data about the antibody used
              for (m in other_keywords) {
                for (k in 1:nrow(sh)) {
                  if (!is.null(sh[k,m]) && !is.na(sh[k,m]) && !trimws(sh[k,m]) %in% c("NA", "", "-")) {
                    flowCore::keyword(ff)[paste0(sh[k,"channel.name"],"_",m)] <- sh[k,m]
                  }
                }
              }

              # save fcs file (overwrite original)
              flowCore::write.FCS(ff, j)
              message(j)
            }
          }
        }
      }
    }
    return(NULL)
  })
}


conjugate_to_channel <- function(conjugates,
                                 channels = NULL,
                                 machine = NULL,
                                 channel_conjugate_match_file = system.file("extdata", "channel_conjugate_matches.xlsx", package = "fcexpr")) {


  ccm <- .check.and.get.ccm(ccm = channel_conjugate_match_file)
  if (!is.null(machine)) {
    if (length(machine) > 1) {
      stop("Please provide only one machine name.")
    }
    if (!machine %in% unique(ccm$machine)) {
      stop("Machine not found in ccm matches. Please choose from: ", paste(unique(ccm$machine), collapse = ", "))
    }
    ccm <- ccm[which(ccm$machine == machine),]
  }
  if (is.null(channels) && !is.null(machine)) {
    channels <- unique(ccm$channel)
  }

  matches <- ccm[intersect(which(ccm[,"channel"] %in% channels), which(tolower(make.names(ccm[,"Conjugate"])) %in% tolower(make.names(conjugates)))), ]

  matches_grouped <- dplyr::group_by(matches, Conjugate) %>% dplyr::summarise(n_ch = nlevels(as.factor(channel)))

  if (any(matches_grouped$n_ch > 1)) {
    print(matches)
    stop("Conjugate could not be matched to channel unambiguously. Try to provide the correct name of the machine: ", paste(unique(matches$machine), collapse = ", "), ".")
  }

  matches <- unique(matches[,c("Conjugate","channel")])
  matches <- stats::setNames(matches[,"channel"], nm = matches[,"Conjugate"])
  names(matches) <- unique(unlist(sapply(names(matches), function(x) grep(paste0("^",x,"$"), conjugates, value = T, ignore.case = T), simplify = T, USE.NAMES = F)))
  if (any(duplicated(names(matches)))) {
    warning("Duplicate matches in ccm. Did you pass channels from multiple Flow Cytometers at once?")
  }
  matches <- matches[order(match(names(matches), conjugates))]
  return(matches)
}

.check.and.get.ccm <- function(ccm) {

  if (is.character(ccm)) {
    if (length(ccm) != 1) {
      stop("ccm has to be length 1 (only one file path).")
    }
    if (!file.exists(ccm)) {
      stop(paste0(ccm, " not found."))
    }
    if ("matches" %in% openxlsx::getSheetNames(ccm)) {
      ccm <- openxlsx::read.xlsx(ccm, sheet = "matches")
    } else {
      ccm <- openxlsx::read.xlsx(ccm)
    }
  } else if (is.data.frame(ccm)) {
    ## tibble or data.frame; is.data.frame(tibble) is TRUE
    ccm <- as.data.frame(ccm)
  } else {
    stop("ccm has to be a file.path or data frame.")
  }

  if (!all(grepl(c("conjugate|channel|machine"), names(ccm), ignore.case = T))) {
    stop("Column names of the ccm have to be conjugte, channel, machine.")
  }
  for (i in c("Conjugte", "channel", "machine")) {
    names(ccm)[which(grepl(i, names(ccm), ignore.case = T))] <- i
  }
  return(ccm)
}

.check_comp_mat <- function(ff) {

  spill_ind <- which(flowCore:::.spillover_pattern %in% names(flowCore::keyword(ff)))
  if (length(spill_ind) == 0) {
    message("No compensation matrix found in flowFrame. Will try create one.")
    spill_key <- flowCore:::.spillover_pattern[1]
  } else if (length(spill_ind) > 1) {
    message("More than one spillover keyword detected. Will use the first one: ", flowCore:::.spillover_pattern[spill_ind])
    spill_key <- flowCore:::.spillover_pattern[min(which(flowCore:::.spillover_pattern %in% names(flowCore::keyword(ff))))]
  } else {
    spill_key <- flowCore:::.spillover_pattern[min(which(flowCore:::.spillover_pattern %in% names(flowCore::keyword(ff))))]
  }

  if (is.null(flowCore::keyword(ff)[[spill_key]]) || all(dim(flowCore::keyword(ff)[[spill_key]]) == 0)) {
    comp_channels <- flowCore::pData(flowCore::parameters(ff))$name
    comp_channels <- comp_channels[which(!grepl("time|hdr|fsc|ssc", comp_channels, ignore.case = T))]
    comp_mat <- diag(x=length(comp_channels))
    colnames(comp_mat) <- comp_channels
    flowCore::keyword(ff)[[spill_key]] <- comp_mat
  }
  return(ff)
}


#' Title
#'
#' @param sampledescription
#' @param FileNames
#' @param FCS.file.folder
#' @param channel_conjugate_match_file
#' @param AbCalcFile_col
#' @param AbCalcSheet_col
#' @param conjugate_to_desc
#' @param other_keywords
#' @param clear_previous
#'
#' @return
#' @export
#'
#' @examples
ab_info_to_fcs <- function (sampledescription,
                            FileNames,
                            FCS.file.folder,
                            channel_conjugate_match_file = system.file("extdata", "channel_conjugate_matches.xlsx", package = "fcexpr"),
                            AbCalcFile_col = "AbCalcFile",
                            AbCalcSheet_col = "AbCalcSheet",
                            conjugate_to_desc = T,
                            other_keywords = c("Isotype", "Clone", "totalDF", "Vendor", "Cat", "Lot"),
                            clear_previous = T) {


  if (!"BiocManager" %in% rownames(utils::installed.packages())) {install.packages("BiocManager")}
  if (!"flowCore" %in% rownames(utils::installed.packages())) {BiocManager::install("flowCore")}

  ## check sd for columns

  sd <- as.data.frame(sampledescription, stringsAsFactors = F)
  if (!AbCalcFile_col %in% names(sd)) {
    stop(paste0(AbCalcFile_col, " not found in sampledescription columns."))
  }
  if (!AbCalcSheet_col %in% names(sd)) {
    stop(paste0(AbCalcSheet_col, " not found in sampledescription columns."))
  }

  if (missing(FCS.file.folder)) {
    stop("FCS.file.folder missing. Please provide.")
  }
  if (missing(FileNames)) {
    stop("FileNames missing. Please provide.")
  }

  ccm <- openxlsx::read.xlsx(channel_conjugate_match_file, sheet = "matches")
  if (!identical(names(ccm), c("Conjugate", "channel", "machine"))) {
    # avoid case sensitivity
    stop("Column names of the channel_conjugate_match_file have to be Conjugte, channel, machine.")
  }
  s <- sd[which(sd[,1] %in% FileNames), ]
  if (length(s) == 0) {
    stop("Non of FileNames found in sampledescription column 1.")
  }

  for (i in 1:nrow(s)) {
    fcs.path <- list.files(FCS.file.folder, pattern = s[i,1], recursive = T, full.names = T)
    ff <- flowCore::read.FCS(fcs.path, truncate_max_range = F, emptyValue = F)
    channels <- flowCore::pData(flowCore::parameters(ff))[,"name"]
    channels.inv <- stats::setNames(names(channels),channels)

    if (!is.null(s[i,AbCalcFile_col]) && !is.na(s[i,AbCalcFile_col])) {
      file <- file.path(basename(FCS.file.folder), "Protocols", s[i,AbCalcFile_col])
      if (!file.exists(file)) {
        stop(paste0(file, " not found. ", s[i,AbCalcFile_col], " is expected to be in a folder named Protocols in ", basename(FCS.file.folder), "."))
      }

      sh <- openxlsx::read.xlsx(file, sheet = s[i,AbCalcSheet_col])
      sh <- sh[which(!is.na(sh[,"Antigen"])),]

      if (!is.na(sh[1,"Live.Dead.Marker"])) {
        sh <- dplyr::bind_rows(sh, data.frame(Conjugate = sh[1,"Live.Dead.Marker"]))
        sh[nrow(sh),which(names(sh) != "Conjugate")] <- NA
      }

      if (length(unique(sh[,"Antigen"])) != length(sh[,"Antigen"])) {
        stop(paste0("Duplicate Antigen found in Ab.calc.sheet (",  s[i,AbCalcSheet_col], "). ",  "Please check."))
      }

      if (!"Conjugate" %in% names(sh)) {
        print("Conjugate column not found in AbCalcSheet")
      }

      if (length(!other_keywords %in% names(sh)) > 0) {
        print(paste0(paste(other_keywords, collapse = ", "), " not found in AbCalcSheet."))
      }

      # extract relevant rows from the reference database (matched channels and fluorochromes)
      # check for matching channels and notify of missing ones
      selection <- ccm[intersect(which(ccm[,"channel"] %in% channels), which(tolower(make.names(ccm[,"Conjugate"])) %in% tolower(make.names(sh[,"Conjugate"])))), ]
      selection <- unique(selection[,c(1,2)])

      # join; create a helper column. not all mistypings will be rescued.
      sh[,"conj.join"] <- tolower(make.names(sh[,"Conjugate"]))
      selection[,"Conjugate"] <- tolower(make.names(selection[,"Conjugate"]))
      sh <- merge(sh, selection, by.x = "conj.join", by.y = "Conjugate")
      sh <- sh[which(!names(sh) %in% "conj.join")]
      sh[,"Antigen.Conjugate"] <- gsub("NA-", "", paste(sh[,"Antigen"], sh[,"Conjugate"], sep = "-"))

      # collapse multiple information for same channels
      collapse.fun <- function(x) {
        paste0(x, collapse = ", ")
      }
      other_keywords <- base::intersect(other_keywords, names(sh))
      sh <-
        sh %>%
        dplyr::group_by(channel) %>%
        dplyr::summarise(dplyr::across(c(Antigen, Conjugate, Antigen.Conjugate, other_keywords), collapse.fun), .groups = "drop")
      sh <- as.data.frame(sh)

      # same order
      sh <- sh[order(match(sh[,"channel"], flowCore::pData(flowCore::parameters(ff))[,"name"])),]
      sh[,"channel.name"] <- stringr::str_extract(channels.inv[sh[,"channel"]], "P[:digit:]{1,}")
      print(sh[,c("Antigen", "Conjugate", "channel")])

      # write into flowframe
      sh.ag.conj <- as.data.frame(sh %>% dplyr::distinct(Antigen, Conjugate, Antigen.Conjugate, channel, channel.name))
      sh.ag.conj <- sh.ag.conj[order(match(sh.ag.conj[,"channel"], flowCore::pData(flowCore::parameters(ff))[,"name"])),]

      if (clear_previous && any(!is.na(flowCore::pData(flowCore::parameters(ff))[,"desc"]))) {
        flowCore::pData(flowCore::parameters(ff))[,"desc"] <- NA
      }
      if (conjugate_to_desc) {
        flowCore::pData(flowCore::parameters(ff))[,"desc"][which(flowCore::pData(flowCore::parameters(ff))[,"name"] %in% selection[,"channel"])] <- sh.ag.conj[,"Antigen.Conjugate"]
      } else {
        flowCore::pData(flowCore::parameters(ff))[,"desc"][which(flowCore::pData(flowCore::parameters(ff))[,"name"] %in% selection[,"channel"])] <- sh.ag.conj[,"Antigen"]
      }

      # other meta data about the antibody used
      for (j in other_keywords) {
        if (j %in% names(sh)) {
          for (k in 1:nrow(sh)) {
            if (!is.null(sh[k,j]) && !is.na(sh[k,j]) && sh[k,j] != "") {
              flowCore::keyword(ff)[paste0(sh[k,"channel.name"],j)] <- sh[k,j]
            }
          }
        }
      }

      # save fcs file (overwrite original)
      flowCore::write.FCS(ff, fcs.path)
      print(fcs.path)

    } else {
      print(paste0("No Ab.calc.file entered for file ", s[i, 1], "."))
    }
  }
}

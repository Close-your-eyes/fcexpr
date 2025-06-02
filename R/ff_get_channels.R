#' Get channel names from flow frame
#'
#' @param ff flow frame
#' @param timeChannel time channel name
#' @param channels channels to return/match, can be channel name, desc or
#' keyword (see examples); if NULL all channels are returned
#' @param part_match allow partial matching of channels identifiers
#' @param rm_scatter remove scatter channels from return
#' @param rm_time remove time channel
#' @param rm_wo_desc remove all channels without desc (empty channels)
#'
#' @return character vector
#' @export
#'
#' @examples
#' \dontrun{
#' ff_get_channels(ff, channels = c("LP502", "PE"))
#' ff_get_channels(ff, channels = c("P4N", "P6N"))
#' ff_get_channels(ff, channels = c("CD3", "CD8"))
#' }
ff_get_channels <- function(ff,
                            timeChannel = NULL,
                            channels = NULL,
                            part_match = T,
                            rm_scatter = T,
                            rm_time = T,
                            rm_wo_desc = F,
                            replace_NA_desc = F,
                            return = c("vector", "data.frame")) {

  stopifnot("only one flow frame" = length(ff) == 1)

  return <- rlang::arg_match(return)

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
    if (rm_time) {
      channels <- channels[which(!channels %in% timeChannel)]
    }
    if (rm_scatter) {
      channels <- channels[which(!grepl("FSC|SSC", channels, ignore.case = T))]
    }
    if (rm_wo_desc) {
      channels <- channels[!is.na(names(channels))]
    }
    if (length(channels) == 0) {
      message("no channels left after filtering.")
    }
  } else {
    channels <- trimws(channels)
    inds123 <- NULL
    if (part_match) {
      inds1 <- unique(unlist(purrr::map(channels, ~which(grepl(.x, flowCore::pData(flowCore::parameters(ff))$name, ignore.case = T)))))
      inds2 <- unique(unlist(purrr::map(channels, ~which(grepl(.x, flowCore::pData(flowCore::parameters(ff))$desc, ignore.case = T)))))
      inds3 <- unique(unlist(purrr::map(channels, ~which(grepl(.x, names(flowCore::pData(flowCore::parameters(ff))$name), ignore.case = T)))))
      inds123 <- unique(c(inds1, inds2, inds3))
      channels <- flowCore::pData(flowCore::parameters(ff))$name[inds123]
    }
    inds <- unique(c(which(names(flowCore::pData(flowCore::parameters(ff))$name) %in% channels),
                     which(flowCore::pData(flowCore::parameters(ff))$name %in% channels),
                     which(flowCore::pData(flowCore::parameters(ff))$desc %in% channels)))
    inds <- unique(c(inds123, inds))
    inds <- sort(inds)
    if (!part_match) {
      notfound <- channels[intersect(which(!channels %in% flowCore::pData(flowCore::parameters(ff))$name),
                                     which(!channels %in% flowCore::pData(flowCore::parameters(ff))$desc))]
      if (length(notfound) > 0) {
        message(paste("These channels were not found in all flowFrames: ", notfound, collapse = ", "), ".")
      }
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
    message("No channels matched to those in the flowFrame.")
  }
  if (replace_NA_desc) {
    names(channels)[which(is.na(names(channels)))] <- channels[which(is.na(names(channels)))]
  }

  if (return == "vector") {
    channels <- stats::setNames(as.character(channels), names(channels))
  } else if (return == "data.frame") {
    channels <- stack(channels)
    channels$ind <- as.character(channels$ind)
    names(channels) <- c("name", "desc")
  }

  return(channels)
}


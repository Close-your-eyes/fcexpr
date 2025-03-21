#' Apply compensation matrix to flowframe
#'
#' In the end, it is flowCore::compensate only. But since channel names from flowframe
#' and column names of comp mats may not match, this function was made to handle that.
#' The problem arises due to conversion of special characters. Channel names may contain
#' a forward slash since this is how filters of a flow cytometer are described.
#' E.g. 760/30 for 760 nm +/- 15 nm. However, in FCS keywords the forward slash gets
#' converted to an underscore. Channel names in a flowframe maintain the slash. So,
#' flowCore::compensate complains.
#'
#' @param ff a flowframe
#' @param compmat a compensation matrix; matrix with column names
#' @param match_channels try to match channels names from ff with column names
#' of compmat
#'
#' @return compensated flowframe
#' @export
#'
#' @examples
#' \dontrun{
#' compmats <- fcexpr::wsx_get_compmats(ws = path/to/flowjo/workspace.wsp)
#' ff <- flowCore::read.FCS("path/to/fcsfile.fcs", truncate_max_range = F)
#' correct_mat <- compmmats |> dplyr::fiilter(FileName == fcsfile.fcs) |> dplyr::pull(mat)
#' ff_compensated_exprs_slot <- ff_apply_compensation(ff, correct_mat)
#' }
ff_apply_compensation <- function(ff,
                                  compmat,
                                  match_channels = T) {
  if (is.list(compmat) && length(compmat) == 1) {
    compmat <- compmat[[1]]
  } else if (is.list(compmat)) {
    stop("When compmat is a list, its length has to be 1.")
  }
  if (!is.matrix(compmat)) {
    stop("compmat has to be a matrix.")
  }
  if (is.null(colnames(compmat))) {
    stop("compmat needs colnames.")
  }
  if (!methods::is(ff, "flowFrame")) {
    stop("ff should be a flowframe.")
  }

  nc <- ncol(compmat)
  channels <- ff@parameters@data[["name"]]
  fluo_channels <- channels[which(!grepl("FSC|SSC|Time", channels, ignore.case = T))]
  nfluo <- length(fluo_channels)
  if (nfluo != nc) {
    message("column number and number of fluorescence channels in ff may not match. correct compmat?")
    message(colnames(compmat))
    message(unname(fluo_channels))
  }

  if (match_channels) {
    for (n in seq_along(colnames(compmat))) {
      if (!colnames(compmat)[n] %in% channels) {
        message("compmat col ", colnames(compmat)[n], " not found in channels of ff.")
        best_match_channel <- channels[which.min(adist(colnames(compmat)[n], channels))]
        message("Will use ", best_match_channel, " as best matching channel.")
        colnames(compmat)[n] <- best_match_channel
      }
    }
  }

  tryCatch(
    {
      ff <- flowCore::compensate(x = ff, spillover = compmat)
    },
    error = function(cond) {
      message("Error in flowCore::compensate. ff unchanged.")
    }
  )
  return(ff)
}

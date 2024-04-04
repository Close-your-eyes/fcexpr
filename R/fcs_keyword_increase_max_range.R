#' Title
#'
#' @param fcs_file_path
#' @param channel_reg_expr
#' @param value
#'
#' @return
#' @export
#'
#' @examples
fcs_keyword_increase_max_range <- function(fcs_file_path,
                                           channel_reg_expr = NULL,
                                           value = 5000) {



  ff <- flowCore::read.FCS(fcs_file_path, truncate_max_range = F, emptyValue = F)

  ## TODO
  ## problem with %P15R names etc. - no matching with grepl
  if (is.null(channel_reg_expr)) {
   stop("channel_reg_expr missing.")
  }
'  if (is.null(channel_reg_expr)) {
    message("no channel_reg_expr provided. Will modify all keyword with this pattern: P[[:digit:]]{1,}R")
    channel_reg_expr <- grep("P[[:digit:]]{1,}R", flowCore::keyword(ff), value = T)
    channel_reg_expr <- channel_reg_expr[which(!grepl("min$", channel_reg_expr))]
    message(paste(channel_reg_expr, collapse = ","))
  }
'

  keys_to_mod <- flowCore::keyword(ff)[which(grepl(channel_reg_expr, names(flowCore::keyword(ff))))]
  keys_to_mod <- keys_to_mod[which(!grepl("min$", names(keys_to_mod)))]

  flowCore::keyword(ff)[names(keys_to_mod)] <- value
  flowCore::write.FCS(ff, fcs_file_path)
  message(fcs_file_path)
}

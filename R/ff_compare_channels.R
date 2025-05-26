#' Compare channel meta data of flow frames for equality
#'
#' @param ff_list list of flow frames
#' @param channels relevant channels
#' @param strict accept differences like different desc
#'
#' @return NULL when there are no issues
#' @export
#'
#' @examples
ff_compare_channels <- function(ff_list,
                                channels = NULL,
                                strict = T) {


  if (strict) {
    out <- purrr::map(.x = ff_list, .f = ~flowCore::parameters(.x)$name)
    out <- purrr::pmap_lgl(out, ~length(unique(.x)) == 1)
    if (!all(out)) {
      warning("Channels of flowFrames do not have the same names. This cannot be handled. Will return data frame(s) of channel names.")
      return(purrr::map(.x = ff_list, .f = ~purrr::map(.x = .x, .f = ~flowCore::parameters(.x)$name)))
    }

    out <- purrr::map(.x = ff_list, .f = ~stats::setNames(flowCore::pData(flowCore::parameters(.x))[,"desc"], flowCore::pData(flowCore::parameters(.x))[,"name"]))
    out <- purrr::pmap_lgl(out, ~length(unique(.x)) == 1)
    if (!all(out)) {
      warning("Channel description are not equal across flowFrames.")
    }
    return(NULL)
  }

  if (!strict) {

    out <- purrr::map(.x = ff_list, .f = ~flowCore::parameters(.x)$name)
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
    out <- purrr::map(.x = ff_list, .f = ~stats::setNames(flowCore::pData(flowCore::parameters(.x))[,"desc"], flowCore::pData(flowCore::parameters(.x))[,"name"]))
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

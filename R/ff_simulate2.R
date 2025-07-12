#' Create a random flow frame
#'
#' @param exprs list of expression matrices (or data frames); if NULL then
#' fcexpr::make_hd_pops is used to create data
#' @param m number of exprs to create with make_hd_pops --> number of flow frames
#' @param ... arguments to make_hd_pops and fcexpr::create_pop
#'
#' @return list of flow frames
#' @export
#'
#' @examples
#' # populations with different properties, pre-made
#' pops1 <- make_hd_pops(npop = 3, range_sd_fct = 20, range_modes = 2, rm_ann_cols = T)
#' pops2 <- make_hd_pops(npop = 3, range_sd_fct = 3:6, range_modes = 3, rm_ann_cols = T)
#' ff_list1 <- ff_simulate2(exprs = dplyr::bind_rows(pops1, pops2))
#' # have populations made within function
#' ff_list <- ff_simulate2(npop = 2, dims = 10, m = 3)
ff_simulate2 <- function(exprs = NULL,
                         m = 1,
                         annotate_channels = F,
                         ...) {

  seed <- 42
  dots <- list(...)
  if (is.null(exprs)) {
    exprs <- replicate(m, do.call(make_hd_pops, args = c(list(rm_ann_cols = T), dots)), simplify = F)
  } else{
    if (!is.list(exprs)) {
      exprs <- list(exprs)
    }
  }

  if (annotate_channels) {
    # make channel names and fluorochromes
    allcol <- unique(unlist(purrr::map(exprs, names)))
    ncols <- length(allcol)
    # system.file("extdata", "channel_conjugate_matches.tsv", package = "fcexpr")
    ccmatch <- read.delim(system.file("extdata", "channel_conjugate_matches.tsv", package = "fcexpr"))
    ccmatchsum <-
      ccmatch |>
      dplyr::distinct(channel, machine) |>
      dplyr::group_by(machine) |>
      dplyr::filter(!grepl("Macs", machine)) |>
      dplyr::tally() |>
      dplyr::filter(n >= ncols)
    if (nrow(ccmatchsum) > 0) {
      # channels from one machine
      machine <-
        dplyr::slice_sample(ccmatchsum, n = 1) |>
        dplyr::pull(machine)
      channels <-
        ccmatch |>
        dplyr::filter(type == "marker") |>
        dplyr::filter(machine == !!machine) |>
        dplyr::group_by(channel) |>
        dplyr::slice_sample(n = 1) |>
        dplyr::ungroup() |>
        dplyr::slice_sample(n = ncols)|>
        dplyr::arrange(channel)
      channels$col <- allcol
    } else {
      # arbitrary channels
      channels <-
        ccmatch |>
        dplyr::filter(type == "marker") |>
        dplyr::filter(machine %in% c("Animal", "Wayne", "Symphony")) |>
        dplyr::group_by(fluorochrome) |>
        dplyr::slice_sample(n = 1) |>
        dplyr::ungroup() |>
        dplyr::distinct(channel, .keep_all = T) |>
        dplyr::slice_sample(n = ncols) |>
        dplyr::arrange(channel)
      channels$col <- allcol
    }
  }

  ff_list <- purrr::map(exprs, function(y) {
    y <- as.matrix(y)

    if (annotate_channels) {
      descs <- paste0(colnames(y), "-", channels[match(channels$col, colnames(y)),"fluorochrome", drop=T])
      colnames(y) <- channels[match(channels$col, colnames(y)),"channel", drop=T]
    }

    BED <- fcexpr:::random_BTIM_ETIM_DATE(seed = seed + round(rnorm(1,100,50),0))
    y <- cbind(y, matrix(seq(runif(1,0.2,0.9), as.numeric(BED[4]), length.out = nrow(y)), ncol = 1,
                         dimnames = list(NULL, "Time")))


    kw_par <- fcexpr:::get_kw_and_pars(exprs = y)
    if (annotate_channels) {
      kw_par[["params"]]@data[["desc"]] <- c(descs, "Time") # time channel
    }

    # generate a few random keywords
    kw_par[["keywrd"]][["$OP"]] <- fcexpr:::random_OP(seed = seed + round(rnorm(1,100,50),0))
    kw_par[["keywrd"]][["$FIL"]] <- fcexpr:::random_FIL(seed = seed + round(rnorm(1,100,50),0))
    kw_par[["keywrd"]][["$BTIM"]] <- BED[1]
    kw_par[["keywrd"]][["$ETIM"]] <- BED[2]
    kw_par[["keywrd"]][["$DATE"]] <- BED[3]
    kw_par[["keywrd"]][["$CYT"]] <- "in silico"
    kw_par[["keywrd"]][["FCSversion"]] <- "3"


    y <- methods::new(
      "flowFrame",
      exprs = y,
      parameters = kw_par[["params"]],
      description = kw_par[["keywrd"]]
    )

    return(y)
  })

  if (is.null(names(ff_list))) {
    datetimes <- purrr::map_chr(ff_list, ~paste(unlist(flowCore::keyword(.x, keyword = c("$DATE", "$BTIM"))), collapse = " "))
    ff_list <- ff_list[order(as.POSIXct(datetimes))]
    nums <- as.character(seq(1,length(ff_list)))
    names(ff_list) <- paste0("sample_", stringr::str_pad(nums, width = max(nchar(nums)), pad = "0"))
  }

  return(ff_list)
}

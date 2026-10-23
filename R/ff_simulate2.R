#' Create simulated flow frames
#'
#' Generate simulated marker measurements with make_hd_pops(), or convert
#' supplied expression data, into a list of flowCore flowFrame objects. Each
#' frame receives an additional Time channel and simulated acquisition metadata.
#'
#' @param exprs NULL to generate measurements internally, a numeric expression
#'   matrix for one sample, or a list of numeric expression matrices or data
#'   frames for multiple samples. Rows represent cells and columns represent
#'   marker measurements. Supply marker columns only, without cell, group,
#'   pop, or an existing Time column. Wrap a single data frame in a list:
#'   a bare data frame is treated as a list of columns by this implementation.
#'   For channel annotation, use a list of data frames with identical column
#'   names in the same order; the annotation code reads names() from each
#'   sample, which does not retrieve column names from ordinary matrices.
#' @param m Positive integer specifying the number of samples to generate when
#'   exprs is NULL. Each sample is generated independently. Ignored when exprs
#'   is supplied; the supplied data determine the number of output frames.
#' @param annotate_channels Logical; if TRUE, randomly assign marker channel
#'   names and fluorochromes using the packaged channel-conjugate lookup table.
#'   Channel descriptions combine the original marker name and fluorochrome.
#'   The selected mapping is shared across samples. If FALSE, retain the input
#'   marker column names. The additional channel is named Time in either case.
#' @param ... Arguments passed to make_hd_pops() and, through it, create_pop()
#'   when exprs is NULL. These include npop, dims, range_n, range_modes,
#'   range_mean, range_sd_fct, and mean_prob_vec_args. Do not supply rm_ann_cols:
#'   it is set to TRUE internally. These arguments are ignored when exprs is
#'   supplied.
#'
#' @details
#' The generation procedure is:
#' \enumerate{
#'   \item If exprs is NULL, call make_hd_pops() once per requested sample.
#'     Each sample contains npop independently generated populations with dims
#'     marker dimensions. Population, cell, and mode-group annotation columns
#'     are removed before conversion. Otherwise, use the supplied measurements.
#'   \item If annotate_channels is TRUE, select marker channel names and
#'     fluorochromes from the packaged lookup table. The code first attempts
#'     to select channels from one instrument, falling back to a pooled set
#'     when no instrument meets its channel-count criterion.
#'   \item Convert each sample to a matrix and append a Time channel with
#'     evenly spaced values over a simulated acquisition duration. Time follows
#'     the existing cell order; rows are not shuffled.
#'   \item Build parameter metadata and simulated keywords for the operator,
#'     file name, acquisition date, and start and end times. The cytometer
#'     keyword is set to "in silico" and the FCS version to "3". Construct
#'     one flowFrame per sample.
#'   \item If the resulting list is unnamed, sort frames by simulated
#'     acquisition date and start time, then assign sequential sample names.
#'     If the supplied list has names, preserve its names and order.
#' }
#'
#' npop counts generated populations per sample, whereas range_modes controls
#' the number of modes within each marker dimension of a population. Cells
#' assigned to the same mode are drawn around a common mean with a common
#' standard deviation, creating clusters. Shared cell blocks coordinate mode
#' membership across dimensions. Larger range_sd_fct values make modes tighter.
#' Modes and populations may overlap, and a single population may contain
#' several visible clusters, so npop does not guarantee a visible cluster count.
#' See create_pop() for the sampling procedure and mean-weighting behavior.
#'
#' range_n controls cells per population rather than cells per sample. dims
#' controls marker dimensions only: Time is always added as another column,
#' independently of annotate_channels. For example, npop set to 5, dims set to
#' 2, m set to 3, and range_n set to 500 generate three frames, each with 2500
#' cells and three columns: two markers and Time. Population labels are not
#' retained in the resulting expression matrices.
#'
#' @return A list of flowCore flowFrame objects, including when only one sample
#'   is requested. Each frame contains the marker measurements, an additional
#'   Time column, parameter metadata, and simulated acquisition keywords.
#'   Named input lists retain their names. Unnamed inputs receive names such
#'   as sample_1, with zero padding when needed for the number of samples.
#' @seealso make_hd_pops(), create_pop(), prob_vec()
#' @export
#'
#' @examples
#' \dontrun{
#' # Three samples with five populations and two marker dimensions each
#' ff_list <- ff_simulate2(npop = 5, dims = 2, m = 3, range_n = 500,
#'                         range_modes = 1, annotate_channels = TRUE)
#' length(ff_list)
#' # Two annotated marker channels plus Time
#' colnames(flowCore::exprs(ff_list[[1]]))
#'
#' # Combine populations with different properties into one supplied sample
#' pops1 <- make_hd_pops(npop = 3, range_sd_fct = 20, range_modes = 2,
#'                       rm_ann_cols = TRUE)
#' pops2 <- make_hd_pops(npop = 3, range_sd_fct = 3:6, range_modes = 3,
#'                       rm_ann_cols = TRUE)
#' supplied <- list(mixed_sample = dplyr::bind_rows(pops1, pops2))
#' ff_list1 <- ff_simulate2(exprs = supplied, annotate_channels = TRUE)
#'
#' # A single numeric matrix can be supplied without an outer list
#' ff_list2 <- ff_simulate2(exprs = as.matrix(pops1))
#' }
ff_simulate2 <- function(exprs = NULL,
                         m = 1,
                         annotate_channels = F,
                         ...) {
  fcexpr:::.ensure_packages(c("flowCore"))

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
    ccmatch <- utils::read.delim(system.file("extdata", "channel_conjugate_matches.tsv", package = "fcexpr"))
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

    BED <- fcexpr:::random_BTIM_ETIM_DATE(seed = seed + round(stats::rnorm(1,100,50),0))
    y <- cbind(y, matrix(seq(stats::runif(1,0.2,0.9), as.numeric(BED[4]), length.out = nrow(y)), ncol = 1,
                         dimnames = list(NULL, "Time")))


    kw_par <- fcexpr:::get_kw_and_pars(exprs = y)
    if (annotate_channels) {
      kw_par[["params"]]@data[["desc"]] <- c(descs, "Time") # time channel
    }

    # generate a few random keywords
    kw_par[["keywrd"]][["$OP"]] <- fcexpr:::random_OP(seed = seed + round(stats::rnorm(1,100,50),0))
    kw_par[["keywrd"]][["$FIL"]] <- fcexpr:::random_FIL(seed = seed + round(stats::rnorm(1,100,50),0))
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

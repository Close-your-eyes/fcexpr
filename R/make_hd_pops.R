#' Create a high dimensional space with multiple populations
#'
#' Generate multiple populations with create_pop() and combine their cells in
#' one data frame. Each population is generated independently using the same
#' supplied parameter settings and receives a population identifier.
#'
#' @param npop Positive integer specifying how many populations to generate.
#'   This counts calls to create_pop(), not necessarily visibly separated
#'   clusters. A single population can contain several modes, and different
#'   populations can overlap.
#' @param popindexstart Integer label for the first population. Subsequent
#'   populations receive consecutive labels, stored as character values in pop.
#'   Use a new starting label when appending independently generated populations
#'   to an existing data frame.
#' @param rm_ann_cols Logical; if TRUE, retain only columns whose names contain
#'   protein, removing the cell, group, and pop annotation columns. The default
#'   FALSE retains these annotations.
#' @param ... Arguments passed to create_pop() for every population, including
#'   range_n, dims, range_mean, range_sd_fct, range_modes, and
#'   mean_prob_vec_args. See create_pop() for their sampling behavior.
#'
#' @details
#' The generation procedure is:
#' \enumerate{
#'   \item Call create_pop() independently once for each requested population.
#'     Each call samples its own cell count from range_n and generates dims
#'     marker dimensions with the requested candidate mode counts, means,
#'     and spread factors.
#'   \item Assign consecutive population labels starting at popindexstart and
#'     combine all populations by rows. Store the labels in pop.
#'   \item If rm_ann_cols is TRUE, remove annotation columns and return only
#'     the simulated protein columns.
#' }
#'
#' range_n controls cells per population, not the total number of cells in the
#' combined data. For example, npop set to 5 and range_n set to 500 request
#' five populations of 500 cells each, giving 2500 cells in total.
#'
#' dims controls the number of marker columns. When annotations are retained,
#' cell, group, and pop are additional columns. Cell identifiers are local to
#' each population; use pop together with cell to identify a cell in the
#' combined data. group encodes mode membership across dimensions within a
#' population and is not a population identifier.
#'
#' Setting range_modes to 1 requests one normal mode per dimension within each
#' population. Allowing larger mode counts can produce several visible clusters
#' within one population. Modes and independently generated populations can
#' also overlap, so npop does not guarantee the number of visible clusters.
#' Mode means are sampled separately for each dimension, and the default
#' U-shaped sampling weights favor both low and high candidate means.
#' Consequently, some dimensions may appear close to zero on a wide axis.
#'
#' This function creates one combined data frame and does not add a Time
#' channel. When ff_simulate2() generates data internally, it calls this
#' function once per sample, removes the annotations, and adds a Time channel.
#' Its m argument controls the number of samples; each resulting flow frame
#' contains the requested marker dimensions plus Time.
#'
#' @return A data frame containing the combined cells from all populations.
#'   Marker columns are named protein1, protein2, and so on. Unless rm_ann_cols
#'   is TRUE, the result also contains cell identifiers, group labels encoding
#'   mode memberships, and the population labels in pop.
#' @seealso create_pop(), ff_simulate2(), prob_vec()
#' @export
#'
#' @examples
#' # Five populations, each with 500 cells and one mode per marker dimension
#' df <- make_hd_pops(npop = 5, range_n = 500, dims = 2, range_modes = 1)
#' table(df$pop)
#' # Return only the two marker columns
#' proteins <- make_hd_pops(npop = 1, dims = 2, range_modes = 1,
#'                          rm_ann_cols = TRUE)
#' # default settings
#' df1 <- make_hd_pops()
#' # plot with ggplot
#' ggplot(df1, aes(protein1, protein2)) +
#'   geom_point(aes(color = pop)) +
#'   scale_x_log10() +
#'   scale_y_log10() +
#'   brathering::scale_color_custom()
#' ggplot(df1, aes(protein3, protein4)) +
#'   geom_point(aes(color = pop)) +
#'   scale_x_log10() +
#'   scale_y_log10() +
#'   brathering::scale_color_custom()
#' # dimension reduction with umap
#' df11 <- cbind(df1, as.data.frame(uwot::umap(as.matrix(df1[,1:5]), verbose = T)))
#' brathering::plot2(df11[,c("V1", "V2", "pop", "group")], cex = 3, color = "pop", legend = NULL)
#' brathering::plot3d(df1, transform = "log10")
#' #add pops with other properties
#' df2 <- dplyr::bind_rows(df1, make_hd_pops(popindexstart = 6,
#'                                           range_sd_fct = 2,
#'                                           range_modes = 1))
#' ggplot(df2, aes(protein1, protein2)) +
#'   geom_point(aes(color = pop)) +
#'   scale_x_log10() +
#'   scale_y_log10() +
#'   brathering::scale_color_custom()
make_hd_pops <- function(npop = 5,
                         popindexstart = 1,
                         rm_ann_cols = F,
                         ...) {
  dots <- list(...)
  df <- replicate(npop, do.call(create_pop, args = dots), simplify = F)
  names(df) <- popindexstart:(npop+popindexstart-1)
  df <- dplyr::bind_rows(df, .id = "pop")
  df <- dplyr::relocate(df, pop, .after = names(df)[ncol(df)])
  if (rm_ann_cols) {
    df <- df[,which(grepl("protein", colnames(df)))]
  }
  return(df)
}

#' Create a random population in high dimensional space
#'
#' Generate one population with one or more normally distributed modes in each
#' marker dimension. Mode counts, means, and spreads are sampled independently
#' for each dimension; cell membership is coordinated across dimensions.
#'
#' @param range_n Integer vector of candidate cell counts for this population.
#'   One entry is sampled uniformly. A single value fixes the cell count;
#'   c(500, 1000) selects either count, whereas 500:1000 allows every integer
#'   count in that interval.
#' @param dims Number of marker dimensions, named protein1, protein2, and so on.
#'   Annotation columns are additional to these dimensions.
#' @param range_mean Numeric vector of candidate mode means. Each mode mean is
#'   sampled with replacement using weights from prob_vec(). A single value
#'   fixes all mode means. The default candidates span 1 to 100000.
#' @param range_sd_fct Positive numeric vector of candidate spread factors.
#'   One entry is sampled uniformly with replacement for each mode. The mode's
#'   standard deviation is its mean divided by the sampled factor. Larger
#'   factors produce tighter modes; smaller factors produce wider modes.
#'   A single value fixes the factor for all modes.
#' @param range_modes Positive integer vector of candidate mode counts. One
#'   entry is sampled uniformly for each dimension. Set to 1 for one mode per
#'   dimension, to 3 for three modes per dimension, or to 1:3 to choose randomly
#'   among one, two, and three modes. Mode counts do not specify the number of
#'   populations or guarantee the number of visibly separated clusters.
#' @param mean_prob_vec_args Named list of additional arguments passed to
#'   prob_vec() to control the weights used to sample mode means. Weights follow
#'   the order of entries in range_mean. By default, the weights are U-shaped,
#'   favoring both ends of the ordered default candidate vector.
#' @param ... Additional arguments; currently unused by create_pop().
#'
#' @details
#' The generation procedure is:
#' \enumerate{
#'   \item Sample one cell count from range_n for the whole population.
#'   \item For each marker dimension, sample a mode count from range_modes,
#'     then sample a mean and a spread factor for each mode.
#'   \item Split the cell indices into shared blocks using
#'     random_split_with_bounds(), with as many blocks as the largest mode
#'     count across dimensions. For dimensions with fewer modes, randomly
#'     combine blocks using collapse_vectors(). This coordinates mode
#'     membership across dimensions instead of assigning each cell independently
#'     in every dimension.
#'   \item For each mode, generate values from a normal distribution with the
#'     selected mean and standard deviation. Combine the dimensions by cell
#'     identity and retain each cell's combination of mode memberships in group.
#' }
#'
#' One generated population can therefore contain several visible clusters.
#' Modes can also overlap, so the number of visible clusters need not match
#' either the mode counts or the number of populations requested elsewhere.
#' With range_modes set to 1, each population has one normal mode per dimension,
#' but separate populations may still overlap.
#'
#' Means are sampled separately for each dimension. With the default U-shaped
#' weights, a population may have several dimensions with very low means and
#' others with high means. Low values can appear to be zero on a wide axis.
#' If range_mean includes zero and it is selected, that mode has zero standard
#' deviation and all its generated values are exactly zero. Normal draws are
#' not truncated, so positive mode means can also produce negative values.
#'
#' make_hd_pops() calls create_pop() once per requested population; its npop
#' argument controls the population count. ff_simulate2() generates m samples,
#' each with independently generated populations when exprs is not supplied.
#' It removes the annotation columns, appends a Time channel, and optionally
#' renames and annotates marker channels when annotate_channels is TRUE.
#' Thus, dims specifies marker dimensions; each resulting flow frame also
#' contains the additional Time channel.
#'
#' @return A data frame with one row per cell, marker columns named protein1,
#'   protein2, and so on, a cell identifier column, and a group column encoding
#'   the combination of mode memberships across dimensions. This function does
#'   not add a population identifier or a Time channel.
#' @seealso make_hd_pops(), ff_simulate2(), prob_vec()
#' @export
#'
#' @examples
#' \dontrun{
#' # One mode per dimension
#' df <- create_pop(range_modes = 1, dims = 3)
#' brathering::plot3d(df, color = "group")
#' # Randomly choose one to three modes separately for each dimension
#' df <- create_pop(range_modes = 1:3, dims = 3)
#' brathering::plot3d(df, color = "group")
#' # Tighter modes using a larger spread factor
#' df <- create_pop(range_modes = 1:3, dims = 3, range_sd_fct = 20)
#' brathering::plot3d(df, color = "group")
#' # Fix the cell count, mode count, means, and spread factor
#' df <- create_pop(range_n = 500, dims = 3, range_modes = 1,
#'                  range_mean = 10000, range_sd_fct = 20)
#' #umap
#' tt <- ff_simulate2(
#'   npop = 6,
#'   range_n = 1e3:1e4,
#'   dims = 6,
#'   m = 2,
#'   annotate_channels = T,
#'   range_modes = 3
#' )
#' um <- fcexpr::ff_calc_umap_tsne(tt[[1]])
#' brathering::plot2(um)
#' }
create_pop <- function(range_n = c(5e2:1e4),
                       dims = 5,
                       range_mean = c(1:1e5),
                       range_sd_fct = c(2:10),
                       range_modes = c(1:3),
                       mean_prob_vec_args = list(),
                       ...) {

  fcexpr:::.ensure_packages(c("brathering"))
  n <- range_n[sample.int(length(range_n), 1)]
  # one list entry of m modalities for each dim

  pars <- purrr::map(stats::setNames(1:dims, paste0("protein", 1:dims)), function(x) {
    m <- range_modes[sample.int(length(range_modes), 1)]
    means <- range_mean[
      sample.int(
        n = length(range_mean),
        size = m,
        replace = TRUE,
        prob = do.call(
          prob_vec,
          args = c(list(n = length(range_mean)), mean_prob_vec_args)
        )
      )
    ]

    sd_fct <- range_sd_fct[sample.int(length(range_sd_fct), m, replace = TRUE)]
    sds <- means / sd_fct
    return(data.frame(means, sds))
  })

  # avoid that cells have random belongings to modalities in different dimensions
  nn <- random_split_with_bounds(total = n, splits = max(purrr::map_int(pars, nrow)))
  nn <- c(1,cumsum(nn))
  inds <- brathering::seq2(nn[-length(nn)], nn[-1])
  inds[-1] <- purrr::map(inds[-1], ~.x[-1])
  pars <- purrr::map(pars, function(x) {
    x[["inds"]] <- collapse_vectors(inds, nrow(x))
    return(x)
  })

  out <-  purrr::map_dfr(pars, function(x) {
    purrr::pmap_dfr(x, function(means,sds,inds) {
      data.frame(cell = inds,
                 FI = stats::rnorm(n = length(inds), mean = means, sd = sds))
    }, .id = "channel_mode")
  }, .id = "channel")
  out[["cell"]] <- stringr::str_pad(out[["cell"]], max(nchar(as.character(out[["cell"]]))), pad = "0")

  #modeinfo <- tidyr::nest(out, data = c(cell, FI))
  #attr(out, "modeinfo") <- modeinfo
  out <-
    out |>
    dplyr::mutate(pm = paste0(gsub("protein", "p", channel), "m", channel_mode)) |>
    dplyr::group_by(cell) |>
    dplyr::mutate(group = paste(unique(pm), collapse = "")) |>
    dplyr::select(-pm) |>
    dplyr::select(-channel_mode) |>
    tidyr::pivot_wider(names_from = channel, values_from = FI)
  out <- dplyr::relocate(out, cell, group, .after = names(out)[ncol(out)])

  # um <- fcexpr::ff_calc_umap_tsne(exprs = as.matrix(out[,1:(ncol(out)-2)]))
  # brathering::plot2(cbind(um, out[,ncol(out)]), color = "group")

  if (anyNA(out)) {
    message("NA found")
  }
  return(out)
}

#' Create a probability vector
#'
#' @param n number of elements
#' @param type type, roughly shape
#' @param power for M-type, exponent to initial probs, higher power --> higher peaks
#' @param peakshift for M-type, peakshift towards middle
#' @param beta_ab_range for U-type, range of ab values to dbeta, can also
#' be a vector of one or two values only
#'
#' @return numeric vector which sums to 1
#' @export
#'
#' @examples
#' library(ggplot2)
#' # type U: low beta_ab_range --> steep ends, high beta_ab_range --> close to unif
#' probs <- purrr::map_dfr(setNames(c(0.3,0.5,0.75,0.95),c(0.3,0.5,0.75,0.95)),
#'                         ~data.frame(x = prob_vec(type = "U", beta_ab = .x)) |>
#'                           dplyr::mutate(index = dplyr::row_number()),
#'                         .id = "beta_ab")
#' ggplot(probs, aes(index,x)) +
#'   geom_line() +
#'   facet_wrap(vars(beta_ab))
#' # type M: higher power to increase peaks
#' probs <- purrr::map_dfr(setNames(c(2,4,8,12),c(2,4,8,12)),
#'                         ~data.frame(x = prob_vec(type = "M", power = .x)) |>
#'                           dplyr::mutate(index = dplyr::row_number()),
#'                         .id = "power")
#' ggplot(probs, aes(index,x)) +
#'   geom_line() +
#'   facet_wrap(vars(power))
prob_vec <- function(n = 1e3,
                     type = c("U", "M"),
                     power = 5,
                     peakshift = c(0.2, 0.2),
                     beta_ab = c(0.2, 0.4)) {
  type <- rlang::arg_match(type)
  if (type == "U") {
    # value closer to 1 --> slopes less steep
    if (length(beta_ab) == 1) {
      beta_ab <- c(beta_ab, beta_ab)
    }

    prob <- stats::dbeta(seq(0,1, length.out = n+2), beta_ab[1], beta_ab[2])
    prob <- prob[-c(1,length(prob))] # rm Inf

  }

  if (type == "M") {
    #peakpos <- sample(seq(0.05, 0.2, length.out = 100), 1)
    if (length(peakshift) == 1) {
      peakshift <- c(peakshift, peakshift)
    }
    peak1 <- stats::dnorm(seq(0, 1, length.out = n), mean = peakshift[1], sd = 0.08)
    peak2 <- stats::dnorm(seq(0, 1, length.out = n), mean = 1-peakshift[2], sd = 0.08)
    prob <- peak1 + peak2
    prob <- prob^power
  }

  prob <- prob / sum(prob)
  return(prob)
}


random_split_with_bounds <- function(total = 100,
                                     splits = 3,
                                     min_frac = 0.1,
                                     max_frac = 0.5,
                                     max_iter = 200) {

  min_val <- ceiling(min_frac * total)
  max_val <- floor(max_frac * total)
  x <- rep(min_val, splits)
  iter <- 1
  while((any(x<min_val) | any(x>max_val) | sum(x)<total) && iter<=max_iter) {
    x <- sample(min_val:max_val, splits-1)
    x <- c(x, total-sum(x))
    iter <- iter + 1
  }
  return(x)
}


collapse_vectors <- function(vec_list, m) {
  n <- length(vec_list)
  if (m > n) stop("m must be less than or equal to n")
  if (m == n) {
    return(vec_list)
  }
  # Randomly assign each vector to one of m groups
  group_ids1 <- sample(1:m)
  group_ids2 <- sample(1:m, n-length(group_ids1), replace = TRUE)
  group_ids <- c(group_ids1, group_ids2)

  # Merge vectors by group
  collapsed <- lapply(unique(group_ids), function(i) {
    # Get all vectors assigned to group i and concatenate them
    do.call(c, vec_list[group_ids == i])
  })

  return(collapsed)
}

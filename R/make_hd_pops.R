#' Create a high dimensional space with multiple populations
#'
#' @param npop number of population to generate
#' @param popindexstart identity of first population
#' @param ... arguments to create_pop
#' @param rm_ann_cols remove annotation columns like cell and pop;
#' keep protein cols only
#'
#' @return data frame
#' @export
#'
#' @examples
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
#' @param range_n range of random total size
#' @param dims number of dimensions (columns)
#' @param range_mean range of random population means
#' @param range_sd_fct population spread: mean/sample(range_sd_fct,1)
#' @param range_modes random number peaks (modes, modalities) per dim,
#' set 1 to have unimodality only, range_n is split by modes randomly
#' with fcexpr:::random_split_with_bounds
#' @param mean_prob_vec_args create a probability vector with
#' fcexpr::prob_vec and provide its arguments here
#'
#' @return data frame
#' @export
#'
#' @examples
#' # unimodal only
#' df <- create_pop(range_modes = 1, dims = 3)
#' brathering::plot3d(df, color = "group")
#' # trimodal
#' df <- create_pop(range_modes = 1:3, dims = 3)
#' brathering::plot3d(df, color = "group")
#' # tight modes
#' df <- create_pop(range_modes = 1:3, dims = 3, range_sd_fct = 2)
#' brathering::plot3d(df, color = "group")
create_pop <- function(range_n = c(5e1:1e3),
                       dims = 5,
                       range_mean = c(1:1e5),
                       range_sd_fct = c(2:10),
                       range_modes = c(1:3),
                       mean_prob_vec_args = list(),
                       ...) {
  if (!requireNamespace("brathering", quietly = T)) {
    pak::pak("Close-your-eyes/brathering")
  }
  n <- sample(range_n, 1)

  # one list entry of m modalities for each dim
  pars <- purrr::map(stats::setNames(1:dims, paste0("protein", 1:dims)), function(x) {
    m <- sample(range_modes, 1)

    means <- sample(
      x = range_mean,
      size = m,
      replace = T,
      prob = do.call(prob_vec, args = c(list(n = length(range_mean)), mean_prob_vec_args))
    )
    sd_fct <- sample(range_sd_fct, m, replace = T)
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
    purrr::pmap_dfr(asplit(x,2), function(means,sds,inds) {

      data.frame(cell = inds,
                 FI = rnorm(n = length(inds), mean = means, sd = sds))
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
                     type = c("U",
                              "M"),
                     power = 5,
                     peakshift = c(0.2, 0.2),
                     beta_ab = c(0.2, 0.4)) {
  type <- rlang::arg_match(type)
  if (type == "U") {
    # value closer to 1 --> slopes less steep
    if (length(beta_ab) == 1) {
      beta_ab <- c(beta_ab, beta_ab)
    }

    prob <- dbeta(seq(0,1, length.out = n+2), beta_ab[1], beta_ab[2])
    prob <- prob[-c(1,length(prob))] # rm Inf

  }

  if (type == "M") {
    #peakpos <- sample(seq(0.05, 0.2, length.out = 100), 1)
    if (length(peakshift) == 1) {
      peakshift <- c(peakshift, peakshift)
    }
    peak1 <- dnorm(seq(0, 1, length.out = n), mean = peakshift[1], sd = 0.08)
    peak2 <- dnorm(seq(0, 1, length.out = n), mean = 1-peakshift[2], sd = 0.08)
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

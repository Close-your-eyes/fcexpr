library(ggplot2)
create_pop <- function(range_n = c(5e1:1e3),
                       dims = 5,
                       range_mean = c(1:1e5),
                       range_sd_fct = c(2:10),
                       range_modes = c(1:3),
                       mean_prob_vec_args = list(type = "U",
                                                 power = 5)) {
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

  if (anyNA(out)) {
    message("NA found")
  }
  return(out)
}

prob_vec <- function(n = 1e3,
                     type = c("U",
                              "M",
                              "unif"),
                     power = 5) {
  type <- rlang::arg_match(type)
  if (type == "U") {
    # value closer to 1 --> slopes less steep
    ab_vals <- sample(seq(0.5, 0.95, length.out = 100), 1)
    prob <- dbeta(seq(0,1, length.out = n+2), ab_vals, ab_vals)
    prob <- prob[-c(1,length(prob))] # rm Inf

  }

  if (type == "M") {
    peakpos <- sample(seq(0.05, 0.2, length.out = 100), 1)
    peak1 <- dnorm(seq(0, 1, length.out = n), mean = peakpos, sd = 0.08)
    peak2 <- dnorm(seq(0, 1, length.out = n), mean = 1-peakpos, sd = 0.08)
    prob <- peak1 + peak2
  }

  if (type == "unif"){
    prob <- 1
  }
  prob <- prob^power

  prob <- prob / sum(prob)
  return(prob)
}


create_sample <- function(npop = 5, popindexstart = 1, ...) {

  dots <- list(...)
  df <- replicate(npop, do.call(create_pop, args = dots), simplify = F)
  names(df) <- popindexstart:(npop+popindexstart-1)
  df <- dplyr::bind_rows(df, .id = "pop")
  df <- dplyr::relocate(df, pop, cell, group, .after = names(df)[ncol(df)])
  return(df)
}
df <- create_sample(range_modes = 1:3, mean_prob_vec_args = list(type = "U"))
attributes(df)
plot(df[,c(1,2)])
plot(df[,c(1,3)])

ggplot(df, aes(protein1, protein2)) +
  geom_point(aes(color = pop)) +
  scale_x_log10() +
  scale_y_log10() +
  brathering::scale_color_custom()

um <- uwot::umap(as.matrix(df[,1:5]), verbose = T)
df2 <- cbind(df, as.data.frame(um))
plot2(df2[,c("V1", "V2", "pop", "group")], cex = 3, color = "group", legend = NULL)

plot(prob_vec(n = 1000, type = "U"))
plot(prob_vec(n = 1000, type = "M"))
brathering::plot2(as.data.frame(df[,c(2,3)]))
plot3d(df)



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

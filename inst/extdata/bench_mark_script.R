install.packages("bench")

n<-1000

?bench::press
?bench::mark
pca_bench <- bench::press(
  n = c(1e3, 1e4, 1e5),

  {
    X <- as.matrix(data.frame(a = rnorm(n, 0, 1),
                              b = rnorm(n, 10, 2),
                              c = rnorm(n, 20, 0.5),
                              d = rnorm(n, 4, 8),
                              e = rnorm(n, 100, 4),
                              f = rnorm(n, 1, 2)))
    #browser()
    bench::mark(
      check = FALSE,
      time_unit = 's',
      max_iterations = 5,
      "stats::prcomp()" = {
        prcomp(X, center = TRUE, scale = TRUE, rank. = 20)
      },
      "irlba::prcomp_irlba()" = {
        irlba::prcomp_irlba(x = X, n = 4, center = TRUE, scale. = TRUE)
      },
      "irlba::irlba()" = {
        X_center <- colMeans(X)
        X_scale <- proxyC::colSds(X)
        suppressWarnings({
          retval <- irlba::irlba(A = X, nv = 4, center = X_center, scale = X_scale)
        })
        retval$x <- retval$u %*% diag(retval$d)
        retval
      }
    )
  }
)
pca_bench$method = attr(pca_bench$expression, "description")
d <- as_tibble(pca_bench) %>% select(-expression, -result, -memory, -time, -gc)

library(tidyverse)
remotes::install_github("AllanCameron/geomtextpath")


# https://github.com/slowkow/slowkow.com/blob/master/content/notes/pca-benchmark/index.Rmd
ggplot(d) +
  aes(n, median, color = method, label = reorder(method, n / median)) +
  geomtextpath::geom_textline(aes(label = method), size = 5.5, text_smoothing = 30) +
  guides(
    color = guide_legend(
      reverse = FALSE,
      override.aes = list(
        label = levels(fct_reorder(d$method, d$n / d$median)),
        size = 5
      ),
      label = FALSE
    )
  ) +
  ggrepel::geom_text_repel(
    data = d %>% filter(n == max(n)),
    mapping = aes(label = glue::glue("{signif(median, 3)}s")),
    size = 6, nudge_x = 0.1, box.padding = 0.1,
    direction = "y", hjust = 0, segment.size = NA
  ) +
  scale_x_log10(
    labels = scales:: label_number_si(), expand = expansion(mult = c(0, 0.3))
  ) +
  scale_y_log10(
    labels = function(x) (x), name = "Seconds",
    expand = expansion(mult = c(0, 0.2), add = 0)
  ) +
  theme(
    # legend.box.margin = margin(0, 0, 0, 0, "mm"),
    # legend.background = element_rect(fill = NA),
    # legend.position = c(0.9, 0.1),
    # legend.justification = c(1, 0),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.3, color = "grey90")
  ) +
  annotation_logticks(side = "lb") +
  labs(
    title = "Seconds elapsed",
    x = "Cells"
  )


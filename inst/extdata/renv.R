library(renv)
install.packages("renv")
res <- dependencies() |>
  dplyr::filter(grepl("fcexpr/R", Source)) |>
  dplyr::distinct(Package)


res$Package

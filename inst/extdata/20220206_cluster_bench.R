## ---- setup ----------
library(tidyverse)
library(fcexpr)
library(ClusterR)
wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
drfz <- "/Volumes/AG_Hiepe/Christopher.Skopnik/202010-202111_Experiments/20210602_DVD_PC2/FJ.workspaces"

wsps <- list.files(drfz, pattern = "\\.wsp$", recursive = T, full.names = T)
wsp <- wsps[1]
ffs <- fcexpr::wsp_get_ff(wsp, groups = "Pack 8", population = "PCsub_IgA_TACI", FCS.file.folder = file.path(dirname(drfz), "FCS_files"),
                          lapply_fun = parallel::mclapply, mc.cores = 3)
ffs <- ffs[[1]]


channels <- setNames(ffs[["logicle"]][[1]]@parameters@data[["name"]], ffs[["logicle"]][[1]]@parameters@data[["desc"]])
channels <- channels[which(!is.na(names(channels)))]
expr <- do.call(rbind, lapply(ffs[["logicle"]], function(x) flowCore::exprs(x)))[,channels]
nrow(expr)
# inspired by:
# https://slowkow.com/notes/pca-benchmark/

# https://github.com/r-lib/bench
# https://allancameron.github.io/geomtextpath/

clust_bench <- bench::press(
  n = c(1e3, 5e3, 1e4, 5e4, 1e5), #n = c(1e3, 5e3, 1e4, 5e4, 1e5, 2.5e5, 5e5, 1e6),

  {
    X <- expr[sample(1:nrow(expr), n, replace = n > nrow(expr)),]
    rownames(X) <- 1:nrow(X)
    #browser()
    bench::mark(
      check = FALSE,
      time_unit = 's',
      max_iterations = 3,

      "stats::kmeans" = {
        stats::kmeans(X, centers = 10)
      },
      "ClusterR::KMeans_arma" = {
        ClusterR::KMeans_arma(X, clusters = 10)
      },
      "ClusterR::MiniBatchKmeans" = {
        predict_MBatchKMeans(X, MiniBatchKmeans(X, clusters = 10, batch_size = 20, num_init = 5, max_iters = 100,
                                                   init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,
                                                   verbose = F)$centroids)
      },
      "louvain_seurat" = {
        Seurat::FindClusters(Seurat::FindNeighbors(X, annoy.metric = "cosine")$snn)
      },
      "leiden" = {
        leiden::leiden(Seurat::FindNeighbors(X, annoy.metric = "cosine")$snn)
      },
      "hclust" = {
        stats::cutree(stats::hclust(stats::dist(X)), k = 10)
      },
      "flowClust" = {
        flowClust::flowClust(X, K = 10)
      },
      "MUDAN_igraph::cluster_walktrap" = {
        message("cluster_walktrap")
        MUDAN::getComMembership(X, k = 30, method = igraph::cluster_walktrap)
      },
      "MUDAN_igraph::cluster_louvain" = {
        message("cluster_louvain")
        MUDAN::getComMembership(X, k = 30, method = igraph::cluster_louvain)
      },
      "MUDAN_igraph::cluster_infomap" = {
        message("cluster_infomap")
        MUDAN::getComMembership(X, k = 30, method = igraph::cluster_infomap)
      }
    )
  }
)



clust_bench$method = attr(clust_bench$expression, "description")
d <- as_tibble(clust_bench) %>% select(-expression, -result, -memory, -time, -gc)
ggplot2::autoplot(clust_bench)

# https://github.com/slowkow/slowkow.com/blob/master/content/notes/pca-benchmark/index.Rmd
p1 <- ggplot(d, aes(n, median, color = method, label = reorder(method, n / median))) +
  geomtextpath::geom_textline(aes(label = method), size = 5.5, text_smoothing = 30) +
  guides(color = guide_legend(reverse = FALSE,
                              override.aes = list(label = levels(fct_reorder(d$method, d$n / d$median)), size = 5),
                              label = FALSE)) +
  ggrepel::geom_text_repel(data = d %>% filter(n == max(n)),
                           mapping = aes(label = glue::glue("{signif(median, 3)}s")),
                           size = 6,
                           nudge_x = 0.1,
                           box.padding = 0.1,
                           direction = "y", hjust = 0,
                           segment.size = NA) +
  scale_x_log10(labels = scales::label_number_si(), expand = expansion(mult = c(0, 0.3))) +
  scale_y_log10(labels = function(x) (x), name = "Seconds", expand = expansion(mult = c(0, 0.2), add = 0)) +
  theme(legend.position = "none",panel.grid.major = element_line(size = 0.3, color = "grey90")) +
  annotation_logticks(side = "lb") +
  labs(title = "Seconds elapsed", x = "Cells")

p2 <- ggplot(d, aes(n, mem_alloc, color = method, label = reorder(method, -n / median))) +
  geomtextpath::geom_textline(size = 5.5, text_smoothing = 30) +
  guides(color = guide_legend(reverse = TRUE, override.aes = list(label = "", size = 3))) +
  ggrepel::geom_text_repel(data = d %>% filter(n == max(n)),
                           mapping = aes(label = as.character(mem_alloc)),
                           size = 6,
                           nudge_x = 0.1,
                           box.padding = 0.1,
                           hjust = 0,
                           direction = "y") +
  scale_x_log10(labels = scales::label_number_si(), expand = expansion(mult = c(0, 0.3))) +
  scale_y_log10(labels = function(x) (x / 1e9), name = "GB", expand = expansion(mult = c(0, 0.2), add = 0)) +
  theme(legend.position = "none", panel.grid.major = element_line(size = 0.3, color = "grey90")) +
  annotation_logticks(side = "lb") +
  labs(title = "GB allocated", x = "Cells")




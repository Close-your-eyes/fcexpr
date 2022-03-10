## ---- setup ----------
library(tidyverse)
library(fcexpr)
library(ClusterR)
wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
drfz <- "/Volumes/AG_Hiepe/Christopher.Skopnik/202010-202111_Experiments/20210602_DVD_PC2/FJ.workspaces"

wsps <- list.files(drfz, pattern = "\\.wsp$", recursive = T, full.names = T)
wsp <- wsps[1]





ffs <- fcexpr::wsp_get_ff(wsp, groups = c("Pack 1", "Pack 2", "Pack 3", "Pack 4", "Pack 5", "Pack 6", "Pack 7", "Pack 8", "Pack 9"), population = "PCsub_IgA_TACI", FCS.file.folder = file.path(dirname(drfz), "FCS_files"),
                          lapply_fun = parallel::mclapply, mc.cores = 6)
ffs <- ffs[[1]]


channels <- setNames(ffs[["logicle"]][[1]]@parameters@data[["name"]], ffs[["logicle"]][[1]]@parameters@data[["desc"]])
channels <- channels[which(!is.na(names(channels)))]
expr <- do.call(rbind, lapply(ffs[["logicle"]], function(x) flowCore::exprs(x)))[,channels]
nrow(expr)
# inspired by:
# https://slowkow.com/notes/pca-benchmark/

# https://github.com/r-lib/bench
# https://allancameron.github.io/geomtextpath/


#### round 1

# cluster_walktrap failed at 50k:
# Error in graph.adjacency.dense(adjmatrix, mode = mode, weighted = weighted,  : long vectors not supported yet: /Volumes/Builds/R4/R-4.1.0/src/include/Rinlinedfuns.h:537
# cluster louvain failed at 50k cells (memory exhausted)
# hclust failed at 100k cells


clust_bench <- bench::press(
  n = c(1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 1e7), #n = c(1e3, 5e3, 1e4, 5e4, 1e5, 2.5e5, 5e5, 1e6),

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
      "ClusterR::KMeans_rcpp" = {
        ClusterR::KMeans_rcpp(X, clusters = 10)
      },
      "ClusterR::MiniBatchKmeans" = {
        predict_MBatchKMeans(X, MiniBatchKmeans(X, clusters = 10, batch_size = 20, num_init = 5, max_iters = 100,
                                                   init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,
                                                   verbose = F)$centroids)
      },
      "louvain_seurat" = {
        if (nrow(X) <= 5e5) {
          Seurat::FindClusters(Seurat::FindNeighbors(X, annoy.metric = "cosine")$snn)
        } else {
          NULL
        }
      },
      "leiden" = {
        if (nrow(X) <= 1e4) {
          leiden::leiden(Seurat::FindNeighbors(X, annoy.metric = "cosine")$snn)
        } else {
          NULL
        }
      },
      "hclust" = {
        if (nrow(X) <= 5e4) {
          stats::cutree(stats::hclust(stats::dist(X)), k = 10)
        } else {
          NULL
        }
      },
      "flowClust" = {
        if (nrow(X) <= 1e5) {
          flowClust::flowClust(X, K = 10)
        } else {
          NULL
        }
      },
      "MUDAN_igraph::cluster_walktrap" = {
        if (nrow(X) <= 1e4) {
          MUDAN::getComMembership(X, k = 30, method = igraph::cluster_walktrap)
        } else {
          NULL
        }
      },
      "MUDAN_igraph::cluster_louvain" = {
        if (nrow(X) <= 1e4) {
          MUDAN::getComMembership(X, k = 30, method = igraph::cluster_louvain)
        } else {
          NULL
        }
      },
      "MUDAN_igraph::cluster_infomap" = {
        if (nrow(X) <= 1e4) {
          MUDAN::getComMembership(X, k = 30, method = igraph::cluster_infomap)
        } else {
          NULL
        }
      }
    )
  }
)

clust_bench$method = attr(clust_bench$expression, "description")
d <- as_tibble(clust_bench) %>% select(-expression, -result, -memory, -time, -gc)
d <- d[which(d$mem_alloc != "0B"),]
saveRDS(d, file.path(wd, "20220309_bench_result.RDS"))

# https://github.com/slowkow/slowkow.com/blob/master/content/notes/pca-benchmark/index.Rmd
p1 <- ggplot(d, aes(n, median, color = method, label = reorder(method, n / median))) +
  geomtextpath::geom_textline(aes(label = method), size = 5.5, text_smoothing = 30) +
  guides(color = guide_legend(reverse = FALSE,
                              override.aes = list(label = levels(fct_reorder(d$method, d$n / d$median)), size = 5),
                              label = FALSE)) +
  ggrepel::geom_text_repel(data = d %>% dplyr::group_by(method) %>% filter(n == max(n)),
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
  ggrepel::geom_text_repel(data = d %>% dplyr::group_by(method) %>% filter(n == max(n)),
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

p_grid <- cowplot::plot_grid(p1, p2)
ggsave(p_grid, filename = paste0("1k_to_10M_cells.png"), device = "png", path = wd, dpi = "retina", width = 18, height = 8)




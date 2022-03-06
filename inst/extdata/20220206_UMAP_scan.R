## ---- setup ----------
library(tidyverse)
library(fcexpr)
library(ClusterR)
wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

# https://towardsdatascience.com/umap-dimensionality-reduction-an-incredibly-robust-machine-learning-algorithm-b5acb01de568

wsps <- list.files(wd, pattern = "\\.wsp$", recursive = T, full.names = T)
wsps <- wsps[1]
ffs <- fcexpr::wsp_get_ff(wsps, groups = "Pack 8", population = "PCsub_IgA_TACI", FCS.file.folder = file.path(wd, "FCS_files"),
                          lapply_fun = parallel::mclapply, mc.cores = 3)
ffs <- ffs[[1]]


channels <- setNames(ffs[["logicle"]][[1]]@parameters@data[["name"]], ffs[["logicle"]][[1]]@parameters@data[["desc"]])
channels <- channels[which(!is.na(names(channels)))]
expr <- do.call(rbind, lapply(ffs[["logicle"]], function(x) flowCore::exprs(x)))[,channels]
nrow(expr)


## --- variation of umap pars ---------
expr_sub <- expr[sample(1:nrow(expr), 5e3),]
rownames(expr_sub) <- 1:nrow(expr_sub)
clust <- Seurat::FindClusters(Seurat::FindNeighbors(expr_sub, annoy.metric = "cosine")$snn)


var <- expand.grid(#n_neighbors = c(10,20,30,50,100),
  spread = c(0.1, 0.5, 1, 2, 3, 5, 10),
  min_dist = c(0.001, 0.01, 0.05, 0.1, 0.5))

res <- lapply(split(var, 1:nrow(var)), function (x) as.data.frame(uwot::umap(expr_sub, spread = x[["spread"]], min_dist = x[["min_dist"]], verbose = F)) %>% dplyr::mutate(spread = x[["spread"]], min_dist = x[["min_dist"]]))


## --- variation of splitting data for umap + umap_transform ---------
rownames(expr) <- 1:nrow(expr)
clust <- Seurat::FindClusters(Seurat::FindNeighbors(expr, annoy.metric = "cosine")$snn)
model_frac <- c(0.05, 0.1, 0.2, 0.3, 0.5)
model_rows <- 1:round(0.1*nrow(expr))
transform_rows <- (max(model_rows)+1):nrow(expr)

umap_model <- uwot::umap(expr[model_rows,], verbose = T, ret_model = T)
umap_embed <- rbind(umap_model[["embedding"]],uwot::umap_transform(expr[transform_rows,], model = umap_model, verbose = T))


# benchmark
umap_bench <- bench::press(
  model_frac = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1),

  {
    model_rows <- 1:round(0.1*nrow(expr))
    transform_rows <- (max(model_rows)+1):nrow(expr)
    bench::mark(
      check = FALSE,
      time_unit = 's',
      max_iterations = 5,
      "uwot::umap" = {
        umap_model <- uwot::umap(expr[model_rows,], verbose = T, ret_model = T)
        umap_embed <- rbind(umap_model[["embedding"]],uwot::umap_transform(expr[transform_rows,], model = umap_model, verbose = T))
      }
    )
  }
)


umap_bench$method = attr(umap_bench$expression, "description")
d <- as_tibble(umap_bench) %>% select(-expression, -result, -memory, -time, -gc)
# https://github.com/slowkow/slowkow.com/blob/master/content/notes/pca-benchmark/index.Rmd
p1 <- ggplot(d, aes(model_frac, median, color = method, label = reorder(method, n / median))) +
  geomtextpath::geom_textline(aes(label = method), size = 5.5, text_smoothing = 30) +
  guides(color = guide_legend(reverse = FALSE,
                              override.aes = list(label = levels(fct_reorder(d$method, d$model_frac / d$median)), size = 5),
                              label = FALSE)) +
  ggrepel::geom_text_repel(data = d %>% filter(n == max(model_frac)),
                           mapping = aes(label = glue::glue("{signif(median, 3)}s")),
                           size = 6,
                           nudge_x = 0.1,
                           box.padding = 0.1,
                           direction = "y", hjust = 0,
                           segment.size = NA) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Seconds elapsed", y = "seconds")

p2 <- ggplot(d, aes(model_frac, mem_alloc, color = method, label = reorder(method, -n / median))) +
  geomtextpath::geom_textline(size = 5.5, text_smoothing = 30) +
  guides(color = guide_legend(reverse = TRUE, override.aes = list(label = "", size = 3))) +
  ggrepel::geom_text_repel(data = d %>% filter(n == max(model_frac)),
                           mapping = aes(label = as.character(mem_alloc)),
                           size = 6,
                           nudge_x = 0.1,
                           box.padding = 0.1,
                           hjust = 0,
                           direction = "y") +
  scale_y_log10(labels = function(x) (x / 1e8), expand = expansion(mult = c(0, 0.2), add = 0)) +
  theme_bw() +
  theme(legend.position = "none") +
  annotation_logticks(side = "lb") +
  labs(title = "GB allocated")






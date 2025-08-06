library(ggplot2)

ff <- read.FCS("/Volumes/CMS_SSD_2TB/Experiment_data/20230222_TNF_IFNg_CD107a_induction_by_IL15_IFNb/FCS_files/Exp_part_1/0006_-_CMS_IL15.fcs")
colnames(ff@exprs)
hist(ff@exprs[,18])


?MASS::mvrnorm
Sigma <- matrix(c(10,3,3,2),2,2)
Sigma
mv <- MASS::mvrnorm(n = 1000, rep(0, 2), Sigma)
hist(mv[,1])
hist(mv[,2])

library(mclust)
data <- ff@exprs[sample(1:nrow(ff@exprs), 50000),]
ncol(data)
um <- uwot::umap(data)
plot(um)

?mclust::mclustModel()

gmm <- mclust::mclustBIC(data)
datamodel <- mclustModel(data, gmm)
datasim <- sim(modelName = datamodel$modelName,
               parameters = datamodel$parameters,
               n = 10000)
umsim <- uwot::umap(datasim[,-1])
plot(umsim)


plot(datasim[,c(2,5)])
plot(data[,c(1,4)])


datasim_list <-  lapply(1:10, function(X) {
  sim(modelName = datamodel$modelName,
      parameters = datamodel$parameters,
      n = 10000)[,-1]
})

um_list <- lapply(datasim_list, uwot::umap, verbose = T)
names(um_list) <- 1:length(um_list)

um_df <- purrr::map_dfr(um_list, as.data.frame, .id = "simulation")

ggplot(um_df, aes(V1,V2)) +
  geom_point() +
  facet_wrap(vars(simulation))

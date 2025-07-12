#library(keras)
library(keras3)
library(tensorflow)
library(fcexpr)
library(brathering)
ff <- flowCore::read.FCS("/Volumes/CMS_SSD_2TB/Experiment_data/20230222_TNF_IFNg_CD107a_induction_by_IL15_IFNb/FCS_files/Exp_part_1/0006_-_CMS_IL15.fcs")
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
source("/Users/vonskopnik/Documents/R/2024_neural_nets/R/vae_model.R")
# Example numeric matrix input
channels <- ff_get_channels(ff, rm_wo_desc = T)
rowinds <- brathering::sample_leverage(ff@exprs[,channels], size = 10000, replace = T, n = 1000)

x_data <- lapply(rowinds, function(x) scale(ff@exprs[x,channels]))
dim(x_data[[1]])
x_data <- abind::abind(x_data, along = 3)
x_data <- aperm(x_data, c(3,1,2))
dim(x_data)
class(x_data)
typeof(x_data)
dims <- 10000*length(channels)
x_data <- reticulate::array_reshape(x_data, c(1000, dims), order = "F")

rspheres <- replicate(5, hdos:::random_sphere(), simplify = F)
rtorus <- replicate(5, hdos:::random_torus(), simplify = F)
rcuboid <- replicate(5, hdos:::random_cuboid(), simplify = F)
landscape <- hdos::bind_objects(c(rspheres, rtorus, rcuboid), missing_dim_fill = "runif")
landscape <- hdos::add_noise_points(landscape, density = 0.005)
landscape2 <- as.matrix(landscape[["coord"]])
rowinds <- brathering::sample_leverage(landscape2, size = 10000, replace = T, n = 1000)
x_data <- lapply(rowinds, function(x) landscape2[x,])
x_data <- abind::abind(x_data, along = 3)
x_data <- aperm(x_data, c(3,1,2))
dims <- 10000*ncol(landscape2)
x_data <- reticulate::array_reshape(x_data, c(1000, dims), order = "F")

orig_dim = dims
lat_dim = 2
int_dim <- 64
batch_size = 20
epochs = 10

## encoder
encoder_inputs <- keras3::layer_input(shape = orig_dim)
int_layer <-
  encoder_inputs |>
  keras3::layer_dense(units = 512, activation = 'relu') |>
  keras3::layer_batch_normalization() |>
  keras3::layer_dropout(rate = 0.3) |>
  keras3::layer_dense(units = 256, activation = 'relu') |>
  keras3::layer_batch_normalization() |>
  keras3::layer_dropout(rate = 0.3) |>
  keras3::layer_dense(units = 128, activation = 'relu') |>
  keras3::layer_batch_normalization() |>
  keras3::layer_dropout(rate = 0.3)


#int_layer <- keras3::layer_dense(encoder_inputs, units = int_dim, activation = "relu")
z_mean <- keras3::layer_dense(int_layer, units = lat_dim, name = "z_mean")
z_log_var <- keras3::layer_dense(int_layer, units = lat_dim, name = "z_log_var")
encoder <- keras3::keras_model(inputs = encoder_inputs,
                               outputs = list(z_mean, z_log_var),
                               name = "encoder")

## decoder
latent_inputs <- keras3::layer_input(shape = c(lat_dim))
decoder_outputs <-
  latent_inputs |>
  keras3::layer_dense(int_dim, activation = "relu") |>
  keras3::layer_dense(units = 128, activation = "relu") |>
  keras3::layer_dense(units = 256, activation = "relu") |>
  keras3::layer_dense(units = 512, activation = "relu") |>
  keras3::layer_dense(units = orig_dim, activation = "linear")
decoder <- keras3::keras_model(inputs = latent_inputs,
                               outputs = decoder_outputs,
                               name = "decoder")

callback_list <- list(
  keras3::callback_early_stopping(monitor = "kl_loss", patience = 5, restore_best_weights = TRUE, mode = "min")
)

## compile and train the vae
vae <- do.call(model_vae, args = list(encoder, decoder, loss_fn = keras3::loss_huber)) # loss_fn <- loss_huber(delta = 1.0)
keras3::compile(vae, optimizer = keras3::optimizer_adam())
keras3::fit(
  object = vae,
  x = x_data,
  epochs = epochs,
  batch_size = batch_size,
  shuffle = TRUE,
  callbacks = callback_list
)


um <- uwot::umap(lapply(rowinds[1], function(x) scale(ff@exprs[x,channels]))[[1]], verbose = T)
plot2(um, cex = 2)


p1 <- predict(vae$decoder, matrix(c(1,-4), nrow = 1)) |> matrix(ncol = 7, byrow = T)
colnames(p1) <- 1:7
plot3d(p1)
plot3d(p1[,3:5])
plot3d(landscape[[1]][,2:4][sample(1:nrow(landscape[[1]]), 10000),])
um2 <- uwot::umap(p1, verbose = T)
clusts <- fcexpr::get_louvain_cluster(um2)
um2 <- cbind(um2, clusts)
plot2(um2, color = "res.0.1", cex = 2)


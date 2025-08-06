library(keras3)
library(tensorflow)
library(fcexpr)
ff <- flowCore::read.FCS("/Volumes/CMS_SSD_2TB/Experiment_data/20230222_TNF_IFNg_CD107a_induction_by_IL15_IFNb/FCS_files/Exp_part_1/0006_-_CMS_IL15.fcs")
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
source(file.path(wd, "vae_model.r"))
# https://divingintogeneticsandgenomics.com/post/how-to-code-a-variational-autoencoder-vae-in-r-using-mnist-dataset/

# Example: assume `your_data` is a numeric matrix
# x_data <- scale(as.matrix(your_data))
# input_dim <- ncol(x_data)
# latent_dim <- 2  # you can adjust this
channels <- ff_get_channels(ff, rm_wo_desc = T)
rowinds <- sample_leverage(ff@exprs[,channels], size = 10000, replace = T, n = 100)
dim(rowinds)

x_train <- apply(rowinds, 2, function(x) ff@exprs[x,channels], simplify = F)
dim(x_train[[1]])
x_train <- abind::abind(x_train, along = 3)
x_train <- aperm(x_train, c(3,1,2))
dim(x_train)
class(x_train)
typeof(x_train)
#x_train <- as.double(x_train)
x_train <- reticulate::array_reshape(x_train, c(100, 10000*length(channels)), order = "F")

original_dim <- as.integer(10000*length(channels))
latent_dim <- 2L
intermediate_dim <- 512
batch_size<- 20

encoder_inputs <- keras3::layer_input(shape = 10000 * length(channels))
x <- keras3::layer_dense(encoder_inputs, intermediate_dim, activation = "relu")

z_mean    <- keras3::layer_dense(x, latent_dim, name = "z_mean")
z_log_var <- keras3::layer_dense(x, latent_dim, name = "z_log_var")
encoder <- keras3::keras_model(encoder_inputs, list(z_mean, z_log_var), name = "encoder")
encoder


layer_sampler <- keras3::new_layer_class(
  classname = "Sampler",
  call = function(z_mean, z_log_var) {
    epsilon <- tf$random$normal(shape = tf$shape(z_mean))
    z_mean + exp(0.5 * z_log_var) * epsilon }
)

latent_inputs <- keras3::layer_input(shape = c(latent_dim))

decoder_outputs <-
  latent_inputs |>
  keras3::layer_dense(intermediate_dim, activation = "relu") |>
  keras3::layer_dense(original_dim, activation = "sigmoid")
decoder <- keras3::keras_model(latent_inputs, decoder_outputs, name = "decoder")
decoder

vae <- model_vae(encoder, decoder)
keras3::compile(vae, optimizer = keras3::optimizer_adam())
keras3::fit(vae, x_train, epochs = 5, shuffle = TRUE)

#library(keras)
library(keras3)
library(tensorflow)
library(fcexpr)
ff <- flowCore::read.FCS("/Volumes/CMS_SSD_2TB/Experiment_data/20230222_TNF_IFNg_CD107a_induction_by_IL15_IFNb/FCS_files/Exp_part_1/0006_-_CMS_IL15.fcs")
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
source(file.path(wd, "vae_model.r"))
# Example numeric matrix input
channels <- ff_get_channels(ff, rm_wo_desc = T)
rowinds <- sample_leverage(ff@exprs[,channels], size = 10000, replace = T, n = 100)
#x_data <- scale(as.matrix(ff@exprs[rowinds[,1],channels]))
input_dim <- ncol(x_data)
latent_dim <- 2

x_data <- apply(rowinds, 2, function(x) ff@exprs[x,channels], simplify = F)
dim(x_data[[1]])
x_data <- abind::abind(x_data, along = 3)
x_data <- aperm(x_data, c(3,1,2))
dim(x_data)
class(x_data)
typeof(x_data)
x_data <- reticulate::array_reshape(x_data, c(100, 10000*length(channels)), order = "F")


# Encoder layers
encoder_input <- layer_input(shape = input_dim)
encoder_dense <- layer_dense(units = 64, activation = "relu")(encoder_input)
z_mean <- layer_dense(units = latent_dim)(encoder_dense)
z_log_var <- layer_dense(units = latent_dim)(encoder_dense)

sampling <- function(args) {
  z_mean <- args[[1]]
  z_log_var <- args[[2]]
  epsilon <- keras::k_random_normal(shape = keras::k_shape(z_mean))
  z_mean + keras::k_exp(0.5 * z_log_var) * epsilon
}

z <- layer_lambda(f = sampling)(list(z_mean, z_log_var))

# Decoder layers
decoder_input <- layer_input(shape = latent_dim)
decoder_dense <- layer_dense(units = 64, activation = "relu")
decoder_output <- layer_dense(units = input_dim, activation = "linear")
h_dec <- decoder_dense(z)
reconstructed <- decoder_output(h_dec)

# Define the VAE model
vae <- keras_model(encoder_input, reconstructed)

# Define custom loss function
vae_loss_fn <- function(x, x_decoded) {
  recon_loss <- keras::k_sum(keras::k_square(x - x_decoded), axis = -1)
  kl_loss <- -0.5 * keras::k_sum(1 + z_log_var - keras::k_square(z_mean) - keras::k_exp(z_log_var), axis = -1)
  keras::k_mean(recon_loss + kl_loss)
}

vae %>% compile(optimizer = "adam", loss = vae_loss_fn)

# Train the model
vae %>% fit(
  x = x_data,
  epochs = 5,
  batch_size = 20)

# Build decoder separately for simulation
decoder_input_layer <- layer_input(shape = latent_dim)
decoder_h <- decoder_dense(decoder_input_layer)
decoder_out <- decoder_output(decoder_h)
decoder <- keras_model(decoder_input_layer, decoder_out)

# Simulate new data from the learned latent structure
n_samples <- 1000
z_samples <- matrix(rnorm(n_samples * latent_dim), ncol = latent_dim)
generated_data <- predict(decoder, z_samples)

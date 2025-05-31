ff <- flowCore::read.FCS("/Volumes/CMS_SSD_2TB/Experiment_data/20230222_TNF_IFNg_CD107a_induction_by_IL15_IFNb/FCS_files/Exp_part_1/0006_-_CMS_IL15.fcs")
colnames(ff@exprs)
hist(ff@exprs[,18])
your_data <- ff@exprs[sample(1:nrow(ff@exprs), 50000),]

# https://divingintogeneticsandgenomics.com/post/how-to-code-a-variational-autoencoder-vae-in-r-using-mnist-dataset/

# Install and load keras
#install.packages("keras")
library(keras)
#install_keras()
# pip3 install tensorflow

# Example: assume `your_data` is a numeric matrix
x_data <- scale(as.matrix(your_data))
input_dim <- ncol(x_data)
latent_dim <- 2  # you can adjust this

# Encoder
inputs <- layer_input(shape = input_dim)
h <- layer_dense(inputs, units = 64, activation = "relu")
z_mean <- layer_dense(h, units = latent_dim)
z_log_var <- layer_dense(h, units = latent_dim)

sampling <- function(args) {
  z_mean <- args[[1]]
  z_log_var <- args[[2]]
  epsilon <- k_random_normal(shape = c(k_shape(z_mean)[[1]], latent_dim))
  z_mean + k_exp(z_log_var / 2) * epsilon
}

z <- layer_lambda(list(z_mean, z_log_var), sampling)

# Decoder
decoder_h <- layer_dense(units = 64, activation = "relu")
decoder_mean <- layer_dense(units = input_dim, activation = "linear")
h_decoded <- decoder_h(z)
x_decoded_mean <- decoder_mean(h_decoded)

# VAE model
vae <- keras_model(inputs, x_decoded_mean)

# VAE loss
vae_loss <- function(x, x_decoded_mean) {
  reconstruction_loss <- loss_mean_squared_error(x, x_decoded_mean)
  kl_loss <- -0.5 * k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var))
  reconstruction_loss + kl_loss
}

vae %>% compile(optimizer = "adam", loss = vae_loss)

# Train the model
vae %>% fit(x_data, x_data,
            epochs = 50,
            batch_size = 32,
            validation_split = 0.1)

# Build the decoder model for simulation
decoder_input <- layer_input(shape = latent_dim)
h_dec <- decoder_h(decoder_input)
x_dec <- decoder_mean(h_dec)
decoder <- keras_model(decoder_input, x_dec)

# Simulate new data
n_samples <- 1000
z_samples <- matrix(rnorm(n_samples * latent_dim), ncol = latent_dim)
generated_data <- predict(decoder, z_samples)

# Reverse scaling (optional, if you stored mean and sd)
# generated_data <- sweep(generated_data, 2, attr(x_data, "scaled:scale"), "*")
# generated_data <- sweep(generated_data, 2, attr(x_data, "scaled:center"), "+")

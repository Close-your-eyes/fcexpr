
model_vae <- new_model_class(
  classname = "VAE",

  initialize = function(encoder, decoder, ...) {
    super$initialize(...)
    self$encoder <- encoder
    self$decoder <- decoder
    self$sampler <- layer_sampler()
    self$total_loss_tracker <-
      metric_mean(name = "total_loss")
    self$reconstruction_loss_tracker <-
      metric_mean(name = "reconstruction_loss")
    self$kl_loss_tracker <-
      metric_mean(name = "kl_loss")
  },

  metrics = mark_active(function() {
    list(
      self$total_loss_tracker,
      self$reconstruction_loss_tracker,
      self$kl_loss_tracker
    )
  }),

  train_step = function(data) {
    with(tf$GradientTape() %as% tape, {

      c(z_mean, z_log_var) %<-% self$encoder(data)
      z <- self$sampler(z_mean, z_log_var)

      reconstruction <- decoder(z)
      reconstruction_loss <-
        loss_binary_crossentropy(data, reconstruction) |>
        sum(axis = c(1)) |>
        mean()

      kl_loss <- -0.5 * (1 + z_log_var - z_mean^2 - exp(z_log_var))
      total_loss <- reconstruction_loss + mean(kl_loss)
    })

    grads <- tape$gradient(total_loss, self$trainable_weights)
    self$optimizer$apply_gradients(zip_lists(grads, self$trainable_weights))

    self$total_loss_tracker$update_state(total_loss)
    self$reconstruction_loss_tracker$update_state(reconstruction_loss)
    self$kl_loss_tracker$update_state(kl_loss)

    list(total_loss = self$total_loss_tracker$result(),
         reconstruction_loss = self$reconstruction_loss_tracker$result(),
         kl_loss = self$kl_loss_tracker$result())
  }
)

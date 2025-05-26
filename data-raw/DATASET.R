## code to prepare `DATASET` dataset goes here

non_fluo_channels <- c("Time", "HDR-T", "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W",
                       "Center", "Offset", "Width", "Residual", "Event_length")
usethis::use_data(non_fluo_channels, internal = T, overwrite = T)

## code to prepare `DATASET` dataset goes here

non_fluo_channels <- c("Time", "HDR-T", "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W",
                       "Center", "Offset", "Width", "Residual", "Event_length")

random_operators <- c("Gangolf Eierschmalz",
                      "Walter Frosch",
                      "Fled Nanders",
                      "Dean Norm",
                      "Francesco Rosinetti",
                      "Ali Sweidel",
                      "Frau Kepetry",
                      "Ken Guru",
                      "Fresh Dumbledore",
                      "Mam Bagera",
                      "Enrico Pallazzo",
                      "Friedrich Quecksilber",
                      "Frank Drebin",
                      "Coward Harpendale",
                      "Houg Daffernan",
                      "Rustin Cohle",
                      "Jens Bloedermann",
                      "Lasse Samenstroem",
                      "Baracko Bama",
                      "Grave Dohl",
                      "Roland Habeck",
                      "Man Jarsalek",
                      "Mark David Chapman",
                      "Pepe Mujica",
                      "Aibert Huwanger")

ff_params <- readRDS(system.file("extdata", "ff_params_empty.rds", package = "fcexpr"))


usethis::use_data(non_fluo_channels, random_operators, ff_params, internal = T, overwrite = T)

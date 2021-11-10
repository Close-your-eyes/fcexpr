# xlsx file names have to match whats gonna be written to the $CYT keyword in FCS files
wd <- file.path(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)), "ignore_folder")
tab <- lapply(list.files(wd, "\\.xlsx$", full.names = T)[which(!grepl("\\~", list.files(wd, "\\.xlsx$")))], function(x) {
  tab <- do.call(rbind, lapply(openxlsx::getSheetNames(x)[-1], function(y) {
    tab <- openxlsx::read.xlsx(x, sheet = y, rows = c(9,14))
    tab$channel <- y
    return(tab)
    #read_xlsx(file.path(wd, "Results_Fortessa.xlsx"), sheet = x, range = "A9:N15") %>% dplyr::rename("variable" = `...1`) %>% tidyr::gather(key = "voltage", value = "value", -c(variable)) %>% dplyr::filter(str_detect(variable, "K-factor")) %>% dplyr::mutate(channel = x)
  }))
  tab <- tidyr::pivot_longer(tab, cols = -channel, names_to = "volt", values_to = "k")
  tab$volt <- as.numeric(tab$volt)
  tab <- tab[which(!is.na(tab$k)),]
  return(as.data.frame(tab))
})
names(tab) <- gsub("\\.xlsx$", "", basename(list.files(wd, "\\.xlsx$", full.names = T)[which(!grepl("\\~", list.files(wd, "\\.xlsx$")))]))


tab <- lapply(tab, function(y) {
  dplyr::bind_rows(lapply(unique(y$channel), function(x) {
    out <- approx(x = y[which(y$channel == x), "volt"], y = y[which(y$channel == x), "k"], n = max(y[which(y$channel == x), "volt"])-min(y[which(y$channel == x), "volt"]) + 1)
    data.frame(volt = out$x,
               k = out$y,
               channel = x)
  }))
})
saveRDS(tab, file.path(dirname(wd), "inst", "extdata", "k_factors.rds"))


ggplot(tab[[2]], aes(x = volt, y = k, group = channel)) +
  geom_point(size = 0.2) +
  geom_line() +
  theme_bw() +
  facet_wrap(vars(channel), scales = "free_y")



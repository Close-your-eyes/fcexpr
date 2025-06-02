library(fcexpr)
library(mclust)
library(ggplot2)
ws <- "/Volumes/CMS_SSD_2TB/Experiment_data/20230222_TNF_IFNg_CD107a_induction_by_IL15_IFNb/FJ_workspaces/Exp_part_1_2_3_4_5_6.wsp"

# check samples, groups, populations
groups <- fcexpr::wsx_get_groups(ws)
paths <- fcexpr::wsx_get_fcs_paths(ws, filter_AllSamples = T)
samplenames <- basename(grep("AnKr", paths$Exp_part_1, value = T))
pops <- fcexpr::wsx_get_poppaths(ws, groups = "Exp_part_1")


# read flow frames
fflist <- fcexpr::wsp_get_ff(wsp = ws,
                             FCS.file.folder = "/Volumes/CMS_SSD_2TB/Experiment_data/20230222_TNF_IFNg_CD107a_induction_by_IL15_IFNb/FCS_files",
                             population = "Single Cells",
                             samples = samplenames[c(1,4)])
ff <- purrr::list_flatten(fflist[["flowframes"]], name_spec = "{outer}")
ff <- purrr::list_flatten(ff, name_spec = "{outer}")

modellist <- ff_model_GMM(ff = ff[[2]],
                          sample_n = 50000,
                          mclustBIC_args = list(modelNames = "VVV", G = 10),
                          source_file = names(ff)[2])
model <- modellist[[2]]
# saveRDS(model[which(names(model) %in%  c("original_data", "z")], file.path("/Users/vonskopnik/Documents/R_packages/fcexpr/inst", paste0(names(ff)[2], "_mclustmodel.rds")))


# compare simulation to original wit umaps
ff_list <- ff_simulate(model = model, m = 9, n = 3e4)
orig_umap_coords <- ff_calc_umap_tsne(exprs = model[["original_data"]][,ff_get_channels(ff_list[[1]], rm_wo_desc = T)])
umap_coords <- purrr::map(ff_list, ff_calc_umap_tsne, channels = ff_get_channels(ff_list[[1]], rm_wo_desc = T))
umap_coords2 <- dplyr::bind_rows(c(list(orig = as.data.frame(orig_umap_coords)),
                                   purrr::map(umap_coords, as.data.frame)), .id = "simulation")

ggplot(umap_coords2, aes(UMAP_1, UMAP_2)) +
  geom_point() +
  facet_wrap(vars(simulation))


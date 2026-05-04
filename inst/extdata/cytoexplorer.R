library(devtools)
install_github("DillonHammill/openCyto", force = TRUE)
install_github("federicomarini/shinyBS")
install_github("DillonHammill/DataEditR")
install_github("DillonHammill/HeatmapR")
# CytoExploreRData
#pak::pak("DillonHammill/CytoExploreRData")
# CytoExploreR
pak::pak("DillonHammill/CytoExploreR")


fs <- flowCore::read.flowSet(c("/Volumes/CMS_SSD_2TB/Experiment_data/20210707_IL15_NKG2D_MICAB_target_cell_killing/FCS_files/20220831_Exp.part.21/1141_-_IL2_IL7_IL15_IFNb_quadruple_cond_4_DuDa.fcs",
                               "/Volumes/CMS_SSD_2TB/Experiment_data/20210707_IL15_NKG2D_MICAB_target_cell_killing/FCS_files/20220831_Exp.part.21/1125_-_IL2_IL7_IL15_IFNb_quadruple_cond_3_DuDa.fcs"))
CytoExploreR::cyto_plot(fs, channels = c("FSC-A", "SSC-A"))

scattergate <- CytoExploreR::cyto_gate_draw(
  fs,
  channels = c("FSC-A", "SSC-A"),
  alias = "scatter"
)

trans <- CytoExploreR::cyto_transformer_logicle(fs)
fs <- CytoExploreR::cyto_transform(fs,
                                   trans = trans)

CytoExploreR::cyto_plot(fs, channels = c("FSC-A", "SSC-A"), gate = scattergate)
CytoExploreR::cyto_plot(fs, channels = c("FSC-A", "SSC-A"), gate = scattergate, parent = "scatter")
CytoExploreR::cyto_fluor_channels(fs)
CytoExploreR::cyto_markers_extract(fs, channels = CytoExploreR::cyto_fluor_channels(fs))
CytoExploreR::cyto_channels(fs)


# Extract gated cells as a new flowFrame
ff_gated <- flowCore::Subset(fs, scattergate[[1]][[1]])
CytoExploreR::cyto_plot(ff_gated, channels = c("FSC-A", "SSC-A"), gate = scattergate)
CytoExploreR::cyto_plot(ff_gated, channels = c("v-LP600 610/20-D-A", "v-LP505 525/50-E-A"), axes_trans = trans)

# ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Exp_part_20_21.wsp"
# ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/ExpPart_6_for_pub.wsp"
# ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/20231005_FJ_exp_wsp.wsp" # error: ids are wrong, multiple ids for AndOrNodes - fixed!

ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/no_group_noOrAndGates.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_more_gate_types.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_different_gating_trees.wsp"

ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_1D_gates_range_and_2Sector.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_1D_gates_range_and_2Sector_with_NotGate.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_1D_gates_range_and_2Sector_with_NotGate_differentGatingTrees.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Multiple_OrNodes_AndNodes_sameDims_sameGatingTrees.wsp" # one full path missing - check why; same id assigned to different AndNodes?!; check what happens if one of the multiple nodes is removed
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Multiple_OrNodes_AndNodes_sameDims_differentGatingTrees.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Multiple_OrNodes_AndNodes_NotNode_on_OrAndNodes_sameDims_sameGatingTrees.wsp" # NotNodes of OrNodes or AndNodes are different! - check if code lines are compatible with ordinary NotNodes
ws <- "/Users/vonskopnik/Desktop/example_workspaces/Multiple_OrNodes_AndNodes_NotNode_on_OrAndNodes_sameDims_sameGatingTrees.wsp"
ws <- "/Users/vonskopnik/Desktop/example_workspaces/Multiple_OrNodes_AndNodes_NotNode_on_OrAndNodes_with_children_sameDims_sameGatingTrees.wsp"
ws <- "/Users/vonskopnik/Desktop/example_workspaces/OrAndNodes_from_different_Dims_sameGatingTrees.wsp" # some gate types missing
ws <- "/Users/vonskopnik/Desktop/example_workspaces/Exp_part_20_21.wsp" # large ws
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Complicated_OrAndGates_OrGate_at_diff_hierachies_sameGatingTree.wsp" # error

library(fcexpr)
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
# wslist <- list.files(file.path(wd, "example_wsp"), full.names = T)
# purrr::map(wslist, R.utils::gunzip, remove = F)
wslist <- list.files(file.path(wd, "example_wsp"), full.names = T)

## run new and old method sequentially and have results compared
# errors: index 2
# very long: index 3
# diffs: index 6
out <- wsx_get_popstats(wslist[6])

wsx_get_popstats2(wslist[2], strip_data = T)
### manual comparison of new and old method
system.time(out1 <- wsx_get_popstats_legacy(wslist[6], strip_data = T))
system.time(out2 <- wsx_get_popstats2(wslist[6], strip_data = T))

df1 <- out1[["counts"]]
df2 <- out2[["counts"]]
#fix columns of df2 a bit

comcols <- intersect(colnames(df1), colnames(df2))
df1 <- dplyr::arrange(df1, FileName, Count)
df2 <- dplyr::arrange(df2, FileName, Count)

df1_sub <- df1[,which(names(df1) %in% c("FileName", "PopulationFullPath", "Population", "Count"))]#, "gate_id", "parentgate_id"))]
df2_sub <- df2[,which(names(df2) %in% c("FileName", "PopulationFullPath", "Population", "Count"))]# "id", "parent_id"))]

waldo::compare(df1_sub, df2_sub)



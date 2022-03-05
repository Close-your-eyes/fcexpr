.onLoad <- function(libname = find.package("scexpr"), pkgname = "scexpr") {

  # https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when

  # CRAN Note avoidance
  if(getRversion() >= "2.15.1") {
    utils::globalVariables(c('$BEGINDATA','$FIL','$TOT','Antigen','Antigen.Conjugate','Antigen.x','Antigen.y',
                           'Box','Box.x','Box.y','Conjugate','Conjugate.x','Conjugate.y','FileName','FilePath','Lot','Lot.x',
                           'Lot.y','Population','PopulationFullPath','antigen_case','channel','cluster','cluster_1','cluster_2',
                           'config','conjugate_case','diff_sign','diptest_not_p','diptest_p','diptest_p_1','diptest_p_2',
                           'feature','gate.level','gate.path.full','group','mean_1','mean_2','mean_diff','mean_not','pval','pvalue','row.num','x','y','ws', 'sp'))
  }

}



#string <- c("$BEGINDATA $FIL $TOT Antigen Antigen.Conjugate Antigen.x Antigen.y Box Box.x Box.y Conjugate Conjugate.x Conjugate.y FileName FilePath Lot Lot.x Lot.y Population PopulationFullPath all_of antigen_case channel cluster cluster_1 cluster_2 combn config conjugate_case diff_sign diptest_not_p diptest_p diptest_p_1 diptest_p_2 feature gate.level gate.path.full group is_logical mean_1 mean_2 mean_diff mean_not pval pvalue row.num setNames sp stack ws x y")
#string2 <- strsplit(string, " ")[[1]]
#cat(paste(string2, collapse = "','"))

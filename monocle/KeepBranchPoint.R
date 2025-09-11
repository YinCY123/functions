KeepBranchPoint <- function(obj, branches){
  obj@auxOrderingData[[obj@dim_reduce_type]]$branch_points <- obj@auxOrderingData[[obj@dim_reduce_type]]$branch_points[branches]
  return(obj)
}
outdeg <- function(obj){
  out <- sapply(obj, function(obj)sum(obj$outdeg))
  return(out)
}

indeg <- function(obj){
  out <- sapply(obj, function(obj)sum(obj$indeg))
  return(out)
}
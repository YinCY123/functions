sparse2h5 <- function(samples, 
                      sample.names,
                      dir, 
                      type = "sparse", 
                      compressed = TRUE,
                      row.names = "symbol", 
                      col.names = TRUE, 
                      ...){
  # loading required packages
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(DropletUtils))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(fs))
  
  # does the dir exist
  if(!dir_exists(dir)){
    dir_create(path = dir)
    message("dir does not exist, created!!")
  }
  
  for(i in seq_along(samples)){
    tmp <- read10xCounts(samples = samples[[1]], 
                         sample.names = sample.names[[i]], 
                         col.names = col.names, 
                         row.names = row.names, 
                         type = type, 
                         compressed = compressed)
    
    path <- paste0(dir, ifelse(str_detect(dir, "\\/$"), "", "/"), sample.names[[i]], ".h5")
    write10xCounts(path = path, 
                   x = counts(tmp), 
                   gene.id = rowData(tmp)$ID, 
                   gene.symbol = rownames(tmp), 
                   type = "HDF5")
  }
}
sparse2h5 <- function(samples, 
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
  
  for(sample in samples){
    tmp <- read10xCounts(samples = sample, 
                         sample.names = basename(sample), 
                         col.names = col.names, 
                         row.names = row.names, 
                         type = type, 
                         compressed = compressed)
    
    path <- paste0(dir, ifelse(str_detect(dir, "\\/$"), "", "/"), basename(sample), ".h5")
    write10xCounts(path = path, 
                   x = counts(tmp), 
                   gene.id = rowData(tmp)$ID, 
                   gene.symbol = rownames(tmp), 
                   type = "HDF5")
  }
}
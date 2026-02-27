geth5 <- function(file, 
                  data = "/X/data", 
                  col_ids = "/X/indices", 
                  feature = "/var/_index", 
                  indptr = "/X/indptr",
                  barcode = "/obs/cell_ID_mask", 
                  cell_meta = "/obs", 
                  feature_meta = "/var", 
                  spatial_loc = "/obsm/spatial", 
                  ...){
  suppressPackageStartupMessages(library(magrittr))
  suppressPackageStartupMessages(library(rhdf5))
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(readr))
  
  message(paste0("contructing expressiong matrix ..."))
  values <- h5read(file = file, name = data) %>% as.vector()
  j_indices <- as.vector(h5read(file = file, name = col_ids)) + 1
  indptr <- h5read(file = file, name = indptr) %>% as.vector()
  
  i_indices <- length(indptr) - 1
  i_indices <- rep(seq_len(i_indices), times = diff(indptr))
  
  features <- h5read(file = file, name = feature) %>% as.character()
  barcodes <- h5read(file = file, name = barcode) %>% as.character()

  mtx <- sparseMatrix(
    i = i_indices, 
    j = j_indices, 
    x = values,
    dims = c(length(features), length(barcodes)), 
    dimnames = list(features, barcodes)
  )

  # get cell metadata
  message(paste0("get cell metadata ..."))
  vars <- h5read(file, name = cell_meta) %>% names
  cell_meta_info <- vector(mode = "list")
  for(var in vars){
    cell_meta_info[[var]] <- h5read(file, paste0(cell_meta, "/", var))
  }
  cell_meta_info <- do.call(cbind, cell_meta_info) %>% as.data.frame() %>% 
    readr::type_convert()
  
  # get feature metadata
  message(paste0("get feature metadata ..."))
  vars <- h5read(file, feature_meta) %>% names
  feature_meta_info <- vector(mode = "list")
  for(var in vars){
    feature_meta_info[[var]] <- h5read(file, paste0(feature_meta, "/", var))
  }
  
  feature_meta_info <- do.call(cbind, feature_meta_info) %>% as.data.frame() %>% readr::type_convert()

  # get spatial information
  spatial_loc <- h5read(file, name = spatial_loc) %>% t %>% as.matrix
  
  tmp <- list(
    data = mtx, 
    cell_meta = cell_meta_info, 
    feature_meta = feature_meta_info, 
    spatial = spatial_loc
  )
  return(tmp)
}

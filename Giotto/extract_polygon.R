extract_polygon <- function(image, is_cosmx = FALSE, fov = NULL){
  # loading packages
  library(EBImage)
  library(sf)
  
  mask <- readImage(image)
  mask <- round(mask * (1/min(mask[mask != 0])) * 1.05)
  
  contours <- ocontour(mask)
  
  # ensure polygon closed
  contours <- lapply(contours, function(x){
    x <- as.matrix(x)
    if(!identical(x[1,], x[nrow(x), ])){
      x <- rbind(x, x[1,,drop = FALSE])
    }
    return(x)
  })
  
  poly_list <- lapply(seq_along(contours), function(i){
    tmp <- data.frame(
      cell_ID = i, 
      vertex_x = contours[[i]][, 1],
      vertex_y = contours[[i]][, 2]
    )
    
    if(is_cosmx){
      tmp$fov = fov
    }
    
    return(tmp)
  })
  
  poly_df <- do.call(rbind, poly_list)
  
  return(poly_df)
}

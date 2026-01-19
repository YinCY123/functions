spatial_coord_norm <- function(spes_list, factor = 1000, spatial_name = "raw", spatial_unit = "cell", ...){
    library(magrittr)
    library(Giotto)

    scales <- function(x){
        ((x - min(x))/(max(x) - min(x))) * factor
    }
    distance <- sapply(spes_list, function(x){
        x %>% getSpatialLocations(output = "data.table") %>% as.data.frame %>% dplyr::pull(1) %>% diff %>% max
    })
    factors <- max(distance)/distance
    names(factors) <- names(spes_list)

    # normalized each location
    for(x in names(spes_list)){
        loc <- spes_list[[x]] %>% getSpatialLocations(output = "data.table") %>% as.data.frame
        loc[[1]] <- scales(loc[[1]] * factors[[x]])
        loc[[2]] <- scales(loc[[2]] * factors[[x]])

        spes_list[[x]] <- setSpatialLocations(gobject = spes_list[[x]], 
                x = createSpatLocsObj(coordinates = loc, name = spatial_name, spat_unit = spatial_unit), 
                name = spatial_name, 
                spat_unit = spatial_unit)
    }
    return(spes_list)
}
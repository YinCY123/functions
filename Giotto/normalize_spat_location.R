normalize_spat_location <- function(gobject, factor = 1000){
    library(magrittr)
    
    scales <- function(x){
        ((x - min(x)) / (max(x) - min(x))) * factor
    }

    tmp <- gobject %>% getSpatialLocations(output = "data.table") %>% as.data.frame()
    tmp[[1]] <- scales(tmp[[1]])
    tmp[[2]] <- scales(tmp[[2]])
    # tmp <- apply(tmp, 2, scales)

    # convert to spatLocsObj
    spatLoc <- createSpatLocsObj(tmp, name = "raw", spat_unit = "cell")

    gobject <- setSpatialLocations(gobject, x = spatLoc)
    
    return(gobject)
}
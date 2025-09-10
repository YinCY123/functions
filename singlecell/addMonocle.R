addMonocle <- function(sce, cds, ...){
    # get pseudotime coordinates
    reducedDim(sce, "monocle") <- cds@reducedDimS %>% t

    # get pseudotime and state
    vars <- c("Pseudotime", "State")
    for(var in vars){
        colData(sce)[[var]] <- pData(cds)[[var]]
    }
    return(sce)
}
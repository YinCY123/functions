cellchat_process <- function(cellchat, 
    type = "triMean", 
    trim = 0.1,
    min_cell = 10, 
    contact.range = NULL,
    ...){
        # loading required packages
        suppressPackageStartupMessages(library(magrittr))
        suppressPackageStartupMessages(library(CellChat))

        cellchat <- cellchat %>% 
            subsetData() %>% 
            identifyOverExpressedGenes(do.fast = TRUE, min.cells = min_cell) %>% 
            identifyOverExpressedInteractions() %>% 
            computeCommunProb(type = type, trim = trim, contact.range = contact.range) %>% 
            filterCommunication(min.cells = min_cell) %>% 
            computeCommunProbPathway() %>% 
            aggregateNet() %>% 
            netAnalysis_computeCentrality()
        return(cellchat)
}
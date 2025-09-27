cellchat_process <- function(cellchat, 
    type = "triMean", 
    trim = 0.1,
    min_cell = 10, 
    ...){
        # loading required packages
        suppressPackageStartupMessages(require(magrittr))
        suppressPackageStartupMessages(require(CellChat))

        cellchat <- cellchat %>% 
            subsetData() %>% 
            identifyOverExpressedGenes(do.fast = TRUE, min.cells = min_cell) %>% 
            identifyOverExpressedInteractions() %>% 
            computeCommunProb(type = type, trim = trim) %>% 
            filterCommunication(min.cells = min_cell) %>% 
            computeCommunProbPathway() %>% 
            aggregateNet() %>% 
            netAnalysis_computeCentrality()
        return(cellchat)
}

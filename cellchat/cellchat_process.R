cellchat_process <- function(cellchat, 
    type = "triMean", 
    trim = 0.1,
    min_cell = 10, 
    ...){
        # loading required packages
        suppressPackageStartupMessages(library(magrittr))
        suppressPackageStartupMessages(library(CellChat))

        # cellchat <- cellchat %>% 
        #     subsetData() %>% 
        #     identifyOverExpressedGenes(do.fast = TRUE, min.cells = min_cell) %>% 
        #     identifyOverExpressedInteractions() %>% 
        #     computeCommunProb(type = type, trim = trim) %>% 
        #     filterCommunication(min.cells = min_cell) %>% 
        #     computeCommunProbPathway() %>% 
        #     aggregateNet() %>% 
        #     netAnalysis_computeCentrality()

        message("subsetData ...")
        cellchat <- subsetData(cellchat)

        message("identifyOverExpressedGenes ...")
        cellchat <- identifyOverExpressedGenes(cellchat, do.faster = TRUE, min.cells = min_cell)

        message("identifyOverExpressedInteractions ...")
        cellchat <- identifyOverExpressedInteractions(cellchat)

        message("computeCommunProb ...")
        cellchat <- computeCommunProb(type = type, trim = trim)

        message("filterCommunication ...")
        cellchat <- filterCommunication(cellchat, min.cells = min_cell)

        message("computeCommunProbPathway ...")
        cellchat <- computeCommunProbPathway(cellchat)

        message("aggregateNet ...")
        cellchat <- aggregateNet(cellchat)

        message("netAnalysis_computeCentrality ...")
        cellchat <- netAnalysis_computeCentrality(cellchat)

        return(cellchat)
}

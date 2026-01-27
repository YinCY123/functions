scTenifoldKnk <- function (countMatrix, qc = TRUE, gKO = NULL, qc_mtThreshold = 0.1, 
    qc_minLSize = 1000, nc_lambda = 0, nc_nNet = 10, nc_nCells = 500, 
    nc_nComp = 3, nc_scaleScores = TRUE, nc_symmetric = FALSE, 
    nc_q = 0.9, td_K = 3, td_maxIter = 1000, td_maxError = 1e-05, 
    td_nDecimal = 3, ma_nDim = 2, nCores = parallel::detectCores()) 
{
    if (isTRUE(qc)) {
        countMatrix <- scQC(countMatrix, mtThreshold = qc_mtThreshold, 
            minLSize = qc_minLSize)
    }
    # if (ncol(countMatrix) > 500) {
    #     countMatrix <- countMatrix[rowMeans(countMatrix != 0) >=  0.05, ]
    # }
    # else {
    #     countMatrix[rowSums(countMatrix != 0) >= 25, ]
    # }
    WT <- scTenifoldNet::makeNetworks(X = countMatrix, q = nc_q, 
        nNet = nc_nNet, nCells = nc_nCells, scaleScores = nc_scaleScores, 
        symmetric = nc_symmetric, nComp = nc_nComp, nCores = nCores)
    WT <- scTenifoldNet::tensorDecomposition(xList = WT, K = td_K, 
        maxError = td_maxError, maxIter = td_maxIter, nDecimal = td_nDecimal)
    WT <- WT$X
    WT <- strictDirection(WT, lambda = nc_lambda)
    WT <- as.matrix(WT)
    diag(WT) <- 0
    WT <- t(WT)
    KO <- WT
    KO[gKO, ] <- 0
    MA <- manifoldAlignment(WT, KO, d = ma_nDim, nCores = nCores)
    DR <- dRegulation(MA, gKO)
    outputList <- list()
    outputList$tensorNetworks$WT <- Matrix(WT)
    outputList$tensorNetworks$KO <- Matrix(KO)
    outputList$manifoldAlignment <- MA
    outputList$diffRegulation <- DR
    return(outputList)
}

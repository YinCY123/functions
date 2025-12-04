scTenifoldKnk <- function (countMatrix, qc = TRUE, gKO = NULL, 
    qc_mtThreshold = 0.1, qc_feature = 0.05, seq_deeps = 500,
    qc_minLSize = 1000, nc_lambda = 0, nc_nNet = 10, nc_nCells = 500, 
    nc_nComp = 3, nc_scaleScores = TRUE, nc_symmetric = FALSE, 
    nc_q = 0.9, td_K = 3, td_maxIter = 1000, td_maxError = 1e-05, 
    td_nDecimal = 3, ma_nDim = 2, nCores = parallel::detectCores()) 
{
    if (isTRUE(qc)) {
        countMatrix <- scQC(countMatrix, mtThreshold = qc_mtThreshold, 
            minLSize = qc_minLSize)
    }
    if (ncol(countMatrix) > seq_deeps) {
        countMatrix <- countMatrix[rowMeans(countMatrix != 0) >= qc_feature, ]
    }
    else {
        countMatrix <- countMatrix[rowSums(countMatrix != 0) >= 25, ]
    }
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
    # validate and translate gKO to numeric row indices to avoid subscript out of bounds
    if (!is.null(gKO)) {
        if (is.character(gKO)) {
            if (is.null(rownames(KO))) {
                stop("gKO provided as character vector but KO has no row names")
            }
            idx <- match(gKO, rownames(KO))
            if (any(is.na(idx))) {
                stop("Some gKO names not found in KO rownames: ", paste(gKO[is.na(idx)], collapse = ", "))
            }
        } else if (is.numeric(gKO)) {
            idx <- as.integer(gKO)
            if (any(idx < 1 | idx > nrow(KO))) {
                stop("Numeric gKO indices out of range")
            }
        } else {
            stop("gKO must be NULL, a character vector of row names, or integer indices")
        }
        KO[idx, ] <- 0
    }
    MA <- manifoldAlignment(WT, KO, d = ma_nDim, nCores = nCores)
    DR <- dRegulation(MA, gKO)
    outputList <- list()
    outputList$tensorNetworks$WT <- Matrix(WT)
    outputList$tensorNetworks$KO <- Matrix(KO)
    outputList$manifoldAlignment <- MA
    outputList$diffRegulation <- DR
    return(outputList)
}


# helper functions
scQC <- function(X, mtThreshold = 0.1, minLSize = 1000){
  if(class(X) == 'Seurat'){
    countMatrix <- X@assays$RNA@counts
  } else {
    countMatrix <- X
  }
  librarySize <- colSums(countMatrix)
  countMatrix <- countMatrix[,librarySize >= minLSize]
  librarySize <- colSums(countMatrix)
  mtGenes <- grep('^MT-',toupper(rownames(countMatrix)))
  nGenes <- colSums(countMatrix != 0)

  genesLM <- lm(nGenes~librarySize)
  genesLM <- as.data.frame(predict(genesLM, data.frame(librarySize), interval = 'prediction'))

  if(isTRUE(length(mtGenes) > 0)){
    mtCounts <- colSums(countMatrix[grep('^MT-',toupper(rownames(countMatrix))),])
    mtProportion <- mtCounts/librarySize
    mtLM <- lm(mtCounts~librarySize)
    mtLM <- as.data.frame(predict(mtLM, data.frame(librarySize), interval = 'prediction'))
    selectedCells <- mtCounts > mtLM$lwr & mtCounts < mtLM$upr & nGenes > genesLM$lwr & nGenes < genesLM$upr & mtProportion <= mtThreshold & librarySize < 2 * mean(librarySize)
  } else {
    selectedCells <- nGenes > genesLM$lwr & nGenes < genesLM$upr & librarySize < 2 * mean(librarySize)
  }
  selectedCells <- colnames(countMatrix)[selectedCells]
  if(class(X) == 'Seurat'){
    X <- subset(X, cells = selectedCells)
  } else {
    X <- countMatrix[,selectedCells]
  }
  return(X)
}

strictDirection <- function(X, lambda = 1){
  S <- as.matrix(X)
  S[abs(S) < abs(t(S))] <- 0
  O <- (((1-lambda) * X) + (lambda * S))
  O <- Matrix::Matrix(O)
  return(O)
}
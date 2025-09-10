addModuleScore <- function(obj, 
                           param = c("ssgsea", "gsva", "plage", "zscore"),                            
                           # common parameter
                           genesets = NULL, 
                           minSize = 10, 
                           maxSize = 5000,
                           ncore = 6,
                           
                           # gsva
                           kcdf = "Gaussian",
                           tau = 1, 
                           maxDiff = TRUE,
                           
                           # ssgsea
                           alpha = 0.25,
                           normalize = TRUE,
                           ...){
  # loading required packages
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(GSVA))
  suppressPackageStartupMessages(require(SingleCellExperiment))
  suppressPackageStartupMessages(require(BiocParallel))
  suppressPackageStartupMessages(require(stringr))
  
  # param
  param <- match.arg(param)

  # convert to sparse matrix
  mtx <- logcounts(obj)
  mtx <- as(mtx, Class = "CsparseMatrix")
  
  if(str_to_lower(param) == "gsva"){
    p = gsvaParam(exprData = mtx, 
                  geneSets = genesets, 
                  minSize = minSize, 
                  maxSize = maxSize, 
                  kcdf = kcdf, 
                  tau = tau, 
                  maxDiff = maxDiff)
    gsva_results <- gsva(param = p, 
                         verbose = TRUE, 
                         BPPARAM = SnowParam(workers = ncore))
    
    # add enrichment score to obj
    for(name in rownames(gsva_results)){
      colData(obj)[[name]] <- gsva_results[name, ,drop = T][rownames(colData(obj))]
    }
    return(obj)
  }else if(param == "plage"){
    p = plageParam(exprData = mtx, 
                   geneSets = genesets, 
                   minSize = minSize,
                   maxSize = maxSize)
    plage_results <- gsva(param = p, 
                          verbose = TRUE, 
                          BPPARAM = SnowParam(workers = ncore))
    
    # add enrichment scoreto obj
    for(name in rownames(plage_results)){
      colData(obj)[[name]] <- plage_results[name, ,drop = T][rownames(colData(obj))]
    }
    return(obj)
    
  }else if(str_to_lower(param) == "zscore"){
    p = zscoreParam(exprData = mtx, 
                    geneSets = genesets, 
                    minSize = minSize, 
                    maxSize = maxSize)
    zscore_results <- gsva(param = p, 
                           verbose = TRUE, 
                           BPPARAM = SnowParam(workers = ncore))
    
    # add enrichment score to obj
    for(name in rownames(zscore_results)){
      colData(obj)[[name]] <- zscore_results[name, ,drop = TRUE][rownames(colData(obj))]
    }
    return(obj)
  }else if(str_to_lower(param) == "ssgsea"){
    p = ssgseaParam(exprData = mtx, 
                    geneSets = genesets, 
                    minSize = minSize, 
                    maxSize = maxSize, 
                    alpha = alpha, 
                    annotation = GSEABase::SymbolIdentifier(),
                    normalize = normalize)
    ssgsea_results <- gsva(param = p, 
                           verbose = TRUE, 
                           BPPARAM = SnowParam(workers = ncore))
    
    # add enrichment score to obj
    for(name in rownames(ssgsea_results)){
      colData(obj)[[name]] <- ssgsea_results[name, , drop = TRUE][rownames(colData(obj))]
    }
    return(obj)
  }else{
    message(paste0("This method: ", param, " currently not supported..."))
  }
}

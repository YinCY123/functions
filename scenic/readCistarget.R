readCistarget <- function(x, ...){
    library(stringr)
    library(magrittr)
    
    ct <- read.table(x, sep = ",", header = TRUE)
    
    r1 <- ct[1, , drop = TRUE] %>% unlist(use.names = F)
    r2 <- ct[2, , drop = TRUE] %>% unlist(use.names = F)
    nms <- ifelse(nchar(r2) < 1, r1, r2)
    ct <- ct[-c(1, 2), ]
    colnames(ct) <- nms

    ct <- ct %>% dplyr::select(TF, MotifID, AUC, NES, MotifSimilarityQvalue, OrthologousIdentity, Annotation, Context, TargetGenes) %>% 
        tibble::as_tibble()

    df <- ct %>% 
        dplyr::mutate(TargetGenes = str_extract_all(TargetGenes, "\\(.+\\)")[[1]]) %>% 
        tidyr::separate_longer_delim(TargetGenes, delim = "),") %>% 
        dplyr::mutate(TargetGenes = str_trim(TargetGenes, side = "both"), 
            TargetGenes = str_remove_all(TargetGenes, "'|\\(|\\)")) %>% 
        tidyr::separate(col = TargetGenes, sep = ",", into = c("target", "weight")) %>% 
        dplyr::mutate(weight = str_trim(weight, side = "both"), 
            weight = as.numeric(weight), 
            AUC = as.numeric(AUC), 
            NES = as.numeric(NES), 
            MotifSimilarityQvalue = as.numeric(MotifSimilarityQvalue), 
            OrthologousIdentity = as.numeric(OrthologousIdentity))
    
    return(df)
}

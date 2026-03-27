readCistarget <- function(x, ...){
    library(stringr)
    library(magrittr)
    
    ct <- read.table(x, sep = ",", header = TRUE)
    
    r1 <- ct[1, , drop = TRUE] %>% unlist(use.names = F)
    r2 <- ct[2, , drop = TRUE] %>% unlist(use.names = F)
    nms <- ifelse(nchar(r2) < 1, r1, r2)
    ct <- ct[-c(1, 2), ]
    colnames(ct) <- nms
    print(nms)

    ct <- ct %>% dplyr::select(dplyr::any_of(nms)) %>% 
        tibble::as_tibble()

    df <- ct %>% 
        dplyr::mutate(TargetGenes = str_extract_all(TargetGenes, "\\(.+\\)")) %>% 
        dplyr::mutate(TargetGenes = unlist(TargetGenes)) %>% 
        tidyr::separate_longer_delim(cols = TargetGenes, delim = "),") %>% 
        dplyr::mutate(TargetGenes = str_remove_all(TargetGenes, "'|\\(|\\)")) %>% 
        tidyr::separate(col = TargetGenes, sep = ",", into = c("target", "weight")) %>% 
        dplyr::mutate(target = str_trim(target, side = "both"), 
            weight = str_trim(weight, side = "both"), 
            weight = as.numeric(weight), 
            weight = round(weight, 4)) %>% 
        dplyr::mutate(
            AUC = as.numeric(AUC), 
            NES = as.numeric(NES), 
            MotifSimilarityQvalue = as.numeric(MotifSimilarityQvalue), 
            OrthologousIdentity = as.numeric(OrthologousIdentity))
    
    return(df)
}

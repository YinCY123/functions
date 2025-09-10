addCytotrace <- function(sce, 
    cytotrace, 
    ...){
    vars <- cytotrace %>% colnames
    for(var in vars){
        colData(sce)[[var]] <- setNames(cytotrace[[var]], rownames(cytotrace))[rownames(colData(sce))]
    }
    return(sce)
}
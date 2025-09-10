get_ppi_hub <- function(g, top = 50, ...){
    # loading required package
    suppressPackageStartupMessages(require(tidygraph))
    
    df <- g %>% activate(nodes) %>% as.data.frame()
    methods <- df %>% colnames() %>% .[-1]
    
    hub_list <- lapply(methods, function(x){
        df %>% dplyr::arrange(!!sym(x)) %>% head(top) %>% dplyr::pull(name)
    })
    names(hub_list) <- methods
    
    return(hub_list)
}
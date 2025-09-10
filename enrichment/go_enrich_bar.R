go_enrich_bar <- function(obj, 
                          ont = "ALL", 
                          x_axis = "pvalue", 
                          y_axis = "Description", 
                          fill = "pvalue", 
                          ylab = NULL, 
                          ncol = 1,
                          str_width = 50, 
                          top = 10, 
                          file = NULL, 
                          width = 7, 
                          height = 7, 
                          scale = 1, 
                          low_color = "blue", 
                          high_color = "red", 
                          ...){
    # loading package
    suppressPackageStartupMessages(require(magrittr))
    suppressPackageStartupMessages(require(stringr))
    suppressPackageStartupMessages(require(ggplot2))
    
    # tbl <- obj@result
    tbl <- obj
    
    xlab <- switch(x_axis, 
        "pvalue" = "-log10(p value)", 
        "p.adjust" = "-log10(p.adjust)", 
        "qvalue" = "-log10(q value)", 
        x_axis)

    fill_lab <- switch(fill, 
        "pvalue" = "-log10(p value)", 
        "p.adjust" = "-log10(p.adjust)", 
        "qvalue" = "-log10(q value)", 
        fill)

    if(fill %in% c("pvalue", "p.adjust", "qvalue")){
        tbl <- tbl %>% 
            dplyr::mutate(fill = -log10(!!sym(fill)))
    }else{
        tbl <- tbl %>% 
            dplyr::mutate(fill = !!sym(fill))
    }
    
    logx <- ifelse(x_axis %in% c("pvalue", "p.adjust", "qvalue"), TRUE, FALSE)

    tbl <- tbl %>% 
        dplyr::mutate(y = str_wrap(!!sym(y_axis), width = str_width))
    if(logx){
        tbl <- tbl %>% dplyr::mutate(x = -log10(!!sym(x_axis)))
    }else{
        tbl <- tbl %>% dplyr::mutate(x = !!sym(x_axis))
    }

    if(ont == "ALL"){
        if(logx){
            tbl = tbl %>% 
                dplyr::group_by(ONTOLOGY) %>% 
                dplyr::arrange(!!sym(x_axis)) %>% 
                dplyr::slice_head(n = top)
        }else{
            tbl = tbl %>% 
                dplyr::group_by(ONTOLOGY) %>% 
                dplyr::arrange(desc(!!sym(x_axis))) %>% 
                dplyr::slice_head(n = top)
        }
    }else{
        if(logx){
            tbl = tbl %>% 
                dplyr::arrange(!!sym(x_axis)) %>% 
                dplyr::slice_head(n = top)
        }else{
            tbl = tbl %>% 
                dplyr::arrange(desc(!!sym(x_axis))) %>% 
                dplyr::slice_head(n = top)
        }
    }


    
    p <- tbl %>% 
        ggplot(aes(x, reorder(y, x, decreasing = F))) +
        geom_bar(aes(fill = fill), stat = "identity") +
        scale_x_continuous(name = xlab) +
        scale_y_discrete(name = ylab) +
        # scale_fill_continuous(name = fill_lab) +
        # theme_bw() +
        theme(strip.background = element_blank(), 
            panel.background = element_blank(),
            axis.line = element_blank(), 
            panel.border = element_rect(fill = NA, linewidth = 1), 
            panel.grid.major = element_line(linetype = "longdash", linewidth = 0.1, color = "grey")) +
        scale_fill_gradient(name = fill_lab, low = low_color, high = high_color)
    
    if(ont == "ALL"){
        p <- p + facet_wrap(vars(ONTOLOGY), ncol = ncol, scales = "free_y")
    }
    
    if(is.null(file)){
        return(p)
    }else{
        ggsave(filename = file, width = width, height = height, scale = scale, limitsize = FALSE)
    }
}
go_enrich_dot <- function(obj,
                          ont = "ALL",
                          ncol = 3, 
                          x_axis = "FoldEnrichment",
                          y_axis = "Description",
                          title = NULL,
                        #   logx = FALSE, 
                          color = "p.adjust", 
                        #   xlab = NULL, 
                          ylab = NULL,
                          point_size = "Count",
                        #   color_label = "-log10(p value)",
                          log_x = FALSE,
                          str_width = 40,
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
    logx <- ifelse(x_axis %in% c("pvalue", "p.adjust", "qvalue"), TRUE, FALSE)
    
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
    
    # wrap description length
    tbl <- tbl %>% 
        dplyr::mutate(y = str_wrap(!!sym(y_axis), width = str_width))
    
    # prepare x axis
    if(logx){
        tbl <- tbl %>% dplyr::mutate(x = -log10(!!sym(x_axis)))
    }else{
        tbl <- tbl %>% dplyr::mutate(x = !!sym(x_axis))
    }

    # 
    if(color %in% c("pvalue", "p.adjust", "qvalue")){
        tbl <- tbl %>% dplyr::mutate(logp = -log10(!!sym(color)))
    }

    # color label
    color_label <- switch(color, 
        "pvalue" = "-log10(p value)", 
        "p.adjust" = "-log10(p adjust)", 
        "qvalue" = "-log10(q value)", 
        x)
    
    # x lab
    if(is.null(xlab)){
        xlab = x_axis
    }
        
    p <- tbl %>% 
        ggplot(aes(x, reorder(y, x, decreasing = F))) +
        geom_point(aes(color = logp, size = !!sym(point_size))) +
        scale_x_continuous(name = xlab) +
        scale_y_discrete(name = ylab) +
        # scale_fill_continuous(name = fill_label) +
        theme_bw() +
        scale_color_gradient(name = color_label, low = low_color, high = high_color) +
        theme(panel.grid = element_line(linetype = 2, linewidth = 0.1, color = "grey"), 
            panel.border = element_rect(linewidth = 1.2, fill = NA), 
            panel.background = element_rect(fill = "white"), 
            strip.background = element_blank(), 
            strip.text.x = element_text(size = 10, face = "bold"), 
            axis.text = element_text(size = 10, face = "bold"), 
            axis.title = element_text(size = 12, face = "bold")) +
        ggtitle(label = title)

    if(ont == "ALL"){
        p <- p + facet_wrap(vars(ONTOLOGY), ncol = ncol, scales = "free_y")
    }
    
    if(is.null(file)){
        return(p)
    }else{
        ggsave(filename = file, width = width, height = height, scale = scale, limitsize = FALSE)
    }
}
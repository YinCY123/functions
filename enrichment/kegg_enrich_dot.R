kegg_enrich_dot <- function(obj,
    x = "pvalue",
    y = "Description",
    ylab = NULL,
    title = NULL,
    logx = TRUE,
    top_n = 20,
    str_width = 40,
    file = NULL,
    width = 7,
    height = 7,
    scale = 1,
    ...){
        # loading required packages
        require(magrittr)
        require(ggplot2)
        require(dplyr)
        require(stringr)

        xlab <- switch(x, 
            "pvalue" = "-log10(p value)", 
            "p.adjust" = "-log10(p.adjust)", 
            "qvalue" = "-log10(q value)", 
            x)

        # tmp <- obj@result
        tmp <- obj
        tmp <- tmp %>% dplyr::arrange(pvalue) %>% 
            dplyr::slice_head(n = top_n) %>% 
            dplyr::mutate(Description = str_wrap(!!sym(y), width = str_width))
        
        if(logx){
            tmp <- tmp %>% dplyr::mutate(xx = -log10(!!sym(x)))
        }else{
            tmp <- tmp %>% dplyr::mutate(xx = !!sym(x))
        }
        
        p <- tmp %>% 
            ggplot(aes(xx, reorder(!!sym(y), xx, decreasing = FALSE))) +
            geom_point(aes(color = xx, size = Count)) +
            scale_x_continuous(name = xlab) +
            scale_y_discrete(name = ylab) +
            scale_color_continuous(name = xlab, low = "blue", high = "red") +
            theme(panel.background = element_blank(),
                panel.border = element_rect(fill = NA, linewidth = 0.6, color = "black"), 
                panel.grid.major = element_line(linetype = 2, linewidth = 0.1, color = "grey"), 
                axis.text = element_text(size = 10, face = "bold"), 
                axis.title = element_text(size = 12, face = "bold")) +
            ggtitle(label = title)

        if(is.null(file)){
            return(p)
        }else{
            ggsave(plot = p, file = file, width = width, height = height, scale = scale)
        }
}

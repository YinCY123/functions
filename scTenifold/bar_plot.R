bar_plot <- function(x, 
    top_n = 10,
    fill_color = c("steelblue", "tomato"), 
    fill_title = "group",
    ...){

    # loading required packages
    library(magrittr)
    library(ggplot2)

    # selection plot data
    tmp <- rbind(
        x %>% dplyr::filter(Z > 0) %>% 
            dplyr::arrange(desc(Z)) %>% 
            dplyr::slice_head(n = top_n) %>% 
            dplyr::mutate(group = "up"), 
        x %>% dplyr::filter(Z < 0) %>% 
            dplyr::arrange(Z) %>% 
            dplyr::slice_head(n = top_n) %>% 
            dplyr::mutate(group = "down")
    )

    # x-axis limit
    lim <- c(-max(ceiling(abs(tmp$Z))), max(ceiling(abs(tmp$Z))))

    # x-axis breaks
    if(diff(lim)/5 > 3){
        x_breaks <- c(seq(floor(lim[1]), 0, 5), seq(5, ceiling(lim[2]), 5))
    }else{
        x_breaks <- seq(floor(lim[1]), ceiling(lim[2]), 1)
    }

    # fill color
    if(is.null(fill_color)){
        fill_color = c("steelblue", "tomato")
        names(fill_color) <- c("down", "up")
    }else{
        names(fill_color) <- c("down", "up")
    }
    

    p <- tmp %>% 
        ggplot(aes(Z, reorder(Gene, Z, decreasing = F))) +
        geom_bar(aes(fill = group), stat = "identity") +
        scale_x_continuous(name = "Z score", limit = lim, breaks = x_breaks) +
        scale_y_discrete(name = NULL) +
        scale_fill_manual(name = fill_title, values = fill_color) +
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA), 
            panel.grid.major.y = element_line(linetype = 2, color = "grey", linewidth = 0.2), 
            ...)
    
    return(p)
}

scatter_plot <- function(df, 
    x_name = "Z score", 
    y_name = "-log10(p value)",
    upper = 1, 
    lower = -1,
    point_colors = c("#663171", "grey", "#ea7428"), 
    point_size = 1, 
    linecolor = "blue",
    text_size = 3, 
    legend.key.size = 3,
    linewidth = 0.5,
    ...){
    # loading required packages
    library(magrittr)
    library(ggplot2)
    library(MetBrewer)
    library(ggrepel)


    # prepare plotting data
    df <- df %>% 
        dplyr::mutate(direction = ifelse(Z > 0, "Z > 0", "Z <= 0"), 
        sig = ifelse(Z > upper, "up", ifelse(Z < lower, "down", "ns")), 
        y = -log10(p.value))
        # dplyr::filter(is.finite(y))
    
    # anno
    anno_up <- df %>% dplyr::filter(Z>0) %>% dplyr::arrange(desc(Z)) %>% dplyr::slice_head(n = 5)
    print(anno_up)

    anno_down <- df %>% dplyr::filter(Z<0) %>% dplyr::arrange(Z) %>% dplyr::slice_head(n = 5)
    print(anno_down)

    anno <- rbind(anno_up, anno_down)

    if(!is.null(point_colors)){
        names(point_colors) = c("down", "ns", "up")
    }else{
        point_color = met.brewer("Java", 2)
        point_colors <- c(point_color[1], "grey", point_color[2])
        names(point_color) = c("down", "ns", "up")
    }
    p <- df %>% 
        dplyr::arrange(Z) %>% 
        ggplot(aes(Z, y)) +
        geom_point(aes(color = sig), size = point_size) +
        geom_text_repel(data = anno, aes(label = Gene), size = text_size) +
        geom_smooth(aes(x = Z, y = y), 
            se = F, 
            linetype = 2,
            linewidth = linewidth, 
            color = linecolor,
            # formula = "y ~ poly(x, 2)",
            method = "loess") +
        # facet_wrap(vars(direction), nrow = 1) +
        scale_x_continuous(name = x_name) +
        scale_y_continuous(name = y_name, expand = expansion(mult = c(0.01, 0.1))) +
        scale_color_manual(name = NULL, values = point_colors) +
        guides(color = guide_legend(override.aes = list(size = legend.key.size))) +
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA), 
            strip.background = element_blank(), 
            ...)

    return(p)
}

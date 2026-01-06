dotPlot <- function(sces, markers, group, 
    x = "celltype", 
    y = "markers",
    scale = TRUE, 
    colors = c("blue", "white", "red"),
    palette = "Spectral",
    color_title = "Mean Expression", 
    size_title = "Percent Expression",
    size_range = c(0, 10),
    midpoint = 0,
    ...){
    # loading required packages
    suppressPackageStartupMessages(library(scater))
    suppressPackageStartupMessages(library(SingleCellExperiment))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(magrittr))

    tmp <- markers %>% unlist
    df <- data.frame(
        markers = tmp, 
        celltype = str_remove(names(tmp), "[:digit:]{1}$"), 
        row.names = NULL
    )

    mtx <- scater::plotDots(
        object = sces, 
        features = intersect(unlist(markers), rownames(sces)), 
        group = group, 
        scale = scale, 
        center = scale
    )$data

    colnames(mtx) <- c("markers", "celltype", "percent", "avg")
    mtx <- mtx %>% dplyr::mutate(
        markers = factor(markers, levels = unique(df$markers)), 
        celltype = factor(celltype, levels = unique(df$celltype))
    )

    # color
    if(is.null(colors)){
        colors <- c("blue", "white", "red")
    }

    p <- mtx %>% 
        # dplyr::filter(!is.na(celltype)) %>% 
        ggplot(aes(!!sym(x), !!sym(y))) +
        geom_point(aes(size = percent, color = avg)) +
        scale_x_discrete(name = NULL) +
        scale_y_discrete(name = NULL) +
        scale_color_gradient2(name = color_title, 
            low = colors[[1]], 
            mid = colors[[2]], 
            high = colors[[3]], 
            midpoint = midpoint) +
        scale_radius(name = size_title, range = size_range, label = scales::percent) +
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA), 
            panel.grid.major = element_line(linetype = 2, color = "grey", linewidth = 0.1),
            ...)
    
    return(p)
}
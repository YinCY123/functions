plotDots <- function(x, features, 
    group_by = "celltype", 
    group_levels = NULL, 
    colors = c("blue", "white", "red"), 
    midpoint = 0, 
    point_range = c(0, 6), 
    ...){

    # loading required packages
    suppressPackageStartupMessages(library(SingleCellExperiment))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(tidyr))
    suppressPackageStartupMessages(library(dplyr))

    # filter features
    features <- unlist(features, use.name = F)
    features <- intersect(features, rownames(x))

    df <- makePerCellDF(x, features = features, use.coldata = TRUE, use.dimred = FALSE) %>% 
        tidyr::pivot_longer(cols = dplyr::any_of(features), names_to = "symbol", values_to = "expr")

    input <- df %>% dplyr::group_by(symbol, !!sym(group_by)) %>% 
        dplyr::summarise(avg = mean(expr), 
            per = mean(expr > 0))
    
    input <- input %>% 
        dplyr::group_by(symbol) %>% 
        dplyr::mutate(avg = scale(avg, scale = FALSE, center = TRUE))

    # # grps <- x$group_by %>% unique()
    # mtx <- x[features, ] %>% logcounts %>% as.matrix
    # mtx <- apply(mtx, 1, function(x){
    #     x - mean(x)
    # }) %>% t
    
    # input <- mtx %>% as.data.frame() %>% 
    #     tibble::rowid_to_column("symbol") %>% 
    #     tidyr::pivot_longer(cols = -symbol, names_to = "barcodes", values_to = "expr")
    # meta <- x %>% colData %>% 
    #     as.data.frame() %>% 
    #     tibble::rownames_to_column("barcodes")

    # input <- dplyr::left_join(input, meta, by = "barcodes")
    # input <- input %>% dplyr::group_by(symbol, !!sym(group_by)) %>% 
    #     dplyr::summarise(avg = mean(expr), 
    #         per = mean(expr > 0))
    
    if(!is.null(group_levels) & length(group_levels) > 0){
        input[[group_by]] <- factor(input[[group_by]], levels = group_levels)
    }

    p <- input %>% 
        dplyr::filter(!is.na(celltype)) %>% 
        dplyr::arrange(celltype, desc(avg)) %>% 
        ggplot(aes(symbol, !!sym(group_by))) +
        geom_point(aes(size = per, color = avg)) +
        scale_color_gradient2(name = "Scaled Expression", low = colors[[1]], mid = colors[[2]], high = cols[[3]], midpoint = midpoint) +
        scale_x_discrete(name = NULL) +
        scale_y_discrete(name = NULL) +
        scale_radius(name = "Percent Expression", range = point_range, label = scales::percent)  +
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA), 
            ...)

    return(p)
}

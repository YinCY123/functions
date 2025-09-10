heatmap <- function(x, 
    scale = TRUE, 
    colors = NULL,
    ncolors = 100,
    annotation_col = NULL, 
    annotation_row = NULL, 
    cluster_col = FALSE,
    cluster_row = TRUE,
    fontsize = 10, 
    fontsize_row = 10, 
    fontsize_col = 10,
    angle_col = 315, 
    show_rownames = TRUE, 
    show_colnames = FALSE,
    clustering_method = "ward.D2",
    main = NULL, 
    width = 7, 
    height = 6,
    expand = 1,
    file = NA,
    ...){
    # loading packages
    suppressPackageStartupMessages(require(pheatmap))
    suppressPackageStartupMessages(require(RColorBrewer))

    # color
    if(is.null(colors)){
        colors <- rev(colorRampPalette(colors = brewer.pal(11, "Spectral"))(ncolors))
    }else{
        colors <- colorRampPalette(colors = colors)(ncolors)
    }

    # adjust width and height
    width = width * expand
    height = height * expand

    # scale
    scale = ifelse(scale, "row", "none")

    pheatmap(mat = x,
        color = colors,
        scale = scale,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        cluster_cols = cluster_col,
        cluster_rows = cluster_row, 
        angle_col = angle_col,
        show_rownames = show_rownames,
        show_colnames = show_colnames, 
        fontsize = fontsize,
        fontsize_row = fontsize_row, 
        fontsize_col = fontsize_col, 
        clustering_method = clustering_method,
        height = height,
        width = width,
        file = file,
        ...)
}

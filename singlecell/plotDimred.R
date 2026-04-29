plotDimred <- function(sces, 
    dimred = "UMAP", 
    features = NULL,
    group_by = "celltype",
    arrow_length = 3,
    text_by = "celltype",
    text_size = 2.5, 
    max_overlaps = 20,
    label_size = 2.5,
    scattermore = TRUE, 
    point_size = 1,
    colors = NULL,
    color_title = NULL,
    legend_ncol = 1,
    seed = 101,
    x_nudge = 1,
    y_nudge = -0.5,
    ...){

        # loading required packages
        suppressPackageStartupMessages(library(ggplot2))
        suppressPackageStartupMessages(library(ggarrow))
        suppressPackageStartupMessages(library(scran))
        suppressPackageStartupMessages(library(stringr))
        suppressPackageStartupMessages(library(ggrepel))
        suppressPackageStartupMessages(library(scran))
        suppressPackageStartupMessages(library(dplyr))
        suppressPackageStartupMessages(library(tidyr))
        suppressPackageStartupMessages(library(scattermore))

        # extract data from single cell experiment object
        df <- makePerCellDF(x = sces, use.coldata = TRUE, use.dimred = TRUE, features = features)

        if(!is.null(features) & length(features) > 1){
            df <- df %>% tidyr::pivot_longer(cols = dplyr::any_of(features), names_to = "symbol", values_to = "expr")
        }

        # color
        cols <- c("#cb2426", "#ea3b21", "#aa1a7d", "#d83890", "#ed6b9f", 
            "#d64b23", "#f08c44", "#ca9424", "#f6be2a", "#ffd4af", 
            "#4c54a0", "#68559d", "#7f7cb6", "#9f9ac4", "#bebed8", 
            "#2573b4", "#0277b9", "#4392c4", "#6baed5", "#9bc9dd", 
            "#035830", "#148843", "#3bab5a", "#76c277", "#a2d59b")
        if(is.null(colors) & length(unique(df[[group_by]])) <= length(cols)){
            set.seed(seed)
            colors <- sample(cols, length(unique(df[[group_by]])), replace = FALSE)
        }else if(is.null(colors) & length(unique(df[[group_by]])) > length(cols)){
            set.seed(seed)
            colors <- sample(cols, length(unique(df[[group_by]])), replace = TRUE)        
        }

        # create cell and text location data
        if(str_to_upper(dimred) == "UMAP"){
            cell_loc <- df %>% 
                dplyr::group_by(!!sym(group_by)) %>% 
                dplyr::summarise(x = median(UMAP.1), y = median(UMAP.2))

            # create arrow needs data
            rg1 <- range(df$UMAP.1) %>% diff
            rg2 <- range(df$UMAP.2) %>% diff
            arrow_df <- data.frame(
                x = c(seq(min(df$UMAP.1), min(df$UMAP.1) + arrow_length, length.out = 100), 
                    rep(min(df$UMAP.1), 100)), 
                y = c(rep(min(df$UMAP.2), 100), 
                    seq(min(df$UMAP.2), min(df$UMAP.2) + arrow_length, length.out = 100)), 
                group = rep(c(1, 2), each = 100))
            arrow_txt <- data.frame(
                x = c(min(df$UMAP.1) + x_nudge, min(df$UMAP.1) + y_nudge), 
                y = c(min(df$UMAP.2) + y_nudge, min(df$UMAP.2) + x_nudge), 
                label = c("UMAP 1", "UMAP 2"), 
                angle = c(0, 90)
            )

            # visualization
            p <- df %>% 
                ggplot(aes(UMAP.1, UMAP.2)) +
                geom_arrow(data = arrow_df, aes(x = x, y = y, group = group)) +
                geom_text(data = arrow_txt, aes(x, y, label = label, angle = angle), size = label_size) +
                scale_color_manual(name = color_title, values = colors, na.value = "grey90") +
                guides(color = guide_legend(override.aes = list(size = 3), ncol = legend_ncol)) +
                coord_fixed() +
                theme(panel.background = element_blank(), 
                    panel.border = element_blank(), 
                    axis.title = element_blank(), 
                    axis.text = element_blank(), 
                    axis.ticks = element_blank(), 
                    ...)
            
            if(scattermore){
                p <- p + geom_scattermore(aes(color = !!sym(group_by)), size = point_size)
            }else{
                p <- p + geom_point(aes(color = !!sym(group_by)), size = point_size)
            }

            if(!is.null(text_by)){
                cell_loc <- df %>% dplyr::group_by(!!sym(text_by)) %>% dplyr::summarise(x = median(UMAP.1), y = median(UMAP.2))

                p <- p +
                    geom_text_repel(data = cell_loc, aes(x, y, label = !!sym(text_by)), size = text_size)
            }
                
        }else if(str_to_upper(dimred) == "TSNE"){
            cell_loc <- df %>% 
                dplyr::group_by(!!sym(group_by)) %>% 
                dplyr::summarise(x = median(TSNE.1), y = median(TSNE.2))

            # create arrow needs data
            rg1 <- range(df$TSNE.1) %>% diff
            rg2 <- range(df$TSNE.2) %>% diff
            arrow_df <- data.frame(
                x = c(seq(min(df$TSNE.1), min(df$TSNE.1) + arrow_length, length.out = 100), 
                    rep(min(df$TSNE.1), 100)), 
                y = c(rep(min(df$TSNE.2), 100), 
                    seq(min(df$TSNE.2), min(df$TSNE.2) + arrow_length, length.out = 100)), 
                group = rep(c(1, 2), each = 100))
            arrow_txt <- data.frame(
                x = c(min(df$TSNE.1) + x_nudge, min(df$TSNE.1) + y_nudge), 
                y = c(min(df$TSNE.2) + y_nudge, min(df$TSNE.2) + x_nudge), 
                label = c("TSNE 1", "TSNE 2"), 
                angle = c(0, 90)
            )

            # visualization
            p <- df %>% 
                ggplot(aes(TSNE.1, TSNE.2)) +
                # geom_scattermore(aes(color = !!sym(group_by)), size = point_size) +
                # geom_point(aes(color = !!sym(group_by)), size = point_size) +
                geom_arrow(data = arrow_df, aes(x = x, y = y, group = group)) +
                geom_text(data = arrow_txt, aes(x, y, label = label, angle = angle), size = label_size) +
                geom_text(data = cell_loc, aes(x, y, label = !!sym(group_by)), size = text_size) +
                scale_color_manual(name = color_title, values = colors, na.value = "grey90") +
                guides(color = guide_legend(override.aes = list(size = 3), ncol = legend_ncol)) +
                coord_fixed() +
                theme(panel.background = element_blank(), 
                    panel.border = element_blank(), 
                    axis.title = element_blank(), 
                    axis.text = element_blank(), 
                    axis.ticks = element_blank(), 
                    ...)

            if(scattermore){
                p <- p + geom_scattermore(aes(color = !!sym(group_by)), size = point_size)
            }else{
                p <- p + geom_point(aes(color = !!sym(group_by)), size = point_size)
            }

            if(!is.null(text_by)){
                cell_loc <- df %>% dplyr::group_by(!!sym(text_by)) %>% dplyr::summarise(x = median(UMAP.1), y = median(UMAP.2))

                p <- p +
                    geom_text_repel(data = cell_loc, aes(x, y, label = !!sym(text_by)), size = text_size)
            }
        }else{
            message("The dimensions you provide is not supported...")
        }
        return(p)
}

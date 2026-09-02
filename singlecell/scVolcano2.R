scVolcano <- function(tbl, 
    group_by = "celltype", 
    facet_by = "celltype", 
    ncol = 5,
    x = "summary.logFC", 
    y = "y", 
    sign = "summary.logFC", 
    log10 = "p.value",
    colors = c("#7F3B08", "#2D004B"),
    size_rng = c(0, 3),
    add_anno = TRUE,
    top_anno = 5,
    text_size = 2,
    legend.position = "none", 
    ...){
        # loading required packages
        library(magrittr)
        library(ggplot2)
        library(ggrepel)
        library(circlize)

        # prepare plot data
        tbl <- tbl %>% 
            dplyr::mutate(y = sign(!!sym(sign)) * -log10(!!sym(log10)))

        # prepare annotation data
        anno <- rbind(
            tbl %>% 
                dplyr::group_by(!!sym(group_by)) %>% 
                dplyr::filter(y > 0) %>% 
                dplyr::arrange(desc(y)) %>% 
                dplyr::slice_head(n = top_anno),
            tbl %>% 
                dplyr::filter(y < 0) %>% 
                dplyr::group_by(!!sym(group_by)) %>% 
                dplyr::arrange(y) %>% 
                dplyr::slice_head(n = top_anno)
        )

        # prepare point color vector
        breaks <- c(min(tbl$y)/2, 0, max(tbl$y)/2)
        col_fun <- colorRamp2(breaks, colors = c(colors[2], "white", colors[1]))
        cols <- col_fun(seq(min(tbl$y), max(tbl$y), length.out = 100))

        p <- tbl %>% 
            ggplot(aes(summary.logFC, y)) +
            geom_point(aes(color = y, size = abs(summary.logFC))) +
            geom_hline(yintercept = 0) +
            facet_wrap(vars(celltype), ncol = ncol, scale = "free") +
            scale_color_gradientn(colors = cols, limits = c(min(tbl$y)/2, max(tbl$y)/2)) +
            scale_x_continuous(name = "log2 (Fold Change)") +
            scale_y_continuous(name = "sign(log2FC) * -log10(p value)") +
            scale_radius(range = size_rng) +
            theme(panel.background = element_blank(), 
                panel.border = element_rect(fill = NA), 
                legend.position = legend.position, 
                ...)

        if(add_anno){
            p <- p +
                geom_text_repel(data = anno, aes(label = symbol), size = text_size)
        }

        return(p)
    }

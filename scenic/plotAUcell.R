plotAUcell <- function(aucell, 
    sces, 
    dir, 
    width = 10, 
    height = 6, 
    group, 
    group_order, 
    trim = 0.1,
    prob = 0.5, 
    top_anno = 5, 
    size_range = c(0, 6),
    legend.position = "top", 
    n_col = 11,
    palette = "Spectral",
    ...){
    # loading required packages
    suppressPackageStartupMessages(library(magrittr))
    suppressPackageStartupMessages(library(rlang))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(ggrepel))

    # read aucell data
    if(is.character(aucell)){
        aucell = read.table(aucell, sep = ",", header = TRUE, check.names = FALSE) %>% 
            tidyr::pivot_longer(cols = -Cell, names_to = "regulon", values_to = "activity")
    }else{
        aucell <- aucell %>% tidyr::pivot_longer(cols = -Cell, names_to = "regulon", values_to = "activity")
    }

    cell_meta <- sces %>% colData %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("barcode") %>% 
        dplyr::select(barcode, !!sym(group), celltype)

    # merge aucell and cell meta
    input <- dplyr::left_join(aucell, cell_meta, by = c("Cell" = "barcode"))
    
    # summary
    input <- input %>% 
        dplyr::group_by(celltype, !!sym(group), regulon) %>% 
        dplyr::summarise(avg = mean(activity, trim = trim), 
            percent = mean(activity > quantile(input$activity, probs = prob), trim = trim))
    
    # anno
    anno <- input %>% 
        dplyr::group_by(!!sym(group)) %>% 
        dplyr::arrange(desc(avg)) %>% 
        dplyr::slice_head(n = top_anno)
    
    # color
    cols <- colorRampPalette(RColorBrewer::brewer.pal(n_col, palette))(100) %>% rev()

    # plot
    p1 <- input %>% 
        dplyr::mutate(group = factor(!!sym(group), levels = group_order)) %>% 
        ggplot(aes(reorder(regulon, avg, decreasing = TRUE), avg)) +
        geom_point(aes(size = percent, color = avg)) +
        geom_text_repel(data = anno, aes(regulon, avg, label = regulon), size = 3) +
        facet_wrap(vars(group), nrow = 1) +
        scale_color_gradientn(name = "Mean Activity", colors = cols) +
        scale_radius(name = "Percent Activity", range = size_range, label = scales::percent) +
        scale_x_discrete(name = NULL, expand = c(0.03, 0.03)) +
        scale_y_continuous(name = "Regulon Activity") +
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA), 
            panel.grid.major.y = element_line(linetype = 2, color = "grey", linewidth = 0.1), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            strip.background = element_blank(), 
            strip.text = element_text(size = 12, face = "bold"), 
            legend.position = legend.position)

    dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
    file1 <- paste0(dir, "regulon_activity_bubble_01.pdf")
    ggsave(plot = p1, file = file1, width = width, height = height)

    p2 <- input %>% 
        ggplot(aes(reorder(regulon, avg), !!sym(group))) +
        geom_point(aes(size = percent, color = avg)) +
        scale_x_discrete(name = NULL) +
        scale_y_discrete(name = NULL) +
        scale_color_gradientn(name = "Mean Activity", colors = cols) +
        scale_radius(name = "Percent Activity", range = c(0,6), label = scales::percent) +
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA), 
            panel.grid.major = element_line(linetype = 2, color = "grey", linewidth = 0.2), 
            legend.position = "top", 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6), 
            ...)
    
    file2 <- paste0(dir, "regulon_activity_bubble_02.pdf")
    ggsave(plot = p2, file = file2, width = width, height = height)
}

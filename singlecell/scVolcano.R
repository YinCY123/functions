scVolcano <- function(fmarkers,
    top_n = 5, 
    sort_by = "summary.logFC",
    vars = c("symbol", "p.value", "FDR", "summary.logFC"),
    group_by = "direction",
    add_gene = TRUE, 
    text_size = 3, 
    label_size = 3,
    colors = c("tomato", "steelblue"),
    y_lab = "log2 (Fold Change)",
    ...){

    suppressPackageStartupMessages(library(magrittr))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(ggrepel))

    # process data
    vars <- append(vars, "direction")
    fmarkers <- lapply(fmarkers, function(y){
        y %>% as.data.frame() %>% 
            tibble::rownames_to_column("symbol") %>% 
            tibble::as_tibble() %>% 
            dplyr::mutate(direction = ifelse(!!sym(sort_by) > 0, "up", "down")) %>% 
            dplyr::select(dplyr::all_of(vars))
    })

    # construct plot data
    cells <- names(fmarkers)
    fmarkers <- lapply(cells, function(y){
        fmarkers[[y]] %>% dplyr::mutate(celltype = y)
    })
    df <- do.call(rbind, fmarkers)

    # construct annotation data
    anno <- df %>% 
        dplyr::group_by(celltype, direction) %>% 
        dplyr::arrange(desc(abs(!!sym(sort_by)))) %>% 
        dplyr::slice_head(n = 5)

    # construct rect data
    rf <- data.frame(
        celltype = cells, 
        y = 0
    )

    p <- df %>% 
        ggplot(aes(celltype, !!sym(sort_by))) +
        geom_hline(yintercept = 0, linetype = 1, color = "black") +
        geom_jitter(aes(size = abs(!!sym(sort_by)), color = direction)) +
        geom_label(data = rf, aes(x = celltype, y = y, label = celltype), size = text_size) +
        geom_text_repel(data = anno, aes(celltype, !!sym(sort_by), label = symbol), size = label_size) +
        scale_x_discrete(name = NULL) +
        scale_y_continuous(name = y_lab) +
        scale_color_manual(values = c("up" = colors[1], "down" = colors[2])) +
        scale_radius(range = c(0.1, 3)) +
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            legend.position = "none",
            ...)

    # extract x-axis break point information
    
    
    return(p)
}

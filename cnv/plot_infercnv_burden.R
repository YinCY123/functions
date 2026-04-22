plot_infercnv_burden <- function(obj, 
    sces,
    group_by = "celltype", 
    text_size = 3, 
    text_color = "red",
    viridis_option = "B",
    x = "UMAP.1", 
    y = "UMAP.2", 
    legend_breaks = NULL, 
    legend_labels = NULL,
    legend_position = "top",
    summary_fun = "quantile", 
    bins = 100,
    fun_args = list(probs = 0.8),
    ...){

    # loading required packages
    library(magrittr)
    library(qs)
    library(SingleCellExperiment)
    library(ggplot2)
    library(viridis)
    library(rlang)

    # extract cnv matrix
    cnv_mtx <- slot(obj, "expr.data")

    # calculate per cell cnv burden
    cnv_burden <- apply(cnv_mtx, 2, function(x){
        sum(abs(x - 1))
    })

    # get cell metadata
    cell_meta <- makePerCellDF(sces, use.coldata = TRUE, use.dimred = TRUE) %>% tibble::rownames_to_column("barcodes")
    cell_meta <- cell_meta %>% dplyr::filter(barcodes %in% names(cnv_burden))

    cell_loc <- cell_meta %>% 
        dplyr::group_by(!!sym(group_by)) %>% 
        dplyr::summarise(x = median(!!sym(x)), y = median(!!sym(y)))
    cell_meta <- cell_meta %>% dplyr::mutate(cnv_burden = cnv_burden[barcodes])

    p <- cell_meta %>% 
        ggplot(aes(!!sym(x), !!sym(y))) +
        stat_summary_hex(aes(z = cnv_burden), 
            bins = bins, 
            fun = summary_fun, 
            fun.args = fun_args, 
            linewidth = 0.1, 
            show.legend = TRUE) +
        geom_text(data = cell_loc, aes(x, y, label = !!sym(group_by)), size = text_size, color = text_color) +
        viridis::scale_fill_viridis(name = "CNV", option = viridis_option, 
            # breaks = legend_breaks, 
            # labels = legend_labels
            ) +
        scale_x_continuous(name = x) +
        scale_y_continuous(name = y) +
        theme(panel.background = element_rect(fill = NA, color = "black"), 
            legend.position = legend_position, 
            ...)
    return(p)
}

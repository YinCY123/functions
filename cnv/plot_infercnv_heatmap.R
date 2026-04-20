plot_infercnv_heatmap <- function(x, 
    sces, 
    ref_cells, 
    observed_cells,
    out_dir = NULL,
    height = 8, 
    width = 12,
    height_ref = 3, 
    height_cells = 9, 
    order_cell_by = "celltype",
    cell_annotation_color = NULL,
    title_padding = 4,
    heatmap_colors = c("navy", "white", "firebrick"), 
    heatmap_breaks = c(0, 1, 2),
    chromes = c(1:22, "X", "Y"),
    fontsize_row = 8, 
    fontsize_col = 10,
    ...){
        # loading required heatmap
        library(magrittr)
        library(infercnv)
        library(circlize)
        library(ComplexHeatmap)
        library(dplyr)
        library(tidyr)
        library(tibble)

        # heatmap title padding
        ht_opt$TITLE_PADDING = unit(title_padding, "mm")

        # extract gene order
        gene_order <- slot(x, "gene_order")
        gene_order <- gene_order %>% 
            dplyr::mutate(chr = factor(chr, levels = chromes)) %>% 
            dplyr::filter(!is.na(chr)) %>% 
            tibble::rownames_to_column("symbol") %>% 
            dplyr::arrange(chr, start)

        # extract cnv matrix
        cnv_mtx <- slot(x, "expr.data")

        # extract cell metadata
        cell_meta <- sces %>% colData %>% 
            as.data.frame %>% 
            tibble::rownames_to_column("barcodes") %>% 
            dplyr::select(barcodes, !!sym(order_cell_by)) %>% 
            dplyr::filter(barcodes %in% colnames(cnv_mtx)) %>% 
            dplyr::arrange(!!sym(order_cell_by))

        ref_cell_barcodes <- cell_meta %>% dplyr::filter(!!sym(order_cell_by) %in% ref_cells) %>% dplyr::pull(barcodes)
        observed_cell_barcodes <- cell_meta %>% dplyr::filter(!!sym(order_cell_by) %in% observed_cells) %>% dplyr::pull(barcodes)

        # construct ref amd observed cnv matrix
        cnv_mtx_ref <- cnv_mtx[gene_order$symbol, ref_cell_barcodes]
        cnv_mtx_observed <- cnv_mtx[gene_order$symbol, observed_cell_barcodes]

        if(is.null(cell_annotation_color)){
            cell_annotation_color <- RColorBrewer::brewer.pal(9, "Set1")
        }

        # ref heatmap annotation
        ht_anno_ref <- HeatmapAnnotation(
            df = cell_meta %>% dplyr::filter(!!sym(order_cell_by) %in% ref_cells) %>% tibble::column_to_rownames("barcodes"), 
            col = list(
                assign(order_cell_by, setNames(cell_annotation_color[1:length(unique(ref_cells))], unique(ref_cells))
            ), 
            name = " ", 
            show_legend = F, 
            which = "row"
        ))

        ht_anno_cells <- HeatmapAnnotation(
            df = cell_meta %>% dplyr::filter(!!sym(order_cell_by) %in% observed_cells) %>% tibble::column_to_rownames("barcodes"), 
            col = list(
                assign(order_cell_by, setNames(cell_annotation_color[1:length(unique(observed_cells))], unique(observed_cells))
            ), 
            name = " ", 
            show_legend = TRUE, 
            which = "row"
        ))

        col_fun <- colorRamp2(breaks = heatmap_breaks, colors = heatmap_colors)
        
        # reference heatmap
        p1 <- Heatmap(t(cnv_mtx_ref), 
            col = col_fun,
            name = "CNV",
            cluster_rows = TRUE, 
            cluster_columns = F, 
            use_raster = TRUE, 
            show_row_names = F, 
            show_column_names = F,
            column_split = gene_order$chr, 
            column_gap = unit(0, "mm"), 
            border = TRUE, 
            border_gp = gpar(col = "black", lwd = 0.2), 
            column_title_gp = gpar(fontsize = fontsize_column), 
            show_row_dend = F, 
            row_title = paste0("Reference (", str_c(unique(ref_cells), ";"), ")"), 
            row_title_side = "right", 
            row_title_gp = gpar(fontsize = fontsize_row, face = "bold"), 
            left_annotation = ht_anno_ref, 
            heatmap_legend_param = list(
                title_position = "lefttop", 
                legend_direction = "horizontal", 
                nrow = 1
            ), 
            raster_by_magick = F, 
            height = unit(height_ref, "cm")
        )

        # observed heatmap
        p1 <- Heatmap(t(cnv_mtx_cell), 
            col = col_fun, 
            name = "CNV", 
            cluster_rows = F, 
            cluster_columns = F, 
            use_raster = TRUE, 
            show_row_names = F, 
            show_column_names = F, 
            column_split = gene_order$chr, 
            column_gap = unit(0, "mm"), 
            border = TRUE, 
            border_gp = gpar(col = "black", lwd = 0.2), 
            column_title_gp = gpar(fontsize = fontsize_col), 
            show_row_dend = F, 
            row_title = "Observed Cell", 
            row_title_side = "right", 
            row_title_gp = gpar(fontsize = fontsize_row, face = "bold"), 
            left_annotation = ht_anno_cell, 
            heatmap_legend_param = list(
                title_position = "lefttop", 
                legend_direction = "horizontal", 
                nrow = 1
            ), 
            raster_by_magick = FALSE, 
            height = unit(height_cells, "cm")
        )


        # generate Heatmap list object
        p <- p1 %v% p2

        if(is.null(out_dir)){
            out_dir = getwd()
        }
        out_dir <- ifelse(grepl("/$", out_dir), out_dir, paste0(out_dir, "/"))

        file <- paste0(out_dir, "infercnv_heatmap.pdf")
        pdf(file = file, width = width, height = height)
        draw(p, merge_legend = TRUE, heatmap_legend_side = "top", ...)
        dev.off()
    }
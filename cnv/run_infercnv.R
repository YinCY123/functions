run_infercnv <- function(sces, 
    celltype_col = "celltype", 
    cells = NULL, 
    reference_cells = NULL, 
    out_dir = NULL, 
    species = "human",
    bm_dataset = "hsapiens_gene_ensembl", 
    ...){

    # loading required packages
    library(magrittr)
    library(SingleCellExperiment)
    library(biomaRt)
    library(stringr)
    library(infercnv)


    # subset SingleCellExperiment
    if(is.null(cells)){
        cells <- sces[[celltype_col]] %>% unique
        sub <- sces
    }else(
        sub <- sces[, sces[[celltype_col]] %in% cells]
    )

    # check reference cells
    if(is.null(reference_cells)){
        message("Please set reference cells...")
    }

    # extract matrix
    mtx <- sub %>% counts

    # annotation data
    anno <- sub %>% colData %>% 
        as.data.frame %>% 
        # tibble::rownames_to_column("barcode")
        dplyr::select(!!sym(celltype_col))

    # gene order
    message("Using biomaRt to retrive gene position...")

    hmart <- useEnsembl("ensembl", 
        dataset = bm_dataset, 
        mirror = "asia")

    # get all symbols
    if(str_to_lower(species) == "human"){
        library(org.Hs.eg.db)
        symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    }else if(str_to_lower(species) == "mouse"){
        library(org.Mm.eg.db)
        symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    }else(
        message("The species is not supported...")
    )

    gene_order <- getBM(
        attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"), 
        filters = "hgnc_symbol", 
        values = symbols, 
        mart = hmart
    )

    gene_order <- gene_order %>% 
        dplyr::filter(!duplicated(hgnc_symbol)) %>% 
        tibble::column_to_rownames("hgnc_symbol")
    
    obj <- CreateInfercnvObject(
        raw_counts_matrix = mtx, 
        gene_order_file = gene_order, 
        annotations_file = anno, 
        ref_group_names = reference_cells
    )

    if(is.na(out_dir)){
        out_dir = getwd()
    }

    tmp <- infercnv::run(
        infercnv_obj = obj, 
        cutoff = 0.1, 
        out_dir = out_dir, 
        cluster_by_groups = TRUE, 
        plot_steps = FALSE, 
        denoise = TRUE, 
        HMM = FALSE, 
        no_prelim_plot = TRUE, 
        png_res = 300
    )
}

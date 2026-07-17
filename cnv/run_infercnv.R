run_infercnv <- function(sces, 
    celltype_col = "celltype", 
    cells = NULL, 
    assay.type = "counts",
    out_dir = NULL, 

    num_cells = 50000, 

    # gene orders
    gtf = NULL, 

    # biomaRt args
    species = "human", 
    bm_dataset = "hsapiens_gene_ensembl", 
    # host = "https://may2025.archive.ensembl.org", 
    mirror = "asia", 

    # infercnv args
    reference_cells = NULL, 
    max_cells_per_group = NULL, 
    min_max_counts_per_cell = c(100, +Inf), 
    chr_exclude = c("chrX", "chrY", "chrM"), 
    
    cutoff = 0.1, 
    cluster_by_groups = TRUE, 
    plot_steps = FALSE, 
    denoise = TRUE, 
    HMM = FALSE, 
    no_prelim_plot = TRUE, 
    png_res = 300, 
    ...){

    # loading required packages
    suppressPackageStartupMessages(library(magrittr))
    suppressPackageStartupMessages(library(rlang))
    suppressPackageStartupMessages(library(SingleCellExperiment))
    suppressPackageStartupMessages(library(biomaRt))
    suppressPackageStartupMessages(library(qs2))
    # suppressPackageStartupMessages(library(rtracklayer))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(infercnv))

    # subset SingleCellExperiment
    if(is.null(cells)){
        cells <- sces[[celltype_col]] %>% unique
        sub <- sces
    }else{
        sub <- sces[, sces[[celltype_col]] %in% cells]
    }
    message(paste0("number of genes: ", nrow(sub), "; number of cells: ", ncol(sub), " after subset celltype."))


    # sample down if number of cells more than 50,000
    set.seed(101)
    if(ncol(sub) > num_cells){
        factors <- 1/(ncol(sub)/num_cells)
        sub <- sub[, sample(ncol(sub), ncol(sub)*factors)]
    }
    message(paste0("number of genes: ", nrow(sub), "; number of cells: ", ncol(sub), " after down sample."))

    # check reference cells
    if(is.null(reference_cells)){
        message("Please set reference cells...")
    }

    if(is.null(out_dir)){
        out_dir = getwd()
    }

    out_dir <- ifelse(grepl("/$", out_dir), out_dir, paste0(out_dir, "/"))


    # extract matrix
    mtx <- assay(sub, assay.type)
    message(paste0("number of genes: ", nrow(mtx), " and number of cells: ", ncol(mtx), " used for downstream analysis..."))

    # annotation data
    anno <- sub %>% colData %>% 
        as.data.frame %>% 
        # tibble::rownames_to_column("barcode")
        dplyr::select(!!sym(celltype_col))

    # gene order
    message("Using biomaRt to retrive gene position...")
    hmart <- useEnsembl("ensembl", 
        dataset = bm_dataset, 
        mirror = mirror)
    # hmart <- useMart(
    #     biomart = "ENSEMBL_MART_ENSEMBL", 
    #     dataset = bm_dataset, 
    #     host = host
    # )

    # get all symbols
    if(str_to_lower(species) == "human"){
        library(org.Hs.eg.db)
        symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    }else if(str_to_lower(species) == "mouse"){
        library(org.Mm.eg.db)
        symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")
    }else{
        message("The species is not supported...")
    }

    gene_order <- getBM(
        attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"), 
        filters = "hgnc_symbol", 
        values = symbols, 
        mart = hmart
    )

    gene_order <- gene_order %>% 
        dplyr::filter(!duplicated(hgnc_symbol)) %>% 
        tibble::column_to_rownames("hgnc_symbol")
    qs_save(gene_order, file = paste0(out_dir, "gene_order.qs2"))
    print(head(gene_order))

    message("create infercnv object...")
    obj <- CreateInfercnvObject(
        raw_counts_matrix = mtx, 
        gene_order_file = gene_order, 
        annotations_file = anno, 
        ref_group_names = reference_cells
    )

    message("run infercnv...")
    tmp <- infercnv::run(
        infercnv_obj = obj, 
        cutoff = cutoff, 
        out_dir = out_dir, 
        cluster_by_groups = cluster_by_groups, 
        plot_steps = plot_steps, 
        denoise = denoise, 
        HMM = HMM, 
        no_prelim_plot = no_prelim_plot, 
        png_res = png_res, 
        ...
    )
}

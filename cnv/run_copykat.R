
run_copykat <- function(sces, 
    celltype_col = "celltype", 
    cells = NULL, 
    ref_cells = NULL,
    sample_col = "Sample",
    output_dir = NULL, 
    id.type = "S", 
    cell.line = "no", 
    ngene.chr = 5, 
    LOW.DR = 0.05, 
    UP.DR = 0.1, 
    win.size = 25, 
    norm.cell.names = NULL, 
    KS.cut = 0.1, 
    sam.name = NULL, 
    distance = "euclidean", 
    output.seg = FALSE, 
    plot.genes = TRUE, 
    genome = "hg20", 
    n.cores = 20, 
    ...){
    
    # loading packages
    library(magrittr)
    library(SingleCellExperiment)
    library(copykat)

    if(is.null(cells)){
        message("cells is NULL, please enter which cells to anaysis...")
        break
    }

    # subset sces to insure using about 50,000 cells to run copykat
    sces <- sces[, sces[[celltype_col]] %in% cells]

    if(ncol(sces) > 50000){
        factors <- 1/(ncol(sces)/50000)

        set.seed(101)
        sces <- sces[, sample(ncol(sces), ncol(sces)*factors)]
        ref_barcodes <- sces[, sces[[celltype_col]] %in% ref_cells] %>% colnames
    }else{
        ref_barcodes <- sces[, sces[[celltype_col]] %in% ref_cells] %>% colnames
    }

    samples <- sces[[sample_col]] %>% unique

    res <- copykat(
        rawmat = as.matrix(counts(sces)), 
        id.type = id.type, 
        cell.line = cell.line, 
        ngene.chr = ngene.chr, 
        LOW.DR = LOW.DR, 
        UP.DR = UP.DR, 
        win.size = win.size, 
        norm.cell.names = ref_barcodes, 
        KS.cut = KS.cut, 
        sam.name = sam.name, 
        distance = distance, 
        output.seg = output.seg, 
        plot.genes = plot.genes, 
        genome = genome, 
        n.cores = n.cores
    )

    if(is.null(output_dir)){
        output_dir <- getwd()
    }else(
        output_dir <- ifelse(grepl("/$", output_dir), output_dir, paste0(output_dir, "/"))
    )

    output_file <- paste0(output_dir, "/", "copykat_outputs.qs")
    qsave(res, file = output_file)
}
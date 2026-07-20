
run_copykat <- function(sces, 
    celltype_col = "celltype", 
    cells = NULL, 
    number_cells = 50000,
    assay.type = "counts",
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
    library(qs2)

    if(is.null(cells)){
        stop("cells is NULL; please supply a vector of cells or metadata group values to analyze.")
    }

    if(all(cells %in% colnames(sces))){
        sces <- sces[, cells]
    } else if(celltype_col %in% colnames(colData(sces))){
        sces <- sces[, colData(sces)[[celltype_col]] %in% cells]
    } else {
        stop("Unable to subset cells: `cells` are not cell names and `celltype_col` is not present in sces colData.")
    }

    if(ncol(sces) == 0){
        stop("No cells selected. Check `cells` values and `celltype_col`.")
    }

    if(ncol(sces) > number_cells){
        factors <- 1/(ncol(sces)/number_cells)
        set.seed(101)
        sces <- sces[, sample(ncol(sces), floor(ncol(sces)*factors))]
    }

    if(!is.null(ref_cells)){
        if(all(ref_cells %in% colnames(sces))){
            ref_barcodes <- ref_cells
        } else if(celltype_col %in% colnames(colData(sces))){
            ref_barcodes <- colnames(sces)[colData(sces)[[celltype_col]] %in% ref_cells]
        } else {
            stop("Unable to identify reference cells: `ref_cells` are not cell names and `celltype_col` is not present in sces colData.")
        }
    } else {
        ref_barcodes <- NULL
    }

    if(length(ref_barcodes) == 0){
        stop("No reference cells found in the selected dataset. Check `ref_cells` values and whether they match `celltype_col` metadata.")
    }

    message(paste0("there are ", ncol(sces), " cells and ", nrow(sces), " genes retained..."))
    print(table(colData(sces)[[celltype_col]]))
    message(paste0("selected ", length(ref_barcodes), " reference barcodes"))
    print(str(ref_barcodes))

    samples <- unique(colData(sces)[[sample_col]])

    if(is.null(output_dir)){
        output_dir <- getwd()
    } else {
        output_dir <- ifelse(grepl("/$", output_dir), output_dir, paste0(output_dir, "/"))
    }
    if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    res <- copykat(
        rawmat = as.matrix(assay(sces, assay.type)), 
        id.type = id.type, 
        cell.line = cell.line, 
        ngene.chr = ngene.chr, 
        LOW.DR = LOW.DR, 
        UP.DR = UP.DR, 
        win.size = win.size, 
        norm.cell.names = "", 
        KS.cut = KS.cut, 
        sam.name = paste0(output_dir, sam.name), 
        distance = distance, 
        output.seg = output.seg, 
        plot.genes = plot.genes, 
        genome = genome, 
        n.cores = n.cores, 
        ...
    )

    output_file <- paste0(output_dir, "copykat_outputs.qs2")
    qs_save(res, file = output_file)
}
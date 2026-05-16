mtx2sparse <- function(mtx,
    dir = NULL,
    sample = NULL,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    org = NULL,
    ...){

    # loading required packages
    suppressPackageStartupMessages(require(magrittr))
    suppressPackageStartupMessages(require(Matrix))
    suppressPackageStartupMessages(require(fs))
    suppressPackageStartupMessages(require(R.utils))
    # require(org) is wrong: substitute() captures the formal name "org" and
    # tries to load package "org". Resolve character package name or OrgDb.
    if (is.null(org)) {
        stop("'org' must be an OrgDb object or the package name (e.g. \"org.Xl.eg.db\").",
            call. = FALSE)
    }
    if (is.character(org)) {
        pkg <- org[[1L]]
        suppressPackageStartupMessages(require(pkg, character.only = TRUE))
        org_db <- get(pkg, envir = as.environment(paste0("package:", pkg)))
    } else {
        org_db <- org
    }
    suppressPackageStartupMessages(require(AnnotationDbi))

    message("add feature data...")
    features <- data.frame(ensembl = rownames(mtx),
            symbol = mapIds(org_db,
                keys = rownames(mtx),
                keytype = keytype,
                column = column),
            type = mapIds(org_db,
                keys = rownames(mtx),
                keytype = keytype,
                column = "GENETYPE"))
        if (tolower(keytype) == "ensembl") {
            features <- features %>% dplyr::select(ensembl, symbol, type)
        }else{
            features <- features %>% dplyr::select(symbol, ensembl, type)
        }

    # if(species == "human"){
    #     suppressPackageStartupMessages(require(org.Hs.eg.db))
    #     features <- data.frame(ensembl = rownames(mtx),
    #         symbol = mapIds(org.Hs.eg.db,
    #             keys = rownames(mtx),
    #             keytype = keytype,
    #             column = column),
    #         type = mapIds(org.Hs.eg.db,
    #             keys = rownames(mtx),
    #             keytype = keytype,
    #             column = "GENETYPE"))
    #     if(keytype == "ensembl"){
    #         features <- features %>% dplyr::select(ensembl, symbol, type)
    #     }else{
    #         features <- features %>% dplyr::select(symbol, ensembl, type)
    #     }

    # }else if(species == "mouse"){
    #     suppressPackageStartupMessages(require(org.Mm.eg.db))
    #     features <- data.frame(ensembl = rownames(mtx),
    #         symbol = mapIds(org.Mm.eg.db,
    #             keys = rownames(mtx),
    #             keytype = keytype,
    #             column = column),
    #         type = mapIds(org.Mm.eg.db,
    #             keys = rownames(mtx),
    #             keytype = keytype,
    #             column = "GENETYPE"))
    #     if(keytype == "ensembl"){
    #         features <- features %>% dplyr::select(ensembl, symbol, type)
    #     }else{
    #         features <- features %>% dplyr::select(symbol, ensembl, type)
    #     }
    # }

    message("get barcode...")
    barcodes <- colnames(mtx)

    # convert to sparse matrix
    if (!inherits(mtx, "CsparseMatrix")) {
        mtx <- as(mtx, "CsparseMatrix")
    }
    
    # create directory
    dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
    dir1 <- paste0(dir, sample, "/")
    dir_create(path = dir1)

    # save file
    message("writing to the disk and compression...")
    write.table(barcodes, 
        file = paste0(dir1, "barcodes.tsv"),
        row.names = F,
        col.names = F,
        sep = "\t",
        quote = F)
    compressFile(filename = paste0(dir1, "barcodes.tsv"), 
        destname = paste0(dir1, "barcodes.tsv.gz"), 
        ext = "gz", 
        FUN = gzfile)
    
    write.table(features,
        file = paste0(dir1, "features.tsv"),
        row.names = F,
        col.names = F,
        quote = F,
        sep = "\t")
    compressFile(filename = paste0(dir1, "features.tsv"), 
        destname = paste0(dir1, "features.tsv.gz"), 
        ext = "gz", 
        FUN = gzfile)

    writeMM(mtx, file = paste0(dir1, "matrix.mtx"))
    compressFile(filename = paste0(dir1, "matrix.mtx"), 
        destname = paste0(dir1, "matrix.mtx.gz"), 
        ext = "gz", 
        FUN = gzfile)
    message("finished...")
}
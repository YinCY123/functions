mtx2sparse <- function(mtx,
    dir = NULL,
    sample = NULL,
    species = "human",
    keytype = "ENSEMBL",
    column = "SYMBOL",
    ...){

    # loading required packages
    suppressPackageStartupMessages(require(magrittr))
    suppressPackageStartupMessages(require(Matrix))
    suppressPackageStartupMessages(require(fs))
    suppressPackageStartupMessages(require(R.utils))

    message("add feature data...")
    if(species == "human"){
        suppressPackageStartupMessages(require(org.Hs.eg.db))
        features <- data.frame(ensembl = rownames(mtx),
            symbol = mapIds(org.Hs.eg.db,
                keys = rownames(mtx),
                keytype = keytype,
                column = column),
            type = mapIds(org.Hs.eg.db,
                keys = rownames(mtx),
                keytype = keytype,
                column = "GENETYPE"))
        if(keytype == "ensembl"){
            features <- features %>% dplyr::select(ensembl, symbol, type)
        }else{
            features <- features %>% dplyr::select(symbol, ensembl, type)
        }

    }else if(species == "mouse"){
        suppressPackageStartupMessages(require(org.Mm.eg.db))
        features <- data.frame(ensembl = rownames(mtx),
            symbol = mapIds(org.Mm.eg.db,
                keys = rownames(mtx),
                keytype = keytype,
                column = column),
            type = mapIds(org.Mm.eg.db,
                keys = rownames(mtx),
                keytype = keytype,
                column = "GENETYPE"))
        if(keytype == "ensembl"){
            features <- features %>% dplyr::select(ensembl, symbol, type)
        }else{
            features <- features %>% dplyr::select(symbol, ensembl, type)
        }
    }

    message("get barcode...")
    barcodes <- colnames(mtx)

    # convert to sparse matrix
    # mtx <- as.matrix(mtx)
    message("convert to sparse...")
    mtx <- as(mtx, "CsparseMatrix")
    
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
splitMtx <- function(mtx, barcode_list, ...) {
    mtx_list <- vector(mode = "list")
    for (name in names(barcode_list)) {
        mtx_list[[name]] <- mtx[, barcode_list[[name]]]
    }
    return(mtx_list)
}
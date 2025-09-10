exportTable <- function(x, 
    file = NULL, 
    row.names = FALSE,
    col.names = FALSE,
    ...){
    write.table(x = x, file = file, 
                row.names = row.names, 
                col.names = col.names, 
                sep = "\t", 
                quote = F)
}
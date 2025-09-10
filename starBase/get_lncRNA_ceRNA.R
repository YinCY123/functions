get_lncRNA_ceRNA <- function(ceRNA_name, 
                                  dir = NULL,
                                  assembly = "hg38",
                                  miRNA_number = 2, 
                                  p_value = 0.05, 
                                  fdr = 0.05,
                                  pancancer = 0,
                                  save = FALSE,
                                  ...){
    base_url <- "https://rnasysu.com/encori/moduleDownload.php?source=ceRNA&type=txt&value=%s;lncRNA;%s;%d;%0.2f;%0.2f;%d"
    url <- sprintf(base_url, assembly, ceRNA_name, miRNA_number, p_value, fdr, pancancer)
    
    if(is.null(dir)){
        dir = getwd()
        message("The file will be saved at current working directory.\n\n")
    }
    
    if(save){
        dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
        file <- paste0(dir, ceRNA_name, ".txt")
        download.file(url = url, destfile = file)
    }else{
        read.table(url, header = T, skip = 1)
    }
}
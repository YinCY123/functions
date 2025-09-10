get_mRNA_ceRNA <- function(ceRNA_name, 
                                dir = NULL, 
                                miRNA_number = 2,
                                p_value = 0.05, 
                                fdr = 0.05, 
                                pancancer = 0,
                                assembly = "hg38", # mm10 for mouse
                                save = FALSE,
                                ...){
    
    if(is.null(dir)){
        dir <- getwd()
        message("The file will be saved at curret working dictionary.\n\n")
    }
    
    base_url <- "https://rnasysu.com/encori/moduleDownload.php?source=ceRNA&type=txt&value=%s;mRNA;%s;%d;%0.2f;%0.2f;%d"
    url <- sprintf(fmt = base_url, assembly, ceRNA_name, miRNA_number, p_value, fdr, pancancer)
    

    if(save){
        dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
        file = paste0(dir, ceRNA_name, ".txt")
        download.file(url = url, destfile = file)
    }else{
        read.table(url, header = TRUE, skip = 1)
    }

}
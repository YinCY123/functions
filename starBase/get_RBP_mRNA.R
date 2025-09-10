get_RBP_mRNA <- function(target_gene = NULL, 
                              RBP = "all", 
                              assembly = "hg38",
                              pancancer = 0,
                              clip_data = 1,
                              dir = NULL,
                              ...){
    base_url <- "https://rnasysu.com/encori/moduleDownload.php?source=rbpClipRNA&type=txt&value=%s;mRNA;%s;%d;%d;%s"
    url <- sprintf(base_url, assembly, RBP, clip_data, pancancer, target_gene)
    
    if(is.null(dir)){
        dir <- getwd()
        message("The file will be saved at current working directory.\n")
    }
    dir = ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
    file = paste0(dir, target_gene, ".txt")
    download.file(url = url, destfile = file)
}

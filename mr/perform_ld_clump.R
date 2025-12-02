perform_ld_clump <- function(x, 
    clump_r2 = 0.01, 
    clump_kb = 10000, 
    bfile = "/home/yincy/BioHome/datasets/mr/EUR", 
    plink_bin = "/home/yincy/BioHome/tools/plink/plink-1.07-x86_64/plink", 
    ...){
        lx <- lapply(unique(x$snp), function(y){
            tmp <- x %>% dplyr::filter(snp == y) %>% dplyr::mutate(rsid = snp) 
            tmp <- ld_clump(dat = tmp, 
                    clump_r2 = clump_r2, 
                    clump_kb = clump_kb, 
                    bfile = bfile, 
                    plink_bin = plink_bin) 
        })
        lx <- do.call(rbind, lx) %>% 
                as.data.frame() %>% 
                dplyr::select(-rsid, -pval, -id)
        return(lx)
}

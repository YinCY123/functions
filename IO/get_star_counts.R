# combine multiple STAR out sample's count into a data frame

get_star_counts <- function(samples,
    count_id = 2,
    ...){
        # loading package
        require(magrittr)
        require(stringr)

        mtx_list <- vector(mode = "list")
        for(sample in samples){

            # parse sample name
            ids <- gregexec(pattern = "/", text = sample)[[1]]
            start = ids[length(ids)-1] + 1
            stop = ids[length(ids)] - 1
            sample_name <- substr(sample, start = start, stop = stop)

            tmp <- read.table(sample, sep = "\t", head = F, skip = 4) %>% 
                dplyr::select(1, all_of(count_id)) %>% 
                magrittr::set_colnames(value = c("ensembl", sample_name)) %>% 
                dplyr::mutate(ensembl = str_remove(ensembl, "\\.[:digit:]{1,}$"))
            
            if(any(duplicated(tmp$ensembl))){
                tmp <- aggregate(. ~ ensembl, FUN = mean, data = tmp) %>% 
                    tibble::column_to_rownames("ensembl") %>% 
                    round(digit = 0)
            }else{
                tmp <- tmp %>% 
                    tibble::column_to_rownames("ensembl")
            }
            mtx_list[[sample_name]] <- tmp
        }
        
        rids <- Reduce(intersect, lapply(mtx_list, rownames))
        mtx_list <- lapply(mtx_list, function(x){x[rids, , drop = F]})
        df <- do.call(cbind, mtx_list) %>% as.data.frame()
        return(df)
    }
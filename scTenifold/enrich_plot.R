enrich_plot <- function(x, 
    variable = "Z",
    fun_path = "../../github/functions/enrichment/enrichbar.R", 
    gene_sets = NULL,
    minGSSize = 10, 
    maxGSSize = 1000,
    ...){
        # loading required packages
        library(magrittr)
        library(clusterProfiler)

        # gene rank
        grk <- setNames(x[[variable]], x$Gene) %>% sort(decreasing = TRUE)

        gsea_res <- GSEA(
            geneList = grk, 
            minGSSize = minGSSize, 
            maxGSSize = maxGSSize, 
            eps = 0, 
            pvalueCutoff = 1, 
            TERM2GENE = gene_sets, 
            seed = TRUE
        )

        # loading function
        source(fun_path)
        p <- enrichbar(res = gsea_res, ...)
        return(list(p, gsea_res))
}

TCGA_analysis <- function(x, gene, dir, gsea = TRUE, kegg = FALSE, go = FALSE, logFC = 0.5, pvalue = 0.05, ...){
    # loading required packages
    suppressPackageStartupMessages(require(SummarizedExperiment))
    suppressPackageStartupMessages(require(magrittr))
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(tidyr))
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(org.Hs.eg.db))
    suppressPackageStartupMessages(require(tibble))
    suppressPackageStartupMessages(require(DESeq2))
    suppressPackageStartupMessages(require(clusterProfiler))
    suppressPackageStartupMessages(require(msigdbr))
    suppressPackageStartupMessages(require(openxlsx))
    suppressPackageStartupMessages(require(ggplot2))
    suppressPackageStartupMessages(require(GseaVis))
    suppressPackageStartupMessages(require(stringr))
    source("/home/yincy/git/bior/functions/enrichbar.R")
    
    # if(class(x) %in% c("SummarizedExperiment", "RangedSummarizedExperiment")){
    #     mtx = assay(x, "unstranded")
    # }else if(class(as.matrix(x)) == "matrix"){
    #     mtx = as.matrix(x)
    # }
    mtx <- assay(x, "unstranded")    
    # preprocessing
    mtx <- mtx %>% as.data.frame %>% 
        tibble::rownames_to_column("ensembl") %>% 
        mutate(ensembl = str_remove(ensembl, pattern = "\\.[0-9]{1,}$"), 
               symbol = mapIds(org.Hs.eg.db, keys = ensembl, keytype = "ENSEMBL", column = "SYMBOL")) %>% 
        relocate(symbol) %>% 
        select(-ensembl) %>% 
        filter(!is.na(symbol))
    
    # remove duplicated genes
    mtx <- aggregate(. ~ symbol, data = mtx, FUN = "mean") %>% 
        column_to_rownames("symbol") %>% 
        round(0)
    print("generated the expression matrix...")
    
    # create sample information
    df <- data.frame(sample_type_id = substr(colnames(mtx), 14, 15), 
                     row.names = colnames(mtx)) %>% 
        dplyr::filter(sample_type_id == "01") 
    
    # select tumor samples
    mtx <- mtx[, rownames(df)]
    print(mtx[1:3, 1:3])
    print("sample selected...")
    
    # create gene expression matrix
    group_idx <- mtx[gene, , drop = T] %>% unlist()
    group_idx <- ifelse(group_idx > median(group_idx), "high", "low")
    print(head(group_idx))
    
    # add group info
    df <- df %>% mutate(group = group_idx[rownames(df)])
    
    print(table(df$group))
    # differential expression analysis
    print("performing differential expression analysis...")
    dds <- DESeqDataSetFromMatrix(countData = mtx,
                                  colData = df, 
                                  design = ~group)
    dds <- DESeq(object = dds)
    print("after differential expression test...")
    tables <- results(dds, contrast = c("group", "high", "low")) %>% 
        as.data.frame %>% rownames_to_column("symbol")
    print(head(tables))
    
    # GSEA
    if(gsea){
        print("GSEA analysis...")
        gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")
        gene_sets <- gene_sets %>% dplyr::filter(gene_symbol %in% tables$symbol)
        
        grk <- setNames(-log10(tables$pvalue) * sign(tables$log2FoldChange), tables$symbol) %>% 
            sort(decreasing = TRUE)
        gsea_res <- GSEA(geneList = grk, 
                         minGSSize = 10, 
                         maxGSSize = 500, 
                         eps = 0, 
                         pvalueCutoff = 1, 
                         pAdjustMethod = "BH", 
                         TERM2GENE = gene_sets[, c("gs_name", "gene_symbol")])
        
        # save GSEA result
        write.xlsx(gsea_res@result, file = paste0(dir, "/gsea_res.xlsx"))
        saveRDS(gsea_res, file = paste0(dir, "/", "gsea_res.rds"))
        
        # visulize GSEA result
        print("Visualize top 3 terms of either direction...")
        enrichbar(res = gsea_res@result, 
                  x = "pvalue", 
                  y = "Description", 
                  top = 20, 
                  text_size = 2,
                  file = paste0(dir, "/GSEA_result.pdf"), 
                  width = 10, 
                  height = 10, 
                  step = 5, 
                  scale = 1)
        
        top3_term <- c(gsea_res@result %>%
                           dplyr::filter(NES > 0) %>%
                           dplyr::arrange(pvalue) %>%
                           dplyr::slice_head(n = 3) %>% dplyr::pull(Description),
                       gsea_res@result %>%
                           dplyr::filter(NES < 0) %>%
                           dplyr::arrange(pvalue) %>%
                           dplyr::slice_head(n = 3) %>% dplyr::pull(Description))
        
        for(id in top3_term){
            gseaNb(object = gsea_res,
                   geneSetID = id,
                   addPval = TRUE,
                   pvalY = 0.85)
            ggsave(file = paste0(dir, "/", id, ".pdf"), width = 7, height = 7)
        }
    }
    
    if(any(kegg, go)){
        degs <- tables %>% 
            dplyr::filter(abs(log2FoldChange) > logFC, pvalue < pvalue) %>% 
            dplyr::pull(symbol)
        print(str(degs))
        if(length(degs) <= 0){
            print("there are no DEGs between high and low group...")
            break
        }
        degs <- mapIds(org.Hs.eg.db, keys = degs, keytype = "SYMBOL", column = "ENTREZID") %>% na.omit
        print("Converted to entrez id...")
        print(str(degs))
    }
    
    # KEGG
    if(kegg){
        print("KEGG enrichment analysis...")
        kegg_enrich <- enrichKEGG(gene = degs, 
                                  organism = "hsa", 
                                  pvalueCutoff = 1, 
                                  qvalueCutoff = 1)
        
        print("writing result to file...")
        write.xlsx(list(kegg_enrich@result), file = paste0(dir, "/kegg_enrich.xlsx"))
        top_kegg <- kegg_enrich@result %>% 
            dplyr::arrange(pvalue) %>% 
            dplyr::slice_head(n = 10) %>% 
            dplyr::mutate(logp = -log10(pvalue))
        
        print("visualize top KEGG...")
        top_kegg %>% 
            ggplot(aes(logp, reorder(Description, logp))) +
            geom_bar(stat = "identity", aes(fill = logp)) +
            scale_fill_continuous(name = "-log10(P Value)") +
            scale_x_continuous(name = "-log10(P Value)") +
            scale_y_discrete(name = NULL) +
            theme_classic()
        
        ggsave(filename = paste0(dir, "/kegg_enrich.pdf"), 
               width = 7, height = 7)
        
    }
    
    # GO
    if(go){
        print("GO enrichment analysis...")
        go_enrich <- enrichGO(gene = degs, 
                              OrgDb = "org.Hs.eg.db", 
                              ont = "ALL", 
                              pvalueCutoff = 1, 
                              qvalueCutoff = 1)
        write.xlsx(list(go_enrich@result), file = paste0(dir, "/go_enrich.xlsx"))
        top_go <- go_enrich@result %>% 
            dplyr::group_by(ONTOLOGY) %>% 
            dplyr::arrange(pvalue) %>% 
            dplyr::slice_head(n = 5) %>% 
            dplyr::mutate(logp = -log10(pvalue))
        
        top_go %>% 
            ggplot(aes(logp, reorder(Description, logp))) +
            geom_bar(aes(fill = ONTOLOGY), stat = "identity") +
            facet_wrap(vars(ONTOLOGY), scale = "free", ncol = 1) +
            scale_x_continuous(name = "-log10(P Value)") +
            scale_y_discrete(name = NULL) +
            theme_classic()  
        ggsave(file = paste0(dir, "/go_enrich.pdf"), 
               width = 7, height = 7)
    }
    print("All the analysis is done...")
}

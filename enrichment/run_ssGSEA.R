run_ssGSEA <- function(x, 
    gene = NULL,
    type = c("cor", "hl"),
    cor_type = c("pearson", "spearman", "kendall"),
    hl_method = c("median", "mean"),
    geneset = c("KEGG", "GO"), 
    kegg_df = NULL,
    species = "human",
    seq_tech = c("RNAseq", "microarray"),
    gsea_x = "pvalue", 
    gsea_y = "Description",
    pval_x = 0.9,
    pval_y = 0.9, 
    top_n = 15,
    left = -5, 
    right = 5,
    width = 7, 
    height = 7,
    scale = 1,
    dir = NULL,
    ...){
        # loading required packages
        require(magrittr)
        require(clusterProfiler)
        require(GseaVis)

        # parameter parse
        type <- match.arg(type)
        cor_type <- match.arg(cor_type)
        geneset <- match.arg(geneset)
        # species <- match.arg(species)
        seq_tech <- match.arg(seq_tech)
        hl_method <- match.arg(hl_method)

        # convert x to matrix
        message("convert to matrix...")
        x <- as.matrix(x)

        # check directory existence
        if(is.null(dir)){
            dir <- getwd()
        }

        # calculate gene rank
        if(type == "cor"){
            message("perform ssGSEA with correlation method...")
            features <- x %>% rownames
            features <- setdiff(features, gene)

            cor_list <- vector(mode = "list")
            # calculate correlation
            for(feature in features){
                tmp <- cor.test(x[gene, , drop = T], x[feature, , drop = T], 
                    method = cor_type)
                df <- data.frame(
                    cor = tmp$estimate %>% unname(), 
                    pval = tmp$p.value,
                    zscore = tmp$statistic
                )
                cor_list[[feature]] <- df
            }

            # convert list to dataframe
            cor_df <- do.call(rbind, cor_list) %>% 
                as.data.frame() %>% 
                tibble::rownames_to_column("symbol")

            # save correlation data
            dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
            file <- paste0(dir, "correlation_df.rds")
            saveRDS(cor_df, file)

            # adjust multiple test
            cor_df$p.adjust <- p.adjust(cor_df$pval, method = "BH")
            grk <- setNames(-log10(cor_df$pval) * sign(cor_df$cor), cor_df$symbol) %>% 
                sort(decreasing = TRUE)
        }else{
            message("perform ssGSEA with high low method...")
            # divde samples into high low groups
            df <- data.frame(samples = x[gene, ] %>% names, 
                value = x[gene, ])
            fun <- ifelse(hl_method == "median", median, mean)
            df <- df %>% 
                dplyr::mutate(group = ifelse(value > fun(value), "high", "low"), 
                    group = factor(group, levels = c("high", "low")))

            if(seq_tech == "RNAseq"){
                message("perform DEG analysis with DESeq2...")
                # load required package
                require(DESeq2)
                dds <- DESeqDataSetFromMatrix(
                    countData = x, 
                    colData = df,
                    design = ~group
                )
                dds <- DESeq(dds)
                tables <- results(dds, contrast = c("group", "high", "low")) %>% 
                    as.data.frame() %>% 
                    tibble::rownames_to_column("symbol")
                # save DEG table
                dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
                file <- paste0(dir, "DEG_tables_hl_DESeq2.rds")
                # saveRDS(tables, file)
                grk <- setNames(-log10(tables$pvalue) * sign(tables$log2FoldChange), tables$symbol) %>% 
                    sort(decreasing = TRUE)
            }else{
                message("perform DEG analysis with limma...")
                # loading required packages
                require(limma)
                design <- model.matrix( ~ 0 + df$group)
                colnames(design) <- c("high", "low")
                cont <- makeContrasts(high-low, levels = design)
                fit <- lmFit(x, design = design) %>% 
                    contrasts.fit(contrasts = cont) %>% 
                    eBayes()
                tables <- topTable(fit, number = Inf) %>% 
                    as.data.frame() %>% 
                    tibble::rownames_to_column("symbol")

                # save DEG table
                dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
                file <- paste0(dir, "DEG_tables_hl_limma.rds")
                # saveRDS(tables, file)
                grk <- setNames(-log10(tables$P.Value) * sign(tables$logFC), tables$symbol) %>% 
                    sort(decreasing = TRUE)
            }
        }
        # choose gene set
        if(geneset == "KEGG" & is.null(kegg_df)){
            source("/home/yincy/git/bior/functions/enrichment/get_kegg_geneset.R")
            s = switch(species, 
                "human" = "hsa", 
                "mouse" = "mmu", 
                species)
            gs <- get_kegg_geneset(species = s) %>% 
                dplyr::select(name, symbol)
        }else{
            gs = kegg_df
        }
        # TODO
        # get gene set from msigdbr


        # perform GSEA analysis
        message("head 3: ", names(head(grk, 3)), "\n\t: ", head(grk, 3))
        message("tail 3: ", names(tail(grk, 3)), "\n\t: ", tail(grk, 3))
        print(head(gs, 3))
        gsea_result <- GSEA(
                geneList = grk, 
                minGSSize = 10, 
                maxGSSize = 1000, 
                eps = 0, 
                pvalueCutoff = 1, 
                TERM2GENE = gs, 
                seed = TRUE)
            # save GSEA result
        dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
        file <- paste0(dir, gene, "_GSEA_results.qs")
        qsave(gsea_result, file)

        # visulization
        source("/home/yincy/git/functions/enrichment/enrichbar.R")
        args(enrichbar)
        # save GSEA bar plot
        dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
        file <- paste0(dir, gene, "_GSEA_bar.pdf")
        enrichbar(gsea_result, 
                x = gsea_x, 
                y = gsea_y, 
                top = top_n, 
                left = left, 
                right = right, 
                file = file,
                width = width, 
                height = height, 
                scale = scale)

        # get top 3 terms
        top_terms <- c(
            gsea_result@result %>% 
                dplyr::filter(NES > 0) %>% 
                    dplyr::arrange(pvalue) %>% 
                    dplyr::slice_head(n = 3) %>% 
                    dplyr::pull(Description),
            gsea_result@result %>% 
                dplyr::filter(NES < 0) %>% 
                dplyr::arrange(pvalue) %>% 
                dplyr::slice_head(n = 3) %>% 
                dplyr::pull(Description)
        )

        # save results
        dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
        for(term in top_terms){
            gseaNb(gsea_result, 
                geneSetID = term,
                addPval = TRUE, 
                pvalY = pval_y, 
                pvalX = pval_x)
            file <- paste0(dir, gene, "_", term, ".pdf")
            ggsave(file, width = 8, height = 6)
        }
    }
sces_pseudotime <- function(sces, 
    cell_col = "celltype", 
    cells = NULL,
    sample_cell = FALSE,
    prop = 1,
    seed = 101,
    feature_threshold = 0.01,
    gene_id = "symbol",
    fullModelFormulaStr = "~group + celltype",
    reducedModelFormulaStr = "~1", 
    order_genes = NULL, 
    use_dispersion = FALSE,
    top_n = 1000,
    ncores = 10,
    reduction_method = "DDRTree",
    p_val = 0.05,
    ...){
        # loading required packages
        suppressPackageStartupMessages(require(magrittr))
        suppressPackageStartupMessages(require(monocle))
        suppressPackageStartupMessages(require(SingleCellExperiment))
        suppressPackageStartupMessages(require(stringr))
        suppressPackageStartupMessages(require(qs))

        # filter cells
        message("filtering cells...")
        sces <- sces[, colData(sces)[, cell_col] %in% cells]
        if(sample_cell){
            set.seed(seed)
            sces <- sces[, sample(ncol(sces), ncol(sces) * prop)]
        }
        message(paste0("Do pseudotime analysis on: ", str_c(unique(colData(sces)[, cell_col]), collapse = ", ")))
        
        # build CellDataSet object
        message("building CellDataSet object...")
        pd <- new("AnnotatedDataFrame", data = sces %>% colData %>% as.data.frame())
        fd <- new("AnnotatedDataFrame", data = sces %>% rowData %>% as.data.frame())
        cds <- newCellDataSet(cellData = logcounts(sces) %>% as(Class = "sparseMatrix"),
            phenoData = pd,
            featureData = fd,
            lowerDetectionLimit = 0.1,
            expressionFamily = VGAM::negbinomial.size())

        
        # estimate size factors and dispersions
        message("estimating size factors and dispersions...")
        cds <- estimateSizeFactors(cds)
        cds <- estimateDispersions(cds)
        
        # filtering low-expressed genes
        message("filtering out low-expressed genes...")
        cds <- detectGenes(cds, min_expr = 0.1)
        fData(cds)$percent_cell_expressed <- fData(cds)$num_cells_expressed / ncol(cds)
        ids <- fData(cds)$percent_cell_expressed > feature_threshold
        cds <- cds[ids, ]

        if(!is.null(order_genes)){
            message("using user provided order genes...")
            ogs <- order_genes
            qsave(ogs, "order_genes.qs")
        }else if(use_dispersion){
            message("using dispersion for order genes...")
            disp_table <- dispersionTable(cds)
            qsave(disp_table, "dispersion_table.qs")
            ogs <- disp_table %>% 
                dplyr::filter(dispersion_empirical > dispersion_fit & mean_expression > 0.1) %>% 
                dplyr::pull(gene_id)
            message(paste0("There total: ", length(ogs), " dispersion genes..."))
            if(length(ogs) > top_n){
                ogs <- ogs[1:top_n]
            }
            qsave(ogs, "order_genes.qs")
        }else{
            message("using ordering genes from differentialGeneTest...")
            diff_test_results <- differentialGeneTest(cds = cds, 
                fullModelFormulaStr = fullModelFormulaStr, 
                reducedModelFormulaStr = reducedModelFormulaStr, 
                cores = ncores)
            qsave(diff_test_results, "differentialGeneTest_table.qs")
            ogs <- diff_test_results %>% 
                dplyr::filter(pval < p_val) %>% 
                dplyr::arrange(pval) %>% 
                dplyr::pull(gene_id) %>% 
                unique()
            message(paste0("There are total: ", length(ogs), " ordering genes."))
            # print(str_c(rep("-", 20), collapse = ""))
            # message(paste0(str_c(head(ogs, 5), collapse = ", "), "..."))
            if(length(ogs) > top_n){
                ogs <- ogs[1:top_n]
            }
            qsave(ogs, "order_genes.qs")
        }
        cds <- setOrderingFilter(cds, ogs)

        # reduce dimentionality
        message("reducing dimentionality...")
        cds <- reduceDimension(cds = cds,
            reduction_method = reduction_method,
            max_components = 2)
        
        # order cell
        message("ordering cells...")
        cds <- orderCells(cds = cds)

        message("analysis is done...")
        return(cds)
    }
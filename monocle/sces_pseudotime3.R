sces_pseudotime3 <- function(sces,
    celltype_col = "celltype",
    cells = NULL,
    cell_meta = NULL,
    feat_meta = NULL,
    batch = "Sample",
    sample_cells = FALSE,
    prop = 1,
    order_genes = NULL,
    npcas = 30,
    reduction_method = "UMAP",
    ...){
        # Check required packages
        if(!requireNamespace("monocle3", quietly = TRUE)){
            stop("Package 'monocle3' is required but not installed.")
        }
        if(!requireNamespace("SingleCellExperiment", quietly = TRUE)){
            stop("Package 'SingleCellExperiment' is required but not installed.")
        }

        # Input validation
        if(sample_cells && (prop <= 0 || prop > 1)){
            stop("prop must be between 0 and 1 when sample_cells is TRUE")
        }

        if(npcas <= 0 || !is.numeric(npcas)){
            stop("npcas must be a positive integer")
        }

        # Create cell dataset from different input types
        if(is.matrix(sces) || is(sces, "dgCMatrix")){
            if(is.null(cell_meta) || is.null(feat_meta)){
                stop("cell_meta and feat_meta cannot be NULL when sces is a matrix.")
            }

            # Create cell dataset
            cds <- monocle3::new_cell_data_set(
                expression_data = sces,
                cell_metadata = cell_meta,
                gene_metadata = feat_meta
            )
        } else if(inherits(sces, "SingleCellExperiment")){
            # Filter by cell type if specified
            if(!is.null(cells)){
                if(!celltype_col %in% names(SummarizedExperiment::colData(sces))){
                    stop(paste0("Column '", celltype_col, "' not found in colData"))
                }
                sces_sub <- sces[, SummarizedExperiment::colData(sces)[[celltype_col]] %in% cells]
            } else {
                sces_sub <- sces
            }

            # Sample cells if requested
            if(sample_cells){
                n_cells <- ncol(sces_sub)
                sample_size <- floor(n_cells * prop)
                if(sample_size < 1){
                    stop("sample_size results in less than 1 cell. Increase prop or set sample_cells = FALSE")
                }
                set.seed(123)  # For reproducibility
                sces_sub <- sces_sub[, sample(n_cells, sample_size)]
            }

            # Create cell dataset
            cds <- monocle3::new_cell_data_set(
                expression_data = SingleCellExperiment::counts(sces_sub),
                cell_metadata = as.data.frame(SummarizedExperiment::colData(sces_sub)),
                gene_metadata = as.data.frame(SummarizedExperiment::rowData(sces_sub))
            )
        } else {
            stop("sces must be a matrix, sparse matrix, or SingleCellExperiment object")
        }

        # Normalize and pre-process
        message("Normalize and pre-process...")
        if(!is.null(order_genes)){
            cds <- monocle3::preprocess_cds(cds,
                num_dim = npcas,
                method = "PCA",
                use_genes = order_genes)
        } else {
            cds <- monocle3::preprocess_cds(cds,
                num_dim = npcas,
                method = "PCA")
        }

        # Remove batch effects if specified
        if(!is.null(batch)){
            message("Remove batch effects...")
            if(!batch %in% names(SummarizedExperiment::colData(cds))){
                warning(paste0("Batch column '", batch, "' not found in colData. Skipping batch correction."))
            } else {
                cds <- monocle3::align_cds(cds, alignment_group = batch)
            }
        }

        # Reduce dimensionality
        message("Reduce dimensionality...")
        cds <- monocle3::reduce_dimension(cds, reduction_method = reduction_method)

        # Cluster cells
        message("Clustering cells...")
        cds <- monocle3::cluster_cells(cds, reduction_method = reduction_method)

        # Learn trajectory graph
        message("Learn graph...")
        cds <- monocle3::learn_graph(cds)

        # Order cells along pseudotime
        message("Order cells...")
        cds <- monocle3::order_cells(cds, reduction_method = reduction_method)

        return(cds)
    }

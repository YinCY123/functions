sces_process <- function(sces, 
                         sample = "Sample", 
                         qc_dir = NULL, 
                         qc_prefix = NULL, 
                         nrow = 3,
                         test_empty = FALSE, 
                         lower = 500, 
                         n_batch_dim = 50, 
                         mito_qc = TRUE, 
                         cell_qc = TRUE, 
                         log = TRUE,
                         feature_qc = TRUE, 
                         nmads = 3, 
                         detected_threshold = 0.01, 
                         ncores = 10, 
                         top_n = Inf, 
                         nn = seq(50, 100, 50), 
                         seed = 101, 
                         ...){
  
  # Input validation
  if (!is(sces, "SingleCellExperiment")) {
    stop("Input must be a SingleCellExperiment object")
  }
  
  if (!sample %in% colnames(colData(sces))) {
    stop("Sample column not found in colData")
  }
  
  # Load packages using requireNamespace instead of require
  suppressPackageStartupMessages(library("magrittr"))
  suppressPackageStartupMessages(library("scran"))
  suppressPackageStartupMessages(library("SingleCellExperiment"))
  suppressPackageStartupMessages(library("scuttle"))
  suppressPackageStartupMessages(library("scater"))
  suppressPackageStartupMessages(library("bluster"))
  suppressPackageStartupMessages(library("qs"))
  suppressPackageStartupMessages(library("batchelor"))
  suppressPackageStartupMessages(library("BiocParallel"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("rlang"))
  if(test_empty){
    suppressPackageStartupMessages(library("DropletUtils"))
  }
  
  # create paralell object
  # bp_param <- BiocParallel::MulticoreParam(workers = ncores, progressbar = TRUE)
  bp_param <- BiocParallel::SnowParam(workers = ncores, progressbar = TRUE)

  # safe vars used for QC plotting (defined regardless of cell_qc)
  vars <- if (mito_qc) c("sum", "detected", "subsets_mito_percent") else c("sum", "detected")

  # remove empty drops
  if(test_empty){
    message("remove empty drops...")
    set.seed(seed)
    emp <- emptyDrops(m = counts(sces), lower = lower, BPPARAM = bp_param)
    qc_dir <- ifelse(is.null(qc_dir), getwd(), qc_dir)
    qc_dir <- ifelse(grepl("/$", qc_dir), qc_dir, paste0(qc_dir, "/"))
    qsave(emp, paste0(qc_dir, "empty_drops.qs"))
    
    is_cell <- which(emp$FDR < 0.001)
    message("There ", length(is_cell), " cells retained...")
    sces <- sces[, is_cell]
    rm(emp); gc()
  }

  # QC
    if (cell_qc) {
      message("Performing Cell QC...")
      is_mito <- grepl("^mt-", rownames(sces), ignore.case = TRUE)
      sces <- addPerCellQCMetrics(
        sces,
        subsets = list(mito = is_mito),
        BPPARAM = bp_param
      )

      # before qc (vars defined earlier)
      df <- makePerCellDF(sces, use.coldata = TRUE, use.dimred = F) %>% 
        tidyr::pivot_longer(cols = dplyr::any_of(vars), 
                            names_to = "vv", 
                            values_to = "value")

        p <- df %>% ggplot2::ggplot(ggplot2::aes_string(x = sample, y = "value")) + 
          geom_violin(mapping = ggplot2::aes_string(fill = sample), scale = "width", width = 0.8) + 
          geom_jitter(size = 0.5, width = 0.4) +
          facet_wrap(vars(vv), nrow = nrow, scale = "free") + 
          scale_x_discrete(name = NULL) + 
          scale_y_continuous(name = NULL) + 
          theme(legend.position = "none", 
            panel.background = element_blank(), 
            panel.border = element_rect(fill = NA), 
            strip.background = element_blank(), 
            strip.text = element_text(size = 14, face = "bold"), 
            panel.grid.major = element_line(linetype = 2, color = "grey", linewidth = 0.2), 
            axis.text = element_text(size = 12, face = "bold"))
      qc_dir <- ifelse(is.null(qc_dir), getwd(), qc_dir)
      qc_dir <- ifelse(grepl("/$", qc_dir), qc_dir, paste0(qc_dir, "/"))

      file <- ifelse(is.null(qc_prefix), paste0(qc_dir, "before_qc.pdf"), paste0(qc_dir, qc_prefix, "_before_qc", ".pdf"))
      ggsave(plot = p, file = file, 
             width = length(unique(colData(sces)[[sample]])) * 3, 
             height = 8, scale = 0.8, limitsize = FALSE)

      # Vectorized QC operations
      qc <- tryCatch({
        qc_metrics <- list(
          sum = sces$sum,
          detected = sces$detected
        )
        if (mito_qc) {
          qc_metrics$mito <- sces$subsets_mito_percent
        }

        type <- if (mito_qc) c("lower", "lower", "higher") else c("lower", "lower")
        logs <- if (log) {
          if (mito_qc) c(TRUE, TRUE, TRUE) else c(TRUE, TRUE)
        } else {
          if (mito_qc) c(FALSE, FALSE, FALSE) else c(FALSE, FALSE)
        }
        qc_results <- mapply(FUN = isOutlier, 
          metric = qc_metrics, 
          type = type, 
          log = logs, 
          MoreArgs = list(nmads = nmads, batch = sces[[sample]]))
        
        rowSums(qc_results) > 0

      }, error = function(e) {
        warning("QC calculation failed: ", e$message)
        rep(FALSE, ncol(sces))
      })
      
      message("Number of cells removed after cell QC: ", sum(qc), ", ", paste0(round(sum(qc)/ncol(sces), 4) * 100, "%"))
      sces <- sces[, !qc]

    rm(qc); gc()
    }
    
    if(feature_qc){
      message("Performing Feature QC...")
      sces <- addPerFeatureQC(sces)
      ids <- rowData(sces)$detected > detected_threshold
      message("Number of features removed after feature QC: ", sum(!ids), ", ", paste0(round(sum(!ids)/nrow(sces), 4) * 100, "%"))
      sces <- sces[ids, ]

      rm(ids); gc()
    }

    # after qc
    is_mito <- grepl("^mt-", rownames(sces), ignore.case = TRUE)
    sces <- addPerCellQCMetrics(
      sces,
      subsets = list(mito = is_mito),
      BPPARAM = bp_param
    )
<<<<<<< HEAD
    df <- makePerCellDF(sces, use.coldata = TRUE, use.dimred = F) %>% 
        tidyr::pivot_longer(cols = dplyr::any_of(vars), 
                            names_to = "vv", 
                            values_to = "value")

    p <- df %>% ggplot(aes(!!sym(sample), value)) + 
          geom_violin(aes(fill = !!sym(sample)), scale = "width", width = 0.8) + 
          geom_jitter(width = 0.4, size = 0.5) + 
          facet_wrap(vars(vv), nrow = nrow, scale = "free") + 
          scale_x_discrete(name = NULL) + 
          scale_y_continuous(name = NULL) + 
          theme(legend.position = "none", 
              panel.background = element_blank(), 
              panel.border = element_rect(fill = NA), 
              strip.background = element_blank(), 
              strip.text = element_text(size = 14, face = "bold"), 
              panel.grid.major = element_line(linetype = 2, color = "grey", linewidth = 0.2), 
              axis.text = element_text(size = 12, face = "bold"))
=======
    df <- makePerCellDF(sces, use.coldata = TRUE, use.dimred = FALSE) %>% 
      tidyr::pivot_longer(cols = dplyr::any_of(vars), 
                names_to = "vv", 
                values_to = "value")

    p <- df %>% ggplot2::ggplot(ggplot2::aes_string(x = sample, y = "value")) + 
        geom_violin(mapping = ggplot2::aes_string(fill = sample), scale = "width", width = 0.8) + 
        geom_jitter(width = 0.4, size = 0.5) + 
        facet_wrap(vars(vv), nrow = nrow, scale = "free") + 
        scale_x_discrete(name = NULL) + 
        scale_y_continuous(name = NULL) + 
        theme(legend.position = "none", 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 14, face = "bold"), 
          panel.grid.major = element_line(linetype = 2, color = "grey", linewidth = 0.2), 
          axis.text = element_text(size = 12, face = "bold"))
>>>>>>> 74c6b11495aa3e3e5e418105bb2691f73cf44bba

    file <- ifelse(is.null(qc_prefix), paste0(qc_dir, "after_qc.pdf"), paste0(qc_dir, qc_prefix, "_after_qc", ".pdf"))
    ggsave(plot = p, file = file, 
             width = length(unique(colData(sces)[[sample]])) * 3, 
             height = 8, scale = 0.8, limitsize = FALSE)
    rm(df, p); gc()

    # normalization
    message("normalization...")
    set.seed(101)
    qcs <- quickCluster(sces, BPPARAM = bp_param, min.size = 100)
    sces <- computePooledFactors(sces, clusters = qcs, BPPARAM = bp_param, positive = TRUE)
    sces <- logNormCounts(sces, size.factors = sizeFactors(sces), BPPARAM = bp_param)
    rm(qcs); gc()

    message("model gene variance...")
    vars <- modelGeneVar(sces, block = colData(sces)[[sample]], BPPARAM = bp_param)
    hvgs <- getTopHVGs(vars, var.threshold = 0)
    if(length(hvgs) > top_n){
      hvgs <- hvgs[1:top_n]
    }
    message("number of hvgs used: ", length(hvgs))
    
    # remove batch effect
    message("remove batch effect...")
    set.seed(seed)
    # reducedDim(sces, "corrected") <- reducedDim(
    #   batchCorrect(sces,
    #                batch = colData(sces)[[sample]],
    #                subset.row = hvgs,
    #                correct.all = TRUE,
    #                PARAM = FastMnnParam(auto.merge = TRUE,
    #                                     BPPARAM = bp_param,
    #                                     d = n_batch_dim)),
    #   "corrected")
    # sces <- correctExperiments(
    #     sces,
    #     batch = as.factor(colData(sces)[[sample]]),  # Ensure batch matches the number of cells
    #     subset.row = hvgs,
    #     correct.all = TRUE,
    #     PARAM = FastMnnParam(
    #         auto.merge = TRUE,
    #         BPPARAM = bp_param,
    #         d = n_batch_dim
    #     )
    # )
    tmp <- fastMNN(
        sces, 
        batch = sces[[sample]], 
        BPPARAM = bp_param, 
        subset.row = hvgs, 
        correct.all = TRUE, 
        auto.merge = TRUE, 
        d = n_batch_dim
    )
    reducedDim(sces, "corrected") <- reducedDim(tmp, "corrected")
    rm(tmp);gc()


    # dimensional reduction
    message("dimensional reduction...")
    set.seed(101)
    sces <- runTSNE(sces, dimred = "corrected", BPPARAM = bp_param, num_threads = ncores) %>% 
      runUMAP(dimred = "corrected", BPPARAM = bp_param, n_threads = ncores)
    
    cluster_parallel <- function(dim_red, k) {
      BiocParallel::bplapply(
        X = list(reducedDim(sces, dim_red)),
        FUN = function(x) clusterRows(x, NNGraphParam(k = k, num.threads = ncores)), 
        BPPARAM = bp_param
      )[[1]]
    }

    # Perform clustering
    message("Performing clustering...")
    for (i in nn) {
      message(sprintf("Processing k = %d", i))
      colData(sces)[[paste0("tsne_nn_", i)]] <- cluster_parallel("TSNE", i)
      colData(sces)[[paste0("umap_nn_", i)]] <- cluster_parallel("UMAP", i)
    }
    return(sces)
}
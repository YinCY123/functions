# Convert to DelayedArray for large matrix operations with proper gene alignment
delayed_merge <- function(sce_list, assay_name = "counts") {
  library(DelayedArray)
  library(HDF5Array)
  library(SingleCellExperiment)
  
  # Get union of all gene names
  all_genes <- Reduce(union, lapply(sce_list, rownames))
  cat("Total genes in union:", length(all_genes), "\n")
  
  # Calculate total cells
  total_cells <- sum(sapply(sce_list, ncol))
  cat("Total cells:", total_cells, "\n")
  
  # Check if we can handle this size with DelayedArray
  matrix_size <- as.numeric(length(all_genes)) * as.numeric(total_cells)
  cat("Matrix size:", format(matrix_size, scientific = TRUE), "elements\n")
  
  # Create aligned DelayedArrays
  delayed_list <- lapply(names(sce_list), function(dataset_name) {
    sce <- sce_list[[dataset_name]]
    cat("Processing dataset:", dataset_name, "\n")
    
    # Get current matrix
    mat <- assay(sce, assay_name)
    
    # Create aligned matrix with union of genes
    aligned_mat <- matrix(0, nrow = length(all_genes), ncol = ncol(mat))
    rownames(aligned_mat) <- all_genes
    colnames(aligned_mat) <- colnames(mat)
    
    # Fill in existing data
    common_genes <- intersect(rownames(mat), all_genes)
    aligned_mat[common_genes, ] <- as.matrix(mat[common_genes, ])
    
    # Convert to DelayedArray
    delayed_mat <- DelayedArray(aligned_mat)
    
    cat("  Genes:", nrow(mat), "->", nrow(delayed_mat), "\n")
    cat("  Cells:", ncol(delayed_mat), "\n")
    
    return(delayed_mat)
  })
  
  # Merge all DelayedArrays
  cat("Merging DelayedArrays...\n")
  merged_delayed <- do.call(cbind, delayed_list)
  
  cat("Final dimensions:", nrow(merged_delayed), "x", ncol(merged_delayed), "\n")
  
  # Create merged SingleCellExperiment object
  cat("Creating merged SingleCellExperiment object...\n")
  
  # Create the SCE object with the merged assay
  merged_sce <- SingleCellExperiment(
    assays = list(counts = merged_delayed)
  )
  
  # Set row and column names
  rownames(merged_sce) <- rownames(merged_delayed)
  colnames(merged_sce) <- colnames(merged_delayed)
  
  # Add dataset origin information to colData
  dataset_origin <- rep(names(delayed_list), sapply(delayed_list, ncol))
  if (length(dataset_origin) != ncol(merged_sce)) {
    warning("Length of dataset_origin (", length(dataset_origin), ") does not match ncol(merged_sce) (", ncol(merged_sce), "). Adjusting to match.")
    dataset_origin <- dataset_origin[seq_len(ncol(merged_sce))]
  }
  colData(merged_sce)$dataset <- dataset_origin
  
  # Merge rowData if available
  if (all(sapply(sce_list, function(x) nrow(rowData(x)) > 0))) {
    # Get all unique rowData columns
    all_rowdata_cols <- unique(unlist(lapply(sce_list, function(x) names(rowData(x)))))
    
    # Create merged rowData
    merged_rowdata <- DataFrame(row.names = all_genes)
    for (col_name in all_rowdata_cols) {
      merged_values <- rep(NA, length(all_genes))
      names(merged_values) <- all_genes
      
      for (sce in sce_list) {
        if (col_name %in% names(rowData(sce))) {
          sce_genes <- rownames(sce)
          merged_values[sce_genes] <- rowData(sce)[[col_name]]
        }
      }
      merged_rowdata[[col_name]] <- merged_values
    }
    rowData(merged_sce) <- merged_rowdata
  }
  
  # Merge colData if available
  if (all(sapply(sce_list, function(x) ncol(colData(x)) > 0))) {
    # Get all unique colData columns
    all_coldata_cols <- unique(unlist(lapply(sce_list, function(x) names(colData(x)))))
    
    # Create merged colData
    merged_coldata <- DataFrame(row.names = colnames(merged_sce))
    merged_coldata$dataset <- dataset_origin
    
    for (col_name in all_coldata_cols) {
      if (col_name != "dataset") {  # Avoid duplicate dataset column
        merged_values <- rep(NA, ncol(merged_sce))
        names(merged_values) <- colnames(merged_sce)
        
        start_idx <- 1
        for (i in seq_along(sce_list)) {
          sce <- sce_list[[i]]
          if (col_name %in% names(colData(sce))) {
            n_cells <- ncol(sce)
            end_idx <- start_idx + n_cells - 1
            merged_values[start_idx:end_idx] <- colData(sce)[[col_name]]
          }
          start_idx <- start_idx + ncol(sce)
        }
        merged_coldata[[col_name]] <- merged_values
      }
    }
    colData(merged_sce) <- merged_coldata
  }
  
  cat("Merged SingleCellExperiment object created successfully!\n")
  return(merged_sce)
}
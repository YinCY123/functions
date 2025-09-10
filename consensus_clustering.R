# Consensus Clustering Function using ConsensusClusterPlus
#
# Parameters:
#   data: Input data matrix (features x samples)
#   normalize: Whether to z-score normalize the data (default: TRUE)
#   maxK: Maximum number of clusters to evaluate
#   reps: Number of subsamples
#   pItem: Proportion of items to sample
#   pFeature: Proportion of features to sample
#   clusterAlg: Clustering algorithm ("hc" for hierarchical, "km" for k-means)
#   distance: Distance metric ("pearson", "spearman", "euclidean", etc.)
#   seed: Random seed for reproducibility
#   plot_format: Output format for plots ("pdf" or "png", default: "pdf")
#   plot_width: Plot width in inches (default: 8)
#   plot_height: Plot height in inches (default: 6)
#   plot_dpi: Resolution for PNG output (default: 300)
#
# Returns: ConsensusClusterPlus object containing clustering results

consensus_cluster <- function(data, 
    normalize = TRUE, 
    maxK = 10, 
    reps = 100, 
    pItem = 0.8, 
    pFeature = 1,
    clusterAlg = "hc", 
    distance = "pearson", 
    finalLinkage = "ward.D2",
    seed = 101,
    dir = NULL, 
    plot = "png", 
    writeTable = TRUE, 
    title = NULL, 
    ...) {
  # Normalize data if requested
  if (normalize) {
    data <- t(scale(t(data)))  # Z-score normalization (feature-wise)
  }
  
  # Check if ConsensusClusterPlus is installed
  if (!require("ConsensusClusterPlus", quietly = TRUE)) {
    stop("Package 'ConsensusClusterPlus' is required but not installed. Please install it using: install.packages('ConsensusClusterPlus')")
  }
  
  
  # Perform consensus clustering
  results <- ConsensusClusterPlus::ConsensusClusterPlus(
    d = data,
    maxK = maxK,
    reps = reps,
    pItem = pItem,
    pFeature = pFeature,
    clusterAlg = clusterAlg,
    distance = distance,
    finalLinkage = finalLinkage,
    seed = seed,
    plot = plot,  # Disable automatic plotting
    writeTable = writeTable, 
    title = title
  )
 
  # Return clustering results
  return(results)
}

# Example usage:
# data <- read.csv("your_data.csv", row.names=1)
# results <- consensus_cluster(as.matrix(data), plot_format="png")
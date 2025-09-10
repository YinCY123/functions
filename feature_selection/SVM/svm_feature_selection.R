#' SVM-based Feature Selection (User-Friendly)
#' 
#' @param df Data frame with a 'group' column (can be any position)
#' @param k Number of folds for mSVM-RFE (default 1 = standard SVM-RFE)
#' @param halve.above Halve features above this number (default 5000)
#' @param nfold Number of cross-validation folds for error curve (default 5)
#' @param max_features Maximum number of top features to evaluate in error curve (default 10)
#' @return List with ranked features and error curve
#' @examples
#' result = svm_feature_selection(df, k=5, nfold=5, max_features=10)
#' print(result$ranking)
#' plot(result$error_curve, type='b')
svm_feature_selection <- function(df, 
    k = 1, 
    halve.above = 5000, 
    nfold = 5, 
    max_features = 10, 
    dir = NULL) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("Package 'e1071' is required. Please install it.")
  }
  if (!"group" %in% colnames(df)) {
    stop("Input data frame must have a 'group' column.")
  }
  # Move 'group' column to first position
  group_col <- which(colnames(df) == "group")
  X <- df[, c(group_col, setdiff(seq_along(df), group_col))]
  # Ensure group is a factor or numeric
  if (!is.numeric(X[[1]]) && !is.factor(X[[1]])) {
    X[[1]] <- as.factor(X[[1]])
  }
  # Source the SVM-RFE code if not already loaded
  if (!exists("svmRFE")) {
    source("/home/yincy/git/bior/functions/feature_selection/SVM/msvmRFE.R")
  }
  nrows <- nrow(X)
  folds <- rep(1:nfold, len=nrows)[sample(nrows)]
  folds <- lapply(1:nfold, function(x) which(folds == x))
  # Perform feature ranking on all training sets
  results <- lapply(folds, svmRFE.wrap, X, k=k, halve.above=halve.above)
  # Obtain top features across ALL folds
  top.features <- WriteFeatures(results, X, save=FALSE)
  # Estimate generalization error using a varying number of top features
  maxf <- min(max_features, ncol(X)-1)
  featsweep <- lapply(1:maxf, FeatSweep.wrap, results, X)
  errors <- sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
  # Plot error curve
  no.info <- min(prop.table(table(X[,1])))
  dir <- ifelse(grepl("/$", dir, paste0(dir, "/")))
  pdf(paste0(dir, "svm_error.pdf"), width = 7, height = 7)
  PlotErrors(errors, no.info=no.info, xlab='Number of Features', ylab=paste(nfold, 'fold CV Error'))
  dev.off()
  # Return both ranking and error curve
  list(ranking=top.features, error_curve=errors)
} 
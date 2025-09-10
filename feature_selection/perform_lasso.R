perform_lasso <- function(data, 
                                    features = NULL,
                                    response = NULL,          # Can be column name or Surv object
                                    survival_time = NULL,     # For survival analysis
                                    survival_status = NULL,   # For survival analysis
                                    method = c("lasso", "elastic_net"),
                                    alpha = 1,                # 1 for LASSO, between 0 and 1 for elastic net
                                    family = c("gaussian", "binomial", "cox"),
                                    nfolds = 10,
                                    seed = 101,
                                    plot = TRUE,
                                    plot_dir = NULL) {
    # Load required packages
    require(glmnet)
    require(dplyr)
    require(survival)
    
    # Input validation
    method <- match.arg(method)
    family <- match.arg(family)
    
    # Validate input parameters based on family
    if (family == "cox") {
        if (is.null(survival_time) || is.null(survival_status)) {
            if (!inherits(response, "Surv")) {
                stop("For Cox regression, either provide survival_time and survival_status, or a Surv object as response")
            }
        }
    } else {
        if (is.null(response)) {
            stop("Response variable must be specified for non-survival analysis")
        }
    }
    
    # Set features if not specified
    if (is.null(features)) {
        features <- setdiff(colnames(data), 
                          c(response, survival_time, survival_status))
    } else {
        features <- intersect(features, colnames(data))
        if (length(features) == 0) {
            stop("No valid features found in data")
        }
    }
    
    # Prepare X matrix
    X <- as.matrix(data[, features])
    
    # Prepare y based on family
    if (family == "cox") {
        y <- if (inherits(response, "Surv")) {
            response
        } else {
            Surv(data[[survival_time]], data[[survival_status]])
        }
    } else if (family == "binomial") {
        if (is.factor(data[[response]])) {
            y <- data[[response]]
        } else {
            y <- factor(data[[response]])
        }
    } else {
        y <- as.numeric(data[[response]])
    }
    
    # Set seed for reproducibility
    set.seed(seed)
    
    # Fit cross-validated model
    cv_fit <- tryCatch({
        cv.glmnet(x = X, 
                  y = y, 
                  alpha = alpha, 
                  family = family,
                  nfolds = nfolds,
                  standardize = TRUE)
    }, error = function(e) {
        stop(paste("Error in cv.glmnet:", e$message))
    })
    
    # Get coefficients at optimal lambda
    coef_min <- coef(cv_fit, s = "lambda.min")
    coef_1se <- coef(cv_fit, s = "lambda.1se")
    
    # Create feature importance dataframes
    get_nonzero_coef <- function(coef_matrix) {
        data.frame(
            feature = rownames(coef_matrix),
            coefficient = as.vector(coef_matrix)
        ) %>%
            filter(coefficient != 0, feature != "(Intercept)") %>%
            arrange(desc(abs(coefficient)))
    }
    
    selected_features_min <- get_nonzero_coef(coef_min)
    selected_features_1se <- get_nonzero_coef(coef_1se)
    
    # Generate plots if requested
    if (plot) {
        if (!is.null(plot_dir)) {
            dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
            
            # Plot 1: Cross-validation
            pdf(file.path(plot_dir, "lasso_cv_plot.pdf"))
            plot(cv_fit)
            dev.off()
            
            # Plot 2: Coefficient path
            pdf(file.path(plot_dir, "lasso_coefficient_path.pdf"))
            plot(cv_fit$glmnet.fit, xvar = "lambda")
            abline(v = log(cv_fit$lambda.min), col = "red", lty = 2)
            abline(v = log(cv_fit$lambda.1se), col = "blue", lty = 2)
            dev.off()
        } else {
            # Display plots in R
            par(mfrow = c(1, 2))
            plot(cv_fit)
            plot(cv_fit$glmnet.fit, xvar = "lambda")
            abline(v = log(cv_fit$lambda.min), col = "red", lty = 2)
            abline(v = log(cv_fit$lambda.1se), col = "blue", lty = 2)
            par(mfrow = c(1, 1))
        }
    }
    
    # Return results
    return(list(
        model = cv_fit,
        selected_features_min = selected_features_min,
        selected_features_1se = selected_features_1se,
        lambda_min = cv_fit$lambda.min,
        lambda_1se = cv_fit$lambda.1se,
        mse_min = min(cv_fit$cvm),
        features_min_n = nrow(selected_features_min),
        features_1se_n = nrow(selected_features_1se)
    ))
}
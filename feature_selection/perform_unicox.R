# Function to perform univariate Cox regression
perform_unicox <- function(data, 
    time_col = "OS_time", 
    event_col = "OS", 
    response = NULL,  # New parameter for non-survival analysis
    survival = TRUE,  # New parameter to indicate if it's survival analysis
    covariates = NULL, 
    p = 0.05, 
    clean_colnames = TRUE,
    ...) {

    # clean colnames
    if(clean_colnames){
        colnames(data) <- stringr::str_replace_all(colnames(data), "-", "_")
    }

    if (is.null(covariates)) {
        if (survival) {
            covariates <- colnames(data)[!(colnames(data) %in% c(time_col, event_col))]
        } else {
            if (is.null(response)) {
                stop("response variable must be specified for non-survival analysis")
                # covariates <- colnames(data)[!(colnames(data) %in% response)]
            }
        }
    }
    
    # Initialize results list
    results <- list()
    
    # Loop through each covariate
    for (var in covariates) {
        print(paste0("Processing variable: ", var, "..."))
        tryCatch({
            if (survival) {
                # Fit univariate Cox model for survival
                formula <- as.formula(paste("survival::Surv(", time_col, ",", event_col, ") ~ ", var))
            } else {
                # Fit univariate Cox model for non-survival
                formula <- as.formula(paste("survival::Surv(", response, ", rep(1, nrow(data))) ~ ", var))
            }
            
            model <- survival::coxph(formula, data = data)
            
            # Extract statistics
            summary_model <- summary(model)
            
            # Store results
            results[[var]] <- data.frame(
                variable = rownames(summary_model$conf.int),
                mean = round(summary_model$conf.int[1], 3),
                lower = round(summary_model$conf.int[3], 3),
                upper = round(summary_model$conf.int[4], 3),
                p_value = round(summary_model$coefficients[5], 3),
                z_score = round(summary_model$coefficients[4], 3), 
                se = round(summary_model$coefficients[3], 3)
            ) %>% 
            dplyr::mutate(ci = paste0("( ", lower, " - ", upper, " )"))
            
        }, error = function(e) {
            message(paste("Error in processing variable:", var))
            message(e$message)
        })
    }
    
    # Combine all results
    results_df <- do.call(rbind, results)
    
    # Add FDR correction
    results_df$p_adj <- round(p.adjust(results_df$p_value, method = "fdr"), 3)
    
    # Sort by p-value
    results_df <- results_df %>% dplyr::filter(p_value < p) %>% dplyr::arrange(p_value)
    
    # restore variable names
    if(clean_colnames){
        results_df <- results_df %>% 
            dplyr::mutate(variable = gsub("_", "-", variable))
        rownames(results_df) <- gsub("_", "-", rownames(results_df))
    }

    return(results_df)
}

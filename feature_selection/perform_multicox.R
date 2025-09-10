perform_multicox <- function(data, 
                           features, 
                           time_col = "OS_time", 
                           event_col = "OS",
                           clean_colnames = TRUE,
                           p = 0.05,
                           response = NULL,
                           survival = TRUE,
                           regression_type = "logistic", # Options: "linear", "logistic", "poisson"
                           ...) {
    # loading package
    suppressPackageStartupMessages(library(survival))
    suppressPackageStartupMessages(library(magrittr))
    suppressPackageStartupMessages(library(stringr))

    # clean colnames
    if(clean_colnames){
        colnames(data) <- str_replace_all(colnames(data), "-", "_")
        features <- str_replace_all(features, "-", "_")
    }

    # Prepare formula
    if(survival) {
        formula_str <- paste("survival::Surv(", time_col, ",", event_col, ") ~ ", 
                            paste(features, collapse = " + "))
    } else {
        if(is.null(response)){
            stop("Response variable must be specified for non-survival analysis")
        }else{
            formula_str <- paste(response, "~", paste(features, collapse = " + "))
        }
    }

    # Fit model based on type
    if(survival) {
        model <- survival::coxph(as.formula(formula_str), data = data)
        model_summary <- summary(model)
        results_df <- data.frame(
            variable = names(model_summary$coefficients[,1]),
            # coef = round(model_summary$coefficients[,1], 2),
            mean = round(model_summary$conf.int[,1], 3),
            lower = round(model_summary$conf.int[,3], 3),
            upper = round(model_summary$conf.int[,4], 3),
            se = round(model_summary$coefficients[,3], 3),
            z_score = round(model_summary$coefficients[,4], 3),
            p_value = round(model_summary$coefficients[,5], 3)
        )
    } else {
        # Choose regression type
        model <- switch(regression_type,
            "linear" = lm(as.formula(formula_str), data = data),
            "logistic" = glm(as.formula(formula_str), data = data, family = binomial),
            "poisson" = glm(as.formula(formula_str), data = data, family = poisson),
            stop("Unsupported regression type")
        )
        model_summary <- summary(model)
        
        # Create results dataframe
        ci <- confint(model)
        results_df <- data.frame(
            variable = names(coef(model))[-1],
            mean = round(coef(model)[-1], 3),
            lower = round(ci[-1, 1], 3),
            upper = round(ci[-1, 2], 3),
            se = round(sqrt(diag(vcov(model)))[-1], 3),
            z_value = round(model_summary$coefficients[-1, 3], 3),
            p_value = round(model_summary$coefficients[-1, 4], 3)
        )
    }
    
    # Sort by p-value and filter
    results_df <- results_df %>%
        dplyr::arrange(p_value) %>% 
        dplyr::filter(p_value < p) %>% 
        dplyr::mutate(ci = paste0("( ", lower, " - ", upper, " )"))
    
    if(clean_colnames){
        results_df$variable <- str_replace_all(results_df$variable, "_", "-")
        rownames(results_df) <- str_replace_all(rownames(results_df), "_", "-")
    }
    
    return(results_df)
}
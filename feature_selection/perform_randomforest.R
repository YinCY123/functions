# Ensure ggplot2 is available
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for plotting. Please install it.")
}

perform_randomforest <- function(x, y,
    ntree = 1000, 
    importance = TRUE, 
    seed = 101, 
    type = c("standard", "survival"), 
    plot_dir = NULL,
    ...) {

    # loading package
    require(ggplot2)
    
    # Input validation
    if (!is.matrix(x) && !is.data.frame(x)) stop("x must be a matrix or data frame")
    type <- match.arg(type)
    
    if (type == "survival") {
        if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
            stop("Package 'randomForestSRC' needed for survival analysis. Please install it.")
        }
        if (!inherits(y, "Surv")) stop("For survival analysis, y must be a Surv object")
        if (nrow(x) != nrow(y)) stop("Number of rows in x must match number of rows in y")
    } else {
        if (!is.factor(y) && !is.numeric(y)) stop("For standard analysis, y must be a factor or numeric vector")
        if (nrow(x) != length(y)) stop("Number of rows in x must match length of y")
    }
    
    # Set random seed for reproducibility
    set.seed(seed)
    if (type == "survival") {
      # Train random forest model for survival
      rf_model <- randomForestSRC::rfsrc(y ~ .,
                                         data = data.frame(y = y, x),
                                         ntree = ntree,
                                         importance = importance)

      # Get feature importance scores
      imp_scores <- randomForestSRC::vimp(rf_model)$importance
      importance_df <- data.frame(
        Feature = names(imp_scores),
        Importance = imp_scores
      )

      # For survival, rf_model$err.rate is a vector of OOB error rates
      error_plot_data <- data.frame(
        Trees = seq_along(rf_model$err.rate),
        `OOB Error` = rf_model$err.rate
      ) |>
        tidyr::pivot_longer(cols = -Trees, names_to = "Error Type", values_to = "error")

    } else {
      # Standard random forest using randomForest package
      require(randomForest)
      rf_model <- randomForest(x = x,
                               y = y,
                               ntree = ntree,
                               importance = importance)

      # Get feature importance scores
      imp_scores <- importance(rf_model)

      # For classification (factor y)
      if (is.factor(y)) {
        importance_df <- data.frame(
          Feature = rownames(imp_scores),
          Importance = imp_scores[, "MeanDecreaseAccuracy"],
          MeanDecreaseGini = imp_scores[, "MeanDecreaseGini"]
        )
      } else {
        # For regression (numeric y)
        importance_df <- data.frame(
          Feature = rownames(imp_scores),
          Importance = imp_scores[, "%IncMSE"],
          IncNodePurity = imp_scores[, "IncNodePurity"]
        )
      }

      # For standard, rf_model$err.rate is a matrix
      error_plot_data <- data.frame(
        Trees = seq_along(rf_model$err.rate[, 1]),
        `OOB Error` = rf_model$err.rate[, "OOB"]/10,
        `Class0` = rf_model$err.rate[, "0"]/10,
        `Class 1` = rf_model$err.rate[, "1"]/10,
        check.names = FALSE
      ) |>
        tidyr::pivot_longer(cols = -Trees, names_to = "Error Type", values_to = "error")
    }

    # Sort by importance score
    importance_df <- importance_df[order(-importance_df$Importance), ]

    # Create error plot
    p <- ggplot(error_plot_data, aes(x = Trees, y = error)) +
      geom_line(aes(color = `Error Type`), linewidth = 0.5) +
      labs(x = "Number of Trees",
           y = "Error Rate") +
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            legend.position = "top",
            panel.grid.major = element_line(linetype = "dashed", color = "grey", linewidth = 0.4))
    
    if(is.null(plot_dir)){
        plot_dir <- getwd()
    }
    plot_dir <- ifelse(grepl("/$", plot_dir), plot_dir, paste0(plot_dir, "/"))

    ggsave(plot = p, file = paste0(plot_dir, "rf_errors.pdf"), 
        width = 8, height = 6)

    # Return results as a list
    return(list(
        model = rf_model,
        importance = importance_df
    ))
}

# Example usage for survival data:
# library(survival)
# surv_obj <- Surv(time = survival_time, event = event_status)
# rf_results <- randomforest_select(x = mtx, y = surv_obj, type = "survival")
# print(rf_results$top_features)

perform_nn <- function(x, y, 
                     hidden_layers = c(64, 32), 
                     top_n = 30,
                     method = c("weights", "permutation", "dropout"),
                     epochs = 100,
                     batch_size = 32,
                     validation_split = 0.2,
                     seed = 42,
                     type = c("standard", "survival")) {
    
    require(keras)
    
    # Input validation
    if (!is.matrix(x) && !is.data.frame(x)) stop("x must be a matrix or data frame")
    method <- match.arg(method)
    type <- match.arg(type)
    
    # Check for required packages
    if (!requireNamespace("keras", quietly = TRUE)) {
        stop("Package 'keras' needed. Please install it.")
    }
    if (type == "survival" && !requireNamespace("tensorflow.surv", quietly = TRUE)) {
        stop("Package 'tensorflow.surv' needed for survival analysis. Please install it.")
    }
    
    # Set random seeds
    set.seed(seed)
    keras::use_session_with_seed(seed)
    
    # Normalize input features
    x_norm <- scale(x)
    
    if (type == "standard") {
        # Prepare target variable
        if (is.factor(y)) {
            num_classes <- length(unique(y))
            y_encoded <- keras::to_categorical(as.numeric(y) - 1)
            output_activation <- "softmax"
            loss <- "categorical_crossentropy"
            output_shape <- num_classes
        } else {
            y_encoded <- y
            output_activation <- "linear"
            loss <- "mse"
            output_shape <- 1
        }
        
        # Build model
        model <- keras::keras_model_sequential() %>%
            keras::layer_input(shape = ncol(x)) %>%
            keras::layer_dense(units = hidden_layers[1], activation = "relu") %>%
            keras::layer_dropout(rate = 0.3)
        
        # Add additional hidden layers
        for (units in hidden_layers[-1]) {
            model <- model %>%
                keras::layer_dense(units = units, activation = "relu") %>%
                keras::layer_dropout(rate = 0.3)
        }
        
        # Add output layer
        model <- model %>%
            keras::layer_dense(units = output_shape, activation = output_activation)
        
    } else {
        # Survival neural network
        if (!inherits(y, "Surv")) stop("For survival analysis, y must be a Surv object")
        
        # Extract time and event from Surv object
        times <- y[, "time"]
        events <- y[, "status"]
        
        # Build survival model using tensorflow.surv
        model <- tensorflow.surv::deepsurv(
            x = x_norm,
            time = times,
            event = events,
            hidden_layers = hidden_layers
        )
    }
    
    # Compile and fit model
    if (type == "standard") {
        model %>% keras::compile(
            optimizer = "adam",
            loss = loss,
            metrics = "accuracy"
        )
        
        history <- model %>% keras::fit(
            x = x_norm,
            y = y_encoded,
            epochs = epochs,
            batch_size = batch_size,
            validation_split = validation_split,
            verbose = 0
        )
    }
    
    # Feature importance calculation based on method
    if (method == "weights") {
        # Use network weights
        if (type == "standard") {
            weights <- keras::get_weights(model)[[1]]
            importance <- apply(abs(weights), 1, sum)
        } else {
            importance <- tensorflow.surv::variable_importance(model)
        }
        
    } else if (method == "permutation") {
        # Permutation importance
        importance <- numeric(ncol(x))
        baseline <- if(type == "standard") {
            keras::evaluate(model, x_norm, y_encoded)[[2]]
        } else {
            tensorflow.surv::concordance_index(model, x_norm, times, events)
        }
        
        for (i in 1:ncol(x)) {
            x_perm <- x_norm
            x_perm[,i] <- sample(x_perm[,i])
            
            score <- if(type == "standard") {
                keras::evaluate(model, x_perm, y_encoded)[[2]]
            } else {
                tensorflow.surv::concordance_index(model, x_perm, times, events)
            }
            
            importance[i] <- baseline - score
        }
        
    } else if (method == "dropout") {
        # Dropout importance
        importance <- numeric(ncol(x))
        baseline <- if(type == "standard") {
            keras::evaluate(model, x_norm, y_encoded)[[2]]
        } else {
            tensorflow.surv::concordance_index(model, x_norm, times, events)
        }
        
        for (i in 1:ncol(x)) {
            x_drop <- x_norm
            x_drop[,i] <- 0
            
            score <- if(type == "standard") {
                keras::evaluate(model, x_drop, y_encoded)[[2]]
            } else {
                tensorflow.surv::concordance_index(model, x_drop, times, events)
            }
            
            importance[i] <- baseline - score
        }
    }
    
    # Create importance dataframe
    importance_df <- data.frame(
        Feature = colnames(x),
        Importance = importance
    )
    
    # Sort by importance
    importance_df <- importance_df[order(-importance_df$Importance), ]
    
    # Select top features
    top_features <- head(importance_df, top_n)
    
    return(list(
        model = model,
        importance = importance_df,
        top_features = top_features,
        selected_features = top_features$Feature,
        history = if(exists("history")) history else NULL
    ))
}

# Example usage:
# library(keras)
# library(tensorflow.surv)  # for survival analysis
# 
# # For standard classification/regression
# groups <- factor(sample_info_gse276391$group)
# nn_results <- nn_select(x = mtx, 
#                        y = groups,
#                        method = "weights")
# 
# # For survival analysis
# library(survival)
# surv_obj <- Surv(time = survival_time, event = event_status)
# nn_results <- nn_select(x = mtx,
#                        y = surv_obj,
#                        type = "survival")

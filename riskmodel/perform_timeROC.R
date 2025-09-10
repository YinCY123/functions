perform_timeROC <- function(data, 
                          pred_col,  # Column name containing prediction scores
                          time_col = "OS_time", 
                          event_col = "OS",
                          times = c(1, 3, 5), 
                        #   dir = NULL, 
                        #   width = 6, 
                        #   height = 6, 
                        #   scale = 1, 
                          ...) {
    
    library(survivalROC)
    
    # Initialize lists to store results
    roc_list <- list()
    auc_list <- numeric(length(times))
    
    # Calculate ROC for each time point
    for (i in seq_along(times)) {
        t <- times[i]
        
        roc <- survivalROC(Stime = data[[time_col]],
                          status = data[[event_col]],
                          marker = data[[pred_col]],
                          predict.time = t,
                          method = "KM")  # Can also use "NNE" method
        
        roc_list[[i]] <- roc
        auc_list[i] <- roc$AUC
    }
    
    # Create plot
    plot_data <- data.frame(
        Time = times,
        AUC = auc_list
    )
    
    # Plot ROC curves
    # dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
    # width = width * scale
    # height = height * scale
    # pdf(file = paste0(dir, "time_ROC.pdf"), widht = width, height = height)
    colors <- rainbow(length(times))
    plot(roc_list[[1]]$FP, 
         roc_list[[1]]$TP, 
         type="l", 
         col=colors[1],
         xlab="False Positive Rate", 
         ylab="True Positive Rate",
         main="Time-dependent ROC curves", 
         lwd = 1.5)
    
    for (i in 2:length(times)) {
        lines(roc_list[[i]]$FP, roc_list[[i]]$TP, col=colors[i], lwd = 1.5)
    }
    
    legend("bottomright", 
           legend=paste0(times, "-year AUC = ", round(auc_list, 3)),
           col=colors, 
           lty=1)
    # dev.off()

    # Return results
    return(list(
        roc_data = roc_list,
        auc_values = data.frame(
            Time = times,
            AUC = auc_list
        )
    ))
}

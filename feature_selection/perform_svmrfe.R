perform_svmrfe <- function(data, 
    k = 10, 
    halve.above = 100, 
    nfold = 10, 
    plot_dir = NULL, 
    width = 10, 
    height = 4,
    seed = 101,
    ...){

    library(e1071)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(parallel)
    source("/home/yincy/disk14/git/functions/feature_selection/SVM/msvmRFE.R")

    # set up folds
    message("set up folds...")
    nfold = nfold
    nrows = nrow(data)

    set.seed(seed)
    folds = rep(1:nfold, len = nrows)[sample(nrows)]
    folds <- lapply(1:nfold, function(x){
        which(folds == x)
    })

    # perform feature ranking
    message("perform feature ranking...")
    set.seed(seed)
    results <- lapply(folds, svmRFE.wrap, data, k = k, halve.above = halve.above)

    # obtain top features across all folds
    top_features <- WriteFeatures(results, data, save = F)

    # estimate generalization error using a varying number of top features
    message("estimate error...")
    featsweep <- lapply(1:(ncol(data) - 1), FeatSweep.wrap, results, data)

    # plot error
    message("ploting error...")
    no.info <- min(prop.table(table(data[, 1])))
    errors <- sapply(featsweep, function(x){
        ifelse(is.null(x), NA, x$error)
    })

    plot_dir <- ifelse(grepl("/$", plot_dir), plot_dir, paste0(plot_dir, "/"))
    df <- data.frame(error = errors, 
        num_feats = seq_along(errors))
    df$min <- ifelse(df$error == min(df$error), paste0(df$num_feats, " - ", round(df$error, 4)), "")
    df$point <- ifelse(df$error == min(df$error), 19, NA)
    df$color <- ifelse(df$error == min(df$error), "red", "white")

    p1 <- df %>% ggplot(aes(num_feats, error)) +
        geom_point(shape = 21, size = 3, aes(color = color)) +
        geom_line(linewidth = 1, color = "grey40") +
        geom_text_repel(aes(label = min), color = "red") +
        scale_color_identity() +
        scale_x_continuous(name = "Number of Features") +
        scale_y_continuous(name = "Error") +
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA))

    # plot accuracy
    message("plotting accuracy...")
    accuracy <- sapply(errors, function(x){
        ifelse(is.na(x), NA, 1 - x)
    })

    df <- data.frame(accuracy = accuracy, 
        num_feats = seq_along(accuracy))
    
    df$max <- ifelse(df$accuracy == max(df$accuracy), paste0(df$num_feats, " - ", round(df$accuracy, 4)), "")
    df$point <- ifelse(df$accuracy == max(df$accuracy), 19, NA)
    df$color <- ifelse(df$accuracy == max(df$accuracy), "red", "white")

    p2 <- df %>% ggplot(aes(num_feats, accuracy)) +
        geom_point(shape = 21, size = 3, aes(color = color)) +
        geom_line(linewidth = 1, color = "grey40") +
        geom_text_repel(aes(label = max), color = "red") +
        scale_color_identity() +
        scale_x_continuous(name = "Number of Features") +
        scale_y_continuous(name = "Accuracy") +
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA))

    p <- p1 + p2 + plot_layout(nrow = 1)
    message("saving the figure...")
    ggsave(plot = p, 
        file = paste0(plot_dir, "SVM-RFE_error.pdf"), 
        width = width, 
        height = height)

    # return SVM features based on error
    ids <- which.min(errors)
    top_features <- top_features[1:ids, ]
    return(top_features)
}
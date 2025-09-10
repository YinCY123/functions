plot_roc <- function(roc_result, 
    dir = NULL, 
    width = 6, 
    height = 6, 
    ...){
    require(pROC)
    file = paste0(dir, "/roc_", gene, ".pdf")

    if(is.null(dir)){
        plot(roc_result, 
            print.auc = TRUE, 
            main = gene, 
            col = "red", 
            lwd = 2, 
            print.auc.cex = 1.5, 
            print.auc.col = "red", 
            grid = TRUE, 
            auc.polygon = F)
    }else{
    pdf(file, width = width, height = height)
        plot(roc_result, 
        print.auc = TRUE, 
        main = gene, 
        col = "red", 
        lwd = 2, 
        print.auc.cex = 1.5, 
        print.auc.col = "red", 
        grid = TRUE, 
        auc.polygon = F)
    dev.off()
}}

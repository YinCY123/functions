myggroc <- function(data, 
  group, 
  var, 
  name, 
  dir, 
  ...){
  suppressPackageStartupMessages(library(pROC))
  suppressPackageStartupMessages(library(ggplot2))
  
  roc_res = roc(data[[group]], data[[var]], smooth = FALSE)
  
  auc <- auc(roc_res)[1]
  auc_text = paste0("AUC = ", round(auc, 4))
  
  line <- data.frame(x = c(1, 0), y = c(0, 1))
  p = ggroc(data = roc_res,
            color = "red",
            size = 1,
            legacy.axes = FALSE) +
    # labs(title = var) +
    theme_classic() +
    geom_segment(data = line, 
      aes(x = 1, y = 0, xend = 0, yend = 1), 
        color = "grey", 
        linetype = 2) +
    annotate(geom = "text",
      x = 0.25,
      y = 0.25,
      label = auc_text,
      color = "red",
      size = 4) +
    ggtitle(name) +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15))
  
  dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
  file <- paste0(dir, name, "_ROC.pdf")
  ggsave(file, plot = p, width = 5, height = 5)
}

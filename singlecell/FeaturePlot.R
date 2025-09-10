FeaturePlot <- function(obj,
                        features = NULL, 
                        dimred = "TSNE", 
                        ncol = 2,
                        low = "grey", 
                        high = "red",
                        point_size = 1,
                        text_by = NULL,
                        text_size = 2,
                        ...){
  # loading packages
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(SingleCellExperiment))
  suppressPackageStartupMessages(require(cowplot))
  suppressPackageStartupMessages(require(scater))
  
  if(is.null(features)){
    message("Please enter at least one features...")
  }

  features <- intersect(features, rownames(obj))
  if(length(features) <= 1){
    tmp <- plotReducedDim(object = obj, 
                          dimred = dimred, 
                          color_by = features, 
                          text_by = text_by, 
                          text_size = text_size,
                          point_size = point_size) +
      scale_color_gradient(name = features, low = low, high = high)
    return(tmp)
    
  }else if(length(features) > 1){
    plot_list <- vector(mode = "list")
    for(feature in features){
      tmp <- plotReducedDim(object = obj, 
                            dimred = dimred, 
                            color_by = feature, 
                            text_by = text_by, 
                            text_size = text_size,
                            point_size = point_size) +
        scale_color_gradient(name = feature, low = low, high = high)
      
      plot_list[[feature]] <- tmp
      
      p <- plot_grid(plotlist = plot_list, ncol = ncol)
    }
    return(p)
  }
}

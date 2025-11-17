FeaturePlot <- function(obj,
                        features = NULL, 
                        dimred = "UMAP", 
                        group_by = "celltype",
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
  # suppressPackageStartupMessages(require(cowplot))
  suppressPackageStartupMessages(require(scater))
  suppressPackageStartupMessages(require(scran))
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(scattermore))


  if(is.null(features)){
    message("Please enter at least one features...")
  }

  features <- intersect(features, rownames(obj))

  df <- makePerCellDF(obj, features = features, use.coldata = TRUE, use.dimred = TRUE) %>% 
    tidyr::pivot_longer(cols = dplyr::any_of(features), names_to = "symbol", values_to = "expr")
  cell_loc <- df %>% dplyr::group_by(!!sym(group_by)) %>% dplyr::summarise(x = median(UMAP.1), y = median(UMAP.2))

  themes <- theme(panel.background = element_blank(), 
    panel.border = element_rect(fill = NA), 
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold"),
    ...)

  if(length(features) < 1){
    message("The input gene(s) does not exist in the datset.")
  }
  if(length(features) < 2){
    if(dimred == "UMAP"){
      p <- df %>% 
        ggplot(aes(UMAP.1, UMAP.2)) +
        geom_scattermore(aes(color = expr), size = point_size) +
        geom_text(data = cell_loc, aes(x, y, label = !!sym(group_by)), size = text_size) +
        scale_x_continuous(name = "UMAP 1") +
        scale_y_continuous(name = "UMAP 2") +
        scale_color_gradient(low = low, high = high) +
        themes
    }
    if(dimred == "TSNE"){
      p <- df %>% 
        ggplot(aes(UMAP.1, UMAP.2)) +
        geom_scattermore(aes(color = expr), size = point_size) +
        geom_text(data = cell_loc, aes(x, y, label = !!sym(group_by)), size = text_size) +
        scale_x_continuous(name = "UMAP 1") +
        scale_y_continuous(name = "UMAP 2") +
        scale_color_gradient(low = low, high = high) +
        themes
    }else{
      message(paste0("The dimension you input is not supported: ", dimred, "..."))
    }
  }else{
    if(dimred == "TSNE"){
      p <- df %>% 
        ggplot(aes(UMAP.1, UMAP.2)) +
        geom_scattermore(aes(color = expr), size = point_size) +
        geom_text(data = cell_loc, aes(x, y, label = !!sym(group_by)), size = text_size) +
        scale_x_continuous(name = "UMAP 1") +
        scale_y_continuous(name = "UMAP 2") +
        scale_color_gradient(low = low, high = high) +
        facet_wrap(vars(symbol), ncol = ncol, scales = "free") +
        themes
    }
    if(dimred == "UMAP"){
      p <- df %>% 
        ggplot(aes(UMAP.1, UMAP.2)) +
        geom_scattermore(aes(color = expr), size = point_size) +
        geom_text(data = cell_loc, aes(x, y, label = !!sym(group_by)), size = text_size) +
        scale_x_continuous(name = "UMAP 1") +
        scale_y_continuous(name = "UMAP 2") +
        facet_wrap(vars(symbol), ncol = ncol, scales = "free") +
        scale_color_gradient(low = low, high = high) +
        themes
    }else{
      message(paste0("The dimension you input is not supported: ", dimred, "..."))
    }
  }
    return(p)
}
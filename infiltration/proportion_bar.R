proportion_bar <- function(data, 
                           x, 
                           y, 
                           fill, 
                           colors = NULL, 
                           file = NULL, 
                           group = NULL, 
                           width = 5, 
                           height = 5,
                           scale = TRUE,
                           ...){
  
  themes <- readRDS("/home/yincy/projects/utils/themes/themes_proportion_bar.rds")
  
  # colors <- readRDS("../utils/colors.rds")
  if(is.null(colors)){
    colors <- c("#4DBBD5", "#E64B35","#00A087", "#3C5488", "#F39B7F")
  }
  
  if(scale){
    data[[y]] <- (data[[y]] - min(data[[y]])) / (max(data[[y]] - min(data[[y]])))
  }
  
  # print(paste0("There are ", length(unique(data[[fill]])), " cell types."))
  
  p <- ggplot(data = data, mapping = aes(x = !!sym(x), y = !!sym(y), fill = !!sym(fill))) +
    # geom_bar(stat = "identity", position = "fill") +
    geom_col(position = "fill", color = NA) +
    facet_wrap(vars(group), scales = "free") +
    labs(fill = "Cell Type", x = "Samples", y = "Estimated Proportion") +
    scale_fill_manual(values = colorRampPalette(colors)(length(unique(data[[fill]])))) +
    themes +
    theme(strip.background = element_blank())
  
  
  if(!is.null(file)){
    pdf(file = file, width = width, height = height)
    print(p)
    dev.off()
  }else{
    print(p)
  }
}

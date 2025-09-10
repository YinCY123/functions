boxplot_enrichment <- function(data, 
                               x, 
                               y, 
                               xlab = NULL, 
                               ylab = NULL,
                               fill = "group", 
                               add_point = FALSE,
                               group = "group", 
                               angle_x = NULL,
                               colors = NULL, 
                               file = NULL, 
                               width = 5, 
                               height = 5,
                               scale = 1,
                               trans = FALSE, 
                               ...){
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(ggpubr))
  suppressPackageStartupMessages(require(ggforce))
  
  if(is.null(colors)){
    colors <- c("#4DBBD5", "#E64B35","#00A087", "#3C5488", "#F39B7F")
  }
  
  enrich_theme <- readRDS("../utils/themes_enrich.rds")
  
  p <- ggplot(data = data, aes(x = !!sym(x), y = !!sym(y), fill = !!sym(fill))) +
    geom_boxplot(outlier.shape = 21,
                 outlier.color = NA,
                 outlier.size = 0.8,
                 alpha = 1,
                 color = "#505050",
                 linewidth = 0.1) +
    # geom_pointrange()
    # geom_sina(aes(color = !!sym(group))) +
    scale_fill_manual(values = c(colors[2], colors[1])) +
    scale_x_discrete(name = xlab) +
    scale_y_continuous(name = ylab) +
    stat_compare_means(aes(group = !!sym(group), label = ..p.signif..), 
                       bracket.size = 0.6, size = 3, 
                       label.y = max(data[[y]]) - 0.02, hide.ns = TRUE, 
                       method = ifelse(length(unique(data[[group]])) == 2, "wilcox.test", "kruskal.test")) +
    enrich_theme +
    theme(axis.text.x = element_text(angle = angle_x, hjust = 0, vjust = 0.5), 
          plot.margin = unit(c(0.2, 1.5, 0.2, 0.2), units = "cm"))
  
  if(add_point){
    p <- p + geom_point(aes(color = !!sym(group)), size = 0.5) +
      scale_color_manual(values = c(colors[2], colors[1]))
  }
  
  if(!is.null(file)){
    ggsave(filename = file, plot = p, width = width, height = height, scale = scale)
  }else{
    print(p)
  }
}

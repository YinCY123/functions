# GO and KEGG enrichment dotplot

dotplot <- function(obj, 
                    x = x, 
                    y = y,
                    sorting = "descending",
                    color = "ONTOLOGY",
                    file = NULL,
                    width = 4, 
                    height = 4,
                    palette = NULL,
                    rotate = TRUE,
                    add = "segments", 
                    xlab = NULL, 
                    ylab = "-log10(P Value)",
                    dot_size = 4.5,
                    ...){
  suppressPackageStartupMessages(library(ggpubr))
  
  # theme
  themes <- theme_pubr() +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.background = element_blank(), 
          plot.margin = margin(t = 2, r = 8, b = 2, l = 2, unit = "pt"),
          #panel.border = element_rect(size = 0.6,fill = NA,colour = "black"),
          axis.ticks = element_line(linewidth = 0.5), 
          axis.line = element_line(linewidth = 0.5),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_text(size = 8, colour = "black"), 
          axis.text.x = element_text(size = 6, colour = "black"),
          axis.title.x = element_text(size = 7, colour = "black"), 
          axis.title.y = element_blank(),
          legend.title = element_blank(), 
          legend.text = element_text(size = 7, colour = "black"),
          legend.background = element_blank(), 
          legend.key = element_blank(),
          legend.position = "top",
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.box.spacing = unit(3, "pt"))
  
  p <- ggdotchart(data = obj, 
             x = x, 
             y = y, 
             xlab = xlab, 
             ylab = ylab,
             sorting = sorting,
             color = color, 
             dot.size = dot_size, 
             rotate = rotate, 
             add = "segments", 
             ggtheme = themes)
  
  if(is.null(file)){
    print(p)
  }else{
    ggsave(filename = file, width = width, height = height, plot = p)
  }
}

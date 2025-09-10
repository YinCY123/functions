# veen plot with ggveen package

vnplot <- function(data, 
                   fill_color = NULL, 
                   set_name_size = 6,
                   stroke_color = "grey",
                   show_percentage = FALSE, 
                   text_size = 4,
                   file = NULL,
                   width = 5,
                   height = 5,
                   ...){
  
  suppressPackageStartupMessages(require(ggvenn))
  
  # if(is.null(file)){
  #   print("The file saved to the curent working directory.")
  # }
  
  if(is.null(fill_color)){
    fill_color <- c("#4DBBD5", "#E64B35","#00A087", "#3C5488", "#F39B7F")
  }
  
  p <- ggvenn(data = data, 
         fill_color = fill_color, 
         stroke_color = stroke_color, 
         show_percentage = FALSE, 
         set_name_size = set_name_size, 
         text_size = text_size)
  
  if(is.null(file)){
    print(p)
  }else{
    ggsave(filename = file, width = width, height = height)
  }
}
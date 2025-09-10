kegg_enrich_bar <- function(obj, 
                            top = 15, 
                            x_axis = "pvalue", 
                            y_axis = "Description", 
                            fill = "pvalue", 
                            log_fill = TRUE, 
                            xlab = "-log10(p value)", 
                            ylab = NULL, 
                            remove_suffix = TRUE, 
                            organism = "mmu", 
                            fill_lab = "-log10(p value)", 
                            low_color = "blue", 
                            high_color = "red", 
                            str_width = 40, 
                            file = NULL, 
                            width = 7, 
                            height = 7, 
                            scale = 1, 
                            ...){
  # loading required packages
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(stringr))
  
  # prepare ploting data
  tmp <- obj@result

  # remove suffix
  if(remove_suffix & organism == "mmu"){
    tmp <- tmp %>% dplyr::mutate(Description = str_remove(Description, " - Mus musculus \\(house mouse\\)"))
  }else if(remove_suffix & organism == "hsa"){
    tmp <- tmp %>% dplyr::mutate(Description = str_remove(Description, " - Homo sapiens \\(human\\)"))
  }

  tmp <- tmp %>% dplyr::mutate(x = -log10(!!sym(x_axis)), 
                               y = str_wrap(!!sym(y_axis), width = str_width)) %>% 
    dplyr::arrange(!!sym(x_axis)) %>% 
    dplyr::slice_head(n = top)

  if(log_fill){
    tmp <- tmp %>% dplyr::mutate(fill = -log10(!!sym(fill)))
  }else{
    tmp <- tmp %>% dplyr::mutate(fill = !!sym(fill))
  }
  
  # generate the plot
  p <- tmp %>% 
    ggplot(aes(x, reorder(y, x, decreasing = F))) +
    geom_bar(stat = "identity", aes(fill = fill)) +
    scale_x_continuous(name = xlab) +
    scale_y_discrete(name = NULL) +
    scale_fill_gradient(name = fill_lab, low = low_color, high = high_color) +
    theme_bw()
  
  # save or not
  if(is.null(file)){
      return(p)
  }else{
      ggsave(filename = file, width = width, height = height, scale = scale, plot = p)
  }
}

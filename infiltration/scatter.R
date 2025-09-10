scatters <- function(expr_data,
                     enrich_data,
                     x, 
                     y, 
                     colors = NULL,
                     dir = NULL, 
                     width = 4.5, 
                     height = 4.5, 
                     units = "cm", 
                     xlab_size = 6, 
                     ylab_size = 6,
                     labelx_prefix = "Expression level of ", 
                     labelx_suffix = NULL, 
                     labely_prefix = NULL, 
                     labely_suffix = "\nlog2(TPM + 1)",
                     cor_type = "pearson",
                     return_cor = TRUE,
                     ...){
  
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(ggside))
  suppressPackageStartupMessages(require(ggpubr))
  
  # convert to matrix
  expr_data <- as.matrix(expr_data)
  enrich_data <- as.matrix(enrich_data)

  df <- data.frame(x_value = if(x %in% rownames(expr_data)){expr_data[x, ]}else{enrich_data[x, ]}, 
                   y_value = if(y %in% rownames(expr_data)){expr_data[y, ]}else{enrich_data[y, ]})
  
  cor <- cor.test(x = df[["x_value"]], 
                  y = df[["y_value"]], 
                  method = cor_type)

  # color
  if(is.null(colors)){
    colors <- c("#4DBBD5", "#E64B35","#00A087", "#3C5488", "#F39B7F")
  }
  
  # define themes
  themes_scaters <- theme(plot.title = element_blank(), 
                          plot.subtitle = element_blank(), 
                          plot.background = element_blank(),
                          plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit="pt"), 
                          panel.border = element_rect(linewidth = 0.45, fill = NA, colour = "black"), 
                          panel.background = element_blank(), 
                          panel.grid = element_blank(), 
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(), 
                          axis.title = element_text(size = 6, colour = "black"), 
                          axis.text.y = element_text(size = 5, colour = "black"))
  
  # plot
  p <- ggplot(data = df, aes(x_value, y_value)) +
    xlim(min(df[["x_value"]]), max(df[["x_value"]])) +
    ylim(min(df[["y_value"]]), max(df[["y_value"]])) +
    labs(x = paste0(labelx_prefix, x, labelx_suffix), y = paste0(labely_prefix, y, labely_suffix)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm", color = ifelse(cor$estimate < 0, colors[1], colors[2])) +
    # geom_smooth(method = "lm", color = "blue") +
    # geom_xsidedensity() +
    # geom_ysidedensity() +
    geom_xsidedensity(color = ifelse(cor$estimate < 0, colors[1], colors[2]),
                      fill = ifelse(cor$estimate < 0, colors[1], colors[2])) +
    geom_ysidedensity(color = ifelse(cor$estimate < 0, colors[1], colors[2]),
                      fill = ifelse(cor$estimate < 0, colors[1], colors[2])) +
    stat_cor(digits = 4, size = 2.5, label.sep = "\n", label.y.npc = ifelse(cor$estimate > 0, "top", "bottom")) +
    themes_scaters
  
  # save or not
  if(cor$p.value < 0.05 & !is.null(dir)){
    ggsave(plot = p, filename = paste0(dir, x, "_", y, "_scatter.pdf"), width = width, height = height, units = units)
  }
  
  # return cor or not
  if(return_cor){
    return(cor)
  }
}
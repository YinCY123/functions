gsea_plot <- function(x, 
                      geneSetID = NULL, 
                      line_color = "green",
                      line_lwd = 1,
                      
                      vline_type = 2, 
                      vline_color = "grey",
                      vline_lwd = 1,
                      
                      hits_color = "black", 
                      hits_lwd = 0.5, 
                      hits_lty = 1,
                      
                      left_color = "red", 
                      right_color = "blue", 
                      
                      ranking_color = "grey",
                      ...){
  
  # loading required packages
  library(ggplot2)
  library(patchwork)
  library(enrichplot)
  
  
  # extract running data
  df <- gseaplot2(x, geneSetID = geneSetID, subplots = 1)$data
  
  # stat anno data
  stat_anno <- x@result %>% dplyr::filter(ID == geneSetID)
  
  # vertical line data
  ids <- which.max(df$runningScore)
  hline <- df[ids, ]
  
  xmid <- nrow(df)/2
  
  p1 <- df %>% 
    ggplot(aes(x, runningScore)) +
    geom_line(color = line_color, linewidth = line_lwd) +
    geom_hline(yintercept = 0) +
    geom_text(data = stat_anno, aes(x = xmid, y = 0.05, 
                                    label = paste0("NES = ", round(stat_anno[,"NES", drop = T], 3), "; ",
                                                   "pvalue = ", round(stat_anno[, "pvalue", drop = T], 3), "; ",  
                                                   "qvalue = ", round(stat_anno[, "qvalue", drop = T], 3)))) +
    scale_x_continuous(name = NULL) +
    scale_y_continuous(name = "Enrichment Score") +
    ggtitle(label = geneSetID) +
    geom_vline(xintercept = hline$x, linetype = vline_type, color = vline_color, linewidth = vline_lwd) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(fill = NA), 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "bold"), 
          ...)
  
  
  # vertical line and color bar
  p21 <- df %>% 
    dplyr::filter(position == 1) %>% 
    ggplot() +
    geom_segment(aes(x = x, xend = x,  y = 0, yend = 1), color = hits_color, linewidth = hits_lwd,  linetype = hits_lty) +
    scale_x_continuous(name = NULL, expand = c(0, 0)) +
    scale_y_continuous(name = NULL, expand = c(0, 0)) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(fill = NA), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  bar_df <- data.frame(x = seq(1, max(df$x)), y = 1)
  p22 <- bar_df %>% 
    ggplot(aes(x, y, fill = x)) +
    geom_raster() +
    scale_fill_gradient(low = left_color, high = right_color) +
    scale_x_continuous(name = NULL, expand = c(0,0)) +
    scale_y_continuous(name = NULL, expand = c(0,0)) +
    theme_void() +
    theme(legend.position = "none")
  
  
  # ranking plot
  rank_df <- data.frame(
    gene = names(gsea@geneList), 
    rank = seq_along(gsea@geneList), 
    logFC = gsea@geneList
  )
  
  
  p3 <- rank_df %>% ggplot() +
    geom_segment(aes(x = rank, xend = rank, y = 0, yend = logFC), color = ranking_color) +
    scale_x_continuous(name = NULL) +
    scale_y_continuous(name = "Ranking Metric") +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA))
  
  p2 <- (p21 / p22) + plot_layout(heights = c(1, 0.25))
  
  p <- (p1 / p2 / p3) + plot_layout(height = c(1, 0.4, 0.3))
  
  return(p)
}

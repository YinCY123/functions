volcano <- function(tbl, 
                    x = "log2FoldChange", 
                    y = "pvalue", 
                    label = "symbol",
                    is_label = TRUE,
                    xlab = NULL, 
                    ylab = NULL, 
                    x_threshold = 1, 
                    y_threshold = 0.05,
                    x_lim = NULL,
                    y_lim = NULL,
                    features = NULL,
                    top = 10,
                    colors = NULL, 
                    left = -5, 
                    right = 5,
                    step = 2,
                    dirction_levels = c("up", "ns", "down"),
                    text_size = 2,
                    point_size = 1,
                    max_overlaps = 10,
                    file = NULL,
                    width = 5, 
                    height = 6, 
                    units = "cm",
                    title = NULL, 
                    scale = 1,
                    ...){
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(ggrepel))
  suppressPackageStartupMessages(require(dplyr))
  
  # remove rows contain NA
  tbl <- na.omit(tbl)
  
  # x-axis limit
  tmp_right = max(abs(tbl[[x]])) %>% ceiling()

  if(right < tmp_right){
    right = tmp_right
  }
  left = -right

  breaks <- c(seq(floor(left/step) * step, 0, by = step) ,
    seq(step, ceiling(right/step) * step, by = step))

  # color
  if(is.null(colors)){
    colors <- c("tomato", "steelblue", "grey")
  }
  
  # add columns
  tbl[["logp"]] = -log10(tbl[[y]])
  tbl <- tbl %>% dplyr::mutate(direction = ifelse(!!sym(x) > x_threshold & !!sym(y) < y_threshold, "up", 
                                           ifelse(!!sym(x) < -x_threshold & !!sym(y) < y_threshold, "down", "ns")), 
                               direction = factor(direction, levels = dirction_levels))

  # annotation data frame
  if(is_label){
    if(!is.null(features)){
      ann_df = tbl %>% dplyr::filter(!!sym(label) %in% features)
    }else{
      ann_df <- rbind(tbl %>% dplyr::filter(direction == "up") %>% dplyr::arrange(!!sym(y)) %>% dplyr::slice_head(n = top), 
                      tbl %>% dplyr::filter(direction == "down") %>% dplyr::arrange(!!sym(y)) %>% dplyr::slice_head(n = top))
    }
  }
  
  p <- ggplot(data = tbl, aes(x = !!sym(x), y = logp)) +
    geom_point(aes(color = direction), size = point_size) +
    geom_vline(xintercept = c(x_threshold, -x_threshold), color = "gray", linetype = 2) +
    geom_hline(yintercept = -log10(y_threshold), color = "gray", linetype = 2) +
    scale_color_manual(name = NULL, values = c("up" = colors[1],  "down" = colors[2], "ns" = colors[3])) +
    scale_x_continuous(name = "log2(Fold Change)",
                       limits = c(left, right),
                       breaks = breaks, 
                       labels = breaks
                       ) +
    scale_y_continuous(name = "-log10 (P Value)") +
    labs(title = title) +
    theme(plot.subtitle = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(linewidth = 0.6, fill = NA, colour = "black"),
          ...) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  if(is_label){
    p <- p + 
          geom_text_repel(data = ann_df, 
                    aes(label = !!sym(label)), 
                    size = text_size, 
                    max.overlaps = max_overlaps)
  }

  # save or not
  if(is.null(file)){
      print(p)
  }else{
      ggsave(filename = file, plot = p, width = width, height = height, scale = scale)
  }
}

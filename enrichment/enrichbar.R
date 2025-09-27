library(fgsea)
library(ggplot2)
library(dplyr)

enrichbar <- function(res, 
                      x, 
                      y, 
                    #   xlab = "-log10(P Value)", 
                    #   ylab = NULL, 
                      filter = "NES", 
                      colors = NULL, 
                      legend.position = "top", 
                      top = 15, 
                      text_size = 3, 
                      left = NULL, 
                      right = NULL, 
                      fill_lab = "group", 
                      step_len = 2, 
                      file = NULL, 
                      width = 7, 
                      height = 7, 
                      scale = 1, 
                      offset = 0.01,
                      ...){
    
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(dplyr))
    
    # add direction info
    up <- res %>% as.data.frame %>% 
        dplyr::filter(!!sym(filter) > 0) %>% 
        dplyr::arrange(!!sym(x)) %>% 
        dplyr::slice_head(n = top)
    up$group <- "up"

    down <- res %>% as.data.frame %>% 
        dplyr::filter(!!sym(filter) < 0) %>% 
        dplyr::arrange(!!sym(x)) %>% 
        dplyr::slice_head(n = top)
    down$group <- "down"
    
    # color
    if(is.null(colors)){
      colors <- c("tomato", "steelblue")
    }
    
    # prepare plot data
    df_to_plot <- rbind(up, down) %>% 
        dplyr::mutate(logp = ifelse(group == "up", -log10(!!sym(x)), log10(!!sym(x)))) %>% 
        dplyr::arrange(desc(logp)) %>% 
        dplyr::mutate(y = nrow(.):1)
    
    # custome x axis
    x_max <- max(abs(df_to_plot$logp))
    x_min <- -x_max

    if(!is.null(left)){
        x_min <- min(x_min, left)
    }
    
    if(!is.null(right)){
        x_max <- max(x_max, right)
    }
    
    # determine step
    # step <- ifelse(((max(df_to_plot$logp) - min(df_to_plot$logp)) %/% 5) > 3, 5, 2)
    if(!is.null(step_len)){
        step = step_len
    }

    xlab = switch(as.character(x), 
        "pvalue" = "-log10(p value)", 
        "p.adjust" = "-log10(p adjust)", 
        "qvalue" = "-log10(q value)", 
        x = x)

    
    # plot
    p <- df_to_plot %>% 
        ggplot(aes(logp, reorder(!!sym(y), logp))) +
        geom_bar(aes(fill = group), stat = "identity") +
        geom_text(data = df_to_plot[df_to_plot$group == "up", ], aes(-offset, !!sym(y), label = !!sym(y), hjust = 1), size = text_size) +
        geom_text(data = df_to_plot[df_to_plot$group == "down", ], aes(offset, !!sym(y), label = !!sym(y), hjust = 0), size = text_size) +
        scale_fill_manual(name = fill_lab, values = c("up" = colors[1], "down" = colors[2])) +
        scale_y_discrete(name = NULL, label = NULL) +
        scale_x_continuous(name = xlab, 
                           limits = c(x_min, x_max), 
                           breaks = seq(x_min, x_max, step_len)) +
        theme_bw() +
        theme(axis.ticks.y = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(), 
              legend.position = legend.position, 
              panel.grid = element_blank())
    
    # save or print
    if(!is.null(file)){
        ggsave(file = file, width = width, height = height, scale = scale, plot = p)
    }else{
      print(p)
    }
}
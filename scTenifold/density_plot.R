density_plot <- function(x, 
    variable, 
    fill_color = NULL,
    x_name = "Z score",

    # density params
    bw = "nrd0", 
    adjust = 1, 
    kernel = "gaussian", 
    weight = NULL, 
    window = "gaussian",

    # geom_ribbon
    lower = -1, 
    upper = 1, 
    scale_y = 1.2,
    ...
    ){

    # loading required packages
    library(ggplot2)
    library(MetBrewer)
    library(magrittr)

    # approximate density
    dens <- density(x = x[[variable]], bw = bw, adjust = adjust, kernel = kernel, weight = weight)

    # model a function
    model <- approxfun(x = dens$x, y = dens$y, yleft = 0, yright = 0, method = "linear")

    lower_df = data.frame(x = dens$x[dens$x < lower], y = model(dens$x[dens$x < lower]))
    lower_df$group <- "down"

    middle <- data.frame(x = dens$x[dens$x >= lower & dens$x <= upper], 
        y = model(dens$x[dens$x >= lower & dens$x <= upper]))
    middle$group <- "ns"

    upper_df = data.frame(x = dens$x[dens$x > upper], y = model(dens$x[dens$x > upper]))
    upper_df$group <- "up"

    anno_df <- Reduce(rbind, list(lower_df, middle, upper_df))

    if(is.null(fill_color)){
        cols <- met.brewer("Java", 2)
        fill_color <- c(cols[1], "grey", cols[2])
        names(fill_color) <- c("down", "ns", "up")
    }

    # up- and down-regulated genes
    num_down <- sum(x[[variable]] < lower)
    num_up <- sum(x[[variable]] > upper)


    x_label <- c(seq(floor(min(anno_df$x)), 0, 1), seq(0, ceiling(max(anno_df$x)), 1))

    p <- anno_df %>% 
        ggplot(aes(x, y)) +
        geom_ribbon(aes(ymin = 0, ymax = y, fill = group)) +
        scale_fill_manual(values = fill_color) +
        scale_x_continuous(name = x_name, breaks = x_label) +
        scale_y_continuous(name = "Density") +
        # geom_text(aes(x = min(anno_df$x), y = max(anno_df$y) * scale_y), )
        annotate(geom = "text", 
            x = min(anno_df$x), 
            y = max(anno_df$y) * scale_y, 
            color = fill_color[c(1,3)],
            label = paste0("number of down-regulated genes: ", num_down, "\nnumber of up-regulated genes: ", num_up), 
            hjust = "inward") +
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA), 
            ...)
    return(p)
}
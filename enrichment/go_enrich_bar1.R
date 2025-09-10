go_enrich_bar1 <- function(go_enrich, 
    x = c("pvalue", "p.adjust", "qvalue", "zScore", "FoldEnrichment", "RichFactor"),
    color = NULL, 
    file = NULL, 
    width = 8, 
    height = 18, 
    ...) {
    library(ggprism)
    library(dplyr)
    library(ggplot2)
    library(tibble)
    library(gground)

    if(is.null(color)){
        pal <- c('#7bc4e2', '#acd372', '#fbb05b')
    }

    x <- match.arg(x, c("pvalue", "p.adjust", "qvalue", "zScore", "FoldEnrichment", "RichFactor"))

    GO <- go_enrich@result
    use_pathway <- group_by(GO, ONTOLOGY) %>%
        dplyr::top_n(10, wt = -p.adjust) %>%
        dplyr::group_by(p.adjust) %>%
        dplyr::top_n(1, wt = Count) %>%
        ungroup() %>%
        dplyr::arrange(ONTOLOGY, p.adjust) %>%
        mutate(Description = factor(Description, levels = Description)) %>%
        tibble::rowid_to_column('index')

    width <- 1
    xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
    rect.data <- group_by(use_pathway, ONTOLOGY) %>%
        reframe(n = n()) %>%
        ungroup() %>%
        mutate(
            xmin = -3 * width,
            xmax = -2 * width,
            ymax = cumsum(n),
            ymin = lag(ymax, default = 0) + 0.6,
            ymax = ymax + 0.4
        )

    p <- use_pathway %>%
        ggplot(aes(x = -log10(p.adjust), y = index, fill = ONTOLOGY)) +
        geom_round_col(
            aes(y = Description), width = 0.6, alpha = 0.8
        ) +
        geom_text(
            aes(x = 0.05, label = Description),
            hjust = 0, size = 5
        ) +
        geom_text(
            aes(x = 0.1, label = geneID, colour = ONTOLOGY), 
            hjust = 0, vjust = 3, size = 3.5, 
            show.legend = FALSE
        ) +
        geom_point(
            aes(x = -width, size = Count),
            shape = 21
        ) +
        geom_text(
            aes(x = -width, label = Count)
        ) +
        scale_size_continuous(name = 'Count', range = c(5, 16)) +
        geom_round_rect(
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = ONTOLOGY),
            data = rect.data,
            radius = unit(2, 'mm'),
            inherit.aes = FALSE
        ) +
        geom_text(
            aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
            data = rect.data,
            inherit.aes = FALSE
        ) +
        geom_segment(
            aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
            linewidth = 1.5,
            inherit.aes = FALSE
        ) +
        labs(y = NULL) +
        scale_fill_manual(name = 'Category', values = pal) +
        scale_colour_manual(values = pal) +
        scale_x_continuous(
            breaks = seq(0, xaxis_max, 2), 
            expand = expansion(c(0, 0))
        ) +
        theme_prism() +
        theme(
            axis.text.y = element_blank(),
            axis.line = element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_text()
        )
    if(is.null(file)){
        print(p)
    }else{
        ggsave(file, plot = p, width = width, height = height, limitsize = FALSE)
    }
}
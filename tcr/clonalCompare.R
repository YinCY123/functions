# Generated from function body. Editing this file has no effect.
clonalCompare <- function (input.data, cloneCall = "strict", chain = "both", samples = NULL, 
    clones = NULL, top.clones = NULL, highlight.clones = NULL, 
    relabel.clones = FALSE, group.by = NULL, order.by = NULL, 
    graph = "alluvial", proportion = TRUE, exportTable = FALSE, 
    palette = "inferno", ...) 
{
    if (!is.null(top.clones) && !is.null(clones)) {
        top.clones <- NULL
    }
    input.data <- .dataWrangle(input.data, group.by, .theCall(input.data, 
        cloneCall, check.df = FALSE, silent = TRUE), chain)
    cloneCall <- .theCall(input.data, cloneCall)
    sco <- .is.seurat.or.se.object(input.data)
    if (!is.null(group.by) && !sco) {
        input.data <- .groupList(input.data, group.by)
    }
    compareColname <- ifelse(proportion, "Proportion", "Count")
    normalizer <- ifelse(proportion, sum, length)
    Con.df <- input.data %>% purrr::imap(function(df, columnNames) {
        tbl <- as.data.frame(table(df[, cloneCall]))
        if (proportion) {
            tbl[, 2] <- tbl[, 2]/normalizer(tbl[, 2])
        }
        colnames(tbl) <- c("clones", compareColname)
        tbl$Sample <- columnNames
        tbl
    }) %>% dplyr::bind_rows()
    if (!is.null(samples)) {
        Con.df <- Con.df[Con.df$Sample %in% samples, ]
    }
    if (!is.null(clones)) {
        Con.df <- Con.df[Con.df$clones %in% clones, ]
    }
    else if (!is.null(top.clones)) {
        top <- Con.df %>% group_by(Sample) %>% slice_max(n = top.clones, 
            order_by = !!sym(compareColname), with_ties = FALSE)
        Con.df <- Con.df[Con.df$clones %in% top$clones, ]
    }
    if (nrow(Con.df) < length(unique(Con.df$Sample)) || nrow(Con.df) == 
        0) {
        stop("Please reasses the filtering strategies here, there are not\n            enough clones to examine.")
    }
    clones.returned <- as.vector(unique(Con.df[order(Con.df[, 
        compareColname], decreasing = TRUE), "clones"]))
    if (relabel.clones) {
        new.clones <- paste0("Clone: ", seq_len(length(clones.returned)))
        names(new.clones) <- clones.returned
        if (!is.null(highlight.clones)) {
            highlight.clones <- unname(new.clones[which(names(new.clones) %in% 
                highlight.clones)])
        }
        Con.df[, "original.clones"] <- Con.df[, "clones"]
        Con.df[, "clones"] <- new.clones[as.vector(Con.df[, "clones"])]
        Con.df[, "clones"] <- factor(Con.df[, "clones"], levels = .alphanumericalSort(unique(Con.df[, 
            "clones"])))
        clones.returned <- as.vector(unique(Con.df[, "clones"]))
    }
    if (!is.null(order.by)) {
        Con.df <- .orderingFunction(vector = order.by, group.by = "Sample", 
            data.frame = Con.df)
    }
    if (exportTable) {
        return(Con.df)
    }
    plot <- ggplot(Con.df, aes(x = .data[["Sample"]], fill = clones, 
        group = clones, stratum = clones, alluvium = clones, 
        y = !!sym(compareColname), label = clones)) + .themeRepertoire() + 
        theme(axis.title.x = element_blank(), legend.text = element_text(size = rel(0.5)), 
            legend.key.size = unit(0.5, "line"), ...)
    if (graph == "alluvial") {
        plot <- plot + geom_stratum() + geom_flow(stat = "alluvium")
    }
    else if (graph == "area") {
        plot <- plot + geom_area(aes(group = clones), color = "black")
    }
    if (!is.null(highlight.clones)) {
        clone.colors <- rep("grey", length(clones.returned))
        pos <- which(clones.returned %in% highlight.clones)
        clone.colors[pos] <- .colorizer(palette, length(pos))
        names(clone.colors) <- clones.returned
        plot <- plot + scale_fill_manual(values = clone.colors)
    }
    else {
        plot <- plot + scale_fill_manual(values = .colorizer(palette, 
            length(unique(Con.df[, "clones"]))))
    }
    return(plot)
}

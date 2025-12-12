# Generated from function body. Editing this file has no effect.
plotDots <- function (object, features, group = NULL, block = NULL, exprs_values = "logcounts", 
    detection_limit = 0, zlim = NULL, colour = color, color = NULL, 
    max_detected = NULL, other_fields = list(), by_exprs_values = exprs_values, 
    swap_rownames = NULL, center = FALSE, scale = FALSE, assay.type = exprs_values, 
    by.assay.type = by_exprs_values) 
{
    if (is.null(group)) {
        group <- rep("all", ncol(object))
    }
    else {
        group <- retrieveCellInfo(object, group, search = "colData")$value
    }
    object <- .swap_rownames(object, swap_rownames)
    features <- .handle_features(features, object)
    group <- factor(group)
    ids <- DataFrame(group = group)
    if (!is.null(block)) {
        ids$block <- retrieveCellInfo(object, block, search = "colData")$value
    }
    summarized <- summarizeAssayByGroup(assay(object, assay.type)[as.character(features), 
        , drop = FALSE], ids = ids, statistics = c("mean", "prop.detected"), 
        threshold = detection_limit)
    ave <- assay(summarized, "mean")
    num <- assay(summarized, "prop.detected")
    group.names <- summarized$group
    if (!is.null(block)) {
        ave <- correctGroupSummary(ave, group = summarized$group, 
            block = summarized$block)
        num <- correctGroupSummary(num, group = summarized$group, 
            block = summarized$block, transform = "logit")
        group.names <- factor(colnames(ave), levels = levels(summarized$group))
    }
    heatmap_scale <- .heatmap_scale(ave, center = center, scale = scale, 
        colour = colour, zlim = zlim)
    evals_long <- data.frame(Feature = rep(features, ncol(num)), 
        Group = rep(group.names, each = nrow(num)), NumDetected = as.numeric(num), 
        Average = as.numeric(heatmap_scale$x))
    if (!is.null(max_detected)) {
        evals_long$NumDetected <- pmin(max_detected, evals_long$NumDetected)
    }
    vis_out <- .incorporate_common_vis_row(evals_long, se = object, 
        colour_by = NULL, shape_by = NULL, size_by = NULL, by.assay.type = by.assay.type, 
        other_fields = other_fields, multiplier = rep(.subset2index(features, 
            object), ncol(num)))
    evals_long <- vis_out$df
    p <- ggplot(evals_long) + geom_point(aes(x = .data$Group, y = .data$Feature, 
        size = .data$NumDetected, col = .data$Average)) + scale_size(limits = c(0, 
        max(evals_long$NumDetected))) + heatmap_scale$colour_scale + 
        theme(panel.background = element_rect(fill = "white"), 
            panel.grid.major = element_line(linewidth = 0.5, 
                colour = "grey80"), panel.grid.minor = element_line(linewidth = 0.25, 
                colour = "grey80"))

    return(list(p, evals_long))
}

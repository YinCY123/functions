netVisual_embeddingPairwise <- function(object, 
          slot.name = "netP", 
          type = c("functional","structural"), 
          comparison = NULL, 
          color.use = NULL, 
          point.shape = NULL, 
          pathway.labeled = NULL, 
          top.label = 1, 
          pathway.remove = NULL, 
          pathway.remove.show = TRUE, 
          dot.size = c(2, 6), 
          label.size = 2.5, 
          dot.alpha = 0.5,
          xlabel = "Dim 1", 
          ylabel = "Dim 2", 
          title = NULL,
          do.label = T, 
          show.legend = T, 
          show.axes = T) {
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("2D visualization of signaling networks from datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  object.names <- setdiff(names(methods::slot(object, slot.name)), "similarity")[comparison]
  prob <- list()
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, slot.name)[[comparison[i]]]
    prob[[i]] = object.net$prob
  }

  if (is.null(point.shape)) {
    point.shape <- c(21, 0, 24, 23, 25, 10, 12)
  }

  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
    # pathway.remove <- sub("--.*", "", pathway.remove)
  }

  if (length(pathway.remove) > 0) {
    for (i in 1:length(prob)) {
      probi <- prob[[i]]
      pathway.remove.idx <- which(paste0(dimnames(probi)[[3]],"--",object.names[i]) %in% pathway.remove)
      #  pathway.remove.idx <- which(dimnames(probi)[[3]] %in% pathway.remove)
      if (length(pathway.remove.idx) > 0) {
        probi <- probi[ , , -pathway.remove.idx]
      }
      prob[[i]] <- probi
    }
  }
  prob_sum.each <- list()
  signalingAll <- c()
  for (i in 1:length(prob)) {
    probi <- prob[[i]]
    prob_sum.each[[i]] <- apply(probi, 3, sum)
    signalingAll <- c(signalingAll, paste0(names(prob_sum.each[[i]]),"--",object.names[i]))
  }
  prob_sum <- unlist(prob_sum.each)
  names(prob_sum) <- signalingAll

  group <- sub(".*--", "", names(prob_sum))
  labels = sub("--.*", "", names(prob_sum))

  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),
                   labels = as.character(labels), clusters = as.factor(clusters), group = factor(group, levels = unique(group)))
  # color dots (light inside color and dark border) based on clustering and no labels
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }
  gg <- ggplot(data = df, aes(x, y)) +
    geom_point(aes(size = Commun.Prob.,fill = clusters, colour = clusters, shape = group)) +
    CellChat_theme_opts() +
    theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) +
    scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE) #+ scale_alpha(group, range = c(0.1, 1))
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  gg <- gg + scale_shape_manual(values = point.shape[1:length(prob)])

  if (do.label) {
    if (is.null(pathway.labeled)) {
      if (top.label < 1) {
        if (length(comparison) == 2) {
          g.t <- rankSimilarity(object, slot.name = slot.name, type = type, comparison1 = comparison)
          pathway.labeled <- as.character(g.t$data$name[(nrow(g.t$data)-ceiling(top.label * nrow(g.t$data))+1):nrow(g.t$data) ])
          data.label <- df[df$labels %in% pathway.labeled, , drop = FALSE]
        }
      } else {
        data.label <- df
      }

    } else {
      data.label <- df[df$labels %in% pathway.labeled, , drop = FALSE]
    }
    gg <- gg + 
      ggrepel::geom_text_repel(mapping = aes(label = labels, colour = clusters, alpha=group), 
          size = label.size, 
          show.legend = F,
          segment.size = 0.2, 
          segment.alpha = 0.5) + 
        scale_alpha_discrete(range = c(1, 0.6))
  }

  if (length(pathway.remove) > 0 & pathway.remove.show) {
    gg <- gg + annotate(geom = 'text', label =  paste("Isolate pathways: ", paste(pathway.remove, collapse = ', ')), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = label.size,fontface="italic")
  }

  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }
  return(list(gg, df))
}
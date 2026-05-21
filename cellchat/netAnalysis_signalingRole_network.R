netAnalysis_signalingRole_network <- function(object, 
    signaling, slot.name = "netP", 
    measure = c("outdeg","indeg","flowbet","info"), 
    measure.name = c("Sender","Receiver","Mediator","Influencer"),
    color.use = NULL, 
    title = NULL,
    rect_gp = NULL,
    color.heatmap = "purple",
    font.size = 8, 
    font.size.title = 10, 
    cluster.rows = FALSE, 
    cluster.cols = FALSE, 
    ...) {
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr[signaling]
  for(i in 1:length(centr)) {
    centr0 <- centr[[i]]
    mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
    mat <- t(mat)
    rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
    if (!is.null(measure)) {
      mat <- mat[measure,,drop = FALSE]
      if (!is.null(measure.name)) {
        if (length(measure.name) != length(measure)) {
          stop("The length of `measure.name` is not the same as that of `measure`! Please modify it! \n")
        }
        rownames(mat) <- measure.name
      }
    }
    mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)

    if (is.null(color.use)) {
      color.use <- scPalette(length(colnames(mat)))
    }
    # color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
    color.heatmap.use <- circlize::colorRamp2(breaks = c(0, 1), colors = c("white", color.heatmap))

    df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
    cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),which = "column",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))

    ht1 = Heatmap(mat, 
      col = color.heatmap.use, 
      na_col = "white", 
      border = TRUE,
      rect_gp = rect_gp,
      name = "Importance",
      bottom_annotation = col_annotation,
      cluster_rows = cluster.rows,
      cluster_columns = cluster.cols,
      row_names_side = "left",
      row_names_rot = 0,
      row_names_gp = gpar(fontsize = font.size),
      column_names_gp = gpar(fontsize = font.size),
      column_title = paste0(names(centr[i]), " signaling pathway network", ifelse(is.null(title), "", paste0(" - ", title))),
      column_title_gp = gpar(fontsize = font.size.title),
      column_names_rot = 45,
      heatmap_legend_param = list(title = "Importance", 
          title_gp = gpar(fontsize = 8, fontface = "plain"),
          title_position = "leftcenter-rot",
          border = NA, 
          at = c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1)),
          legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm")), 
      ...
    )
    return(ht1)
  }
}
netAnalysis_dot <- function(object, slot.name = "netP", pattern = c("outgoing","incoming"), cutoff = NULL, color.use = NULL,
                            pathway.show = NULL, group.show = NULL,
                            shape = 21, dot.size = c(1, 3), dot.alpha = 1, main.title = NULL,
                            font.size = 10, font.size.title = 12,
                            return_data = TRUE, ...){
  pattern <- match.arg(pattern)
  patternSignaling <- methods::slot(object, slot.name)$pattern[[pattern]]
  data1 = patternSignaling$pattern$cell
  data2 = patternSignaling$pattern$signaling
  data = patternSignaling$data
  if (is.null(main.title)) {
    if (pattern == "outgoing") {
      main.title = "Outgoing communication patterns of secreting cells"
    } else if (pattern == "incoming") {
      main.title = "Incoming communication patterns of target cells"
    }
  }
  if (is.null(color.use)) {
    color.use <- scPalette(nlevels(data1$CellGroup))
  }
  if (is.null(cutoff)) {
    cutoff <- 1/length(unique(data1$Pattern))
  }
  options(warn = -1)
  data1$Contribution[data1$Contribution < cutoff] <- 0
  data2$Contribution[data2$Contribution < cutoff] <- 0
  data3 = merge(data1, data2, by.x="Pattern", by.y="Pattern")
  data3$Contribution <- data3$Contribution.x * data3$Contribution.y
  data3 <- data3[,colnames(data3) %in% c("CellGroup","Signaling","Contribution")]
  if (!is.null(pathway.show)) {
    data3 <- data3[data3$Signaling %in% pathway.show, ]
    pathway.add <- pathway.show[which(pathway.show %in% data3$Signaling == 0)]
    if (length(pathway.add) > 1) {
      data.add <- expand.grid(CellGroup = levels(data1$CellGroup), Signaling = pathway.add)
      data.add$Contribution <- 0
      data3 <- rbind(data3, data.add)
    }
    data3$Signaling <- factor(data3$Signaling, levels = pathway.show)
  }
  if (!is.null(group.show)) {
    data3$CellGroup <- as.character(data3$CellGroup)
    data3 <- data3[data3$CellGroup %in% group.show, ]
    data3$CellGroup <- factor(data3$CellGroup, levels = group.show)
  }

  data <- as.data.frame(as.table(data));
  data <- data[data[,3] != 0, ]
  data12 <- paste0(data[,1],data[,2])
  data312 <- paste0(data3[,1],data3[,2])
  idx1 <- which(match(data312, data12, nomatch = 0) ==0)
  data3$Contribution[idx1] <- 0
  data3$id <- data312
  data3 <- data3 %>% group_by(id) %>% top_n(1, Contribution)

  data3$Contribution[which(data3$Contribution == 0)] <- NA

  df <- data3
  gg <- ggplot(data = df, aes(x = Signaling, y = CellGroup)) +
    geom_point(aes(size =  Contribution, fill = CellGroup, colour = CellGroup), shape = shape) +
    scale_size_continuous(range = dot.size) +
    theme_linedraw() +
    scale_x_discrete(position = "bottom") +
    ggtitle(main.title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title, face="plain"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.text.y = element_text(angle = 0, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25)) +
    theme(panel.grid.major = element_line(colour="grey90", size = (0.1)))
  gg <- gg + scale_y_discrete(limits = rev(levels(data3$CellGroup)))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE, na.value = "white")
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE, na.value = "white")
  gg <- gg + guides(colour="none") + guides(fill="none")
  gg <- gg + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8), ...)

  if(return_data){
    return(df)
  }else{
    return(gg)
  }
}

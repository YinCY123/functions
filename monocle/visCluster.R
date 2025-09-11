# Generated from function body. Editing this file has no effect.
visCluster <- function (object = NULL, branches_label = NULL, ht.col.list = list(col_range = c(-2, 
    0, 2), col_color = c("#08519C", "white", "#A50F15")), border = TRUE, 
    plot.type = c("line", "heatmap", "both"), ms.col = c("#0099CC", 
        "grey90", "#CC3333"), line.size = 0.1, line.col = "grey90", 
    add.mline = TRUE, mline.size = 2, mline.col = "#CC3333", 
    ncol = 4, ctAnno.col = NULL, set.md = "median", textbox.pos = c(0.5, 
        0.8), textbox.size = 8, panel.arg = c(2, 0.25, 4, "grey90", 
        NA), ggplot.panel.arg = c(2, 0.25, 4, "grey90", NA), 
    annoTerm.data = NULL, annoTerm.mside = "right", termAnno.arg = c("grey95", 
        "grey50"), add.bar = FALSE, bar.width = 8, textbar.pos = c(0.8, 
        0.8), go.col = NULL, go.size = NULL, by.go = "anno_link", 
    annoKegg.data = NULL, annoKegg.mside = "right", keggAnno.arg = c("grey95", 
        "grey50"), add.kegg.bar = FALSE, kegg.col = NULL, kegg.size = NULL, 
    by.kegg = "anno_link", word_wrap = TRUE, add_new_line = TRUE, 
    add.box = FALSE, boxcol = NULL, box.arg = c(0.1, "grey50"), 
    add.point = FALSE, point.arg = c(19, "orange", "orange", 
        1), add.line = TRUE, line.side = "right", markGenes = NULL, 
    markGenes.side = "right", genes.gp = c("italic", 10, NA), 
    term.text.limit = c(10, 18), mulGroup = NULL, lgd.label = NULL, 
    show_row_names = FALSE, subgroup.anno = NULL, annnoblock.text = TRUE, 
    annnoblock.gp = c("white", 8), add.sampleanno = TRUE, sample.group = NULL, 
    sample.col = NULL, sample.order = NULL, cluster.order = NULL, 
    sample.cell.order = NULL, HeatmapAnnotation = NULL, column.split = NULL, 
    cluster_columns = FALSE, pseudotime_col = NULL, gglist = NULL, 
    row_annotation_obj = NULL, ...) 
{
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        stop("Package 'ComplexHeatmap' is required. Please install it.")
    }
    ComplexHeatmap::ht_opt(message = FALSE)
    if (is.null(ht.col.list[["col_range"]])) {
        col_range = c(-2, 0, 2)
    }
    else {
        col_range = ht.col.list[["col_range"]]
    }
    if (is.null(ht.col.list[["col_color"]])) {
        col_color = c("#08519C", "white", "#A50F15")
    }
    else {
        col_color = ht.col.list[["col_color"]]
    }
    if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required. Please install it.")
    }
    col_fun = circlize::colorRamp2(col_range, col_color)
    plot.type <- match.arg(plot.type)
    if (plot.type == "line") {
        if (object$type %in% c("scRNAdata", "monocle", "wgcna")) {
            data <- data.frame(object$long.res) %>% dplyr::arrange(.data[["cluster"]])
        }
        else {
            if (object$type %in% c("mfuzz", "TCseq")) {
                data <- data.frame(object$long.res) %>% dplyr::arrange(.data[["cluster"]], 
                  .data[["membership"]])
            }
            else {
                data <- data.frame(object$long.res) %>% dplyr::arrange(.data[["cluster"]])
            }
        }
        data$gene <- factor(data$gene, levels = unique(data$gene))
        if (!is.null(sample.order)) {
            data$cell_type <- factor(data$cell_type, levels = sample.order)
        }
        line <- ggplot2::ggplot(data, ggplot2::aes(x = cell_type, 
            y = norm_value))
        if (object$type %in% c("mfuzz", "TCseq")) {
            line <- line + ggplot2::geom_line(ggplot2::aes(color = membership, 
                group = gene), size = line.size) + ggplot2::scale_color_gradient2(low = ms.col[1], 
                mid = ms.col[2], high = ms.col[3], midpoint = 0.5)
        }
        else {
            line <- line + ggplot2::geom_line(ggplot2::aes(group = gene), 
                color = line.col, size = line.size)
        }
        if (add.mline == TRUE) {
            if (object$type == "wgcna") {
                linec <- unique(data$modulecol)
                names(linec) <- linec
                line <- line + ggplot2::geom_line(stat = "summary", 
                  fun = "median", size = mline.size, ggplot2::aes(group = 1, 
                    color = modulecol)) + ggplot2::scale_color_manual(values = linec)
            }
            else {
                line <- line + ggplot2::geom_line(stat = "summary", 
                  fun = "median", colour = mline.col, size = mline.size, 
                  ggplot2::aes(group = 1))
            }
        }
        else {
            line <- line
        }
        line1 <- line + ggplot2::theme_classic(base_size = 14) + 
            ggplot2::ylab("Normalized expression") + ggplot2::xlab("") + 
            ggplot2::theme(axis.ticks.length = ggplot2::unit(0.1, 
                "cm"), axis.text.x = ggplot2::element_text(angle = 45, 
                hjust = 1, color = "black"), strip.background = ggplot2::element_blank()) + 
            ggplot2::facet_wrap(~cluster_name, ncol = ncol, scales = "free")
        return(line1)
    }
    else {
        data <- data.frame(object$wide.res, check.names = FALSE) %>% 
            dplyr::arrange(as.numeric(as.character(cluster)))
        if (object$type %in% c("mfuzz", "TCseq")) {
            mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
                dplyr::select(-gene, -cluster, -membership)
        }
        else if (object$type == "wgcna") {
            mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
                dplyr::select(-gene, -cluster, -modulecol)
        }
        else if (object$type == "scRNAdata") {
            mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
                dplyr::select(-gene, -cluster)
        }
        else if (object$type == "monocle") {
            mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
                dplyr::select(-gene, -cluster)
        }
        else {
            mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
                dplyr::select(-gene, -cluster)
        }
        rownames(mat) <- data$gene
        if (object$geneMode == "all" | ncol(mat) > 20) {
            use_raster = TRUE
        }
        else {
            use_raster = FALSE
        }
        if (!is.null(sample.order)) {
            mat <- mat[, sample.order]
        }
        cl.info <- data.frame(table(data$cluster)) %>% dplyr::mutate(Var1 = as.numeric(as.character(Var1))) %>% 
            dplyr::arrange(Var1)
        cluster.num <- nrow(cl.info)
        subgroup <- lapply(1:nrow(cl.info), function(x) {
            nm <- rep(as.character(cl.info$Var1[x]), cl.info$Freq[x])
            paste("C", nm, sep = "")
        }) %>% unlist()
        if (!is.null(cluster.order)) {
            subgroup <- factor(subgroup, levels = paste("C", 
                cluster.order, sep = ""))
            cluster_row_slices = FALSE
        }
        else {
            cluster_row_slices = TRUE
        }
        if (object$geneMode == "all" & object$type == "scRNAdata") {
            celltype <- sapply(strsplit(colnames(mat), split = "\\|"), 
                "[", 2)
            cell.num.info <- table(celltype)[unique(celltype)]
            if (is.null(sample.cell.order)) {
                column_split = factor(rep(names(cell.num.info), 
                  cell.num.info), levels = unique(celltype))
            }
            else {
                column_split = factor(rep(names(cell.num.info), 
                  cell.num.info), levels = sample.cell.order)
            }
            if (is.null(sample.col)) {
                block.col = 1:length(cell.num.info)
            }
            else {
                block.col = sample.col
            }
        }
        else {
            if (is.null(sample.group)) {
                sample.info = colnames(mat)
                if (is.null(HeatmapAnnotation)) {
                  if (object$geneType == "branched") {
                    if (ncol(mat) == 200) {
                      if(is.null(branches_label)){
                        column_split = rep(c("branch1", "branch2"), 
                          each = 100)
                      }else{
                        column_split <- rep(branches_label, each = 100)
                      }
                    }
                    else {
                      column_split = rep(levels(object$pseudotime), 
                        rep(100, ncol(mat)/100))
                    }
                  }
                  else {
                    column_split = NULL
                  }
                }
                else {
                  column_split = column.split
                }
            }
            else {
                sample.info = sample.group
                column_split = factor(sample.group, levels = unique(sample.group))
            }
            if (object$type != "monocle") {
                sample.info <- factor(sample.info, levels = unique(sample.info))
                if (is.null(sample.col)) {
                  if (!requireNamespace("circlize", quietly = TRUE)) {
                    stop("Package 'circlize' is required. Please install it.")
                  }
                  scol <- circlize::rand_color(n = length(sample.info))
                  names(scol) <- sample.info
                }
                else {
                  scol <- sample.col
                  names(scol) <- sample.info
                }
            }
            else {
                sample.info <- factor(object$pseudotime, levels = unique(object$pseudotime))
                if (is.null(pseudotime_col)) {
                  if (object$geneType == "branched") {
                    if (length(unique(object$pseudotime)) == 
                      3) {
                      pseudotime_col <- c("red", "grey80", "blue")
                    }
                    else {
                      if (!requireNamespace("circlize", quietly = TRUE)) {
                        stop("Package 'circlize' is required. Please install it.")
                      }
                      pseudotime_col <- circlize::rand_color(n = length(unique(object$pseudotime)))
                    }
                  }
                  else {
                    pseudotime_col <- c("blue", "red")
                  }
                }
                else {
                  pseudotime_col <- pseudotime_col
                }
                if (is.null(sample.col)) {
                  if (object$type != "monocle") {
                    if (!requireNamespace("circlize", quietly = TRUE)) {
                      stop("Package 'circlize' is required. Please install it.")
                    }
                    scol <- circlize::rand_color(n = length(sample.info))
                    names(scol) <- sample.info
                  }
                  else {
                    if (object$geneType == "branched") {
                      if (length(unique(object$pseudotime)) == 
                        3) {
                        scol <- rep(pseudotime_col, table(object$pseudotime)[unique(object$pseudotime)])
                        names(scol) <- sample.info
                      }
                      else {
                        scol <- rep(pseudotime_col, table(object$pseudotime)[unique(object$pseudotime)])
                        names(scol) <- sample.info
                      }
                    }
                    else {
                      scol <- (grDevices::colorRampPalette(pseudotime_col))(100)
                      names(scol) <- sample.info
                    }
                  }
                }
                else {
                  scol <- sample.col
                  names(scol) <- sample.info
                }
            }
        }
        if (add.sampleanno == TRUE) {
            if (object$geneMode == "all" & object$type == "scRNAdata") {
                topanno = ComplexHeatmap::HeatmapAnnotation(cluster = ComplexHeatmap::anno_block(gp = grid::gpar(fill = block.col), 
                  labels = NULL), show_annotation_name = FALSE)
            }
            else {
                if (is.null(HeatmapAnnotation)) {
                  topanno = ComplexHeatmap::HeatmapAnnotation(sample = sample.info, 
                    col = list(sample = scol), gp = grid::gpar(col = ifelse(object$type == 
                      "monocle", NA, "white")), show_legend = ifelse(object$type == 
                      "monocle", FALSE, TRUE), show_annotation_name = FALSE)
                }
                else {
                  topanno = HeatmapAnnotation
                }
            }
        }
        else {
            topanno = NULL
        }
        if (is.null(ctAnno.col)) {
            if (!requireNamespace("circlize", quietly = TRUE)) {
                stop("Package 'circlize' is required. Please install it.")
            }
            colanno <- circlize::rand_color(n = cluster.num)
        }
        else {
            colanno <- ctAnno.col
        }
        names(colanno) <- 1:cluster.num
        align_to = split(1:nrow(mat), subgroup)
        anno.block <- ComplexHeatmap::anno_block(align_to = align_to, 
            panel_fun = function(index, nm) {
                npos = as.numeric(unlist(strsplit(nm, split = "C"))[2])
                grid::grid.rect(gp = grid::gpar(fill = colanno[npos], 
                  col = NA))
                if (annnoblock.text == TRUE) {
                  grid::grid.text(label = paste("n:", length(index), 
                    sep = ""), rot = 90, gp = grid::gpar(col = annnoblock.gp[1], 
                    fontsize = as.numeric(annnoblock.gp[2])))
                }
            }, which = "row")
        if (!is.null(markGenes)) {
            rowGene <- rownames(mat)
            annoGene <- markGenes
            gene.col <- data %>% dplyr::select(gene, cluster) %>% 
                dplyr::filter(gene %in% annoGene)
            gene.col <- purrr::map_df(1:cluster.num, function(x) {
                tmp <- gene.col %>% dplyr::filter(as.numeric(cluster) == 
                  x) %>% dplyr::mutate(col = colanno[x])
            })
            gene.col <- gene.col[match(annoGene, gene.col$gene), 
                ]
            if (is.na(genes.gp[3])) {
                gcol = gene.col$col
            }
            else {
                gcol = genes.gp[3]
            }
            index <- match(annoGene, rowGene)
            geneMark = ComplexHeatmap::anno_mark(at = index, 
                labels = annoGene, which = "row", side = markGenes.side, 
                labels_gp = grid::gpar(fontface = genes.gp[1], 
                  fontsize = as.numeric(genes.gp[2]), col = gcol))
        }
        else {
            geneMark = NULL
        }
        right_annotation = ComplexHeatmap::rowAnnotation(gene = geneMark, 
            cluster = anno.block)
        if (object$type == "monocle" | object$geneMode == "all" | 
            ncol(mat) > 20) {
            show_column_names = FALSE
        }
        else {
            show_column_names = TRUE
        }
        if (object$geneType == "non-branched") {
            rg <- range(as.numeric(as.character(sample.info)))
            if (!requireNamespace("circlize", quietly = TRUE)) {
                stop("Package 'circlize' is required. Please install it.")
            }
            col_fun2 = circlize::colorRamp2(c(rg[1], rg[2]), 
                pseudotime_col)
            lgd = ComplexHeatmap::Legend(col_fun = col_fun2, 
                title = "pseudotime")
            lgd_list = list(lgd)
        }
        else if (object$geneType == "branched") {
            if (length(levels(sample.info)) == 3) {
                lgd = ComplexHeatmap::Legend(labels = levels(sample.info), 
                  legend_gp = grid::gpar(fill = pseudotime_col), 
                  title = "branch")
            }
            else {
                lgd = ComplexHeatmap::Legend(labels = levels(sample.info), 
                  legend_gp = grid::gpar(fill = pseudotime_col), 
                  title = "branch")
            }
            lgd_list = list(lgd)
        }
        else {
            lgd_list <- NULL
        }
        if (plot.type == "heatmap") {
            if (!is.null(row_annotation_obj)) {
                left_annotation_ht <- row_annotation_obj
            }
            else {
                left_annotation_ht <- NULL
            }
            htf <- ComplexHeatmap::Heatmap(as.matrix(mat), name = "Z-score", 
                cluster_columns = cluster_columns, show_row_names = show_row_names, 
                border = border, column_split = column_split, 
                row_split = subgroup, cluster_row_slices = cluster_row_slices, 
                column_names_side = "top", show_column_names = show_column_names, 
                top_annotation = topanno, left_annotation = left_annotation_ht, 
                right_annotation = right_annotation, col = col_fun, 
                use_raster = use_raster, ...)
            ComplexHeatmap::draw(htf, merge_legend = TRUE, annotation_legend_list = lgd_list)
        }
        else {
            rg = range(mat)
            if (!is.null(gglist)) {
                anno_ggplot2 = ComplexHeatmap::anno_zoom(align_to = align_to, 
                  which = "row", panel_fun = function(index, 
                    nm) {
                    g <- gglist[[nm]]
                    g <- grid::grid.grabExpr(print(g))
                    grid::pushViewport(grid::viewport())
                    grid::grid.rect()
                    grid::grid.draw(g)
                    grid::popViewport()
                  }, size = grid::unit(as.numeric(ggplot.panel.arg[1]), 
                    "cm"), gap = grid::unit(as.numeric(ggplot.panel.arg[2]), 
                    "cm"), width = grid::unit(as.numeric(ggplot.panel.arg[3]), 
                    "cm"), side = "right", link_gp = grid::gpar(fill = ggplot.panel.arg[4], 
                    col = ggplot.panel.arg[5]))
            }
            else {
                anno_ggplot2 = NULL
            }
            panel_fun = function(index, nm) {
                if (add.box == TRUE & add.line != TRUE) {
                  xscale = c(-0.1, 1.1)
                }
                else {
                  xscale = c(-0.1, 1.1)
                  panel_scale <- c(0.1, 0.9)
                }
                grid::pushViewport(grid::viewport(xscale = xscale, 
                  yscale = c(0, 1)))
                grid::grid.rect()
                if (object$geneMode == "all" & object$type == 
                  "scRNAdata") {
                  mulGroup <- cell.num.info
                  grid::grid.lines(x = c(0, 1), y = rep(0.5, 
                    2), gp = grid::gpar(col = "black", lty = "dashed"))
                  cu <- cumsum(mulGroup)
                  seqn <- data.frame(st = c(1, cu[1:(length(cu) - 
                    1)] + 1), sp = c(cu[1], cu[2:length(cu)]))
                }
                else {
                  if (is.null(mulGroup)) {
                    mulGroup <- ncol(mat)
                    seqn <- data.frame(st = 1, sp = ncol(mat))
                  }
                  else {
                    mulGroup <- mulGroup
                    grid::grid.lines(x = c(0, 1), y = rep(0.5, 
                      2), gp = grid::gpar(col = "black", lty = "dashed"))
                    cu <- cumsum(mulGroup)
                    seqn <- data.frame(st = c(1, cu[1:(length(cu) - 
                      1)] + 1), sp = c(cu[1], cu[2:length(cu)]))
                  }
                }
                if (object$geneMode == "all" & object$type == 
                  "scRNAdata") {
                  cell.ave <- purrr::map_dfr(1:nrow(seqn), function(x) {
                    tmp <- seqn[x, ]
                    tmpmat <- mat[index, c(tmp$st:tmp$sp)]
                    rg <- range(mat[index, ])
                    if (set.md == "mean") {
                      mdia <- mean(rowMeans(tmpmat))
                    }
                    else if (set.md == "median") {
                      mdia <- stats::median(apply(tmpmat, 1, 
                        stats::median))
                    }
                    else {
                      message("supply mean/median !")
                    }
                    res <- data.frame(x = x, val = mdia)
                    return(res)
                  })
                  if (add.line == TRUE) {
                    grid::grid.lines(x = scales::rescale(cell.ave$x, 
                      to = c(0, 1)), y = scales::rescale(cell.ave$val, 
                      to = c(0.1, 0.9)), gp = grid::gpar(lwd = 3, 
                      col = mline.col))
                  }
                }
                else {
                  lapply(1:nrow(seqn), function(x) {
                    tmp <- seqn[x, ]
                    tmpmat <- mat[index, c(tmp$st:tmp$sp)]
                    if (set.md == "mean") {
                      mdia <- colMeans(tmpmat)
                    }
                    else if (set.md == "median") {
                      mdia <- apply(tmpmat, 2, stats::median)
                    }
                    else {
                      message("supply mean/median !")
                    }
                    pos = scales::rescale(1:ncol(tmpmat), to = c(0, 
                      1))
                    if (is.null(boxcol)) {
                      boxcol <- rep("grey90", ncol(tmpmat))
                    }
                    else {
                      boxcol <- boxcol
                    }
                    if (add.box == TRUE) {
                      lapply(1:ncol(tmpmat), function(x) {
                        ComplexHeatmap::grid.boxplot(scales::rescale(tmpmat[, 
                          x], to = c(0, 1), from = c(rg[1] - 
                          0.5, rg[2] + 0.5)), pos = pos[x], direction = "vertical", 
                          box_width = as.numeric(box.arg[1]), 
                          outline = FALSE, gp = grid::gpar(col = box.arg[2], 
                            fill = boxcol[x]))
                      })
                    }
                    if (add.point == TRUE) {
                      grid::grid.points(x = scales::rescale(1:ncol(tmpmat), 
                        to = panel_scale), y = scales::rescale(mdia, 
                        to = c(0, 1), from = c(rg[1] - 0.5, rg[2] + 
                          0.5)), pch = as.numeric(point.arg[1]), 
                        gp = grid::gpar(fill = point.arg[2], 
                          col = point.arg[3]), size = grid::unit(as.numeric(point.arg[4]), 
                          "char"))
                    }
                    if (add.line == TRUE) {
                      grid::grid.lines(x = scales::rescale(1:ncol(tmpmat), 
                        to = panel_scale), y = scales::rescale(mdia, 
                        to = c(0, 1), from = c(rg[1] - 0.5, rg[2] + 
                          0.5)), gp = grid::gpar(lwd = 3, col = mline.col[x]))
                    }
                  })
                }
                grid.textbox <- utils::getFromNamespace("grid.textbox", 
                  "ComplexHeatmap")
                text <- paste("Gene Size:", nrow(mat[index, ]), 
                  sep = " ")
                grid.textbox(text, x = textbox.pos[1], y = textbox.pos[2], 
                  gp = grid::gpar(fontsize = textbox.size, fontface = "italic", 
                    ...))
                grid::popViewport()
            }
            if (!is.null(subgroup.anno)) {
                align_to = split(1:nrow(mat), subgroup)
                align_to = align_to[subgroup.anno]
            }
            else {
                align_to = subgroup
            }
            anno = ComplexHeatmap::anno_link(align_to = align_to, 
                which = "row", panel_fun = panel_fun, size = grid::unit(as.numeric(panel.arg[1]), 
                  "cm"), gap = grid::unit(as.numeric(panel.arg[2]), 
                  "cm"), width = grid::unit(as.numeric(panel.arg[3]), 
                  "cm"), side = line.side, link_gp = grid::gpar(fill = panel.arg[4], 
                  col = panel.arg[5]))
            if (!is.null(annoTerm.data)) {
                termanno <- annoTerm.data
                if (ncol(termanno) == 2) {
                  colnames(termanno) <- c("id", "term")
                }
                else if (ncol(termanno) == 3) {
                  colnames(termanno) <- c("id", "term", "pval")
                }
                else if (ncol(termanno) == 4) {
                  colnames(termanno) <- c("id", "term", "pval", 
                    "ratio")
                }
                else {
                  message("No more than 4 columns!")
                }
                if (is.null(go.col)) {
                  if (!requireNamespace("circlize", quietly = TRUE)) {
                    stop("Package 'circlize' is required. Please install it.")
                  }
                  gocol <- circlize::rand_color(n = nrow(termanno))
                }
                else {
                  gocol <- go.col
                }
                if (is.null(go.size)) {
                  gosize <- rep(12, nrow(termanno))
                }
                else {
                  if (go.size == "pval") {
                    termanno.tmp <- purrr::map_df(unique(termanno$id), 
                      function(x) {
                        tmp <- termanno %>% dplyr::filter(id == 
                          x) %>% dplyr::mutate(size = scales::rescale(-log10(pval), 
                          to = term.text.limit))
                      })
                    gosize <- termanno.tmp$size
                  }
                  else {
                    gosize <- go.size
                  }
                }
                termanno <- termanno %>% dplyr::ungroup() %>% 
                  dplyr::mutate(col = gocol, fontsize = gosize)
                term.list <- lapply(1:length(unique(termanno$id)), 
                  function(x) {
                    tmp = termanno[which(termanno$id == unique(termanno$id)[x]), 
                      ]
                    df <- data.frame(text = tmp$term, col = tmp$col, 
                      fontsize = tmp$fontsize)
                    return(df)
                  })
                names(term.list) <- unique(termanno$id)
                if (!is.null(subgroup.anno)) {
                  align_to2 = split(seq_along(subgroup), subgroup)
                  align_to2 = align_to2[subgroup.anno]
                  term.list = term.list[subgroup.anno]
                }
                else {
                  align_to2 = subgroup
                  term.list = term.list
                }
                textbox = ComplexHeatmap::anno_textbox(align_to2, 
                  term.list, word_wrap = word_wrap, add_new_line = add_new_line, 
                  side = annoTerm.mside, background_gp = grid::gpar(fill = termAnno.arg[1], 
                    col = termAnno.arg[2]), by = by.go)
                if (ncol(termanno) - 2 > 2) {
                  anno_gobar <- function(data = NULL, bar.width = 0.1, 
                    align_to = NULL, panel.arg = panel.arg, ...) {
                    if (ncol(data) - 2 == 3) {
                      data <- data %>% dplyr::mutate(bary = -log10(pval))
                    }
                    else {
                      data <- data %>% dplyr::mutate(bary = ratio)
                    }
                    ComplexHeatmap::anno_zoom(align_to = align_to, 
                      which = "row", panel_fun = function(index, 
                        nm) {
                        grid::pushViewport(grid::viewport(xscale = c(0, 
                          1), yscale = c(0, 1)))
                        grid::grid.rect()
                        tmp <- data %>% dplyr::filter(id == nm)
                        grid::grid.segments(x0 = rep(0, nrow(tmp)), 
                          x1 = scales::rescale(rev(tmp$bary), 
                            to = c(0.1, 0.9)), y0 = scales::rescale(1:nrow(tmp), 
                            to = c(0.1, 0.9)), y1 = scales::rescale(1:nrow(tmp), 
                            to = c(0.1, 0.9)), gp = grid::gpar(lwd = bar.width, 
                            col = rev(tmp$col), lineend = "butt"))
                        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                          "ComplexHeatmap")
                        text <- nm
                        grid.textbox(text, x = textbar.pos[1], 
                          y = textbar.pos[2], gp = grid::gpar(fontsize = textbox.size, 
                            fontface = "italic", col = unique(tmp$col), 
                            ...))
                        grid::popViewport()
                      }, size = grid::unit(as.numeric(panel.arg[1]), 
                        "cm"), gap = grid::unit(as.numeric(panel.arg[2]), 
                        "cm"), width = grid::unit(as.numeric(panel.arg[3]), 
                        "cm"), side = "right", link_gp = grid::gpar(fill = termAnno.arg[1], 
                        col = termAnno.arg[2]), ...)
                  }
                  baranno = anno_gobar(data = termanno, align_to = align_to2, 
                    panel.arg = panel.arg, bar.width = bar.width)
                }
                if (add.bar == TRUE) {
                  baranno
                }
                else {
                  baranno = NULL
                }
            }
            else {
                textbox = NULL
                baranno = NULL
            }
            if (!is.null(annoKegg.data)) {
                termanno <- annoKegg.data
                if (ncol(termanno) == 2) {
                  colnames(termanno) <- c("id", "term")
                }
                else if (ncol(termanno) == 3) {
                  colnames(termanno) <- c("id", "term", "pval")
                }
                else if (ncol(termanno) == 4) {
                  colnames(termanno) <- c("id", "term", "pval", 
                    "ratio")
                }
                else {
                  message("No more than 4 columns!")
                }
                if (is.null(kegg.col)) {
                  if (!requireNamespace("circlize", quietly = TRUE)) {
                    stop("Package 'circlize' is required. Please install it.")
                  }
                  gocol <- circlize::rand_color(n = nrow(termanno))
                }
                else {
                  gocol <- kegg.col
                }
                if (is.null(kegg.size)) {
                  gosize <- rep(12, nrow(termanno))
                }
                else {
                  if (kegg.size == "pval") {
                    termanno.tmp <- purrr::map_df(unique(termanno$id), 
                      function(x) {
                        tmp <- termanno %>% dplyr::filter(id == 
                          x) %>% dplyr::mutate(size = scales::rescale(-log10(pval), 
                          to = term.text.limit))
                      })
                    gosize <- termanno.tmp$size
                  }
                  else {
                    gosize <- kegg.size
                  }
                }
                termanno <- termanno %>% dplyr::ungroup() %>% 
                  dplyr::mutate(col = gocol, fontsize = gosize)
                term.list <- lapply(1:length(unique(termanno$id)), 
                  function(x) {
                    tmp = termanno[which(termanno$id == unique(termanno$id)[x]), 
                      ]
                    df <- data.frame(text = tmp$term, col = tmp$col, 
                      fontsize = tmp$fontsize)
                    return(df)
                  })
                names(term.list) <- unique(termanno$id)
                if (!is.null(subgroup.anno)) {
                  align_to2 = split(seq_along(subgroup), subgroup)
                  align_to2 = align_to2[subgroup.anno]
                  term.list = term.list[subgroup.anno]
                }
                else {
                  align_to2 = subgroup
                  term.list = term.list
                }
                textbox.kegg = ComplexHeatmap::anno_textbox(align_to2, 
                  term.list, word_wrap = word_wrap, add_new_line = add_new_line, 
                  side = annoKegg.mside, background_gp = grid::gpar(fill = keggAnno.arg[1], 
                    col = keggAnno.arg[2]), by = by.kegg)
                if (ncol(termanno) - 2 > 2) {
                  anno_keggbar <- function(data = NULL, bar.width = 0.1, 
                    align_to = NULL, panel.arg = panel.arg, ...) {
                    if (ncol(data) - 2 == 3) {
                      data <- data %>% dplyr::mutate(bary = -log10(pval))
                    }
                    else {
                      data <- data %>% dplyr::mutate(bary = ratio)
                    }
                    ComplexHeatmap::anno_zoom(align_to = align_to, 
                      which = "row", panel_fun = function(index, 
                        nm) {
                        grid::pushViewport(grid::viewport(xscale = c(0, 
                          1), yscale = c(0, 1)))
                        grid::grid.rect()
                        tmp <- data %>% dplyr::filter(id == nm)
                        grid::grid.segments(x0 = rep(0, nrow(tmp)), 
                          x1 = scales::rescale(rev(tmp$bary), 
                            to = c(0.1, 0.9)), y0 = scales::rescale(1:nrow(tmp), 
                            to = c(0.1, 0.9)), y1 = scales::rescale(1:nrow(tmp), 
                            to = c(0.1, 0.9)), gp = grid::gpar(lwd = bar.width, 
                            col = rev(tmp$col), lineend = "butt"))
                        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                          "ComplexHeatmap")
                        text <- nm
                        grid.textbox(text, x = textbar.pos[1], 
                          y = textbar.pos[2], gp = grid::gpar(fontsize = textbox.size, 
                            fontface = "italic", col = unique(tmp$col), 
                            ...))
                        grid::popViewport()
                      }, size = grid::unit(as.numeric(panel.arg[1]), 
                        "cm"), gap = grid::unit(as.numeric(panel.arg[2]), 
                        "cm"), width = grid::unit(as.numeric(panel.arg[3]), 
                        "cm"), side = "right", link_gp = grid::gpar(fill = keggAnno.arg[1], 
                        col = keggAnno.arg[2]), ...)
                  }
                  baranno.kegg = anno_keggbar(data = termanno, 
                    align_to = align_to2, panel.arg = panel.arg, 
                    bar.width = bar.width)
                }
                if (add.kegg.bar == TRUE) {
                  baranno.kegg
                }
                else {
                  baranno.kegg = NULL
                }
            }
            else {
                textbox.kegg = NULL
                baranno.kegg = NULL
            }
            if (line.side == "right") {
                if (markGenes.side == "right") {
                  right_annotation2 = ComplexHeatmap::rowAnnotation(gene = geneMark, 
                    cluster = anno.block, line = anno, anno_ggplot2 = anno_ggplot2, 
                    textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                    baranno.kegg = baranno.kegg)
                  if (!is.null(row_annotation_obj)) {
                    left_annotation <- row_annotation_obj
                  }
                  else {
                    left_annotation <- NULL
                  }
                }
                else {
                  right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block, 
                    line = anno, anno_ggplot2 = anno_ggplot2, 
                    textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                    baranno.kegg = baranno.kegg)
                  if (!is.null(row_annotation_obj)) {
                    left_annotation <- do.call(ComplexHeatmap::rowAnnotation, 
                      modifyList(list(gene = geneMark), row_annotation_obj))
                  }
                  else {
                    left_annotation <- ComplexHeatmap::rowAnnotation(gene = geneMark)
                  }
                }
            }
            else {
                if (markGenes.side == "right") {
                  right_annotation2 = ComplexHeatmap::rowAnnotation(gene = geneMark, 
                    cluster = anno.block, anno_ggplot2 = anno_ggplot2, 
                    textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                    baranno.kegg = baranno.kegg)
                  if (!is.null(row_annotation_obj)) {
                    left_annotation <- do.call(ComplexHeatmap::rowAnnotation, 
                      modifyList(list(line = anno), row_annotation_obj))
                  }
                  else {
                    left_annotation <- ComplexHeatmap::rowAnnotation(line = anno)
                  }
                }
                else {
                  right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block, 
                    anno_ggplot2 = anno_ggplot2, textbox = textbox, 
                    bar = baranno, textbox.kegg = textbox.kegg, 
                    baranno.kegg = baranno.kegg)
                  if (!is.null(row_annotation_obj)) {
                    left_annotation <- do.call(ComplexHeatmap::rowAnnotation, 
                      modifyList(list(gene = geneMark, line = anno), 
                        row_annotation_obj))
                  }
                  else {
                    left_annotation <- ComplexHeatmap::rowAnnotation(line = anno, 
                      gene = geneMark)
                  }
                }
            }
            if (object$type == "monocle" | object$geneMode == 
                "all" | ncol(mat) > 20) {
                show_column_names = FALSE
            }
            else {
                show_column_names = TRUE
            }
            htf <- ComplexHeatmap::Heatmap(as.matrix(mat), name = "Z-score", 
                cluster_columns = cluster_columns, show_row_names = show_row_names, 
                border = border, column_split = column_split, 
                top_annotation = topanno, right_annotation = right_annotation2, 
                left_annotation = left_annotation, column_names_side = "top", 
                show_column_names = show_column_names, row_split = subgroup, 
                cluster_row_slices = cluster_row_slices, col = col_fun, 
                use_raster = use_raster, ...)
            if (is.null(mulGroup)) {
                ComplexHeatmap::draw(htf, merge_legend = TRUE, 
                  annotation_legend_list = lgd_list)
            }
            else {
                if (is.null(lgd.label)) {
                  lgd.label <- paste("group", 1:length(mulGroup), 
                    sep = "")
                }
                else {
                  lgd.label <- lgd.label
                }
                lgd_list2 = ComplexHeatmap::Legend(labels = lgd.label, 
                  type = "lines", legend_gp = grid::gpar(col = mline.col, 
                    lty = 1))
                if (!is.null(lgd_list)) {
                  lgd_list_com <- ComplexHeatmap::packLegend(lgd_list, 
                    lgd_list2)
                }
                else {
                  lgd_list_com <- lgd_list2
                }
                ComplexHeatmap::draw(htf, annotation_legend_list = lgd_list_com, 
                  merge_legend = TRUE)
            }
        }
    }
}

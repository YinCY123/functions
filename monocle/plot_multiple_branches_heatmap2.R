# Generated from function body. Editing this file has no effect.
plot_multiple_branches_heatmap2 <- function (cds = NULL, branches, branches_name = NULL, cluster_rows = TRUE, 
    hclust_method = "ward.D2", num_clusters = 6, hmcols = NULL, 
    add_annotation_row = NULL, add_annotation_col = NULL, show_rownames = FALSE, 
    use_gene_short_name = TRUE, norm_method = c("vstExprs", "log"), 
    scale_max = 3, scale_min = -3, trend_formula = "~sm.ns(Pseudotime, df=3)", 
    return_heatmap = FALSE, cores = 1) 
{
    pseudocount <- 1
    if (!(all(branches %in% Biobase::pData(cds)$State)) & length(branches) == 
        1) {
        stop("This function only allows to make multiple branch plots where branches is included in the pData")
    }
    branch_label <- branches
    if (!is.null(branches_name)) {
        if (length(branches) != length(branches_name)) {
            stop("branches_name should have the same length as branches")
        }
        branch_label <- branches_name
    }
    g <- cds@minSpanningTree
    m <- NULL
    for (branch_in in branches) {
        branches_cells <- row.names(subset(Biobase::pData(cds), 
            State == branch_in))
        root_state <- subset(Biobase::pData(cds), Pseudotime == 
            0)[, "State"]
        root_state_cells <- row.names(subset(Biobase::pData(cds), 
            State == root_state))
        if (cds@dim_reduce_type != "ICA") {
            root_state_cells <- unique(paste("Y_", cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_state_cells, 
                ], sep = ""))
            branches_cells <- unique(paste("Y_", cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[branches_cells, 
                ], sep = ""))
        }
        deg_root <- igraph::degree(g, v = root_state_cells)
        root_cell_all <- root_state_cells[which(deg_root == 1)]
        if (length(root_cell_all) == 0) {
            # fallback: use root state cells with minimal degree
            min_deg_root <- min(deg_root, na.rm = TRUE)
            root_cell_all <- root_state_cells[which(deg_root == min_deg_root)]
            warning("No root cell (degree==1) found for the provided root state; using cell(s) with minimal degree as candidate.")
        }
        deg_tip <- igraph::degree(g, v = branches_cells)
        tip_cell_all <- branches_cells[which(deg_tip == 1)]
        if (length(tip_cell_all) == 0) {
            # fallback: use branch cells with minimal degree (closest to leaves)
            min_deg_tip <- min(deg_tip, na.rm = TRUE)
            tip_cell_all <- branches_cells[which(deg_tip == min_deg_tip)]
            warning("No tip cell (degree==1) found for the provided branch ", branch_in, "; using cell(s) with minimal degree as candidate.")
        }
        # If multiple candidate root/tip cells exist, pick the first one (safe default).
        root_cell <- root_cell_all[1]
        tip_cell <- tip_cell_all[1]
        # Try to use traverseTree; if it fails or returns no path, fall back to igraph shortest path.
        traverse_res <- tryCatch(traverseTree(g, root_cell, tip_cell), error = function(e) NULL)
        path_cells <- NULL
        if (!is.null(traverse_res) && !is.null(traverse_res$shortest_path) && length(traverse_res$shortest_path) >= 1 && !is.null(traverse_res$shortest_path[[1]])) {
            path_cells <- names(traverse_res$shortest_path[[1]])
        } else {
            sp <- tryCatch(igraph::shortest_paths(g, from = root_cell, to = tip_cell, output = "vpath"), error = function(e) NULL)
            if (!is.null(sp) && length(sp$vpath) >= 1 && length(sp$vpath[[1]]) > 0) {
                # igraph::V(g)$name gives vertex names corresponding to the vpath indices
                vnames <- igraph::V(g)$name
                path_cells <- vnames[sp$vpath[[1]]]
            } else {
                stop("No path found between root and tip cells for branch: ", branch_in)
            }
        }
        if (cds@dim_reduce_type != "ICA") {
            pc_ind <- cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex
            path_cells <- row.names(pc_ind)[paste("Y_", pc_ind[, 
                1], sep = "") %in% path_cells]
        }
        cds_subset <- cds[, path_cells]
        newdata <- data.frame(Pseudotime = seq(0, max(Biobase::pData(cds_subset)$Pseudotime), 
            length.out = 100))
        if (requireNamespace("monocle", quietly = TRUE)) {
            tmp <- monocle::genSmoothCurves(cds_subset, cores = cores, 
                trend_formula = trend_formula, relative_expr = T, 
                new_data = newdata)
        }
        else {
            warning("Cannot find monocle 'monocle' is not installed.")
        }
        if (is.null(m)) 
            m <- tmp
        else {
            m <- cbind(m, tmp)
        }
    }
    m = m[!apply(m, 1, sum) == 0, ]
    norm_method <- match.arg(norm_method)
    if (requireNamespace("monocle", quietly = TRUE)) {
        if (norm_method == "vstExprs" && is.null(cds@dispFitInfo[["blind"]]$disp_func) == 
            FALSE) {
            m = monocle::vstExprs(cds, expr_matrix = m)
        }
        else if (norm_method == "log") {
            m = log10(m + pseudocount)
        }
    }
    else {
        warning("Cannot find monocle 'monocle' is not installed.")
    }
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale_max] = scale_max
    m[m < scale_min] = scale_min
    heatmap_matrix <- m
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    if (is.null(hmcols)) {
        bks <- seq(-3.1, 3.1, by = 0.1)
        hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    }
    else {
        bks <- seq(-3.1, 3.1, length.out = length(hmcols))
    }
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        ph <- pheatmap::pheatmap(heatmap_matrix, useRaster = T, 
            cluster_cols = FALSE, cluster_rows = T, show_rownames = F, 
            show_colnames = F, clustering_distance_rows = row_dist, 
            clustering_method = hclust_method, cutree_rows = num_clusters, 
            silent = TRUE, filename = NA, breaks = bks, color = hmcols)
    }
    else {
        warning("Cannot create heatmap. 'pheatmap' is not installed.")
    }
    annotation_col <- data.frame(Branch = factor(rep(rep(branch_label, 
        each = 100))))
    annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
        num_clusters)))
    col_gaps_ind <- c(1:(length(branches) - 1)) * 100
    if (!is.null(add_annotation_row)) {
        old_colnames_length <- ncol(annotation_row)
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
            ])
        colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
    }
    if (use_gene_short_name == TRUE) {
        if (is.null(Biobase::fData(cds)$gene_short_name) == FALSE) {
            feature_label <- as.character(Biobase::fData(cds)[row.names(heatmap_matrix), 
                "gene_short_name"])
            feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
            row_ann_labels <- as.character(Biobase::fData(cds)[row.names(annotation_row), 
                "gene_short_name"])
            row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
        }
        else {
            feature_label <- row.names(heatmap_matrix)
            row_ann_labels <- row.names(annotation_row)
        }
    }
    else {
        feature_label <- row.names(heatmap_matrix)
        row_ann_labels <- row.names(annotation_row)
    }
    row.names(heatmap_matrix) <- feature_label
    row.names(annotation_row) <- row_ann_labels
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    if (!(cluster_rows)) {
        annotation_row <- NA
    }
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        ph_res <- pheatmap::pheatmap(heatmap_matrix[, ], useRaster = T, 
            cluster_cols = FALSE, cluster_rows = cluster_rows, 
            show_rownames = show_rownames, show_colnames = F, 
            clustering_distance_rows = row_dist, clustering_method = hclust_method, 
            cutree_rows = num_clusters, annotation_row = annotation_row, 
            annotation_col = annotation_col, gaps_col = col_gaps_ind, 
            treeheight_row = 20, breaks = bks, fontsize = 12, 
            color = hmcols, silent = TRUE, border_color = NA, 
            filename = NA)
    }
    else {
        warning("Cannot create heatmap. 'pheatmap' is not installed.")
    }
    wide.res <- cbind(heatmap_matrix, annotation_row) %>% data.frame(., 
        check.names = FALSE) %>% dplyr::mutate(gene = rownames(.), 
        .before = 1) %>% dplyr::rename(cluster = Cluster)
    df <- reshape2::melt(wide.res, id.vars = c("cluster", "gene"), 
        variable.name = "cell_type", value.name = "norm_value") %>% 
        dplyr::mutate(cell_type = as.numeric(as.character(cell_type)))
    df$cluster_name <- paste("cluster ", df$cluster, sep = "")
    cl.info <- data.frame(table(wide.res$cluster)) %>% dplyr::mutate(Var1 = as.numeric(as.character(Var1))) %>% 
        dplyr::arrange(Var1)
    id <- unique(df$cluster_name)
    df <- purrr::map_df(seq_along(id), function(x) {
        tmp <- df %>% dplyr::filter(cluster_name == id[x])
        tmp %>% dplyr::mutate(cluster_name = paste(cluster_name, 
            " (", cl.info$Freq[x], ")", sep = ""))
    })
    df$cluster_name <- factor(df$cluster_name, levels = paste("cluster ", 
        cl.info$Var1, " (", cl.info$Freq, ")", sep = ""))
    prepared_data <- list(wide.res = wide.res, long.res = df, 
        type = "monocle", geneMode = "all", geneType = "branched", 
        pseudotime = annotation_col$Branch)
    if (return_heatmap == TRUE) {
        grid::grid.rect(gp = grid::gpar("fill", col = NA))
        grid::grid.draw(ph_res$gtable)
        return(ph_res)
    }
    else {
        return(prepared_data)
    }
}

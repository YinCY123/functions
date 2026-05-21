# ========================================== 1. 参数配置 ==========================================
suppressPackageStartupMessages({
  library(Seurat); library(monocle); library(tidyverse); library(ggridges); library(igraph)
})

seurat_path <- "./newNC_combined.RData"
out_rds     <- "processed_subset.rds" 
out_pdf     <- "Final_Publication_Plots.pdf"
n_cells     <- 20000        
n_genes     <- 150          
start_cl    <- "1-NF"       
custom_cols <- c("NAG"="#4daf4a", "Precancerosis"="#377eb8", "Tumor"="#e41a1c") 

# ============ 2. 数据载入与缓存 ================
if (1) {
  load(seurat_path); seu <- newNC.combined; rm(newNC.combined); gc()
  seu$Disease_Stage <- factor(case_when(
    seu$sample_all %in% "NAG" ~ "NAG",
    seu$sample_all %in% c("CAG", "IM", "Precancerosis") ~ "Precancerosis",
    TRUE ~ "Tumor"
  ), levels = c("NAG", "Precancerosis", "Tumor"))
  set.seed(12345); seu_sub <- subset(seu, cells = sample(Cells(seu), min(n_cells, ncol(seu))))
  Idents(seu_sub) <- "cluster_ordered"
  markers <- FindAllMarkers(seu_sub, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5, max.cells.per.ident = 100)
  top_genes <- markers %>% 
    group_by(cluster) %>% 
    top_n(n = ceiling(n_genes), wt = avg_log2FC) %>% 
    pull(gene) %>% 
    unique()
  saveRDS(list(seu = seu_sub, genes = top_genes), out_rds)
} else {
  cache <- readRDS(out_rds); seu_sub <- cache$seu; top_genes <- cache$genes
}

# ========== 3. Monocle 轨迹构建 ===========
if (inherits(seu_sub[["RNA"]], "Assay5")) seu_sub <- JoinLayers(seu_sub)
norm_data <- as(GetAssayData(seu_sub, layer = "data"), "dgCMatrix")

fd_df <- data.frame(gene_short_name = rownames(norm_data), row.names = rownames(norm_data), stringsAsFactors = FALSE)
fd_df$num_cells_expressed <- Matrix::rowSums(norm_data > 0)
pd <- new("AnnotatedDataFrame", data = seu_sub@meta.data[colnames(norm_data), ])
fd <- new("AnnotatedDataFrame", data = fd_df)

cds <- newCellDataSet(norm_data, phenoData = pd, featureData = fd, expressionFamily = VGAM::uninormal())
ordering_genes <- intersect(top_genes, rownames(norm_data))[1:min(n_genes, length(top_genes))]

# 降维与排序
cds <- setOrderingFilter(cds, ordering_genes) %>% 
       reduceDimension(max_components = 2, method = 'DDRTree', norm_method = 'none') %>% 
       orderCells()

root_state <- names(which.max(table(pData(cds)$State, pData(cds)$cluster_ordered)[, start_cl]))
cds <- orderCells(cds, root_state = root_state)

# ️ 释放内存，防止 system call failed: Cannot allocate memory
rm(seu_sub, norm_data, pd, fd, fd_df, cache, markers); invisible(gc())

# ============ 4. 提取数据自己画 (彻底抛弃老旧画图函数) ==============
# 提取细胞点坐标
plot_df <- pData(cds)
cell_coords <- t(reducedDimS(cds))
plot_df$x <- cell_coords[, 1]
plot_df$y <- cell_coords[, 2]

# 提取骨架树 (MST) 连线坐标
ica_space_df <- as.data.frame(t(reducedDimK(cds)))
edge_list <- as.data.frame(igraph::get.edgelist(minSpanningTree(cds)))
mst_df <- data.frame(
  x = ica_space_df[edge_list$V1, 1], 
  y = ica_space_df[edge_list$V1, 2],
  xend = ica_space_df[edge_list$V2, 1], 
  yend = ica_space_df[edge_list$V2, 2]
)

# 顶刊主题
pub_theme <- theme_classic(base_size = 13) + theme(
  plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
  axis.title = element_text(face = "bold"), legend.position = "right",
  strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12)
)

# 构建底层模板 (黑线骨架 + 点)
base_p <- ggplot() + 
  geom_segment(data = mst_df, aes(x=x, y=y, xend=xend, yend=yend), linewidth = 0.8, color = "black") +
  labs(x = "Component 1", y = "Component 2") + 
  pub_theme

# 生成图表
p1 <- base_p + geom_point(data = plot_df, aes(x, y, color = cluster_ordered), size=1) + ggtitle("P1: Ordered Cluster Trajectory")
p2 <- base_p + geom_point(data = plot_df, aes(x, y, color = Pseudotime), size=1) + scale_color_viridis_c() + ggtitle("P2: Pseudotime")
p3 <- base_p + geom_point(data = plot_df, aes(x, y, color = Pseudotime), size=0.8) + scale_color_viridis_c() + facet_wrap(~cluster_ordered, ncol = 6) + ggtitle("P3: Faceted Cluster") + theme(legend.position = "none")
p4 <- base_p + geom_point(data = plot_df, aes(x, y, color = Disease_Stage), size=1) + scale_color_manual(values = custom_cols) + ggtitle("P4: Disease Stage")
p5 <- base_p + geom_point(data = plot_df, aes(x, y, color = State), size=1) + ggtitle(paste0("P5: States (Root State: ", root_state, ")"))
p6 <- base_p + geom_point(data = plot_df, aes(x, y, color = Pseudotime), size=0.8) + scale_color_viridis_c() + facet_wrap(~Disease_Stage, ncol = 3) + ggtitle("P6: Faceted Disease") + theme(legend.position = "none")

# P7: 翻转排序屋脊图
lvl <- plot_df %>% group_by(cluster_ordered) %>% summarize(m = mean(Pseudotime)) %>% arrange(desc(m)) %>% pull(cluster_ordered)
plot_df$cluster_ordered <- factor(plot_df$cluster_ordered, levels = lvl)
p7 <- ggplot(plot_df, aes(x = Pseudotime, y = cluster_ordered, fill = cluster_ordered)) +
  geom_density_ridges(alpha = 0.85, scale = 1.6, color = "white", linewidth = 0.3) +
  theme_bw(base_size = 13) + theme(plot.title=element_text(face="bold",hjust=0.5), legend.position="none") +
  ggtitle("P7: Pseudotime Trajectory Order") + xlab("Pseudotime") + ylab("Cluster")

# =============== 5. 输出保存 ==================
pdf(out_pdf, width = 16, height = 12)
print(p1); print(p2); print(p3); print(p4); print(p5); print(p6); print(p7)
dev.off()

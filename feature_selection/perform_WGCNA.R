perform_WGCNA <- function(mtx, 
    type = c("unsigned", "signed", "signed hybrid"),
    corType = "pearson",
    trait = NULL, # vector or matrix
    normalize = FALSE,
    filter = FALSE,
    q = 0.25,
    powers = c(1:10, seq(12, 30)),
    minModuleSize = 50,
    TOM_size = 500,
    dir = NULL, 
    ncores = 12,
    ...){
    
    # loading required packages
    require(magrittr)
    require(WGCNA)
    require(qs)

    # enable multi threads 
    enableWGCNAThreads(nThreads = ncores)
    corFun <- ifelse(corType == "pearson", WGCNA::cor, WGCNA::bicor)
    robustY <- ifelse(corType == "pearson", TRUE, FALSE)
    maxPOutliers <- ifelse(corType == "pearson", 1, 0.05)
    type <- match.arg(type)

    # directory parse
    if(is.null(dir)) {
        stop("Output directory 'dir' cannot be NULL")
    }

    dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))
    if(!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
    }

    # filter features
    if(filter){
        mad <- apply(mtx, 1, mad)
        ids <- mad > max(quantile(mad, q), 0.01)
        table(ids)
        mtx <- mtx[ids, ]
    }
    message("after filtering...")

    # transpose, column sample, row feature
    mtx <- mtx %>% t %>% as.data.frame()

    # missing value detection
    gsg <- goodSamplesGenes(mtx, verbose = 3)
    if(!gsg$allOK){
        if(sum(gsg$goodGenes) > 0){
            genes <- colnames(mtx)[!gsg$goodGenes]
            message(paste0("The following genes are removed: ", str_c(genes, collapse = ",")))
        }

        if(sum(gsg$goodSamples) > 0){
            samples <- rownames(mtx)[!gsg$goodSamples]
            message(paste0("The following samples are removed: ", str_c(samples, collapse = ",")))
        }
        mtx <- mtx[gsg$goodSamples, gsg$goodGenes]
    }
    
    nGenes <- mtx %>% ncol
    nSamples <- mtx %>% nrow

    # sample clustering
    message("sample clustering...")
    sampleTree <- hclust(dist(mtx), method = "average")
    
    file <- paste0(dir, "sample_tree.pdf")
    pdf(file = file, width = 8, height = 6)
        plot(sampleTree, 
            main = "Sample Clustering", 
            sub = "", 
            labels = FALSE)
    dev.off()

    # pick softthreshold
    message("pick soft threshold...")
    sft <- pickSoftThreshold(data = mtx, 
        dataIsExpr = TRUE, 
        RsquaredCut = 0.85, 
        powerVector = powers, 
        corFnc = corFun, 
        networkType = type)

    # save object
    file <- paste0(dir, "sft.rds")
    saveRDS(sft, file = file)

    # save figures
    file <- paste0(dir, "picksoftthreshold.pdf")
    pdf(file, width = 10, height = 6)
        par(mfrow = c(1, 2))
        plot(x = sft$fitIndices[, 1], 
            y = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
            xlab = "Soft Threshold (power)", 
            ylab = "Scale Free Topology Model Fit", 
            type = "n", 
            main = "Scale Independence")
        text(x = sft$fitIndices[, 1], 
            y = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
            labels = powers, 
            col = "red")
        abline(h = 0.85, col = "red", lty = 2)

        plot(x = sft$fitIndices[, 1], 
            y = sft$fitIndices[, 5], 
            xlab = "Soft Threshold (power)", 
            ylab = "Mean Connectivity", 
            type = "n", 
            main = "Mean Connectivity")
        text(sft$fitIndices[,1], 
            sft$fitIndices[, 5], 
            labels = powers, 
            cex = 0.9, 
            col = "red")
    dev.off()

    power <- sft$powerEstimate

    # pick empirical power
    if(is.na(power)){
        power = ifelse(nSamples < 20, ifelse(type == "unsigned", 9, 18), 
            ifelse(nSamples < 30, ifelse(type == "unsigned", 8, 16), 
                ifelse(type == "unsigned", 6, 12)))
    }
    message(paste0("The power is: ", power))

    # build network
    message("building network...")
    cor <- WGCNA::cor
    file <- paste0(dir, "WGCNA_TOM")
    net <- blockwiseModules(datExpr = mtx, 
        power = power, 
        maxBlockSize = nGenes, 
        minModuleSize = minModuleSize, 
        TOMType = type, 
        reassignThreshold = 0, 
        mergeCutHeight = 0.25, 
        numericLabels = TRUE, 
        pamRespectsDendro = FALSE, 
        saveTOMs = TRUE, 
        corType = corType, 
        maxPOutliers = maxPOutliers, 
        loadTOM = FALSE, 
        saveTOMFileBase = "WGCNA_TOM", 
        verbose = 3)
    
    # save the object
    message("save network object...")
    file <- paste0(dir, "net.rds")
    saveRDS(net, file = file)

    message("number of modules and size: ")
    print(net$colors %>% table)

    # convert label to color
    moduleLables <- net$colors
    moduleColors <- labels2colors(moduleLables)

    # save module genes
    moduleGenes <- split(names(moduleLables), labels2colors(moduleLables))
    file = paste0(dir, "module_genes.rds")
    saveRDS(moduleGenes, file = file)

    file <- paste0(dir, "dendro_and_color.pdf")
    pdf(file = file, width = 7, height = 7)
        plotDendroAndColors(
            dendro = net$dendrograms[[1]], 
            colors = moduleColors[net$blockGenes[[1]]], 
            groupLabels = "Module Colors", 
            dendroLabels = FALSE, 
            hang = 0.03, 
            addGuide = TRUE, 
            guideHang = 0.05
        )
    dev.off()
    
    MEs <- net$MEs
    MEs_col <- MEs
    colnames(MEs_col) <- paste0("ME", labels2colors(readr::parse_number(colnames(MEs))))
    MEs_col <- orderMEs(MEs_col)
    qsave(MEs_col, paste0(dir, "ME_color.qs"))

    # module cor heatmap
    file = paste0(dir, "module_cor_heatmap.pdf")
    pdf(file = file, width = 7, height = 7)
        plotEigengeneNetworks(
            multiME = MEs_col,
            setLabels = "Eigengene", 
            marDendro = c(3, 3, 2, 4), 
            marHeatmap = c(3 ,4, 2, 2),
            plotDendrograms = T, 
            xLabelsAngle = 90)
    dev.off()

    load(net$TOMFiles[[1]])
    TOM <- TOM %>% as.matrix()
    set.seed(101)
    ids <- sample(ncol(TOM), TOM_size)
    subTOM <- TOM[ids, ids]
    disssubTOM <- 1 - subTOM
    plotTOM <- disssubTOM^power
    hc <- hclust(d = dist(plotTOM), method = "average")
    diag(plotTOM) <- NA

    message("produce TOM heatmap...")
    file <- paste0(dir, "TOM.pdf")
    pdf(file, width = 7, height = 7)
        TOMplot(dissim = plotTOM, 
            dendro = hc, 
            main = "TOM Network heatmap")
    dev.off()

    # module trait correlation
    message("conducte module trait analysis...")
    if(corType == "pearson"){
        modTraitCor <- cor(MEs_col, trait, use = "p")
        modTraitP <- corPvalueStudent(modTraitCor, nSamples)
    }else{
        modTraitCorP <- bicorAndPvalue(MEs_col, trait, robustY = robustY)
        modTraitCor <- modTraitCorP$bicor
        modTraitP <- modTraitCorP$p
    }
    textMtx <- paste0(signif(modTraitCor, 2), 
        "\n(", 
        signif(modTraitP, 1), 
        ")")
    dim(textMtx) <- dim(modTraitCor)
    file = paste0(dir, "module_trait_cor_heatmap.pdf")
    pdf(file, width = 7, height = 7*(nrow(modTraitCor)/25))
        par(mar = c(2,5,2,1))
        labeledHeatmap(Matrix = modTraitCor,
            xLabels = colnames(trait),
            yLabels = colnames(MEs_col),
            ySymbols = colnames(MEs_col),
            cex.lab = 0.5,
            colorLabels = FALSE,
            colors = blueWhiteRed(50),
            setStdMargins = FALSE,
            cex.text = 0.5,
            textMatrix = textMtx,
            # zlim = c(-0.5, 0.5),
            xLabelsAngle = 0,
            xLabelsAdj = 0.5,
            main = "Module Trait Relationships")
    dev.off()
    message("all analysis is done...")
}
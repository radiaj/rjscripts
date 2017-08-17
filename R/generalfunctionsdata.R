#general functions
corrdist <- function(x) as.dist(1-cor(t(x), method="pearson"))
hclust.avl <- function(x) hclust(x, method="ward.D2")

cbind.fill<-function(...){
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

# Functions for miR analysis - DDCT method
getDCtqPCR <- function(x, ei, ti){
    geom <- geometric.mean(x[ei])
    y <- x[ti] - geom
    return(y)
}

# Housekeeping genes
getTargetDCtqPCR <- function(mat, endoctrls=379:381, targets=1:378){
    Dctmat <- apply(mat, 2, getDCtqPCR, ei=endoctrls, ti=targets)
    return(Dctmat)
}


# Get the DDCt for TE and TS sample using dnorm, DctM1, and DctM2
getDDCtqPCRmean <- function(matDct, tumors, calibrators){
    mat <- matDct[, tumors]
    calmat <- matDct[, calibrators]
    calmean <- apply(calmat, 1, mean, na.rm=TRUE)
    print(head(calmat))
    DDCtmat <- sweep(mat, 1, calmean, "-") #- Read more at: http://scl.io/fCPufvMf#gs.t2=lklQ
    return(DDCtmat)
}

getDDCtqPCR <- function(matDct, tumors, calibrators){
    mat <- matDct[, tumors]
    calmat <- matDct[, calibrators]
    calmean <- apply(calmat, 1, median, na.rm=TRUE)
    print(head(calmat))
    DDCtmat <- sweep(mat, 1, calmean, "-") #- Read more at: http://scl.io/fCPufvMf#gs.t2=lklQ
    return(DDCtmat)
}

qPCRDDCtFoldChange <- function(x) 2^c(-x)

# General function for quantile normalization
quantile_normalisation <- function(df){
    df_rank <- apply(df,2,rank,ties.method="min")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)
    
    index_to_mean <- function(my_index, my_mean){
        return(my_mean[my_index])
    }
    
    df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
    rownames(df_final) <- rownames(df)
    return(df_final)
}

getMetaboliteTotalPool <- function(label13c, LABEL="", sampsInds=1:12){
    labInds <- intersect(sampsInds, label13c)
    meta13cGLUC <- MetaTotalPool[labInds[1:2]]
    
    #merged.data.frame = Reduce(function(...) merge(..., all=T, by=c("CodeClass", "Name", "Accession"), suffixes=names(NanoStringData)), NanoStringData)
    merged.data.frame = Reduce(function(...) merge(..., all=T, by=c("Metabolite")), meta13cGLUC)
    
    meta13cGLUC <- MetaTotalPool[labInds[1:2]]
    naive13c = Reduce(function(...) merge(..., all=T, by=c("Metabolite")), meta13cGLUC)
    
    meta13cGLUC <- MetaTotalPool[labInds[3:4]]
    thy1_13c = Reduce(function(...) merge(..., all=T, by=c("Metabolite")), meta13cGLUC)
    
    cd8_13c <- merge(naive13c, thy1_13c, all=T, by=c("Metabolite"), suffixes=c("_naive", "_thy1_1"))
    
    write.csv(cd8_13c, paste("cd8_13c_", LABEL, "_Naive_Thy1_1_TotalPool.csv", sep=""))
    
    # Matrix to look at the data across the samples
    cd8_13cMat <- cd8_13c[, -1]
    rownames(cd8_13cMat) <- cd8_13c[, 1]
    
    M1 <- as.matrix(cd8_13cMat[-1, ])
    
    return(M1)
    
}

getMetaboliteGLMres <- function(y, SampVsCtrl=BvsA, LABEL_13C = "GLUTAMINE_"){
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, contrast=BvsA)
    topTags(lrt)
    
    write.csv(topTags(lrt, n=Inf), file=paste(LABEL_13C, "_glm.csv", sep=""))
    return(as.data.frame(topTags(lrt, n=Inf)))
}

getVolcanoPlot13C_TotalPools <- function(res){
    require("ggplot2")
    require("ggrepel")
    res$P.Value <- sapply(res$PValue, function(x) ifelse(x <1e-291, 1e-291, x))
    res$Significant <- ifelse(res$P.Value < 0.001, "p < 0.001", "Not Sig")
    res$Gene <- rownames(res)
    
    p <- ggplot(res, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = Significant)) +
    scale_color_manual(values = c("grey", "red")) + scale_y_log10() +
    theme_bw(base_size = 16) +
    geom_text_repel(
    data = subset(res, abs(logFC) > 10),
    aes(label = Gene),
    size = 4)
    return(p)
    
}

getVolcanoPlot13CByDayComp <- function(res){
    require("ggplot2")
    require("ggrepel")
    res$P.Value <- sapply(res$PValue, function(x) ifelse(x <1e-291, 1e-291, x))
    res$Significant <- ifelse(res$P.Value < 0.001, "p < 0.001", "Not Sig")
    res$Gene <- rownames(res)
    
    p <- ggplot(res, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = Significant)) +
    scale_color_manual(values = c("grey", "red")) + scale_y_log10() +
    theme_bw(base_size = 16) +
    geom_text_repel(
    data = subset(res, P.Value < 0.001),
    aes(label = Gene),
    size = 4)
    return(p)
    
}

getGAMfile <- function(df, colids, fname){
    gene.de <- df[, colids]
    colnames(gene.de) <- c("ID", "pval", "log2FC", "baseMean")
    print(head(gene.de))
    write.table(gene.de, file=fname, sep="\t", quote=FALSE)
    return(gene.de)
}


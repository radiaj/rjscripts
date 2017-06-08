##################################################################################################################
# R script to run DESEQ2 
# By Radia Johnson, Ph.D
# Date: May. 24th, 2016
##################################################################################################################
# Functions
getCountVector <- function(df, column=3, removerows=1:4){
	x <- df[-removerows, column]
	return(x)
}
#----------------------------------------------------------

library("DESeq2")

setwd("J:/JONES_LAB/COLLABORATIONS/Collaboration_With_Dr_Pollack/pollack_star_counts")

sampleInfo <- read.csv("../Samples - prostate cancer RNAseq.csv", header=FALSE, skip=1, stringsAsFactors=TRUE)
files <- list.files(pattern="ReadsPerGene.out.tab")
countDFlist <- lapply(files, read.delim, header=FALSE, row.names=1)

# Check the dataframes
head(lapply(countDFlist, head))

# rename the list items
names(countDFlist) <- files

countsMatrix <- sapply(countDFlist, getCountVector)

# Add the gene names
identical(rownames(countDFlist[[1]]), rownames(countDFlist[[2]]))

geneNames <- rownames(countDFlist[[1]])[-c(1:4)]
rownames(countsMatrix) <- geneNames

# Organize the columns by group
colnames(countsMatrix)
xx <- colnames(countsMatrix)

sampleNames <- gsub("ReadsPerGene.out.tab", "", xx)

infoNames <- sampleInfo$V1

# Get the indices that match between the GQ names and the counts files 
indsMat <- match(infoNames, sampleNames)
sampleNames[indsMat]

#======================================================================
# Run DESeq2 on the raw counts matrix
#======================================================================
COMPARISON = "VentralvsAnterior"

# Create the ordered counts matrix
cnts <- countsMatrix[, indsMat]
colnames(cnts) <- infoNames[indsMat]
groups = sampleInfo$V3
pairs = as.factor(sampleInfo$V2)

head(cnts)
grp <- factor(groups)
grp
pairs = as.factor(sampleInfo$V2)

coldat=DataFrame(grp=grp, pairs=pairs)

#dds <- DESeqDataSetFromMatrix(cnts, colData=coldat, design = ~batch+ grp)

dds <- DESeqDataSetFromMatrix(cnts, colData=coldat, design = ~ pairs+grp)
dds <- DESeq(dds)

## Extract results from DESeq2 analysis
res <- results(dds)
summary(res)

# Doing the normalization after filtering genes that are not expressed across the samples
dds1 <- dds[rowSums(counts(dds)) > 1, ]
rld1 <- rlogTransformation(dds1, blind=TRUE) # local regression fit was used to transform the data
head(assay(rld1))

write.csv(assay(rld1), paste(COMPARISON,"rlogTransformedCnts.csv", sep=""))
#------------------------------------
# Boxplot of gene expression data
#------------------------------------
M1 <-assay(rld1)
colnames(M1) = sampleInfo$V4

par(mar=c(15,5,2,2))
boxplot(M1, col="blue", las=2, cex.axis=1.5)

save(cnts, dds, dds1, rld1, M1, coldat, countsMatrix, file=paste(COMPARISON, "DESeq2_May12_2017.rds", sep="")) 
write.csv(cnts, paste(COMPARISON, "rawCounts.csv"))

#============================================================================================
# Calculate the sample to sample distances to look at the sample clustering
#============================================================================================
require("RColorBrewer")
require("gplots")

# Color palette
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)

distsRL <- dist(t(assay(rld1)))
mat<- as.matrix(distsRL)
rownames(mat) <- colData(dds1)$grp

#sampleNames <- gsub("ReadsPerGene.out.tab", "", colnames(mat)) # MODIFY THIS LINE ACCORDINGLY
#colnames(mat) <- sampleNames

hc <- hclust(distsRL)
par(mar=c(15,5,2,2))
heatmap.2(mat, Rowv=as.dendrogram(hc),
symm=TRUE, trace='none',
col = rev(hmcol), margin=c(13, 13), key.title="")
#dev.copy(png,'deseq2_heatmap_samplebysampleV3.png')
#dev.off()


#-------------------------------------------------
# Get the list of differentially expressed genes - IOSE (no replicates)
#-------------------------------------------------
# Other functions used when no biological replicate available
getRegLog2FoldChange <- function(rld, ref=1, samp=2){
	res <- data.frame(assay(rld)[, c(ref, samp)], avgLogExpr = ( assay(rld)[,samp] + assay(rld)[,ref] ) / 2, rLogFC = assay(rld)[,samp] - assay(rld)[,ref] )
	return(res)
}
# This ranking put ones genes on top which are strongly downregulated in the second sample compared to the first one (control).
# If you do this with normal log expressions, the weak genes will appear at the top because they are noisiest and hence tend to have exaggerated fold changes. 


getRegLog2FCtable <- function(fname, rf=2, sp=1, rld1){

	resM1 <- getRegLog2FoldChange(rld1, ref=rf, samp=sp)
	#print(head(resM1))

	resM1 <- resM1[ order(resM1$rLogFC), ]
	print(head(resM1))
	print(dim(resM1))

	write.csv(as.data.frame(resM1), file=paste(fname, 'Regularized_log2_FC.csv', sep="_"))
}

colData(dds1)

res <- res[order(res$padj),]

write.csv(as.data.frame(res), file=paste(COMPARISON, "diffExprsGenes.csv", sep=""))

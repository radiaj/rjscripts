#======================================================
# R script to look at PGC1a levels in metastasis
# By Radia Johnson, Ph.D.
# Date: 06/08/2017
#======================================================
# On BRIAREE
#======================================================

# module add R/3.2.1-gcc
# R
library("limma")

setwd("/RQexec/johnsonr/SIEGEL_LAB/R_analysis/LCM")
# write.csv(list.files("2776", full.names=TRUE), "2776_files.csv")
# write.csv(list.files("2792", full.names=TRUE), "2792_files.csv")

list.files()

cl2776targets <- read.csv("cell_line2776targets.csv", row.names=1)
cl2792targets <- read.csv("cell_line2792targets.csv", row.names=1)

# Treat each cell line seperately then combine with the Combat method accounting for the cell line batch effect
# For CELL LINE 2776
targets = cl2776targets
RG <- read.maimages(targets$Filename, source="agilent") # median signal used

# To inspect the results
show(RG)
summary(RG$R)

# For CELL LINE 2792
targets = cl2792targets
RG2 <- read.maimages(targets$Filename, source="agilent") # median signal used
# To inspect the results
show(RG2)
summary(RG2$R)

# Pre-processing and normalizing the data
preProcessAndNormalize <- function(RG) {
    RG.within <- normalizeWithinArrays(RG, method="loess", bc.method="none")
    RG.norm <- normalizeBetweenArrays(RG.within, method="quantile")
    RG.norm <- RG.norm[RG.norm$genes$ControlType == 0, ] # remove the control probes
    agg.M <- aggregate(RG.norm$M,list(RG.norm$genes$ProbeName),function(x) mean(x,na.rm=T))
    agg.A <- aggregate(RG.norm$A,list(RG.norm$genes$ProbeName),function(x) mean(x,na.rm=T))
   	RG.norm <- RG.norm[!duplicated(RG.norm$genes$ProbeName),]
    RG.norm$M <- agg.M[match(RG.norm$genes$ProbeName,agg.M[,1]),2:ncol(agg.M)]
    RG.norm$A <- agg.A[match(RG.norm$genes$ProbeName,agg.A[,1]),2:ncol(agg.M)]

    return(RG.norm)
}


RG1.norm <- preProcessAndNormalize(RG)
RG2.norm <- preProcessAndNormalize(RG2)


# AIMS
#The idea will be to reanalyse old LCM data we have collected (agilent microarrays) and to look at the metabolic aspect.
#Seb had injected 2 cell lines (2776 and 2792) in mammary fat pad and look at the primary tumor or in splenic looking at the liver metastases.
#In the mammary fat pad, we have data for the core and for the margin.
#In the liver metastases, we have data for the tumor cells (core vs margin) and also for the liver (adjacent vs distal) for various time (10days, 2 weeks and 3 weeks).
#I am interested in the glutamine and glucose metabolism and particularly in the lipid metabolism.


preProcessAndNormalize <- function(RG) {
                RG.within <- normalizeWithinArrays(RG, method="loess", bc.method="none")
                RG.norm <- normalizeBetweenArrays(RG.within, method="quantile")
                agg.M <- aggregate(RG.norm$M,list(RG.norm$genes$ProbeName),function(x) mean(x,na.rm=T))
 
                agg.A <- aggregate(RG.norm$A,list(RG.norm$genes$ProbeName),function(x) mean(x,na.rm=T))
   		RG.norm <- RG.norm[!duplicated(RG.norm$genes$ProbeName),]
    		RG.norm$M <- agg.M[match(RG.norm$genes$ProbeName,agg.M[,1]),2:ncol(agg.M)]
    		RG.norm$A <- agg.A[match(RG.norm$genes$ProbeName,agg.A[,1]),2:ncol(agg.M)]
 
#             res = list(RG.norm, agg.M, agg.A)
               return(RG.norm)
}

################################ 
# Extract the gene expression data
# Key points: This expression set has multiple probes per gene along with duplicate probes (dimension= 42405 X 256)
# Duplicate probes are removed by taking the most variant probe per gene
################################ 

edata1 <-sapply(-RG1.norm$M,function(x){x},simplify = T) 
edata2 <-sapply(-RG2.norm$M,function(x){x},simplify = T) 
edata = cbind(edata1, edata2)

# Save the edata objects

metadata <- read.csv("targetsRJ05032017.csv", stringsAsFactors = FALSE)
ccdata <- lapply(colnames(edata), function(x){unlist(strsplit(x, split = "/"))[2]})
ccdata <- unlist(ccdata)

aadata <- unlist(lapply(metadata[,"Filename"], function(x){unlist(strsplit(x, split = "[.]"))[1]}))
aadata <- unlist(lapply(aadata, function(x){unlist(strsplit(x, split = "/"))[2]}))
names(aadata) <- metadata[,"sampleName"]
test <- names(aadata[which(aadata %in% ccdata)])
colnames(edata) <- test

if(identical( RG1.norm$genes$ProbeName, RG2.norm$genes$ProbeName)) {
	rownames(edata) <- RG1.norm$genes$ProbeName
	print("The probe names were added." )
}

head(edata)

# Batch correct the samples by cell line
library("sva")
batchCL = as.factor(metadata$cell_line)
combat_edata = ComBat(edata, batchCL, par.prior=TRUE)

# Save the result objects
save(RG1.norm, RG2.norm, edata, metadata, combat_edata, file="ChristineLCMrenorm06092017.RData")

# Remove the duplicate probes for each gene (if they exist) by choosing max IQR
GeneSymbol = RG1.norm$genes$GeneName
dataM = data.frame(GeneSymbol, combat_edata)

# Functions required 
TopIqrSymbolmatrix <- function(a, sampCols=1:12, geneCol=7){
  ### Collapses a dataset given by data according to most variable probes or by
  ### measure listed in method

  	data <- a[, sampCols]
  
  	a$iqr_value <- apply(data, 1, IQR, na.rm=TRUE)
	
	aa <- a[order( -a$iqr_value ), ] 
	print(head(aa[, 1:3]))
	print(tail(colnames(aa)))

	GC2.exprs <- aa[ !duplicated(aa[, geneCol]), ]    
	print(dim(GC2.exprs))

 	return(GC2.exprs)

}

dataMtopIQR = TopIqrSymbolmatrix(dataM, sampCols=2:65, geneCol=1)
head(dataMtopIQR)

M1 = as.matrix(dataMtopIQR[, 2:65])
rownames(M1) = dataMtopIQR$GeneSymbol
colnames(M1) = gsub("X", "", colnames(M1))
head(M1)

save(M1, dataMtopIQR, file="ChristineLCMmatrix.RData")

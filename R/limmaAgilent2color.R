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

# REACH HERE

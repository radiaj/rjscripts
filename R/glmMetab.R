#=========================================
# R script to analyze the data
# By Radia Johnson, Ph.D.
# Date: January 17th, 2017
#=========================================
source("/Users/radiamariejohnson/rjcode/repos/r-functions/generalfunctionsdata.R")
library("gdata")
list.files()

# Load metabolite data
metData <- read.xls("H1299_GAD1_GLNpulse_ 06.30.16_ New template analysis.xlsx", sheet="data", row.names=1)

group <- rep(c("CTRL", "GAD"), each=9)

timepoint <- rep(rep(c("0.5hrs", "1hrs", "1.5hrs"), each=3), 2)

#--------------------------------
# GLM model for DEG analysis
#--------------------------------
library("edgeR")
gg <- paste(group,  timepoint, sep="_")
grp <- factor(gg, levels=unique(gg)[c(1:3,6, 4:5)])

M1 <- metData
M2 <- apply(M1, c(1, 2), function(x) ifelse(is.na(x), 0, x))
M3 <- apply(M2, c(1, 2), function(x) ifelse(x < 0, 0, x))

y <- DGEList(counts=M3, group=grp)
y <- calcNormFactors(y, method="none")
design <- model.matrix(~0 + grp)
y <- estimateDisp(y,design)


M4 <- cpm(y, log=TRUE)
par(mar=c(12, 3.5, 2, 1))
boxplot(M4, col="cyan", las=2, cex.axis=0.7)

# Tested the TMM normalization ----
y <- DGEList(counts=M3, group=grp)
y <- calcNormFactors(y, method="TMM")
design <- model.matrix(~0 + grp)
y <- estimateDisp(y,design)


M4 <- cpm(y, log=TRUE)
par(mar=c(12, 3.5, 2, 1))
boxplot(M4, col="cyan", las=2, cex.axis=0.7, main="TMM normalization")
#-----------------------------------
# USED no normalization for downstream analysis email Takla the results
#---------------------------------------------------
# Get the differentially expressed metabolites
#---------------------------------------------------

y <- DGEList(counts=M3, group=grp)
y <- calcNormFactors(y, method="none")
design <- model.matrix(~0 + grp)
y <- estimateDisp(y,design)

colnames(design)

# Save the results
RES <- list()
BvsA <- makeContrasts(grpGAD_0.5hrs - grpCTRL_0.5hrs,  levels=design)
yy <- "shGADvsCTRL_30mins"
RES[[yy]]<- getMetaboliteGLMres(y, SampVsCtrl=BvsA, LABEL_13C = yy)

# For the other pairwise comparisons
BvsA <- makeContrasts(grpGAD_1hrs - grpCTRL_1hrs, levels=design)
yy <- "shGADvsCTRL_1hr"
RES[[yy]]<- getMetaboliteGLMres(y, SampVsCtrl=BvsA, LABEL_13C = yy)

BvsA <- makeContrasts(grpGAD_1.5hrs - grpCTRL_1.5hrs,  levels=design)
yy <- "shGADvsCTRL_1.5hrs"
RES[[yy]] <- getMetaboliteGLMres(y, SampVsCtrl=BvsA, LABEL_13C = yy)

lapply(RES, head)

# For the positive ion mode
for(i in c(1:3)){
	pp <- getVolcanoPlot13CByDayComp(RES[[i]])
	ggsave(paste(names(RES)[i],".png", sep=""))
}


# Get the tables for the GAM analysis
library("gdata")
metabData <- read.xls("shGADvsCTRL_1.5hrs_glm.xlsx", stringsAsFactors=FALSE, row.names=1)
# ID, pval, log2FC, baseMean

# Save File for GAM analysis
getGAMfile <- function(df, colids, fname){
                gene.de <- df[, colids]
                colnames(gene.de) <- c("ID", "pval", "log2FC", "baseMean")
                print(head(gene.de))
                write.table(gene.de, file=fname, sep="\t", quote=FALSE)
                return(gene.de)
}

for(i in names(RES)){
	res <- RES[[i]]
	xx <- metabData[rownames(res), 8]
	xx[which(is.na(xx))]="HMDB00208"
	resB <- cbind(gene=xx, res)
	head(resB)
	getGAMfile(resB, c(1, 5, 2, 3), fname=paste(i, "_GAMtable.txt"))
}

# DEG tables

H1299table <- read.csv("H1299_GADsi_vs_H1299_CTL_DEG.csv")
Zzz <- getGAMfile(H1299table, c(8, 6, 3, 2), fname=paste("H1299_GADsi_vs_CTRL", "_GAMtable.txt"))

A549table <- read.csv("A549_HBSS_vs_A549_DMEM_DEG.csv")
Zzz <- getGAMfile(A549table, c(8, 6, 3, 2), fname=paste("A549_HBSS_vs_DMEM", "_GAMtable.txt"))







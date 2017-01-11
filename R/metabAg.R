##################################################################################################################
# R script to analyse Metabolite Data 
# By Radia Johnson, Ph.D
# Date Modified: Jan. 9th, 2017
# INCOMPLETE - DRAFT ONLY
##################################################################################################################
# Functions needed
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
# Running the analysis for 13C GLUCOSE and ACETATE 
setwd("~/Documents/agios data")

# POSITIVE ION MODE
negMaster <- read.delim("A2941_Neg_Master.txt", stringsAsFactors=FALSE)
posMaster <- read.delim("A2941_cells_pos_Correction Process Data.txt", stringsAsFactors=FALSE)

#--------------------------------------------------------
# Create a pool total matrix for all the metabolites
#--------------------------------------------------------
# For negMaster ----
# totalPool <- negMaster[, c(1, 8, 10, 12:14)]
# For posMaster ----
totalPool <- posMaster[, c(1, 8, 10:13)]

Pooldf <- unique(totalPool)

# For negMaster ----
#write.csv(Pooldf, "A2941_Neg_Pool_Total.csv")
# For posMaster ----
write.csv(Pooldf, "A2941_Pos_Pool_Total.csv")


# Create a list object by Tissue
totalPoolist <- list()
tissue <- unique(Pooldf$Tissue)

# For negMaster ----
#for(i in tissue[-2]){
#	totalPoolist[[i]] <- subset(Pooldf, Tissue==i)
#}
# For posMaster ----
for(i in tissue){
	totalPoolist[[i]] <- subset(Pooldf, Tissue==i)
}

# Visualize the data frames
lapply(totalPoolist, head)

# Use dcast from the reshape2 each sheet sperated by label and time point
library("reshape2")
isolabesl <- unique(Pooldf$Isotope.Label)[1:3]
timepoints <- unique(Pooldf$Timepoint)[1:2]

df <- totalPoolist[[1]]
lsDF <- split(df, with(df, interaction(Isotope.Label,Timepoint)), drop = TRUE)

# For each data frame in lsDF get the metabolite total pool
test <- dcast(lsDF[[1]], Metabolite ~ SubjectID, value.var="Pool.Total")

# FOR ALL COMBINATIONS *****************
MetaTotalPool <- list()
for(j in names(totalPoolist)){

	df <- totalPoolist[[j]]
	lsDF <- split(df, with(df, interaction(Isotope.Label,Timepoint)), drop = TRUE)

	for(i in names(lsDF)){
		MetaTotalPool[[paste(j, i, sep="-")]] <- dcast(lsDF[[i]], Metabolite ~ SubjectID, value.var="Pool.Total")
	}
}

# Test the results
lapply(MetaTotalPool, head)

# Save the list object for negMaster
#saveRDS(MetaTotalPool, "MetaTotalPool.rds")

# Save the list object for the posMaster
saveRDS(MetaTotalPool, "posMetaTotalPool.rds")

# 
#=====================================================
# To run DEG analysis on total pool by 13C label
#=====================================================
# Restart from here
# Load the functions above
# setwd("~/Documents/agios data/positiveIon")
# MetaTotalPool <- readRDS("posMetaTotalPool.rds")

library("edgeR")
library("gdata")
SampDetails <- read.delim("../SampleDetails.txt", row.names=1)

# 13C-Glucose, Glutamine, Acetate --------------------------------------------------------
gluc13c <- grep("GLUCOSE", names(MetaTotalPool))
glutamine13c <- grep("GLUTAMINE", names(MetaTotalPool))
acetate13c <- grep("ACETATE", names(MetaTotalPool))
#-------------------------------------------------------------
LAB13C="ACETATE"
xx <- acetate13c # Change according to the 13C label of interest
# CREATE A FUNCTION To RETURN xx

M1 <- getMetaboliteTotalPool(label13c=xx, LABEL=LAB13C)

samptested <- sapply(strsplit(colnames(M1), "_"), function(x) x[[1]])
subsetested <- sapply(strsplit(colnames(M1), "_"), function(x) x[[2]])
timep <- SampDetails[samptested, 2]
gg <- paste(subsetested,  timep, sep="_")
grp <- factor(paste(subsetested,  timep, sep="_"), levels=unique(gg))

M2 <- apply(M1, c(1, 2), function(x) ifelse(is.na(x), 0, x))
M3 <- apply(M2, c(1, 2), function(x) ifelse(x < 0, 0, x))

y <- DGEList(counts=M3, group=grp)
y <- calcNormFactors(y, method="none")
design <- model.matrix(~0 + grp)
y <- estimateDisp(y,design)

M4 <- cpm(y, log=TRUE)
par(mar=c(12, 3.5, 2, 1))
boxplot(M4, col="cyan", las=2, cex.axis=0.7)

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
BvsA <- makeContrasts(grpnaive_6 - grpnaive_2.5,  levels=design)
yy <- paste(LAB13C, "_Naive_Day6_vs_2", sep="")
RES[[yy]]<- getMetaboliteGLMres(y, SampVsCtrl=BvsA, LABEL_13C = yy)
# For the other pairwise comparisons
BvsA <- makeContrasts(grpthy1_2.5 - grpnaive_2.5,  levels=design)
yy <- paste(LAB13C, "_Thy1_1_vs_Naive_Day2", sep="")
RES[[yy]]<- getMetaboliteGLMres(y, SampVsCtrl=BvsA, LABEL_13C = yy)

BvsA <- makeContrasts(grpthy1_6 - grpthy1_2.5,  levels=design)
yy <- paste(LAB13C, "_Thy1_1_Day6_vs_2", sep="")
RES[[yy]] <- getMetaboliteGLMres(y, SampVsCtrl=BvsA, LABEL_13C = yy)

BvsA <- makeContrasts(grpthy1_6 - grpnaive_6,  levels=design)
yy <- paste(LAB13C, "_Thy1_1_vs_Naive_Day6", sep="")
RES[[yy]] <- getMetaboliteGLMres(y, SampVsCtrl=BvsA, LABEL_13C = yy)

# For the Volcano plots
# To run for all 
for(i in c(2, 4)){
	pp <- getVolcanoPlot13C_TotalPools(RES[[i]])
	ggsave(paste(names(RES)[i],".png", sep=""))
}

for(i in c(1, 3)){
	pp <- getVolcanoPlot13CByDayComp(RES[[i]])
	ggsave(paste(names(RES)[i],".png", sep=""))
}

# For the positive ion mode
for(i in c(1:4)){
	pp <- getVolcanoPlot13CByDayComp(RES[[i]])
	ggsave(paste(names(RES)[i],".png", sep=""))
}



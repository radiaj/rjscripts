#--------------------------------------------------------------------
# GAGE Pathway Analysis using the gs.KSTest
#--------------------------------------------------------------------
library("gage")
data(kegg.gs)

# Using proposed method for preprocessing counts 
# Based on https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf 
getGAGEnormalizedCnts <- function(processed.cnts){
	cnts <- processed.cnts
	sel.rn <- rowSums(cnts) !=0
	cnts <- cnts[sel.rn, ]
	libsizes <- colSums(cnts)
	size.factor <- libsizes/exp(mean(log(libsizes)))
	cnts.norm <- t(t(cnts)/size.factor)
	cnts.norm <- log2(cnts.norm+8)
	return(cnts.norm)
}

# Get the GAGE recommended normalized counts
cnts.norm <- getGAGEnormalizedCnts(cnts)

# Custom made function to get and save the significant gene sets as a table using the KSTest instead if the gage test
getSignificantPathwaysGS <- function(cnts.norm, reference, sample, groupsRS, fname, geneset=kegg.gs, ...){
	
	ref.idx <- grep(reference, groupsRS, fixed=TRUE)
	print(ref.idx)
	samp.idx <- grep(sample, groupsRS, fixed=TRUE)
	print(samp.idx)
	filename=paste(fname, sample, "_vs_", reference, sep="")
	cnts.kegg.p <- gage(cnts.norm, gsets=geneset, ref=ref.idx, samp=samp.idx, saaTest=gs.KSTest, compare="as.group", set.size=c(2, 500)) 
	cnts.kegg.sig <- sigGeneSet(cnts.kegg.p, heatmap=F) 
	getTablesKeggSign(filename, cnts.kegg.sig)
	RESULTS <- list(SIGNIFICANT=cnts.kegg.sig, ALL_PATHWAYS=cnts.kegg.p)
	return(RESULTS)
}

getTablesKeggSign <- function(GRPNAME, cnts.kegg.sig){
    filename1 <- paste(GRPNAME, "Signifcant_Pathways.txt", sep="_")
    print(filename1)
    cat("Upregulated Pathways", sep = "\t", file = filename1)
    suppressWarnings(write.table(cnts.kegg.sig$greater, filename1, sep="\t", quote=FALSE, row.names=T, col.names=NA, append=T))
    cat("Downregulated Pathways", sep = "\t", file = filename1, append=T)
    suppressWarnings(write.table(cnts.kegg.sig$less, filename1, sep="\t", quote=FALSE, row.names=T, col.names=NA, append=T))
}


# Create the group combination to test -----------------------------------
grpCombn <- matrix(NA, nrow=2, ncol=2)
grpCombn[1, ] <- c("CTRL", "KO")
grpCombn[2, ] <- c("KO", "CTRL")

M1 <- cnts.norm
sampleCondition = rep(c("KO", "CTRL"), each=4)
colnames(M1) <- as.character(sampleCondition)

# For the human Kegg pathways
for(i in 1:2){
	kegg.sig.g1 <- getSignificantPathwaysGS(M1, fname="Kegg_gs_", geneset=kegg.gs, reference=grpCombn[i, 1], sample=grpCombn[i, 2], groupsRS=colnames(M1))
}


# For the heatmaps of the genes that drive the significance of the pathway
getGeneHeatmapForPathwayRJ <- function(cnts.norm, path, genesets=METFGFR1, reference=ref, sample=samp, ref1=1:67, samp1=68:98, fname="",  groupsRS=colnames(M1), ...){
		outname = paste(fname, path, sep="_")
		gs <- genesets[[path]]

		ref.idx <- grep(reference, groupsRS, fixed=TRUE)
		print(ref.idx)
		samp.idx <- grep(sample, groupsRS, fixed=TRUE)
		print(samp.idx)

		essData <- essGene(gs, cnts.norm, ref =ref.idx, samp =samp.idx, compare='as.group', ...)
		print(head(essData))

		geneData(genes = gs, exprs = essData, ref = ref1, 
    	samp = samp1 , outname = outname, txt = T, heatmap = T,
    	Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = F, margin=c(5, 10), cexRow=0.8, ...)

}


pathwaysToSee <- c("hsa04110 Cell cycle", "hsa03013 RNA transport")
for(j in pathwaysToSee){

	for(i in 1:2){
		XX <- getGeneHeatmapForPathwayRJ(M1, j, genesets=kegg.gs, reference=grpCombn[i, 1], sample=grpCombn[i, 2], ref1=1:4, samp1=5:8, fname="KOvsCTRL")
	}

}


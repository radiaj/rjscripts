
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
colData(dds1)

getRegLog2FCtable <- function(fname, rf=2, sp=1, rld1){

	resM1 <- getRegLog2FoldChange(rld1, ref=rf, samp=sp)
	#print(head(resM1))

	resM1 <- resM1[ order(resM1$rLogFC), ]
	print(head(resM1))
	print(dim(resM1))

	write.csv(as.data.frame(resM1), file=paste(fname, 'Regularized_log2_FC.csv', sep="_"))
}

fname <- "IOSE_shRNA_CDK6-1"
getRegLog2FCtable(fname, rf=1, sp=2, rld1)

fname <- "IOSE_shRNA_CDK6-2"
getRegLog2FCtable(fname, rf=1, sp=3, rld1)

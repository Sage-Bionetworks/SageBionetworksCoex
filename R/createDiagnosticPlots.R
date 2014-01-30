

# This function creates the diagnostic plots for the coexpression analysis
# Inputs:
# fileType  - one of "pdf", "jpeg", "png"
# outputDir - directory or NULL to use system temporary directory
# connectivity - the vector of connectivity (k) values for each gene
#
# Outputs:
# Vector of file paths (directory and file name) to the created files
#
# Files created:
# ScaleFree.<filetype>	"Scale free" plot (P(k) vs. k) for optimum value of beta
# CorHeatmap.<filetype>	heatmap of |cc|^beta, ordered according to dendrogram, with color-coded module membership along the margins
# GeneDendrogram.<filetype>	gene dendrogram with color-coded modules and color names
# TOMHeatmap.<filetype>	heatmap of TOM, rows and columns sorted by gene dendrogram
# GenePCDendrogram.<filetype>	dendrogram of clusters of module representative genes
# ExpressionHeatmap.<filetype>	heatmap of expression data, organized by gene and sample clusters

#
createDiagnosticPlots <-
		function(
				filetype="jpg", 
				outputDir=NULL,
				expressionData,
				dichotCor,
				tomDist,
				connectivity,
				beta,
				geneTree,
				geneModules,
				genePCTree) 
{
	checkVersion(2,13)	
	
	SCALE_FREE = "ScaleFree"
	CORMAP = "CorHeatmap"
	GENE_DENDROGRAM = "GeneDendrogram"
	TOM_HEATMAP = "TOMHeatmap"
	GENE_PC_DENDROGRAM = "GenePCDendrogram"
	## SAMPLE_DENDROGRAM = "SampleDendrogram"
	EXPRESSION_HEATMAP = "ExpressionHeatmap"
	## COEFF_K_IN_MOD = "CoeffKInMod"
	
	
	filetype = tolower(filetype)
	if (!(filetype=="jpeg" || filetype=="jpg" || filetype=="ps" || filetype=="png" || filetype=="pdf")) {
		stop(paste("Expected jpg, jpeg, ps, png, or pdf but found ", filetype, ".", sep=""))
	}
	
	if (is.null(outputDir)) {
		outputDir=tempdir()
	}
	fileList <- list()
	
	# ScaleFree.<filetype>	"Scale free" plot (P(k) vs. k) for optimum value of beta
	filePath=file.path(outputDir, paste(SCALE_FREE, filetype, sep="."))
	fileList$SCALE_FREE <- filePath
	openImgDev(filePath)
	par(mfrow=c(1,1))
	scaleFreePlot(connectivity, nBreaks = 15, truncated = TRUE, main = paste("beta=", beta))
	dev.off()
	
	# Cormap.<filetype>	heatmap of |cc|^beta, ordered according to dendrogram, with color-coded module membership along the margins
	filePath=file.path(outputDir, paste(CORMAP, filetype, sep="."))
	fileList$CORMAP <- filePath
	corHeatmap(filePath, dichotCor, geneTree, geneModules)
	
	# HierModLabeled.<filetype>	gene dendrogram with color-coded modules and color names
	filePath=file.path(outputDir, paste(GENE_DENDROGRAM, filetype, sep="."))
	fileList$GENE_DENDROGRAM <- filePath
	openImgDev(filePath)
	par(mfrow=c(1,1))
	plotDendroAndColors(dendro=geneTree, colors=geneModules)
	dev.off()
	
	# Heatmap.<filetype>	heatmap of TOM, rows and columns sorted by gene dendrogram
	filePath=file.path(outputDir, paste(TOM_HEATMAP, filetype, sep="."))
	fileList$TOM_HEATMAP <- filePath
	tomHeatmap(filePath, tomDist, geneTree, geneModules)
	
	# PCsHierClust.<filetype>	dendrogram of clusters of module representative genes
	filePath=file.path(outputDir, paste(GENE_PC_DENDROGRAM, filetype, sep="."))
	fileList$GENE_PC_DENDROGRAM <- filePath
	openImgDev(filePath)
	par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
	plot(genePCTree, xlab="",ylab="",main="",sub="")
	par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
	dev.off()

	## # SampleHierModules.<filetype>	dendrogram of samples, grouped into two modules
	## filePath=file.path(outputDir, paste(SAMPLE_DENDROGRAM, filetype, sep="."))
	## fileList$SAMPLE_DENDROGRAM <- filePath
	## openImgDev(filePath)
	## par(mfrow=c(1,1))
	## plotDendroAndColors(dendro=sampleTree, colors=sampleModules, main="Sample Dendrogram")
	## dev.off()
	
	# ExprBands.<filetype>	heatmap of expression data, organized by gene and sample clusters
	#Now we create color plots for each module
	#For each module there should be a clear band structure, 
	#i.e. genes #corresponding to a given array (column) should have the same color.
	sampleCluster <-.clusterSamples(expressionData)
	filePath=file.path(outputDir, paste(EXPRESSION_HEATMAP, filetype, sep="."))
	fileList$EXPRESSION_HEATMAP <- filePath
	openImgDev(filePath)
	par(mfrow=c(1,1))
	iExprMatrix = t(expressionData)[geneTree$order, sampleCluster$sampleTree$order]
	heatmap(as.matrix( iExprMatrix ), Rowv=NA, RowSideColors=as.character(geneModules[geneTree$order]),
			Colv=NA, ColSideColors=sampleCluster$module[sampleCluster$sampleTree$order],
			scale="row", revC=F, xlab="", col = .rgcolorsFunc(50))
	dev.off()

	## # CoeffKInMod.<filetype> scatter plot of within-module cluster coefficient vs. node connectivity ('k') for each gene
	## filePath=file.path(outputDir, paste(COEFF_K_IN_MOD, filetype, sep="."))
	## fileList$COEFF_K_IN_MOD <- filePath
	## openImgDev(filePath)
	## par(mfrow=c(1,1))
	## plot(intraModularStatistics$k.in, intraModularStatistics$cc.in, col=as.character(geneModules), xlab="no of connections", ylab="cluster coefficient")
	## dev.off()
	
	fileList
}

# Input:  
# expressionData - a data frame with rows as samples and columns as genes/probes
#  	(i.e. the matrix has the same orientation as that passed to the gene clustering function.
#	This function will transpose the data to use it.)
# nocluster - the number of clusters to create.  The default is 2, e.g. control/treatment
#
# Output:
# the dendrogram of samples
# a vector of memberships, i.e. output[i] is the module to which sample i belongs
.clusterSamples <- 
		function(
				expressionData, 
				nocluster=2) 
{
	# Hierarchical clustering of samples
	#compute correlation coefficient on a per-sample basis
	cormatrix_samples <- cor(t(expressionData), use = "pairwise.complete.obs") 
	
	## print(paste("clusterSamples: dim(cormatrix_samples)[1]", dim(cormatrix_samples)[1], 
	##                 " any NAs?", any(is.na(cormatrix_samples))))
	
	distmatrix_samples = 1-cormatrix_samples
	hierclu_samples <- flashClust(as.dist(distmatrix_samples), method="average") 
	colcolors_samples = .moduleDetectByFixedClusterno(ihcluster=hierclu_samples, nocluster=nocluster)
	results<-NULL
	results$sampleTree<-hierclu_samples
	results$modules<-colcolors_samples
	results
}

# cut dendrogram to get the required number of clusters
#
.moduleDetectByFixedClusterno <- 
		function(ihcluster, nocluster=2) 
{
	maxhei   = max(ihcluster$height)
	minhei   = min(ihcluster$height)
	curheightcutoff = minhei
	curno = 0
	cnt=0
	while(curno != nocluster){
		cnt = cnt + 1
		if (cnt > 6)
			break
		curheightcutoff = (curheightcutoff + maxhei)/2.0
		colcode  = moduleDetectLabel(hiercluster=ihcluster, curheightcutoff, minsize1=1)
		colcode  = as.character(colcode)
		# we need as.character(colcode) to deal with the case of grey is zero
		#table(colcolors)
		#turquoise      blue      grey 
		#12             6         0 
		colr = as.integer(table(colcode))
		curno    = sum(colr>0)
	}
	colcode   
}

# defined in package 'sma', not available in cran
.rgcolorsFunc<-
function(n)
{
	paste("#", as.hexmode(round((0:(n-1))*255/(n-1))), as.hexmode(round(((n-1):0)*255/(n-1))), "00", sep="")
}

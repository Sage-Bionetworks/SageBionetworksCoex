# Compute topological overlap matrix and gene dendrogram given a data frame of normalized gene expression data.
#
# Note:  This is both the most computationally intensive and the most rote part of coexpression analysis.
# Therefore it is broken out as an independent function.  Variations on the downstream steps can then
# be run without repeating this step.
#
# Inputs:
# expression data -  a data frame with rows as samples and columns as genes/probes
# R's limit on array size constrains the number of genes to 46340 (~sqrt(2^31))
#
# beta = NULL causes the function to select the optimal beta between 1-12
#
# RsquaredCut = the threshold for linearity of the 'scale-free' (log(p(k)) vs. log(k)) plot 
#	after for accepting a candidate value of beta.  The function returns the minimum beta
# 	meeting this criterion.
#	TODO: include WGCNA reference
#
# Outputs:
# beta - the 'soft threshold' beta
# sftStatistics - a data frame having the statistics collected while determining beta, or NULL if beta is passed in
# 	TODO:  need reference for the stat's
# dichotCor - |cc|^beta but with zero diagonal
# tomDist - the topological ovelap matrix (TOM) as a distance metric
# geneTree - the gene dendrogram
#
clusterGenes <- 
		function(
				expressionData, 
				beta=NULL, 
				RsquaredCut = 0.9) 
{
	checkVersion(2,13)	
	
	results<-list()
	
	if (is.null(beta)) {
		## this function has a bug
		## sft = pickSoftThreshold.fromSimilarity(corMat, 
		##         RsquaredCut = RsquaredCut,
		##         powerVector = seq(1, 12, by = .5),
		##         blockSize = dim(expressionData)[2]+1) # i.e. use just one block
		
		# compute the scale-free exponent, beta
		sft = pickSoftThreshold(expressionData, 
				RsquaredCut = RsquaredCut, 
				powerVector = seq(1, 12, by = .5), 
				blockSize = dim(expressionData)[2]+1, # i.e. use one block
				corFnc = "cor")
		collectGarbage()
		
		sftStatistics <- sft$fitIndices
		
		# The following is adapted from Bin Zhang and Steve Horvath's 'powercut()' in R-tomfunctions.R
		#
		#in case that the largest R^2 is not bigger than RsquaredCut,
		#we have to reduce it gradually
		corrlcutoff=sft$powerEstimate
		msqrcut = RsquaredCut
		while(is.na(corrlcutoff)){
			ind1   = sftStatistics[,2] > msqrcut
			indcut = NA
			indcut = ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
			if (!is.na(indcut)) corrlcutoff = sftStatistics[indcut,1]
			msqrcut     = msqrcut - 0.05
		}
		
		beta <- corrlcutoff
		
		# no longer need sft
		rm(sft)
	} else {
		sftStatistics <- NULL
	}
	results$beta<-beta
	results$sftStatistics<-sftStatistics
	
	spaceTimeStats<-spaceTime("clusterGenes: 'pickSoftThreshold' done")	
	
	# Compute the gene-gene correlation matrix
	corMat <- cor(expressionData, use = "pairwise.complete.obs") 
	
	spaceTimeStats<-spaceTime("clusterGenes: 'cor' done", spaceTimeStats)
	
	results$dichotCor  <- abs(corMat)^beta
	
	spaceTimeStats<-spaceTime("clusterGenes: 'dichotCor' done", spaceTimeStats)	
	
	# don't need 'corMat' (which may be huge) anymore
	rm(corMat)
	collectGarbage()
	spaceTimeStats<-spaceTime("clusterGenes: rm(corMat)", spaceTimeStats)	
	
	# Compute the topological overlap matrix (TOM)
	# The UCLA library assumes '1's rather than '0's along the diagonal...

    # results$tomDist <- TOMdist(results$dichotCor)
	# the following duplicates the one-line function TOMdist in 'blockwiseModulesC.R'.  
	# We do this simply to avoid the extra function call when passing a large argument
	results$tomDist <- 1-TOMsimilarity(results$dichotCor)
	
	collectGarbage()
	spaceTimeStats<-spaceTime("clusterGenes: 'TOMdist' done", spaceTimeStats)	
	
	# ... but the formal definition calls for zeros along the diagonal
	diag(results$dichotCor)<- 0

	# i.e. 'k'
	results$connectivity<- apply(results$dichotCor,2,sum, na.rm=TRUE) 
	
	collectGarbage()
	
	# Compute the dendrogram
	geneTree <- flashClust(as.dist(results$tomDist),method="average")
	results$geneTree<-geneTree
	
	spaceTimeStats<-spaceTime("clusterGenes: 'flashClust' done (final step)", spaceTimeStats)	
	results$spaceTimeStats<-spaceTimeStats
	
	results
}
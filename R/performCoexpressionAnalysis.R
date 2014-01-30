
# Input:
# expression data -  a data frame with rows as samples and columns as genes/probes
# R's limit on array size constrains the number of genes to 46340 (~sqrt(2^31))
#
# beta = NULL causes the function to select the optimal beta between 1-12
#
# dynamicCutMethod = "tree" or "hybrid".  The former is the 'classic' approach used at Sage,
# while the latter is the default for the UCLA-WGCNA approach.  The details are here:
# http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/BranchCutting/Supplement.pdf
#
# mergeModuleMethod = "conn" or "eigen".  The former is the 'classic' approach used at Sage,
# in which a module's representative gene is the most highly CONNected one.  The latter is
# the default for the UCLA-WGCNA approach, in which the representative gene is the module's
# EIGENvector.
#
# Output:
# beta - the 'soft threshold' beta
# sftStatistics - a data frame having the statistics collected while determining beta, or NULL if beta is passed in
# 	TODO:  need reference for the stat's
# dichotCor - |cc|^beta but with zero diagonal
# tomDist - the topological ovelap matrix (TOM) as a distance metric
# geneTree - the gene dendrogram

#
performCoexpressionAnalysis <-
		function(
			expressionData, 
			beta=NULL, 
			dynamicCutMethod="tree",
			mergeModuleMethod="conn") 
{
	checkVersion(2,13)	
	
	spaceTimeStats<-spaceTime("performCoexpressionAnalysis: Start")
	
	results<-clusterGenes(expressionData, beta)
	
	spaceTimeStats<-rbind(spaceTimeStats, results$spaceTimeStats)
	spaceTimeStats<-spaceTime("performCoexpressionAnalysis: clusterGenes done", spaceTimeStats)	
	
	mfgtResults<-modulesFromGeneTree(
			geneTree=results$geneTree, 
			expressionData=expressionData,
			dichotCor=results$dichotCor)
	
	spaceTimeStats<-spaceTime("performCoexpressionAnalysis: modulesFromGeneTree done", spaceTimeStats)	
	
	results$geneModules<-mfgtResults$geneModules
	results$genePCTree<-mfgtResults$genePCTree
	
	## sampleClusterResults<-clusterSamples(expressionData)
	## results$sampleTree<-sampleClusterResults$sampleTree
	## results$sampleModules<-sampleClusterResults$modules
	
	imsResults<-intraModularStatistics(
			results$dichotCor, results$tomDist, results$geneModules)
	
	results$intraModularStatistics <-imsResults$intraModularStatistics
	spaceTimeStats<-rbind(spaceTimeStats, imsResults$spaceTimeStats)
	
	spaceTimeStats<-spaceTime("performCoexpressionAnalysis: all done", spaceTimeStats)	
	results$spaceTimeStats <- spaceTimeStats
	
	results
}





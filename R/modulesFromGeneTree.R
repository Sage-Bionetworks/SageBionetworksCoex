# Input:
# geneTree - the gene dendrogram
# expressionData - the expression data matrix
# dichotCor - |CC|^beta (with zero diagonal), only required if mergeModuleMethod=="conn"
# tomDist - the topological overlap matrix (TOM) as a distance metric, only required if dynamicCutMethod=="hybrid"
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
# a vector of module memberships, i.e. output[i] is the module to which gene i belongs
# "grey" is a reserved name for genes which belong to no module
#
modulesFromGeneTree <- 
		function(
				geneTree,
				expressionData,
				dichotCor=NULL,
				tomDist=NULL,
				dynamicCutMethod="tree",
				mergeModuleMethod="conn") 
{
	checkVersion(2,13)	
	
	# First pass at cluster definition, using 'dynamic tree cutting'
	heightcutoff    = 0.99
	minModuleSize   = 30 # modules must have this minimum number of genes
	
	if (dynamicCutMethod=="tree") {
		# WGCNA using param's that match the Sage algorithm
		dtcModules <- cutreeDynamic(
				dendro = geneTree, 
				cutHeight = heightcutoff, 
				minClusterSize = minModuleSize, 
				method ="tree")	
	} else if (dynamicCutMethod=="hybrid"){
		dtcModules <- cutreeDynamic(
				dendro = geneTree, 
				cutHeight = heightcutoff, 
				minClusterSize = minModuleSize,
				distM = tomDist,
				deepSplit = 2, 
				pamRespectsDendro = FALSE)		
	} else {
		stop(paste("Expected 'tree' or 'hybrid' but found ", dynamicCutMethod, ".", sep=""))
	}
	# convert labels to colors
	colorModules <- labels2colors(dtcModules)
	
	mergedColors <- NULL
	genePCTree <- NULL
	
	# Merge close modules
	if (mergeModuleMethod=="eigen") {
		MEDissThres = 0.2
		merge = mergeCloseModules(expressionData, colorModules, cutHeight = MEDissThres)
		# The merged module colors
		mergedColors <- merge$colors
		genePCTree <- merge$dendro
		# Eigengenes of the new merged modules: merge$newMEs
	} else if (mergeModuleMethod=="conn") {
		pre_nomodules= length(colorModules)+1
		# we repeat the following until no more modules are combined
		while(pre_nomodules!=length(colorModules)){
			pre_nomodules= length(colorModules)
			
			# take the profile of the most connected gene in each module as PC 
			pcs=moduleHubProfiles(datexpr  =as.matrix(expressionData), 
					adjmatrix=dichotCor, 
					couleur=as.character(colorModules), 
					min_modulesize=10, 
					myheightcutoff=heightcutoff,
					h1row = geneTree, 
					myminModuleSize=minModuleSize,
					mydeepSplit=FALSE)
			
			distCorPCs = 1-abs(cor(pcs[[1]],use="p"))
			distCorPCs = ifelse(is.na(distCorPCs), 0, distCorPCs)
			genePCTree = flashClust(as.dist(distCorPCs),method="a")
			
			# clusters of PCs
			pcheicutoff=0.5
			pccolcode   = moduleDetectLabel(hiercluster=genePCTree, pcheicutoff, minsize1=1)
			
			#check whether only "grey" module is detected
			# we can only merge modules if some PCs have been clustered together
			pcmodules   = names(table(pccolcode))
			noClusters =  (length(pcmodules)==1) & (pcmodules[1]=="grey")
			if(!noClusters){
				colorModules <- mergeClusterByPCA(
						mdendro= geneTree, 
						genecluster = colorModules,
						pccluster= pccolcode,
						pcnames = names(pcs[[1]]))
			}
		} # end of 'while' loop
		mergedColors <- colorModules
	} else {
		stop(paste("Expected 'eigen' or 'conn' but found ", mergeModuleMethod, ".", sep=""))
	}
	
	collectGarbage()
	
	results<-NULL
	results$genePCTree<-genePCTree
	results$geneModules<-mergedColors
	results
}
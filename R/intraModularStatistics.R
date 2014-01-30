# dichotCor |cc|^beta, with zeros on the diagonal
#
# tomDist - the topological overlap matrix encoding of 'dichotCor', converted to a distance metric
#
# modules is a vector of categories, i.e. the category for item i is modules[i]
# 	if modules[i]=="grey" then i is not in any category
#
intraModularStatistics <- 
		function(
				dichotCor, 
				tomDist, 
				modules) 
{
	no.genes<-dim(dichotCor)[1]
	if (no.genes!=length(modules)) {
		stop(paste("dim(dichotCor)[1] (", no.genes, ") must equal length(modules) (", length(modules), ").", sep=""))
	}
	# Network statistics
	# in-module cluster coefficients, per gene
	cmCCResults = computeModuleCC(dichotCor, modules)
	cc.in = cmCCResults$ccin
	spaceTimeStats<-cmCCResults$spaceTimeStats

	spaceTimeStats<-spaceTime("IMS: computeModuleCC done", spaceTimeStats)
	
	k.in  = computeModuleLinks(dichotCor, modules)
	nolinks <- apply(dichotCor,2,sum, na.rm=TRUE)
	k.out = nolinks - k.in
	k.diff= k.in - k.out
	k.all = nolinks
	
	#sum of topological overlaps within modules
	to.in = computeModuleLinks(tomDist, modules, isAdjacency=FALSE)
	
	spaceTimeStats<-spaceTime("IMS: TOM-CML done", spaceTimeStats)	
	
	#sum of topological overlaps 
	to.all       = no.genes-apply(tomDist, 2, sum, na.rm=TRUE)
	to.out = to.all - to.in
	to.diff= to.in  - to.out
	
	
	#-------------------   normalized measures   ---------------------------
	k.in.normed = computeModuleLinks(dichotCor, modules, normalized=TRUE)
	k.all.normed= k.all/length(nolinks)
	
	spaceTimeStats<-spaceTime("IMS: DC-CML done", spaceTimeStats)	
	
	#sum of topological overlaps within modules
	to.in.normed     = computeModuleLinks(tomDist, modules, isAdjacency=FALSE, normalized=TRUE)
	
	spaceTimeStats<-spaceTime("IMS: TOM#2-CML done", spaceTimeStats)	
	
	#sum of topological overlaps 
	to.all.normed = to.all/no.genes
	
	# The following data frame contains all network properties:
	datNetwork=data.frame(k.in,  k.in.normed,  to.in,  to.in.normed,
			k.all, k.all.normed, to.all, to.all.normed,
			k.out, k.diff, to.out, to.diff,
			cc.in)
	
	results<-list()
	results$spaceTimeStats<-spaceTimeStats
	results$intraModularStatistics<-datNetwork
	results
}
require("SageBionetworksCoex") || stop("unable to load SageBionetworksCoex package")

#dyn.load(file.path(paste("../SageBionetworksCoex/libs/SageBionetworksCoex", .Platform$dynlib.ext, sep="")))


require("flashClust")
# test various options for the package
# note, we combine these into one script so that they can share 
# the results of the initial time-consuming steps

# read in the data
inputFile<-file.path("../../SageBionetworksCoex/data", "expressionData.txt")
## temp <- read.delim(file=inputFile, sep="\t");
## allMatrix <- temp[-1];
## rownames(allMatrix) <- t(temp[1]);

allMatrix<-read.table(inputFile, header=T, row.names=1)

# note, we transpose 'allMatrix' so rows are samples and col's are genes
expressionData = t(allMatrix)

# we can use the following for various 'tree cutting' choices
analysisResults<-clusterGenes(expressionData)

# there are four different combinations of 'dynamicCutMethod' and 'mergeModuleMethod'
mfgtResults<-modulesFromGeneTree(geneTree=analysisResults$geneTree, 
		expressionData=expressionData,
		dichotCor=analysisResults$dichotCor, 
		tomDist=analysisResults$tomDist, dynamicCutMethod="tree", mergeModuleMethod="conn")

if (is.null(mfgtResults$genePCTree)) stop("Null genePCTree for dynamicCutMethod=tree and mergeModuleMethod=conn")

# will use these values in plotting
analysisResults$geneModules<-mfgtResults$geneModules
analysisResults$genePCTree<-mfgtResults$genePCTree


# make sure other combinations of parameters work
mfgtResults<-modulesFromGeneTree(geneTree=analysisResults$geneTree, 
        expressionData=expressionData,
        dichotCor=analysisResults$dichotCor, 
        tomDist=analysisResults$tomDist, dynamicCutMethod="tree", mergeModuleMethod="eigen")

if (is.null(mfgtResults$genePCTree)) stop("Null genePCTree for dynamicCutMethod=tree and mergeModuleMethod=eigen")

mfgtResults<-modulesFromGeneTree(geneTree=analysisResults$geneTree, 
        expressionData=expressionData,
        dichotCor=analysisResults$dichotCor, 
        tomDist=analysisResults$tomDist, dynamicCutMethod="hybrid", mergeModuleMethod="conn")

if (is.null(mfgtResults$genePCTree)) stop("Null genePCTree for dynamicCutMethod=hybrid and mergeModuleMethod=conn")

# Note:  This variation in parameters causes the PCTree to have just two branches,
# which causes the tree plotting function to fall over
mfgtResults<-modulesFromGeneTree(geneTree=analysisResults$geneTree, 
        expressionData=expressionData,
        dichotCor=analysisResults$dichotCor, 
        tomDist=analysisResults$tomDist, dynamicCutMethod="hybrid", mergeModuleMethod="eigen")

if (is.null(mfgtResults$genePCTree)) stop("Null genePCTree for dynamicCutMethod=hybrid and mergeModuleMethod=eigen")


print(paste("genePCTree:", mfgtResults$genePCTree))

## sampleClusterResults<-clusterSamples(expressionData)
## analysisResults$sampleTree<-sampleClusterResults$sampleTree
## analysisResults$sampleModules<-sampleClusterResults$modules
## 
analysisResults$intraModularStatistics<-intraModularStatistics(
		analysisResults$dichotCor, analysisResults$tomDist, analysisResults$geneModules)


outputDir<-tempdir()

if(!file.exists(outputDir)){
	dir.create(outputDir)
}

# now try all the types of image files
for (ft in c("jpg", "ps", "png", "pdf")) {
	fileList<-createDiagnosticPlots(
			filetype=ft,
			outputDir=outputDir,
			expressionData,
			analysisResults$dichotCor,
			analysisResults$tomDist,
			analysisResults$connectivity,
			analysisResults$beta,
			analysisResults$geneTree,
			analysisResults$geneModules,
			analysisResults$genePCTree)
}
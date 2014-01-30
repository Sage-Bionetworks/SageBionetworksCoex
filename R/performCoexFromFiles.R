
# 8 images, see 'createDiagnosticsPlots' for details
# powertable.txt	tab-delimited table of the regression parameters for the range of values of beta that were checked for the network's "scale-free" property.
# tommodules.txt	a tab delimted table combining (1) the input data, (2) network statistics, on a per-gene basis, and (3) the genes' module membership.
performCoexFromFiles <- function(
		inputFile, 
		outputDir=NULL,
		imageFileType="jpg",
		beta=NULL) 
{
	checkVersion(2,13)	
	
	spaceTimeStats<-spaceTime("performCoexFromFiles: Start")
	
	## temp <- read.delim(file=inputFile, sep="\t");
	## allMatrix <- temp[-1];
	## rownames(allMatrix) <- t(temp[1]);
	
	allMatrix <- read.table(file=inputFile, header=TRUE, row.names=1)
	
	
	if (nrow(allMatrix)<3) stop(paste("Number of data rows in input matrix is only ",nrow(allMatrix), "."))
	if (ncol(allMatrix)<3) stop(paste("Number of data columns in input matrix is only ",ncol(allMatrix), "."))
	
	if (is.null(outputDir)) {
		outputDir<-tempdir()
	}
	
	if(!file.exists(outputDir)){
		dir.create(outputDir)
	}
	
	spaceTimeStats<-spaceTime("performCoexFromFiles: Data read in", spaceTimeStats)
	
	# note, we transpose 'allMatrix' so rows are samples and col's are genes
	expressionData = t(allMatrix)
	analysisResults <- performCoexpressionAnalysis(expressionData, beta=beta)
	
	## if (class(analysisResults)=='try-error') {
	##     print(traceback()) # doesn't work, I get "No traceback available"
	##     return
	## }
	
	geneTree<-analysisResults$geneTree
	
	spaceTimeStats<-rbind(spaceTimeStats, analysisResults$spaceTimeStats)
	spaceTimeStats<-spaceTime("performCoexFromFiles: Analysis done", spaceTimeStats)
	
	createDiagnosticPlotsResults<-createDiagnosticPlots(
			filetype=imageFileType,
			outputDir=outputDir,
			expressionData,
			analysisResults$dichotCor,
			analysisResults$tomDist,
			analysisResults$connectivity,
			analysisResults$beta,
			analysisResults$geneTree,
			analysisResults$geneModules,
			analysisResults$genePCTree)
	fileList<-createDiagnosticPlotsResults
	
	spaceTimeStats<-spaceTime("performCoexFromFiles: createDiagnosticPlots done", spaceTimeStats)
	
	EXT_NAME<-"txt"
	POWER_TABLE<-"PowerTable"
	NETWORK<-"Network"
	GENE_TREE<-"GeneTree"
	SPACE_TIME_STATS<-"SpaceTimeStats"
	
	filePath<-file.path(outputDir, paste(POWER_TABLE, EXT_NAME, sep="."))
	fileList$POWER_TABLE <- filePath
	write.table(analysisResults$sftStatistics, filePath, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

	filePath<-file.path(outputDir, paste(NETWORK, EXT_NAME, sep="."))
	fileList$NETWORK <- filePath
	network <- createNetwork(expressionData, analysisResults$geneModules, analysisResults$intraModularStatistics)
	write.table(network, filePath, sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
	
	filePath<-file.path(outputDir, paste(GENE_TREE, EXT_NAME, sep="."))
	fileList$GENE_TREE <- filePath
	geneTree<-analysisResults$geneTree
	save(geneTree, file=filePath, ascii=TRUE)
	
	spaceTimeStats<-spaceTime("performCoexFromFiles: all done", spaceTimeStats)
	
	filePath<-file.path(outputDir, paste(SPACE_TIME_STATS, EXT_NAME, sep="."))
	fileList$SPACE_TIME_STATS <- filePath
	save(spaceTimeStats, file=filePath, ascii=TRUE)
	
	results<-list()
	results$fileList<-fileList
	results$spaceTimeStats<-spaceTimeStats
	return(results)
}
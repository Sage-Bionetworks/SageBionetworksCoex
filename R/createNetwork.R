createNetwork <-
		function(expressionData,
				geneModules,
				intraModularStatistics)
{
	checkVersion(2,13)	
	
	sd2 <- apply(expressionData,2, sd, na.rm=T) # stdev columns
	mod.frame <- data.frame(sd2, geneModules)
	colnames(mod.frame) <- c("standard.deviation", "module")
	network <- cbind(t(expressionData), intraModularStatistics, mod.frame)
	
	network
}
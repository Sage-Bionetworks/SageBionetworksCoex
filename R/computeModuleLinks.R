#-------------------------------------------------------------------------
#Function: compute within-module number of connections for each gene
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  in-module number of connections of each gene
#-------------------------------------------------------------------------
computeModuleLinks = function(adjMatrix, colorcodeC, isAdjacency=TRUE, normalized=FALSE, usegreymodule=F)
{
	modnames= names( table(colorcodeC) )
	
	no.genes     = dim(adjMatrix)[1]
	links        = rep(0, no.genes)
	idxseq       = c(1:no.genes)   
	for (each in modnames ){
		if((usegreymodule==F) & (as.character(each)=="grey" ) ){
			next
		}
		
		whichmod    = each==colorcodeC
		module.size = sum(whichmod)
		
		if (module.size==1){
			next
		}
		
		modk <- apply(adjMatrix[whichmod,whichmod],2,sum, na.rm=TRUE) 
		if(!isAdjacency){
			modk <- (module.size -modk)
		}
		
		#normalize against the module size
		if(normalized==TRUE){
			modk = modk/module.size
		}
		
		#indices of genes in the current module
		idxmod = idxseq * whichmod
		
		#put the links's to the buffer      
		links[idxmod] = modk
		
		collectGarbage()
	}
	links
}


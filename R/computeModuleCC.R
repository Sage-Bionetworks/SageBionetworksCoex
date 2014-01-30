#-------------------------------------------------------------------------
#Function: compute within-module cluster coefficient for each gene
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  in-module cluster coefficients
#-------------------------------------------------------------------------
computeModuleCC <-
function(adjMatrix, colorcodeC)
{
	spaceTimeStats<-spaceTime(paste("cMCC:n=", dim(adjMatrix)[1]))
   modnames= names( table(colorcodeC) )

   no.genes     = dim(adjMatrix)[1]
   clustercoeff = rep(0, no.genes)
   idxseq       = c(1:no.genes)   
   for (each in modnames ){
      if(as.character(each)=="grey" ){
         next
      }
      whichmod = each==colorcodeC
      if( sum(whichmod)>1){
		 n<-length(colorcodeC[whichmod])
		 spaceTimeStats<-spaceTime(paste("cMCC:",each,"n=",n), spaceTimeStats)
		 cc<-computeClusterCoefficient(adjMatrix[whichmod,whichmod])
         icc <-cc$cc
		 spaceTimeStats<-rbind(spaceTimeStats, cc$spaceTimeStats)
      }else{
         icc = 1
      }

      #indices of genes in the current module
      idxmod = idxseq * whichmod
      
      #put the cc's to the
      clustercoeff[idxmod] = icc
	  
	  collectGarbage()
   }
   return(list(ccin=clustercoeff, spaceTimeStats=spaceTimeStats))
}


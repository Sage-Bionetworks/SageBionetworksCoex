

tomHeatmap <-
		function(filePath,
				tomDist,
				geneDendrogram,
				geneModules

)
{
	no.genes = dim(tomDist)[1]
	
	options(expressions = 10000)
	openImgDev(filePath)
	par(mfrow=c(1,1))
	#we need reverse the distmatrix and row labels
	orderedDist1  = tomDist[geneDendrogram$order,geneDendrogram$order] 
	labelcol <- as.character(geneModules)
	if(no.genes<=6000) {
		diag(orderedDist1)= 0
		revDist1=apply(orderedDist1, 2, rev)
		revDist1= revDist1*revDist1
		heatmap(revDist1,Rowv=NA,Colv=NA, scale="none",revC=F,
				ColSideColors=as.character(labelcol[geneDendrogram$order]),
				RowSideColors=as.character(rev(labelcol[geneDendrogram$order])) )# , col = rev(rgcolors.func(50)))
	} else{
		fixedSize     = 3000
		fullColorCol  = as.character(labelcol[geneDendrogram$order])
		fold          = as.integer(no.genes/fixedSize)
		sampleSeq     = seq(1, no.genes, fold)
		realSampSize  = length(sampleSeq)
		sampleColorCol= fullColorCol[sampleSeq]
		revDist1      = orderedDist1[sampleSeq, sampleSeq]
		
		diag(revDist1)= 0
		revDist1      = apply(revDist1, 2, rev)    
		
		revDist1= revDist1*revDist1
		heatmap(revDist1,Rowv=NA,Colv=NA, scale="none",revC=F,
				ColSideColors=as.character(sampleColorCol),
				RowSideColors=as.character(rev(sampleColorCol)))# col = rev(rgcolors.func(50)) )
		
	}


	dev.off()
	
	collectGarbage()
}



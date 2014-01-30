


corHeatmap <-
		function(filePath,
				dichotCor,
				geneDendrogram,
				geneModules
				)
{		
	orderedDist1  = abs(dichotCor)[geneDendrogram$order,geneDendrogram$order] 
	no.genes = dim(dichotCor)[1]
	labelcol <- as.character(geneModules)
	
	openImgDev(filePath)
	par(mfrow=c(1,1))
	if(no.genes<=6000) {
		diag(orderedDist1)= 1
		
		revDist1=apply(orderedDist1, 2, rev)
		revDist1=1-as.matrix(revDist1)
		revDist1= revDist1*revDist1
		heatmap(revDist1,Rowv=NA,Colv=NA, scale="none",revC=F,
				ColSideColors=as.character(labelcol[geneDendrogram$order]),
				RowSideColors=as.character(rev(labelcol[geneDendrogram$order]))) # col = rgcolors.func(50) )
		
	} else{
		fixedSize     = 1000
		fullColorCol  = as.character(labelcol[geneDendrogram$order])
		fold          = as.integer(no.genes/fixedSize)
		sampleSeq     = seq(1, no.genes, fold)
		realSampSize  = length(sampleSeq)
		sampleColorCol= fullColorCol[sampleSeq]
		
		revDist1      = orderedDist1[sampleSeq, sampleSeq]
		diag(revDist1)= 1
		
		revDist1      = apply(revDist1, 2, rev)    
		
		revDist1=1-as.matrix(revDist1)
		revDist1= revDist1*revDist1
		heatmap(revDist1,Rowv=NA,Colv=NA, scale="none",revC=F,
				ColSideColors=as.character(sampleColorCol),
				RowSideColors=as.character(rev(sampleColorCol)))
	}
	
	dev.off()
	
	collectGarbage()
}
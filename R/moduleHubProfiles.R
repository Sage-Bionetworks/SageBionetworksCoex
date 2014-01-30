moduleHubProfiles <-
function(datexpr, adjmatrix, couleur, min_modulesize=10, myheightcutoff, h1row, myminModuleSize, mydeepSplit) {

  no.pcs=10

  modlevels= names(table(couleur)) #levels(factor(couleur))

  #the first no.pcs of elements in the list are PCs and the last one is the variations
  listPCs = as.list(rep(NULL, no.pcs+1) )

  for (i in c(1:no.pcs) ){
    listPCs[[i]]     = data.frame(matrix(0,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
    colnames(listPCs[[i]]) <- modlevels
  }

  listPCs[[ no.pcs+1] ]= data.frame(matrix(0,nrow=no.pcs,ncol= length(modlevels))) #varexplained
  colnames(listPCs[[ no.pcs+1] ]) <- modlevels


    # take the profile of the most connected gene in each module as PC
    #
	# modified to use WGCNA version, see below
	## colcode.reduced=cutreeDynamic(hierclust=h1row, deepSplit=mydeepSplit,
	##         maxTreeHeight      =myheightcutoff, 
	##         minModuleSize      =myminModuleSize,
	##         nameallmodules=F, useblackwhite=F)
	
	colcode.reduced <- cutreeDynamic(
			dendro = h1row, 
			cutHeight = myheightcutoff, 
			minClusterSize =myminModuleSize, 
			method ="tree")
	
    colorI = as.character(colcode.reduced)
    kin    = computeModuleLinks(adjmatrix, couleur)

    # find the hub in each module 
    #
    kinIdxed= cbind(kin, c(1:length(couleur)), couleur)
    orderK  = order(-kin)
    kinIdxed= kinIdxed[orderK, ]
    orderK  = order(kinIdxed[,3])
    kinIdxed= kinIdxed[orderK, ]
    
    hubIdx    = rep(0, length(modlevels) )
    for(z in c(1:length(modlevels)) ) {        
        isel      = modlevels[z] == kinIdxed[,3]
        ikinIdxed = kinIdxed[isel, ]

        # extract hubs' profiles
        #
        listPCs[[1] ][,z] = datexpr[,as.integer(ikinIdxed [1,2])]
        hubIdx[z] = ikinIdxed [1,2]
    }

    # extract hubs' profiles
    #
    #listPCs[[1]] = t(datexpr[,as.integer(hubIdx)])
    names(listPCs[[1]]) <- modlevels

    return(listPCs)
}


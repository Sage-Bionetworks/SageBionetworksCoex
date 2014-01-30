computeClusterCoefficient <-
function(adjMatrix, weighted=T) {
		spaceTimeStats<-spaceTime("cCC: start")
        no.genes <- dim(adjMatrix)[1]
        nolinksNeighbors <- c(rep(-666,no.genes))
        total.edge <- c(rep(-666,no.genes))

        #for (i in 1:no.genes){
        #     nolinksNeighbors[i] <-  adjMatrix[i,] %*% adjMatrix %*% adjMatrix[,i]
        #     #total.edge[i] <-  adjMatrix[i,] %*% Pmax %*% adjMatrix[,i]
        #}
	
	
        #nolinksNeighbors <- apply(adjMatrix, 1, .computeLinksInNeighbors, imatrix=adjMatrix)
		nolinksNeighbors <- .ataat(adjMatrix)
		
		spaceTimeStats<-spaceTime("cCC: cLiN", spaceTimeStats)
        plainsum  <- apply(adjMatrix, 1, sum)
        if(weighted) {
           squaresum <- apply(adjMatrix^2, 1, sum)
           total.edge = plainsum^2 - squaresum
        }else{ # for unweighted network, 1^2 = 1
           total.edge = plainsum^2 - plainsum
        }

        # in case of single node, this will not affect the CC computation
        #
        total.edge = ifelse(total.edge==0, 1, total.edge)

        cluster.coef = nolinksNeighbors/total.edge
        cluster.coef = ifelse(total.edge>0,cluster.coef, 0) 

		spaceTimeStats<-spaceTime("cCC: end", spaceTimeStats)
        list(cc=cluster.coef, spaceTimeStats=spaceTimeStats)
}

.computeLinksInNeighbors <- function(x, imatrix)
{
    y= x %*% imatrix %*% x
    y
}

# compute diag(at*a*at), where a is square AND SYMMETRIC
.ataat<-function(a) {
	if (dim(a)[1]!=dim(a)[2]) stop("Matrix is not square.")
	n = dim(a)[1]
	# compute a*at
	# this is dramatically accelerated by calling a compiled function (2 lines below) rather than the following
	# aat <- a %*% t(a)
	aat<-array(0, dim(a))
	result<-.C("aat", as.double(as.matrix(a)), as.integer(n), aat=as.double(as.matrix(aat)))
	aat<-array(result$aat, dim(a))
	# comes back lower triangular, so...
	aat[utIndices(n)]=aat[ltIndices(n)] # could have done this in the 'C' routine...
	# put at and a*t(a) into a 3-D matrix
	z<-array(c(t(a),aat), c(dim(a),2))
	# for each row in a (and the matching one in aat) compute the dot product
	apply(z,2,.dot)
	#diag(z[,,1] %*% z[,,2])
}

# return the lower triangular indices of an nXn matrix (omitting the diagonal)
ltIndices<-function(n) {unlist(apply(as.matrix(1:(n-1)), 1, .ltIntern, n=n))}
.ltIntern<-function(i,n){(n*(i-1)+i+1):(n*i)}

# return the upper triangular indices of an nXn matrix (omitting the diagonal)
utIndices<-function(n) {unlist(apply(as.matrix(1:(n-1)), 1, .utIntern, n=n))}
.utIntern<-function(i,n){seq(i*(n+1), (n-1)*n+i, n)}


# compute the dot product of the two columns.  dim(x)[2] must be 2
.dot<-function(x){
	if (dim(x)[2] != 2) stop(paste("second dim is ", dim(x)[2]))
	x[,1]%*%x[,2]
}





testATAAT<-function(a) {
	ans1<-apply(a, 1, SageBionetworksCoex:::.computeLinksInNeighbors, imatrix=a)
	print(paste("ans1=", paste(ans1, collapse=" ")))
	ans2<-SageBionetworksCoex:::.ataat(a)
	print(paste("ans2=", paste(ans2, collapse=" ")))	
	if (any(ans1!=ans2)) stop(paste("Unequal answers: ", paste(ans1-ans2, collapse=" ")))
}

runAllTests<-function() {
	# symmetric case
	print(testATAAT(array(c(1,2,1,2,2,3,1,3,4), c(3,3))))
	# random, symmetric case
	print(testATAAT(cor(replicate(50, rnorm(100)))))
	# assymetric case
	print(testATAAT(array(1:9, c(3,3))))
}

#dyn.load(file.path(paste("../SageBionetworksCoex/libs/SageBionetworksCoex", .Platform$dynlib.ext, sep="")))

runAllTests()
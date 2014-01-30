require("SageBionetworksCoex") || stop("unable to load SageBionetworksCoex package")

#dyn.load(file.path(paste("../SageBionetworksCoex/libs/SageBionetworksCoex", .Platform$dynlib.ext, sep="")))

# test on 1000 genes by 100 samples

## require(WGCNA)
## require(lattice)
inputFileName<-"expressionData.txt"
inputFile<-file.path("../../SageBionetworksCoex/data", inputFileName)
ans<-performCoexFromFiles(inputFile)
outputFile <- ans$fileList$NETWORK

inputFrame <- read.table(inputFile, header=T, row.names=1)
outputFrame <- read.table(outputFile, header=T, row.names=1)

cat("Input (1st 5 rows, cols):\n")
inputFrame[1:5,1:5]

cat("Output (1st 5 rows, cols):\n")
outputFrame[1:5,1:5]

# check that the first part of 'outputFrame' echos the input
if (!all(colnames(outputFrame)[1:dim(inputFrame)[2]]==colnames(inputFrame))) stop(paste("For test file", inputFileName, "column names of Network file don't match."))
if (!all(row.names(outputFrame)==row.names(inputFrame))) stop(paste("For test file", inputFileName, "row names of Network file don't match."))
if (!all(outputFrame[,1:dim(inputFrame)[2]]==inputFrame)) stop(paste("For test file", inputFileName, "left columns of Network file don't match."))

# try specifying 'beta'
performCoexFromFiles(inputFile, beta=6)




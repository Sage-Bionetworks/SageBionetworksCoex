\name{performCoexFromFiles}
\alias{performCoexFromFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Coexpression From Files
}
\description{
Performs coexpression analysis on normalized expression data from a tabular input file.  Creates diagnostic plots
and output data files for the coexpression network and scale free regression statistics.  The output files go into the 
specified directory or, if a directory is not specified, into the session's temporary directory.
}
\usage{
performCoexFromFiles(inputFile, outputDir = NULL, imageFileType = "jpg", beta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFile}{
A normalized matrix of gene expression data, given as a tab-delimited text file whose rows are genes/probes and 
whose columns are samples.  The file has exactly one header row for sample IDs and one header column for gene/probe IDs, and is 
formatted in as expected for R, i.e. the first row has one fewer columns than the remaining rows.
Note: R's limit on array size constrains the number of genes (rows) to 46340.
}
  \item{outputDir}{
The output directory for the created files, or NULL to use the session's temporary directory.  
Existing files in the specified directory are overwritten.
}
  \item{imageFileType}{
One of "jpg", "pdf", "png", "ps", not case sensitive.  Default is "jpg".
}
  \item{beta}{
  Soft thresholding exponent, typically 5-6.  Setting to NULL (the default) causes the function to select the optimal beta between 1 and 12.
}
}
\details{
This function is a wrapper for performCoexpressionAnalysis() and createDiagnosticPlots.  Additionally it creates
tabular files for (1) the scale-free regression statistics (as produced by WGCNA::pickSoftThreshold()), and
(2) the coexpression network.  The latter is the 'network' file, a column-concatenation of (a) the input data, (b) the
intramodular statistics (as produced by intraModularStatistics()), (c) the gene-wise standard deviation
of the input data, and (d) the gene module assignments.
}
\value{
The file handles for the created plots.
\item{ScaleFree}{File handle to "Scale free" plot (P(k) vs. k) for optimum value of beta}
\item{CorHeatmap}{File handle to heatmap of |cc|^beta, ordered according to dendrogram, with color-coded module membership along the margins}
\item{GeneDendrogram}{File handle to gene dendrogram with color-coded modules and color names}
\item{TOMHeatmap}{File handle to heatmap of TOM, rows and columns sorted by gene dendrogram}
\item{GenePCDendrogram}{File handle to dendrogram of clusters of module representative genes}
\item{ExpressionHeatmap}{File handle to heatmap of expression data, organized by gene and sample clusters}
\item{PowerTable}{Scale-free regression statistics (as produced by WGCNA::pickSoftThreshold())}
\item{Network}{The derived network, as described in the 'Details' section.}
\item{GeneTree}{Gene dendrogram, an R object of the form returned by 'hclust'.}
}
\references{
\url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/}

\href{http://www.bepress.com/sagmb/vol4/iss1/art17/}{ Zhang, B. & Horvath, S. (2005) Statistical Applications in Genetics and Molecular Biology 4, Article 17.}

\href{http://bioinformatics.oxfordjournals.org/content/early/2007/11/16/bioinformatics.btm563.full.pdf+html}{Langfelder P, Zhang B, Horvath S (2007) Defining clusters from a hierarchical cluster tree: 
the Dynamic Tree Cut library for R. Bioinformatics 2008 24(5):719-720}

\url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA/}
}
\author{
Bruce Hoff
}
\note{
Note: R's limit on array size constrains the number of genes to 46340.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Coexpression }

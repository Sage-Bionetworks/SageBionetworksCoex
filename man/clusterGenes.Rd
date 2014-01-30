\name{clusterGenes}
\alias{clusterGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cluster Genes
}
\description{
Compute topological overlap matrix and gene dendrogram given a data frame of normalized gene expression data.
}
\usage{
clusterGenes(expressionData, beta = NULL, RsquaredCut = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expressionData}{
A normalized matrix of gene expression data, given as a data frame with rows as samples and columns as genes/probes.
}
  \item{beta}{
  Soft thresholding exponent, typically 5-6.  Setting to NULL (the default) causes the function to select the optimal beta between 1 and 12.
}
  \item{RsquaredCut}{
the threshold for linearity of the 'scale-free' (log(p(k)) vs. log(k)) plot after for accepting a candidate value of beta.  The function returns the minimum beta meeting this criterion.
Note:  If no value of beta causes the R-squared value of the regression to exceed the given threshold, then the threshold is decremented by 0.05 until a solution is found.
}
}

\value{
\item{beta}{The 'soft threshold' beta, or NULL (the default) if the function is to determine what value to use.}
\item{sftStatistics}{A data frame having the statistics collected while determining beta, as produced by WGCNA::pickSoftThreshold, or NULL if beta is passed in.}
\item{dichotCor}{|cc|^beta (elementwise, where 'cc' is the matrix of inter-gene Pearson correlation coefficients) but with zero diagonal.}
\item{tomDist}{The toplogical overlap matrix, in the form of a distance metric.}
\item{geneTree}{The gene dendrogram created from WGCNA::flashClust using 'tomDist' as the distance metric.}
}
\references{
\url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/}

\href{http://www.bepress.com/sagmb/vol4/iss1/art17/}{ Zhang, B. & Horvath, S. (2005) Statistical Applications in Genetics and Molecular Biology 4, Article 17.}

\url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA/}
}
\author{
Bruce Hoff
}
\note{
This is both the most computationally intensive and the most rote part of coexpression analysis.  Therefore it is broken out as an independent function.  Variations on the downstream steps can then be run without repeating this step.
}

\keyword{Coexpression}

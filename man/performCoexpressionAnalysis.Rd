\name{performCoexpressionAnalysis}
\alias{performCoexpressionAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Coexpression Analysis
}
\description{
Performs the analytic (versus graph or I/O) portion of Coexpression analysis.
}
\usage{
performCoexpressionAnalysis(expressionData, beta = NULL, dynamicCutMethod = "tree", mergeModuleMethod = "conn")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expressionData}{
A normalized matrix of gene expression data, given as a data frame with rows as samples and columns as genes/probes.
}
  \item{beta}{
  Soft thresholding exponent, typically 5-6.  Setting to NULL (the default) causes the function to select the optimal beta between 1 and 12.
}
  \item{dynamicCutMethod}{
"tree" or "hybrid".  The former is the 'classic' approach used at Sage Bionetworks, while 
the latter is the default for the UCLA-WGCNA approach.  The details are here:
http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/BranchCutting/Supplement.pdf}
  \item{mergeModuleMethod}{
"conn" or "eigen".  The former is the 'classic' approach used at Sage Bionetworks, in which a 
module's representative gene is the most highly CONNected one.  The latter is the default for the 
UCLA-WGCNA approach, in which the representative gene is the module's EIGENvector.
}
}

\value{
\item{beta}{The 'soft threshold' beta, or NULL (the default) if the function is to determine what value to use.}
\item{sftStatistics}{A data frame having the statistics collected while determining beta, as produced by WGCNA::pickSoftThreshold, or NULL if beta is passed in.}
\item{dichotCor}{|cc|^beta (elementwise, where 'cc' is the matrix of inter-gene Pearson correlation coefficients) but with zero diagonal.}
\item{tomDist}{The toplogical overlap matrix, in the form of a distance metric.}
\item{geneTree}{The gene dendrogram created from WGCNA::flashClust using 'tomDist' as the distance metric.}
\item{geneModules}{
A vector of module memberships, i.e. geneModules[i] is the module to which gene i belongs. 
"grey" is a reserved name for genes which belong to no module.
}
\item{genePCTree}{the dendrogram of representative vectors of each module}
\item{sampleTree}{the dendrogram of samples}
\item{sampleModules}{A vector of module memberships, i.e. output[i] is the module to which sample i belongs.}
\item{intraModularStatistics}{A data frame of gene-wise connectivity statistics, as returned by intraModularStatistics()}
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Coexpression }


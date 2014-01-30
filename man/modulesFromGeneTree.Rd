\name{modulesFromGeneTree}
\alias{modulesFromGeneTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modules from Gene Tree
}
\description{
Computes discrete groups ("modules") of genes, based on a gene dendrogram.
}
\usage{
modulesFromGeneTree(geneTree, expressionData, dichotCor = NULL, tomDist = NULL, dynamicCutMethod = "tree", mergeModuleMethod = "conn")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{geneTree}{
The gene dendrogram, returned by clusterGenes()
}
\item{expressionData}{
The input data for the analysis, a normalized matrix of gene expression data, given as a data 
frame with rows as samples and columns as genes/probes.
}
\item{dichotCor}{
|cc|^beta, returned by clusterGenes(), only required if mergeModuleMethod=="conn"
}
\item{tomDist}{
Topological Overlap Matrix (TOM) expressed as a distance matrix, returned by clusterGenes(), 
only required if dynamicCutMethod=="hybrid"
}
\item{dynamicCutMethod}{
"tree" or "hybrid".  The former is the 'classic' approach used at Sage Bionetworks, while 
the latter is the default for the UCLA-WGCNA approach.  The details are here:
http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/BranchCutting/Supplement.pdf
}
 \item{mergeModuleMethod}{
"conn" or "eigen".  The former is the 'classic' approach used at Sage Bionetworks, in which a 
module's representative gene is the most highly CONNected one.  The latter is the default for the 
UCLA-WGCNA approach, in which the representative gene is the module's EIGENvector.
}
}
\details{
There are two steps to module finding:  (1) Cut the dendrogram branches to produce subtrees, 
the leaves of which become the modules, (2) Rejoining modules which are too close to be 
considered separate.  Each steps can be done by different algorithms.  For the "classic" version
use dynamicCutMethod="tree" and mergeModuleMethod="conn".  For the default method of UCLA-WGCNA
use dynamicCutMethod="hybrid" and mergeModuleMethod="eigen".  The references for the algorithms
are given below.
}
\value{ 
\item{geneModules}{
A vector of module memberships, i.e. geneModules[i] is the module to which gene i belongs. 
"grey" is a reserved name for genes which belong to no module.
}
\item{genePCTree}{the dendrogram of representative vectors of each module}
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
\keyword{Coexpression}% __ONLY ONE__ keyword per line

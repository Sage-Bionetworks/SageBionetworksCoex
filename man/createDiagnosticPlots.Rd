\name{createDiagnosticPlots}
\alias{createDiagnosticPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create diagnostic Plots
}
\description{
Given the output of Coexpression analysis, this creates a set of diagnostic plots (heatmaps and dendrograms), each in a file.
}
\usage{
createDiagnosticPlots(filetype = "jpg", outputDir = NULL, expressionData, dichotCor, tomDist, connectivity, beta, geneTree, geneModules, genePCTree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filetype}{
One of "jpg", "pdf", "png", "ps", not case sensitive.
}
  \item{outputDir}{
The directory/folder into which the plots are to be written.  If NULL, the session temporary directory is used.
}
  \item{expressionData}{
The input data for the analysis, a normalized matrix of gene expression data, given as a data frame with rows as samples and columns as genes/probes.
}
  \item{dichotCor}{
|cc|^beta, returned by clusterGenes()
}
  \item{tomDist}{
Topological Overlap Matrix (TOM) expressed as a distance matrix, returned by clusterGenes()
}
  \item{connectivity}{
The connectivity value (k) for each gene, as returned by clusterGenes()
}
  \item{beta}{
The scale-free exponent, returned by clusterGenes()
}
  \item{geneTree}{
The gene dendrogram, returned by clusterGenes()
}
  \item{geneModules}{
The module membership of each gene, returned by modulesFromGeneTree()
}
  \item{genePCTree}{
The dendrogram of module 'principal components', returned by modulesFromGeneTree()
}
}
\value{
The file handles for the created plots.
\item{ScaleFree}{File handle to "Scale free" plot (P(k) vs. k) for optimum value of beta}
\item{CorHeatmap}{File handle to heatmap of |cc|^beta, ordered according to dendrogram, with color-coded module membership along the margins}
\item{GeneDendrogram}{File handle to gene dendrogram with color-coded modules and color names}
\item{TOMHeatmap}{File handle to heatmap of TOM, rows and columns sorted by gene dendrogram}
\item{GenePCDendrogram}{File handle to dendrogram of clusters of module representative genes}
\item{ExpressionHeatmap}{File handle to heatmap of expression data, organized by gene and sample clusters}
}
\references{
\href{http://www.bepress.com/sagmb/vol4/iss1/art17/}{ Zhang, B. & Horvath, S. (2005) Statistical Applications in Genetics and Molecular Biology 4, Article 17.}
}
\note{
TOM and Correlation heat maps are not generated correctly in Windows, though they are in Unix.
}
\author{
Bruce Hoff
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Coexpression }

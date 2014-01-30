\name{createNetwork}
\alias{createNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Network
}
\description{
Creates a data frame for the derived network, given the input data and the data objects for the coexpression results.  As such 
this function is purely for data assembly, not for computation.
}
\usage{
createNetwork(expressionData,
				geneModules,
				intraModularStatistics)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expressionData}{
A normalized matrix of gene expression data, given as a data frame with rows as samples and columns as genes/probes.
  }
  \item{geneModules}{
The module membership of each gene, as computed by modulesFromGeneTree()
}
  \item{intraModularStatistics}{
A data frame of intramodule gene statistics, returned by intraModularStatistics()
  }
}

\value{
Returns a data frame with one row per gene and a column-concatenation of (a) the input data, (b) the
intramodular statistics (as produced by intraModularStatistics()), (c) the gene-wise standard deviation
of the input data, and (d) the gene module assignments.
}
\references{
\href{http://www.bepress.com/sagmb/vol4/iss1/art17/}{ Zhang, B. & Horvath, S. (2005) Statistical Applications in Genetics and Molecular Biology 4, Article 17.}

}
\author{
Bin Zhang
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Coexpression }

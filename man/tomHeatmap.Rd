\name{tomHeatmap}
\alias{tomHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tom Heatmap
}
\description{
Creates a heat map of the topological overlap matrix (TOM), ordered according to the given dendrogram, with color-coded module membership along the margins
}
\usage{
tomHeatmap(filePath, tomDist, geneDendrogram, geneModules)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filePath}{
The folder and file into which the graphic output file is placed.  The suffix indicates the type of file to create (jpg, pdf, png, or ps).
}
  \item{tomDist}{
Topological Overlap Matrix (TOM) expressed as a distance matrix, returned by clusterGenes()
}
  \item{geneDendrogram}{
The gene dendrogram, returned by clusterGenes()
}
  \item{geneModules}{
A vector of memberships, i.e. geneModules[i] is the module to which gene i belongs.
}
}
\details{
The graph is created and placed into the file indicated by 'filePath'.    This function works on Unix but not on Windows.
}
\value{
No value is returned.

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
\keyword{ heatmap }

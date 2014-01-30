\name{intraModularStatistics}
\alias{intraModularStatistics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
IntraModular Statistics
}
\description{
Computes gene-wise connectivity statistics, within and between modules.
}
\usage{
intraModularStatistics(dichotCor, tomDist, modules)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dichotCor}{
|cc|^beta, as computed by clusterGenes()
}
  \item{tomDist}{
Topological Overlap Matrix (TOM) expressed as a distance matrix, as computed by clusterGenes()
}
  \item{modules}{
The module membership of each gene, as computed by modulesFromGeneTree()
}
}

\value{
Returns a data frame with one row per gene and the following columns:
%% Bin Zhang needs to review the following definitions
k.in -- in-module number of connections of each gene
k.in.normed -- in-module number of connections of each gene, divided by module size
to.in -- sum of topological overlaps within modules
to.in.normed -- sum of topological overlaps within modules, divided by module size
k.all -- gene connectivity statistic, k.
k.all.normed -- gene connectivity, k, divided by number of genes.
to.all -- sum of topological overlaps 
to.all.normed  -- sum of topological overlaps, divided by number of genes.
k.out -- gene connectivity statistic, k, minus in-module number of connections of each gene
k.diff -- k.in minus k.out
to.out -- to.all minus to.in
to.dif -- to.in  minus to.out
cc.in --  within-module cluster coefficient for each gene
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

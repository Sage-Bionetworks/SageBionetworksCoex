\name{SageBionetworksCoex-package}
\alias{SageBionetworksCoex-package}
\alias{SageBionetworksCoex}
\docType{package}
\title{
Sage Bionetworks Coexpression
}
\description{
Compute gene coexpression based on gene expression data.  Output includes coexpression modules, 
diagnostic plots and supporting statistics.  The package leverages the WGCNA package published to C-RAN by the 
UCLA Array Data Analysis Group, but includes parameter settings and custom functions to replicate the 'flavor' of the algorithm
as executed at Sage Bionetworks.  R version 2.13 or higher is required.
}
\details{
\tabular{ll}{
Package: \tab SageBionetworksCoex\cr
Type: \tab Package\cr
Version: \tab 0.12\cr
Date: \tab 2011-08-27\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
The main entry point for this package is performCoexFromFiles().  Alternatively one may call performCoexpressionAnalysis(), 
and createDiagnosticPlots() directly, to perform the computation and generate the graphic output, respectively.
}
\note{This package requires the 'WGCNA' and (for graphics), 'lattice'.}

\author{
Bruce Hoff
<scicomp@sagebase.org>
}
\references{
\url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/}

\href{http://www.bepress.com/sagmb/vol4/iss1/art17/}{ Zhang, B. & Horvath, S. (2005) Statistical Applications in Genetics and Molecular Biology 4, Article 17.}

\href{http://bioinformatics.oxfordjournals.org/content/early/2007/11/16/bioinformatics.btm563.full.pdf+html}{Langfelder P, Zhang B, Horvath S (2007) Defining clusters from a hierarchical cluster tree: 
the Dynamic Tree Cut library for R. Bioinformatics 2008 24(5):719-720}

\url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA/}
}
\keyword{ package }


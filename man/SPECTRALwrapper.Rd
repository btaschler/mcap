% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPECTRALwrapper.R
\name{SPECTRALwrapper}
\alias{SPECTRALwrapper}
\title{Wrapper to do spectral clustering}
\usage{
SPECTRALwrapper(xx, k, true_labels = NULL, verbose = FALSE)
}
\arguments{
\item{xx}{The data matrix (n x p).}

\item{k}{The number of clusters.}

\item{true_labels}{Vector of true cluster assignments (when provided, it is 
used to compute the Rand index).}

\item{verbose}{Logical, when true: print progress information.}
}
\value{
\item{model_fit}{ Model fit (output of \code{\link{mixglasso}})}.
        \item{aRI}{ Adjusted Rand index (when \code{true_labels} is provided)}.
}
\description{
Spectral clustering using an RBF kernel. 
The wrapper performs multiple restarts in case true labels are provided 
(in order to optimise cluster assignments w.r.t. the Rand index).
}
\examples{
  SPECTRALwrapper(xx = matrix(rnorm(500),50), k=2)
}
\seealso{
\code{\link{MCAPfit}}, 
         \code{\link{GMMwrapper}}, 
         \code{\link[kernlab]{specc}}, 
         \code{\link[mclust]{adjustedRandIndex}}
}
\author{
Bernd Taschler \email{bernd.taschler@dzne.de}
}

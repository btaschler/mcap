% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CentrePerGroup.R
\name{CentrePerGroup}
\alias{CentrePerGroup}
\title{Centre a data matrix to mean zero (per group)}
\usage{
CentrePerGroup(xx, true_labels = numeric(nrow(xx)))
}
\arguments{
\item{xx}{The data matrix (n x p).}

\item{true_labels}{Vector of true cluster assignments. Default: 0 (all 
elements are treated as belonging to one group).}
}
\value{
\item{xx_centred}{ Mean centred data matrix.}
}
\description{
Mean centering of an input matrix. When true labels are provided, the 
groups/clusters are centred individually.
}
\examples{
  ## 10x5 matrix, 1 group:
  CentrePerGroup(xx=matrix(rnorm(50),10,5))
  
  ## 2 groups:
  CentrePerGroup(xx=matrix(rnorm(50),10,5), true_labels=round(runif(10)))
  
  ## 10 groups:
  CentrePerGroup(xx=matrix(rnorm(5000),100,50), true_labels=round(runif(100)*10))
}
\seealso{
\code{\link{colMeans}}, \code{\link{scale}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ComputeSEM.R
\name{ComputeSEM}
\alias{ComputeSEM}
\title{Compute standard error of the mean for an input vector}
\usage{
ComputeSEM(x)
}
\arguments{
\item{x}{Input vector (any NA will be ignored).}
}
\value{
Standard error of the mean (SEM) of \code{x}.
}
\description{
Compute standard error of the mean for an input vector
}
\examples{
  ## random vector:
  ComputeSEM(rnorm(100))
  
  ## with missing entries:
  x <- runif(1000)*5
  x[sample(1000, 50)] <- NA
  ComputeSEM(x)
}

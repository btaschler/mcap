% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rand1.R
\name{Rand1}
\alias{Rand1}
\title{Return a random sign vector}
\usage{
Rand1(n)
}
\arguments{
\item{n}{Number of elements}
}
\value{
Sign vector of length \code{n}.
}
\description{
Create a random sign vector of length \code{n} (iid ~ +1/-1).
}
\examples{
  Rand1(10)
}

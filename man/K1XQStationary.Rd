% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/K1XQStationary.R
\name{K1XQStationary}
\alias{K1XQStationary}
\title{Impose stationarity under the Q-measure}
\usage{
K1XQStationary(StationaryEigenvalues)
}
\arguments{
\item{StationaryEigenvalues}{Binary variable: set "1" if the user whises the largest eigenvalue
to be strictly smaller than 1. Set "0", otherwise}
}
\value{
list
}
\description{
Impose stationarity under the Q-measure
}
\examples{
stat <- 1 # Takes values 1 and 0
K1XQStationary(stat)


}

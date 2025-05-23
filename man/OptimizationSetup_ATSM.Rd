% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Optimization.R
\name{OptimizationSetup_ATSM}
\alias{OptimizationSetup_ATSM}
\title{Optimization routine for the entire selected ATSM}
\usage{
OptimizationSetup_ATSM(AuxVecSet, FFvec, EstType, tol = 1e-04)
}
\arguments{
\item{AuxVecSet}{List containing features for estimation of the risk-neutral parameters.}

\item{FFvec}{Log-likelihood function}

\item{EstType}{Estimation type}

\item{tol}{convergence tolerance (scalar). Default value is set as 1e-4.}
}
\description{
Optimization routine for the entire selected ATSM
}
\keyword{internal}

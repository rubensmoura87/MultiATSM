% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConfidenceBoundsBootstrap.R
\name{ComputeBounds_IRFandGIRF}
\alias{ComputeBounds_IRFandGIRF}
\title{Compute the confidence bounds from the model's numerical outputs}
\usage{
ComputeBounds_IRFandGIRF(
  ModelBootstrap,
  quants,
  FacDim,
  YieDim,
  ModelType,
  Economies,
  ndraws,
  Horiz
)
}
\arguments{
\item{ModelBootstrap}{numerical output set from the bootstrap analysis}

\item{quants}{quantile of the confidence bounds}

\item{FacDim}{dimension of the risk factor set}

\item{YieDim}{dimension of the bond yield set}

\item{ModelType}{Desired model type}

\item{Economies}{Economies that are part of the economic system}

\item{ndraws}{number of draws selected}

\item{Horiz}{horizon of numerical outputs}
}
\description{
Compute the confidence bounds from the model's numerical outputs
}
\keyword{internal}

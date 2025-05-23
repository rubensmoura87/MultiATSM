% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConfidenceBoundsBootstrap.R
\name{YieldQuantile_bs}
\alias{YieldQuantile_bs}
\title{Compute quantiles for model bond yield-related outputs}
\usage{
YieldQuantile_bs(
  DrawSet,
  LabIRF,
  ndraws,
  quants,
  Horiz,
  FacDim,
  YieDim,
  LabelsYies,
  ModelType,
  Ortho = FALSE
)
}
\arguments{
\item{DrawSet}{Draw-specific set}

\item{LabIRF}{vector containing the labels "IRF" and "GIRF"}

\item{ndraws}{number of draws}

\item{quants}{quantile of the confidence bounds}

\item{Horiz}{horizon of numerical outputs}

\item{FacDim}{dimension of the risk factor set}

\item{YieDim}{dimension of the bond yield set}

\item{LabelsYies}{labels of the factor set}

\item{ModelType}{desired model type}

\item{Ortho}{Orthogonolized version for the JLL models. Default is FALSE.}
}
\description{
Compute quantiles for model bond yield-related outputs
}
\keyword{internal}

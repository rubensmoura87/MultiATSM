% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BiasCorrection.R
\name{Get_V_tilde_BC}
\alias{Get_V_tilde_BC}
\title{Compute the variance-covariance matrix after the bias correction procedure}
\usage{
Get_V_tilde_BC(
  Phi_tilde,
  N,
  RiskFactors,
  GVARinputs,
  JLLinputs,
  FactorLabels,
  ModelType
)
}
\arguments{
\item{Phi_tilde}{Feedback matrix resulting from the bias-correction procedure}

\item{N}{number of country-specific spanned factors (scalar)}

\item{RiskFactors}{time series of the risk factors (T x F)}

\item{GVARinputs}{inputs used in the estimation of the GVAR-based models (see "GVAR" function). Default is set to NULL}

\item{JLLinputs}{inputs used in the estimation of the JLL-based models (see "JLL" function). Default is set to NULL}

\item{FactorLabels}{string-list based which contains the labels of all variables present in the model}

\item{ModelType}{string-vector containing the label of the model to be estimated}
}
\description{
Compute the variance-covariance matrix after the bias correction procedure
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bootstrap.R
\name{TimeVarWeights_GVAR}
\alias{TimeVarWeights_GVAR}
\title{Compute the star variables with time-varying weights}
\usage{
TimeVarWeights_GVAR(
  RiskFactors,
  Economies,
  RiskFactors_List,
  ListFactors,
  Wgvar,
  FactorLabels
)
}
\arguments{
\item{RiskFactors}{A matrix of the complete set of risk factors (F x T).}

\item{Economies}{A character vector containing the names of the economies included in the system.}

\item{RiskFactors_List}{List of domestic risk factors (both spanned and unspanned)}

\item{ListFactors}{List of risk factors}

\item{Wgvar}{List of transition matrices}

\item{FactorLabels}{A list of character vectors with labels for all variables in the model.}
}
\description{
Compute the star variables with time-varying weights
}
\keyword{internal}

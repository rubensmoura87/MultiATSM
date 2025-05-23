% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ForecastYields.R
\name{Gen_Forecast_Yields}
\alias{Gen_Forecast_Yields}
\title{Compute the bond yield forecast for any model type}
\usage{
Gen_Forecast_Yields(
  K0Z,
  K1Z,
  A,
  Bfull,
  ZZsubsample,
  C,
  J,
  YieldsLabels,
  ForLabels,
  ForHoriz,
  ModelType
)
}
\arguments{
\item{K0Z}{Intercept from the P-dynamics (F x 1)}

\item{K1Z}{Feedback matrix from the P-dynamics (F x F)}

\item{A}{Intercept of model-implied yields model (J x 1)}

\item{Bfull}{Slope of model-implied yields model (J x N or CJ x CN)}

\item{ZZsubsample}{Sub-sample of risk factors (F x t)}

\item{C}{Number of countries in the economic cohort (scalar)}

\item{J}{Number of country-specific bond yields}

\item{YieldsLabels}{Labels of bond yields}

\item{ForLabels}{Forecast labels (string-based vector)}

\item{ForHoriz}{Forecast horizon (scalar)}

\item{ModelType}{A string-vector containing the label of the model to be estimated}
}
\description{
Compute the bond yield forecast for any model type
}
\keyword{internal}

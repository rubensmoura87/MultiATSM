% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NumOutputs.R
\name{YieldsFitAll}
\alias{YieldsFitAll}
\title{Fit yields for all maturities of interest}
\usage{
YieldsFitAll(MatInt, ModelPara, FactorLabels, ModelType, Economies, YLab)
}
\arguments{
\item{MatInt}{numerical vector containing the fit maturities of interest}

\item{ModelPara}{List of model parameter estimates (See the "Optimization" function)}

\item{FactorLabels}{a string-list based which contains all the labels of all the variables present in the model}

\item{ModelType}{a string-vector containing the label of the model to be estimated}

\item{Economies}{a string-vector containing the names of the economies which are part of the economic system}

\item{YLab}{Label of yields ("Months" or "Yields")}
}
\description{
Fit yields for all maturities of interest
}
\keyword{internal}

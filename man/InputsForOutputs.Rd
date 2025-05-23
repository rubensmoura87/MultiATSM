% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InputsForOutputs.R
\name{InputsForOutputs}
\alias{InputsForOutputs}
\title{Collects the inputs that are used to construct the numerical and the graphical outputs}
\usage{
InputsForOutputs(
  ModelType,
  Horiz,
  ListOutputWished,
  OutputLabel,
  WishStationarityQ,
  DataFrequency,
  WishGraphYields = 0,
  WishGraphRiskFactors = 0,
  WishOrthoJLLgraphs = 0,
  WishForwardPremia = 0,
  LimFP = NULL,
  WishBootstrap = 0,
  ListBoot = NULL,
  WishForecast = 0,
  ListForecast = NULL,
  UnitYields = "Month"
)
}
\arguments{
\item{ModelType}{A character vector indicating the model type to be estimated.}

\item{Horiz}{A numeric scalar specifying the desired analysis horizon for the outputs.}

\item{ListOutputWished}{A list of desired graphical outputs. Available options are: "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremia".}

\item{OutputLabel}{A string for the name of the output label to be stored.}

\item{WishStationarityQ}{A binary variable (1 or 0) indicating whether to impose that the largest eigenvalue under Q is strictly
smaller than 1. Set to 1 to impose the restriction, or 0 otherwise.}

\item{DataFrequency}{A character vector specifying the data frequency. Available options: "Daily All Days", "Daily Business Days",
"Weekly", "Monthly", "Quarterly", "Annually".}

\item{WishGraphYields}{A binary variable (1 or 0) indicating whether the user wishes to generate graphs for yields. Default is 0.}

\item{WishGraphRiskFactors}{A binary variable (1 or 0) indicating whether the user wishes to generate graphs for risk factors. Default is 0.}

\item{WishOrthoJLLgraphs}{A binary variable (1 or 0) indicating whether the user wishes to generate orthogonalized JLL-based graphs.
Default is 0.}

\item{WishForwardPremia}{A binary variable (1 or 0) indicating whether the user wishes to generate forward premia graphs. Default is 0.}

\item{LimFP}{A numeric vector containing the maturities associated with the start and end dates of the loan.}

\item{WishBootstrap}{A binary variable (1 or 0) indicating whether the user wishes to perform bootstrap-based estimation. Default is 0.}

\item{ListBoot}{A List containing the following four elements:
\enumerate{
 \item \code{methodBS}: Desired bootstrap method: (a) 'bs' for standard residual bootstrap, (b) 'wild' for wild bootstrap,
                   or (c) 'block' for block bootstrap.
 \item \code{BlockLength}: If block bootstrap is chosen, specify the block length (numeric scalar).
 \item \code{ndraws}: Number of bootstrap draws.
 \item \code{pctg}: Confidence level expressed in basis points (numeric vector).
}}

\item{WishForecast}{A binary variable (1 or 0) indicating whether the user wishes to generate forecasts. Default is 0.}

\item{ListForecast}{A list containing the following three elements:
\enumerate{
 \item \code{ForHoriz}: forecast horizon;
 \item \code{t0Sample}: Index of the first variable in the information set.
 \item \code{t0Forecast}: Index of the first forecast cut-off date.
 \item \code{ForType}: A string specifying the desired forecast type. Available options are: "Rolling" or "Expanding".
}}

\item{UnitYields}{A character string indicating the maturity unit of yields. Options are: (i) "Month" for yields expressed in months, or (ii) "Year" for yields expressed in years. Default is "Month".}
}
\value{
List of necessary inputs to generate the graphs of the outputs of the desired model
}
\description{
Collects the inputs that are used to construct the numerical and the graphical outputs
}
\examples{

ModelType <- "JPS original"
Horiz <- 100
DesiredOutputGraphs <- c("Fit", "GIRF", "GFEVD")
OutputLabel <- "Test"
WishStationarityQ <- 1
WishGraphRiskFac <- 0
WishGraphYields <- 1

InputsList <- InputsForOutputs(ModelType, Horiz, DesiredOutputGraphs, OutputLabel,
                              WishStationarityQ, WishGraphYields, WishGraphRiskFac)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functionf_vectorized.R
\name{Update_ParaList}
\alias{Update_ParaList}
\title{converts the vectorized auxiliary parameter vector x to the parameters that go directly
into the likelihood function.}
\usage{
Update_ParaList(
  x,
  sizex,
  con,
  FactorLabels,
  Economies,
  JLLinputs = NULL,
  GVARinputs = NULL,
  ListInputSet
)
}
\arguments{
\item{x}{vector containing all the vectorized auxiliary parameters}

\item{sizex}{matrix (6x2) containing the size information of all parameters}

\item{con}{if con = 'concentration', then set the value of the parameter whose name
contains @ to empty}

\item{FactorLabels}{string-list based which contains the labels of all the variables present in the model}

\item{Economies}{string-vector containing the names of the economies which are part of the economic system}

\item{JLLinputs}{Set of necessary inputs used in the estimation of the JLL-based models}

\item{GVARinputs}{Set of necessary inputs used in the estimation of the GVAR-based models}

\item{ListInputSet}{variable inputs used in the optimization (see "Optimization" function)}
}
\value{
same form as the list of parameters, except now the parameters are updated with the values provided by the auxiliary x.
Importantly, by construction, all the constraints on the underlying parameters are satisfied.
}
\description{
converts the vectorized auxiliary parameter vector x to the parameters that go directly
into the likelihood function.
}
\references{
This function is a modified version of the "update_para" function by Le and Singleton (2018). \cr
"A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."\cr
(Euro Area Business Cycle Network Training School - Term Structure Modelling).
Available at: https://cepr.org/40029
}
\keyword{internal}

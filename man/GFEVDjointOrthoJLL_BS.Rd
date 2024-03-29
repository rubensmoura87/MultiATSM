% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NumOutputsBootstrap.R
\name{GFEVDjointOrthoJLL_BS}
\alias{GFEVDjointOrthoJLL_BS}
\title{GFEVDs after bootstrap for JLL-based models}
\usage{
GFEVDjointOrthoJLL_BS(
  ModelType,
  ModelParaBoot,
  GFEVDhoriz,
  FactorLabels,
  Economies
)
}
\arguments{
\item{ModelType}{string-vector containing the label of the model to be estimated}

\item{ModelParaBoot}{list of model parameter estimates (see the "Optimization" function) after a bootstrap draw}

\item{GFEVDhoriz}{single numerical vector conataining the desired horizon of analysis for the GFEVDs}

\item{FactorLabels}{string-list based which contains all the labels of all the variables present in the model}

\item{Economies}{string-vector containing the names of the economies which are part of the economic system}
}
\description{
GFEVDs after bootstrap for JLL-based models
}
\references{
\itemize{
\item This function is a modified and extended version of the "fevd" function from
Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.

\item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
}
}

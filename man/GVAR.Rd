% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GVAR.R
\name{GVAR}
\alias{GVAR}
\title{Estimate a GVAR(1) and a VARX(1,1,1)}
\usage{
GVAR(GVARinputs, N)
}
\arguments{
\item{GVARinputs}{List containing the following necessary inputs for the estimation of the GVAR:
\enumerate{
       \item  Economies:  string-vector containing the names of the economies which are part of the economic system
       \item  'GVARFactors': list of all variables that are used in the estimation of the VARX \cr
               (see e.g. 'CM_Factors_GVAR' file);
       \item 'VARXtype': character-vector containing three possibilities:
 \itemize{
       \item 'unconstrained': model is estimated without any constrained (each equation is estimated individually by OLS);
       \item 'constrained: Spanned Factors': model is estimated taking into account the fact that foreign-pricing-factors
                             do NOT impirge on
                              (i) domestic economic variables and (ii) domestic pricing factors.
                              (equations are estimated by restricted least squares)
       \item 'constrained: ' extended by the name of the risk factor: model is estimated taking into account the fact that
                             the restricted factor is only affected by its own lagged values and the lagged values of its own star variables.
       (equations are estimated by restricted least squares)
         }
         \item 'Wgvar':  GVAR transition matrix (C x C) - see the output from 'Transition_Matrix' function
}}

\item{N}{number of country-specific spanned factors (scalar)}
}
\value{
A list containing
\enumerate{
\item parameters of the country-specific VARX(1,1,1)
\itemize{
\item intercept (M+Nx1);
\item phi_1   (M+N x M+N);
\item phi_1^star (M+N x M+N);
\item phi_g (M+N x M+N);
\item Sigma (M+N x G)
}
\item parameters of the GVAR.
\itemize{
\item F0 (F X 1);
\item F1 (F x F);
\item Sigma_y (F x F)
}
}
}
\description{
Estimate a GVAR(1) and a VARX(1,1,1)
}
\examples{
data(CM_Factors_GVAR)

N <- 3

GVARinputs <- list()
GVARinputs$Economies <- c("China", "Brazil", "Mexico", "Uruguay")
GVARinputs$GVARFactors <- FactorsGVAR
GVARinputs$VARXtype <- "unconstrained"
GVARinputs$Wgvar <- matrix( c(0, 0.83, 0.86, 0.38,
                              0.65, 0, 0.13, 0.55,
                              0.32, 0.12, 0, 0.07,
                              0.03, 0.05, 0.01, 0), nrow = 4, ncol = 4)

GVAR(GVARinputs, N)
}
\references{
Chudik and Pesaran, (2016). "Theory and Practice of GVAR modelling" (Journal of Economic Surveys)
}

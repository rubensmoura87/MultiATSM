% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logdet.R
\name{logdet}
\alias{logdet}
\title{computes the logarithm of determinant of a matrix A}
\usage{
logdet(A, opt)
}
\arguments{
\item{A}{squared matrix}

\item{opt}{"chol" or "NULL" (text)}
}
\description{
computes the logarithm of determinant of a matrix A
}
\details{
Theoretically, this function should be functionally equivalent to log(det(A)). However, it avoids the
overflow/underflow problems that are likely to happen when applying det to large matrices.
The key idea is based on the mathematical fact that the determinant of a triangular matrix equals the
product of its diagonal elements. Hence, the matrix's log-determinant is equal to the sum of their logarithm
values. By keeping all computations in log-scale, the problem of underflow/overflow caused by product of
many numbers can be effectively circumvented. The implementation is based on LU factorization.
}
\references{
This function is based on the "logdet" function by Le and Singleton (2018).\cr
 "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
 (Euro Area Business Cycle Network Training School - Term Structure Modelling).
 Available at: https://cepr.org/40029
}
\keyword{internal}

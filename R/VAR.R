#' Estimates a standard VAR(1)
#'
#' @param RiskFactors A numeric matrix (FTx T) representing the time series of risk factors.
#' @param VARtype String vector with two possible values: 'unconstrained' or 'constrained'.
#' @param Bcon Constraints matrix (F+1 x N), which includes an intercept. If Bcon(i,j) = NA, then B(i,j) is treated as a free parameter. \cr
#'             Default is set to NULL.
#'
#' @return intercept, feedback matrix and the variance-covariance matrix of a VAR(1)
#' @examples
#' data("CM_Factors")
#' # Example 1: unconstrained case
#' VAR(RiskFactors, VARtype= 'unconstrained')
#'
#' # Example 2: constrained case
#' K <- nrow(RiskFactors)
#' Bcon <- matrix(0, nrow = K, ncol = K+1)
#' Bcon[ , 1:3] <- NaN
#' VAR(RiskFactors, VARtype= 'constrained', Bcon)
#'
#' @export

VAR <- function(RiskFactors, VARtype, Bcon = NULL) {
  K <- nrow(RiskFactors)
  T <- ncol(RiskFactors)
  LHS <- RiskFactors[, 2:T]
  RHS <- RiskFactors[, 1:(T-1)]

  if (VARtype == 'unconstrained') {
    RegVAR <- stats::lm(t(LHS) ~ t(RHS)) # VAR(1) under the P.
    K0Z <- t(t(RegVAR$coefficients[1, ]))
    K1Z <- t(RegVAR$coefficients[2:(K+1), ])
    eZ <- RegVAR$residuals
    SSZ <- crossprod(eZ) / (T - 1)
  } else { # i.e. if VARtype == 'constrained'
    intercept <- rep(1, times = T-1)
    RHS <- rbind(intercept, RHS)
    Coeff <- Reg__OLSconstrained(LHS, RHS, Bcon)
    colnames(Coeff) <- rownames(RHS)
    rownames(Coeff) <- rownames(RHS)[-1]

    K0Z <- as.matrix(Coeff[, 1])
    K1Z <- Coeff[, 2:(K+1)]
    eZ <- LHS - Coeff %*% RHS
    SSZ <- crossprod(t(eZ)) / (T - 1)
  }

  return(list(K0Z = K0Z, K1Z = K1Z, SSZ = SSZ))
}

#########################################################################################################
#' Restricted OLS regression
#'
#' @param Y left hand side variables (M x T)
#' @param X regressors (i.e. N-1 variables + the intercept) (N x T)
#' @param Bcon constraints matrix (M x N). If Bcon(i,j) = nan --> B(i,j) is a free parameter
#' @param G weighting matrix (psd) - (M x M). Default is set to be identity
#'
#' @keywords internal
#'
#' @return matrix of coefficient (M x N)
#' @details
#' Estimate of B is obtained by minimizing the objective:
#'   sum_t (Y_t-B X_t)' G^(-1) (Y_t-B*X_t)
#' subject to the constraint that B = Bcon for all non-nan entries of Bcon
#'
#' @references
#' This function is based on the "Reg__OLSconstrained" function by Le and Singleton (2018).
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029
#'

Reg__OLSconstrained <- function(Y, X, Bcon, G = NULL) {
  M <- nrow(Y)
  T <- ncol(Y)
  N <- nrow(X)

  if (is.null(G)) {
    G <- diag(M)
  }

  id <- is.nan(Bcon) # id= TRUE -> free parameter; id= FALSE -> constrained

  B <- matrix(0, nrow = M, ncol = N)
  B[!id] <- Bcon[!id]

  if (any(B != 0)) {
    Y <- Y - B %*% X
  }

  bigX <- kronecker((X%*%t(X))/T, solve(G))
  bigY <- solve(G,Y%*%t(X)) /T


  B[id] <- solve(bigX[id, id], bigY[id])

  return(B)
}

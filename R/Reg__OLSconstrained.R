#' Restricted OLS regression
#"
#'@param Y    left hand side variables (M x T)
#'@param X    regressors (i.e. N-1 variables + the intercept) (N x T)
#'@param Bcon   constraints matrix (M x N). If Bcon(i,j) = nan --> B(i,j) is a free parameter
#'@param G      weighting matrix (psd) - (M x M). Default is set to be identity
#
#
#'@return  matrix of coefficient (M x N)
#'@details
#'# Estimate of B is obtained by minimizing the objective:\cr
#'   sum_t (Y_t-B X_t)' G^{-1} (Y_t-B*X_t) \cr
#' subject to the constraint that B = Bcon for all non-nan entries of Bcon
#'
#'@references
#' This function is based on the "Reg__OLSconstrained" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029



Reg__OLSconstrained <-function(Y, X, Bcon, G){


  M <- dim(Y)[1]
  T <- dim(Y)[2]

  if (!exists('G')||is.null(G)){
    G <- diag(M)
  }

  N <- dim(X)[1]

  id <- is.nan(Bcon) # id= TRUE -> free parameter; id= FALSE -> constrained

  B <- matrix(0, nrow=M, ncol=N)
  B[!id] <- Bcon[!id]

  if (any(B!=0)){
    Y <- Y - B%*%X
  }

  bigX <- kronecker((X%*%t(X))/T, solve(G))
  bigY <- solve(G,Y%*%t(X)) /T


  B[id] <- solve(bigX[id, id], bigY[id])

  return(B)
}

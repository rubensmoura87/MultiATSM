#' Estimates a VAR(1)

#' @param RiskFactors matrix containing all the risk factors (K x T)
#' @param VARtype string-vector which accommodates two possibilities: 'unconstrained' or 'constrained'
#' @param Bcon constraints matrix (K+1 x N) - should contain an intercept. IfBcon(i,j) = nan --> B(i,j) is a free parameter. \cr
#'              Default is set to NULL.

#'@return  intercept, feedback matrix and the variance-covariance matrix of a VAR(1)
#'@examples
#' data("CM_Factors")
#' #Example 1
#' VAR(RiskFactors, VARtype= 'unconstrained')

#' #Example 2
#' K <- nrow(RiskFactors)
#' Bcon <-matrix(0, nrow = K, ncol = K+1)
#' Bcon[,1:3] <- NaN
#' VAR(RiskFactors, VARtype= 'constrained', Bcon)
#'
#'
#'@export

VAR <- function(RiskFactors, VARtype, Bcon = NULL){

  K <- nrow(RiskFactors)
  T <- ncol(RiskFactors)
  LHS <- RiskFactors[,2:T]
  RHS <- RiskFactors[,1:(T-1)]

  # Unconstrained
  if (VARtype == 'unconstrained'){
    RegVAR <-lm( t(LHS)~ t(RHS)) # VAR(1) under the P.
    K0Z <- t(t(RegVAR$coefficients[1,]))
    K1Z <- t(RegVAR$coefficients[2:(K+1),])
    eZ <- RegVAR$residuals
    eZ <- t(eZ)
    SSZ <- eZ%*%t(eZ)/dim(eZ)[2]
  }

  # Constrained
  if (VARtype == 'constrained'){

    intercept <- rep(1, times= T-1)
    RHS <- rbind(intercept, RHS)

    Coeff <- Reg__OLSconstrained(LHS, RHS, Bcon, G=NULL)
    colnames(Coeff) <- rownames(RHS)
    rownames(Coeff) <- rownames(RHS)[-1]

    K0Z <- as.matrix(Coeff[,1])
    K1Z <- Coeff[,2:(K+1)]
    eZ <- LHS - Coeff%*%RHS
    SSZ <- eZ%*%t(eZ)/dim(eZ)[2]
  }

  out <- list(K0Z, K1Z, SSZ)
  names(out) <- c("K0Z", "K1Z", "SSZ")

  return(out)
}

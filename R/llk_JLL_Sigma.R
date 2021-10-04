#' Build the log-likelihood function of the P-dynamics from the JLL-based models
#

#'@param VecPara   vector that contains all the non-zero entries of the lower-triangualar part of the Cholesky factorization
#'@param res       residuals from the VAR of the JLL model (K x T)
#'@param IdxNONzero    vector that contains indexes of the matrix of the non-zero entries of the Cholesky factorization
#'@param K              dimensions of the variance-covariance matrix (scalar)
#'

#
#'@keywords internal
#'@return value of the log-likelihood function (scalar)
#



llk_JLL_Sigma <- function(VecPara, res, IdxNONzero, K){


  Se <- matrix(0, K,K)
  Se[IdxNONzero] <- VecPara # restricted Se matrix
  Sigma_Res <- Se%*%t(Se)
  y <- GaussianDensity(res,Sigma_Res)

  llk <- -mean(y)

  return(llk)
}

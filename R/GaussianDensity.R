#' computes the density function of a gaussian process

#'@param res    matrix of residuals (N x T)
#'@param SS     covariance matrice or array of covariance matrices \cr
#'               (If dim(SS) > 3, then the model has a stochastic volatility) (N x N) or (N x N x T)
#'
#'@param invSS  Inverse of SS (N x N) or (N x N x T) - optional input
#'@param logabsdetSS   log(abs(|SS|)) (1 x T) - optional input
#'
#'@return y   - vector of density (1 x T)
#'
#'@references
#' This function is based on the "Gaussian" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'



GaussianDensity <- function(res,SS, invSS, logabsdetSS){


  nargin <- nargs()


  N <- dim(res)[1]
  T <- dim(res)[2]


  if (is.na(dim(SS)[3])){ # If SS is a matrix and NOT a 3-dimension array.
    if (!missing("invSS")){
      SSres <-  invSS%*%res
    }else{
      SSres <- solve(SS,res, tol = 1e-50)} # SSres <- (SS)^(-1)*res


    if (missing("logabsdetSS")){
      logabsdetSS <- 0.5*logdet(SS%*%t(SS),opt="chol")} # why do we multiply by 0.5? Because in "logdet", we do
                                                        # 2 *sum(log(abs)).

  }else{
    if(nargin<3){
      invSS <- mult__inv(SS, whichoutput= NULL , 2)$inva
      logabsdetSS <- mult__inv(SS, whichoutput= NULL , 2)$logabsdet
    }else if (nargin==3){
      logabsdetSS <- mult__inv(SS, whichoutput= NULL,'logabsdet')$logabsdet
    }

    a <- invSS[[1]]
    b <- array(res, c(N, 1, T))
    SSres <- array(mult__prod(a, b, c=NULL), c(N, T))
    logabsdetSS <- logabsdetSS[[1]]
  }
  y <- -0.5*N*log(2*pi)-0.5*logabsdetSS-0.5*abs(colSums(res*SSres)) # Last term is equal to "t(res)*(SS)^(-1)*res"

  return(y)

}

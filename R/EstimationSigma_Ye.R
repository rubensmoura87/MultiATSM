#' Estimate numerically the Cholesky-factorization from the JLL-based models
#'

#'@param SigmaUnres   unrestricted variance-covariance matrix (K x K)
#'@param res          residuals from the VAR of the JLL model (K x T)
#'@param M          number of domestic unspanned factors per country (scalar)
#'@param G          number of global unspanned factors (scalar)
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'@param DomUnit Name of the economy which is assigned as the dominant unit. \cr
#'               If no dominant unit is assigned, then this variable is defined as "none"
#'
#'@importFrom functional Curry
#'@importFrom neldermead fminsearch optimset
#
#'@keywords internal
#'@return    Cholesky-factorization after the maximization (K x K)
#'


EstimationSigma_Ye <- function(SigmaUnres, res, M, G, Economies, DomUnit){

  # SIGMA_Ye
  K <- nrow(SigmaUnres)
  C <- length(Economies)
  N <- (K - G - M*C)/C

  # Set the constraints in the Sigma matrix
  Se <- CholRestrictionsJLL(SigmaUnres, M, G, N, Economies, DomUnit)
  IdxNONzeroSigmaJLL <- which(Se!=0)
  x <- Se[IdxNONzeroSigmaJLL] # vector containing the initial guesses

  MLfunction <- Curry(llk_JLL_Sigma, res=res, IdxNONzero= IdxNONzeroSigmaJLL, K=K)

  iter <- 'off' # hides the outputs of each iteration. If one wants to display these features then set 'iter'
  options200 <- optimset(MaxFunEvals = 200000*numel(x), Display = iter,
                         MaxIter = 200000, GradObj='off', TolFun= 10^-2, TolX= 10^-2)


  Xmax <- fminsearch(MLfunction, x, options200)$optbase$xopt
  SIGMA_Ye <- matrix(0, K,K)
  SIGMA_Ye[IdxNONzeroSigmaJLL]<- Xmax # Cholesky term (orthogonalized factors)

  #Labels
  rownames(SIGMA_Ye) <- rownames(SigmaUnres)
  colnames(SIGMA_Ye) <- rownames(SigmaUnres)

  return(SIGMA_Ye)
}

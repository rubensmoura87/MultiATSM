#' Compute the cross-section loadings of yields of a canonical A0_N model ("sep Q" models)
#'
#'@param mat       vector of maturities (J x 1). Maturities are in multiples of the discrete interval
#'                 used in the model
#'@param K1XQ      risk neutral feedback matrix (N x N)
#'@param dX        state loadings for the one-period rate (1xN). Default is a vector of ones
#'@param r0        the long run risk neutral mean of the short rate (scalar)
#'@param SSX       the covariance matrix of the errors (N x N)
#'
#'
#'@keywords internal
#'
#'@return
#' List containing:
#' \itemize{
#' \item Intercept (Jx1)
#' \item slope (JxN)
#' \item the betan (JX1, part of the intercepts unrelated to the long run risk neutral mean r0)
#' coefficients of a canonical A_0(N).
#' }
#'@references
#'  \itemize{
#'  \item This function is based on the "A0N__computeBnAn" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'  \item Dai and Singleton (2000). "Specification Analysis of Affine Term Structure Models" (The Journal of Finance)
#'  }
#'




A0N__computeBnAn_sepQ <- function(mat, K1XQ, dX, r0, SSX){


  mat <- round(mat[])
  K <- max(mat) # longest maturity of interest
  N <- dim(K1XQ)[1] # Number of latent factors


  if (!exists("dX", inherits = FALSE) || is.null(dX) ){
    dX <- rep(1, times=N)
  }
  dX <- c(dX)

  BnX <- matrix(0, nrow = K+1, ncol= N)
  AnX <- matrix(0, nrow = K+1, ncol=1)

  # Generate the loadings from a canonical A_0(N) model:
  for ( k in 1:K){
    BnX[k+1,] <- dX + BnX[k,]%*%K1XQ # Note: BnX depends only on K1XQ and dx NOT r0 or SSX.
    if (!exists("r0", inherits = FALSE) || is.null(r0) ){
      AnX <- as.numeric(vector())
    } else {
      AnX[k+1] = r0 + AnX[k] - 0.5*BnX[k,]%*%SSX%*%t(t(BnX[k,]))
    }
  }

  BnX <- BnX[mat+1,]/mat # Adjust the loading for the maturity


  if (!exists("r0", inherits = FALSE) || is.null(r0) ){
    betan <- c()
  } else {
    AnX <- AnX[mat+1]/mat # Adjust the loading for the maturity
    betan <- AnX - r0
  }


  return(list(BnX, AnX, betan))

}


##########################################################################################################
##########################################################################################################

#' Compute the cross-section loadings of yields of a canonical A0_N model ("joint Q" models)
#'
#'@param mat       vector of maturities (J x 1). Maturities are in multiples of the discrete interval
#'                 used in the model
#'@param K1XQ      risk neutral feedback matrix (N x N)
#'@param dX        state loadings for the one-period rate (1xN). Default is a vector of ones
#'@param r0        the long run risk neutral mean of the short rate (scalar)
#'@param SSX       the covariance matrix of the errors (N x N)
#'@param Economies Set of economies that are part of the economic system (vector of text)
#'
#'@keywords internal
#'
#'@return
#' List containing:
#' \itemize{
#' \item Intercept (Jx1)
#' \item slope (JxN)
#' \item the betan (JX1, part of the intercepts unrelated to the long run risk neutral mean r0)
#' coefficients of a canonical A_0(N).
#' }
#'@references
#'  This function is an extended version of the "A0N__computeBnAn" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."\cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029




A0N__computeBnAn_jointQ <- function(mat, K1XQ, dX, r0, SSX, Economies){

  C <- length(Economies)
  mat <- round(mat[])
  K <- max(mat) # longest maturity of interest
  NC <- dim(K1XQ)[1] # Total number of latent factors of the system
  N <- NC/C # number of spanned factors per country
  J <- length(mat)

  if (!exists("dX", inherits = FALSE) || is.null(dX) ){
    dxCS <- rep(1,N)
    dX <- matrix(0, nrow=C, ncol=NC)
    idx0 <-0

    for (j in 1:C){
      idx1<- idx0+N
      dX[j,(idx0+1):idx1] <- dxCS
      idx0<- idx1
    }
  }


  # Generate the loadings from a canonical A_0(N) model for all the countries in the system:
  idx0 <- 0
  for (i in 1:C){
    AnXCS <- matrix(0, nrow = K+1, ncol=1) # Country-specific A
    BnXCS <- matrix(0, nrow = K+1, ncol= NC) # Country-specific B
    for ( k in 1:K){
      idx1 <- idx0+N
      BnXCS[k+1,] <- dX[i,] + BnXCS[k,]%*%K1XQ # Note: BnX depends only on K1XQ and dx NOT r0 or SSX.
      if (!exists("r0", inherits = FALSE) || is.null(r0) ){
        AnXCS <- as.numeric(vector())
      } else {
        AnXCS[k+1] = r0[i] + AnXCS[k] - 0.5*BnXCS[k,]%*%SSX%*%t(t(BnXCS[k,]))
      }
    }
    if (i==1){
      BnX <- BnXCS[mat+1,(idx0+1):idx1]/mat # Adjust the loading for the maturity

      AnXCS <- AnXCS[mat+1]/mat
      AnX <- t(t(AnXCS))

      idx0 <- idx1
    }else{
      BnXCS <- BnXCS[mat+1, (idx0+1):idx1]/mat
      BnX <- magic::adiag(BnX,BnXCS)

      AnXCS <- t(t(AnXCS[mat+1]/mat))
      AnX <- rbind(AnX,AnXCS)

      idx0 <- idx1
    }
  }


  idx0 <- 0
  if (!exists("r0", inherits = FALSE) || is.null(r0) ){
    betan <- as.numeric(vector())
  } else {
    betan <- c()
    for (h in 1:C){
      idx1 <- idx0 +J
      betan[(idx0+1):idx1] <- AnX[(idx0+1):idx1] - t(t(rep(r0[h], times = J)))
      idx0 <- idx1
    }
  }


  return(list(BnX, AnX, betan))

}

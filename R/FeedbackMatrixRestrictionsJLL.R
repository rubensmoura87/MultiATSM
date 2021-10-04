#' Set the zero-restrictions on the feedback matrix of JLL's P-dynamics
#'

#'@param DomUnit Name of the economy which is assigned as the dominant unit. \cr
#'               If no dominant unit is assigned, then this variable is defined as "none"
#'@param K       Total number of risk factors of the economic system (scalar)
#'@param G       Number of global unspanned factors (scalar)
#'@param M       Number of country-specific unspanned factors (scalar)
#'@param N       Number of country-specific spanned factors (scalar)
#
#'@keywords internal
#'@return       matrix containing the restrictions of the feedback matrix (K x K)
#



FeedbackMatrixRestrictionsJLL<- function(DomUnit,K,G,M,N){

  C <- (K-G)/(M+N) # number of countries of the system


  Bcon <- matrix(0, nrow= K, ncol=K+1) # Includes all variables and the intercept
  Bcon[, 1] <- NaN # Intercept
  Bcon[ ,2:(G+1)] <- NaN # Global

  if (DomUnit== "None"){
    IDXrow000 <- G
    IDXcol000 <- G+1

    for (i in 1:C){
      IDXrow111 <- IDXrow000 + M +N
      IDXcol111 <- IDXcol000 + M +N

      Bcon[(IDXrow000+1):IDXrow111,(IDXcol000+1):IDXcol111] <- NaN

      IDXrow000 <- IDXrow111
      IDXcol000 <- IDXcol111
    }
  } else{ # With a dominant unit

    Bcon [, (G+1+1):(G+M+N+1)] <- NaN # Dominant unit

    IDXrow000 <- G + M + N
    IDXcol000 <- G+M+N+1

    for (i in 1:(C-1)){
      IDXrow111 <- IDXrow000 + M +N
      IDXcol111 <- IDXcol000 + M +N

      Bcon[(IDXrow000+1):IDXrow111,(IDXcol000+1):IDXcol111] <- NaN

      IDXrow000 <- IDXrow111
      IDXcol000 <- IDXcol111
    }
  }

  return(Bcon)
}

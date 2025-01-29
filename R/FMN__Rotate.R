#' Performs state rotations
#'
#' @param y0    list  of model parameters as described below
#' @param U1    matrix (N x N)
#' @param U0    vector (N x 1). Optional. Default: vector of zeros.
#'
#'@keywords internal
#'
#'@return      y1 - list of outputs after the transformation, the structure parallels that of y0
#'
#'
#'@details
#' This function performs a rotation from a model with Z as states to one with S = U0 + U1*Z as states. \cr
#' Specifically, each model is characterized by the following inputs organized in a list of variables:\cr
#' (i) K0: intercepts (N x 1);\cr
#' (ii) K1: feedback matrix (N x N*p); \cr
#' (iii) SS: volatility matrices (N x N*(M+1))\cr

#' More specifically, the state Z follows the dynamics: \cr
#' Z_t = N(K0 + K1 [Z_\{t-1\}; Z_\{t-2\}; ...],  SSi[ , , 1] + sum_\{i=1\}^M SSi[ , ,i+1]%*%V_\{i,t\}) \cr
#' where SSi <- array(SS, c(N, N, M+1))

#'@references
#' #' This function is modified version of the "FMN__Rotate" function by Le and Singleton (2018). \cr
#' "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."\cr
#' (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029
#'


FMN__Rotate <- function(y0, U1, U0){


  N <- dim(U1)[1]
  if (missing('U0') || is.null(U0)){ U0 <- matrix(0, nrow = N, ncol = 1) }
  y1 <- list()


  if ("Q" %in% names(y0)){
    y1$Q$K1 <- mrdivide(U1%*%y0$Q$K1,U1) # (U1* K1XQ) K1Z = U1^(-1)
    y1$Q$K0 <- U0+U1%*%y0$Q$K0-y1$Q$K1%*%U0 # K0Z = U0+ U1*K0X - K1Z*U0

    M <- dim(y0$Q$SS)[[2]]/N -1
    SSi <- array(y0$Q$SS, c(N, N, M+1))
    SSi <- mult__prod(mult__prod(U1, drop(SSi), c= NULL), t(U1), c=NULL) # Sigma(Z) = U1*Sigma(X)*(U1)'
    y1$Q$SS <- array(SSi, c(N, N*(M+1)))
  }

  # pricing:  Mapping of parameters X <-> Z.
  if ("B" %in% names(y0)){
    y1$B <- if (N == 1) mrdivide(matrix(y0$B), U1) else mrdivide(y0$B, U1) # BX* BZ = U1 --> BZ = BX*U1^(-1)
    y1$A <- y0$A - y1$B%*%U0 # AX = AZ + BZ*U0 --> AZ = AX - BZ*U0
  }


  # P dynamics:
  if ("P" %in% names(y0)){
    NP <- numel(y0$P$K0) #  NP may be greater than N if there is unspanned variables
    U0 <- rbind(U0, matrix(0, nrow= NP-N, ncol=1))
    U1 <- magic::adiag(U1, diag(NP-N))

    p <- dim(y0$P$K1)[2]/dim(y0$P$K1)[1]

    y1$P$K1 <- array(data=NA, c(size(y0$P$K1)))
    temp <- 0
    for (i in 1:p){
      a <- as.matrix(U1%*%y0$P$K1[, ((i-1)*NP+1):(i*NP)])
      b <-  as.matrix(U1)
      y1$P$K1[, ((i-1)*NP+1):(i*NP)] = mrdivide(a,b)
      temp <- temp + y1$P$K1[ , ((i-1)*NP+1):(i*NP)]
    }
    y1.P.K0 <- U0+U1%*%y0$P$K0-temp%*%U0
    N <- dim(y0$P$SS)[1]
    M <- dim(y0$P$SS)[2]/N-1
    SSi <- array(y0$P$SS, c(N, N, M+1))
    SSi <- mult__prod(mult__prod(as.matrix(U1), drop(SSi), c=NULL), t(U1), c=NULL)
    y1.P.SS <- array(SSi, c(N,N*(M+1)) )
  }


  if ( !("P" %in% names(y0)) && !("Q" %in% names(y0)) && !("B" %in% names(y0)) && !("A" %in% names(y0)) ) {

    NP <- numel(y0$K0)
    if (NP-N !=0){ U0 <- rbind(U0, matrix(0, nrow= NP-N,ncol =1))}
    U1 <- magic::adiag(U1, diag(NP-N))

    p <- if (is.matrix(y0$K1)) dim(y0$K1)[2]/dim(y0$K1)[1] else 1

    if (is.null(dim(y0$K1))) {y1$K1 <- NA } else {y1$K1 <- array(data = NA, dim = dim(y0$K1)) }
    temp <- 0
    for (i in 1:p){
      if (is.matrix(y0$K1)) {
      a <- U1 %*% y0$K1[, ((i-1) * NP + 1):(i * NP)]
      b <- as.matrix(U1)
      y1$K1[, ((i-1)*NP+1):(i*NP)] <- mrdivide(a,b)
      temp <- temp + y1$K1[, ((i-1)*NP+1):(i*NP)]
      } else{
      a <- U1 %*% y0$K1
      b <- as.matrix(U1)
      y1$K1 <- mrdivide(a,b)
      temp <- temp + y1$K1
      }
    }

    y1$K0 <- U0+U1%*%y0$K0-temp%*%U0
    N <-  if (is.null(dim(y0$SS))) 1 else dim(y0$SS)[1]
    M <-  if (is.null(dim(y0$SS))) 0 else dim(y0$SS)[2]/N-1
    SSi <- array(y0$SS, c(N, N, M+1))
    SSi <- mult__prod(mult__prod(U1, drop(SSi), c=NULL), t(U1), c= NULL)
    y1$SS = array(SSi, c(N, N*(M+1)))

  }


  return(y1)
}

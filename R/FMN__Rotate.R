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

FMN__Rotate <- function(y0, U1, U0) {
  N <- dim(U1)[1]
  if ( is.null(U0)) {
    U0 <- matrix(0, nrow = N, ncol = 1)
  }

  y1 <- list()
  U1_inv <- solve(U1)

  # 1) Q transformation
  if ("Q" %in% names(y0)) {
    y1$Q$K1 <- U1 %*% y0$Q$K1 %*% U1_inv # (U1* K1XQ) K1Z = U1^(-1)
    y1$Q$K0 <- U0 + U1 %*% y0$Q$K0 - y1$Q$K1 %*% U0 # K0Z = U0+ U1*K0X - K1Z*U0

    M <- dim(y0$Q$SS)[[2]] / N - 1
    SSi <- array(y0$Q$SS, c(N, N, M + 1))
    SSi <- mult__prod(mult__prod(U1, drop(SSi)), t(U1)) # Sigma(Z) = U1*Sigma(X)*(U1)'
    y1$Q$SS <- array(SSi, c(N, N * (M + 1)))
  }

  # 2) Pricing transformation
  # Mapping of parameters X <-> Z.
  if ("B" %in% names(y0)) {
    y1$B <- if (N == 1) matrix(y0$B) %*% solve(U1) else y0$B %*% U1_inv # BX* BZ = U1 --> BZ = BX*U1^(-1)
    y1$A <- y0$A - y1$B %*% U0 # AX = AZ + BZ*U0 --> AZ = AX - BZ*U0
  }

  # 3) P transformation:
  if ("P" %in% names(y0)) {
    NP <- length(y0$P$K0) #  NP may be greater than N if there are unspanned variables
    U0 <- rbind(U0, matrix(0, nrow = NP - N, ncol = 1))
    U1 <- magic::adiag(U1, diag(NP - N))
    U1_ext <- magic::adiag(U1, diag(NP - N))
    U1_ext_inv <- solve(U1_ext)

    p <- dim(y0$P$K1)[2] / dim(y0$P$K1)[1]

    y1$P <- list()
    y1$P$K1 <- array(data = NA, c(dim(y0$P$K1)))

    temp <- Reduce(`+`, lapply(seq_len(p), function(i) {
      idx <- ((i - 1) * NP + 1):(i * NP)
      y1$P$K1[, idx] <<- U1_ext %*% y0$P$K1[, idx] %*% U1_ext_inv
      y1$P$K1[, idx]
    }))


    y1$P$K0 <- U0 + U1 %*% y0$P$K0 - temp %*% U0
    N <- dim(y0$P$SS)[1]
    M <- dim(y0$P$SS)[2] / N - 1
    SSi <- array(y0$P$SS, c(N, N, M + 1))
    SSi <- mult__prod(mult__prod(as.matrix(U1), drop(SSi)), t(U1))
    y1$P$SS <- array(SSi, c(N, N * (M + 1)))
  }

  # 4) General case (if P, Q, B, A are missing)
  if (!any(names(y0) %in% c("P", "Q", "B", "A"))) {
    NP <- length(y0$K0)
    if (NP - N != 0) {
    U0 <- rbind(U0, matrix(0, nrow = NP - N, ncol = 1))
    }
    U1 <- magic::adiag(U1, diag(NP - N))

    p <- if (is.matrix(y0$K1)) dim(y0$K1)[2] / dim(y0$K1)[1] else 1

    y1$K1 <- if (is.null(dim(y0$K1))) NA else array(NA, dim(y0$K1))

    U1_ext <- magic::adiag(U1, diag(NP - N))
    U1_ext_inv <- solve(U1_ext)

    temp <- Reduce(`+`, lapply(seq_len(p), function(i) {
      if (is.matrix(y0$K1)) {
        idx <- ((i - 1) * NP + 1):(i * NP)
        y1$K1[, idx] <<- U1_ext %*% y0$K1[, idx] %*% U1_ext_inv
        y1$K1[, idx]
      } else {
        y1$K1 <<- U1_ext %*% y0$K1 %*% U1_ext_inv
        y1$K1
      }
    }))

    y1$K0 <- U0 + U1 %*% y0$K0 - temp %*% U0
    N <- if (is.null(dim(y0$SS))) 1 else dim(y0$SS)[1]
    M <- if (is.null(dim(y0$SS))) 0 else dim(y0$SS)[2] / N - 1
    SSi <- array(y0$SS, c(N, N, M + 1))
    SSi <- mult__prod(mult__prod(U1, drop(SSi)), t(U1))
    y1$SS <- array(SSi, c(N, N * (M + 1)))
  }

  return(y1)
}

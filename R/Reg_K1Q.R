#'Estimate the risk-neutral feedbak matrix K1Q using linear regressions
#'
#
#'@param Y        matrix of yields used in estimation  (J x T)
#'@param mat     vector of maturities (in years) of yields used in estimation (J x 1)
#'@param Z        pricing factors (can be yields-based or non-yields/macro variables) (N x T)
#'@param dt       time unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#'@param type     'Jordan' -> K1Q will be of the Jordan type \cr
#                 'NULL' --> no adjustment will be made
#'
#'
#'@return Risk neutral feedback matrix K1Q.
#'@references
#' This function is modified version of the "Reg_K1Q" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal



Reg_K1Q <- function(Y, mat, Z, dt,type){

  # Step 1: Define the boundaries of the maturities that will be interpolated
  M <- round(max(mat)/dt) # longest maturity adjusted for the time units of the model
  m <- round(min(mat)/dt) # shortest maturity adjusted for the time units of the model

  # Step 2: Interpolate yields by spline. Interpolation is made for each time unit of the model.
  # Define target maturities
  target_maturities <- seq(dt, M * dt, by = dt)
  b <- matrix(data=NA, ncol= ncol(Y) , nrow = M )

  # Interpolate using splinefun for each column
  for (i in seq_len(ncol(Y))) {
    spline_func <- tryCatch({
      stats::splinefun(mat, Y[, i], method = "fmm")
    }, error = function(e) {
      stop(sprintf("Error: Spline interpolation failed for column %d: %s", i, e$message))
    })

    interpolated_values <- spline_func(target_maturities)

    if (any(is.na(interpolated_values))) {
      warning(sprintf("Warning: NA values encountered during interpolation for column %d.", i))
    }

    b[, i] <- interpolated_values
  }

  R <- b

  # Step 3: regress by OLS interpolated yields onto pricing factors to obtain B: Y_t = A +B*Z_t
  Bn <- stats::lm( Jmisc::demean(t(R))~ Jmisc::demean(t(Z))-1)$coefficients
  Bn <- t(Bn)

  # Step 4: set h such that b_h will be obtained without the need of extrapolation beyond the current maturity range
  h <- m

  # Step 5: Construct the LHS from some equation
  n <- as.matrix((2*m):M, nrow =(2*m):M , ncol=1)

  d <- matrix(data=NA, ncol= ncol(Bn), nrow=nrow(n) )
  for (j in 1:ncol(Bn)){
    d[,j] <- (h/n)*Bn[h,j]
  }
  LHS <- Bn[n,] - d

  # Step 6: Construct the RHS from some equation
  RHS <- matrix(data=NA, ncol= ncol(Bn), nrow=nrow(n) )
  for (j in 1:ncol(Bn)){
    RHS[,j] <- (1-h/n)*Bn[n-h,j]
  }

  K1Qh <- qr.solve(RHS, LHS)

  # Get K^(1/h)
  eigendecomp <- eigen(K1Qh)

  # Calculate the matrix power using eigen decomposition
  eigenvalues <- eigendecomp$values
  eigenvectors <- eigendecomp$vectors

  # Take the power of the absolute values of the eigenvalues
  powered_eigenvalues <- Mod(eigenvalues)^(1/h) * exp(1i * Arg(eigenvalues)/h)
  if (ncol(Bn)==1){powered_eigenvalues <- matrix(powered_eigenvalues)} # useful for the case of one spanned factor

  # Construct K^(1/h)
  K1Q <- Re(eigenvectors %*% diag(powered_eigenvalues) %*% solve(eigendecomp$vectors))
  # Check numerical issues in final matrix
  if (any(is.nan(K1Q)) || any(is.infinite(K1Q))) {
    stop("Numerical instability detected in the K1Q matrix")
  }


  #Step 7: Convert  K1Q into Jordan form
  if (type == "Jordan"){K1Q <- Convert2JordanForm(K1Q)  }

    return(K1Q)
}

################################################################################################
#' Convert a generic matrix to its Jordan form
#
#'@param K1XQ  squared matrix in non-Jordan form
#'
#'@importFrom pracma numel
#'
#'@return squared matrix in Jordam form
#'@details
#'this Jordan matrix form handles real, complex and repeated eigenvalues.
#'
#'@references
#' This function is based on the "Convert2JordanForm" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal



Convert2JordanForm <- function(K1XQ) {


  x <- t(eigen(K1XQ)$values) # Extract the eigenvalues of matrix K1XQ.

  idx <- numeric(length(x))
  for (j in 1:length(x)){ # Exclude eigenvalues with imaginary values
    if (Im(x[j])==0){
      idx [j] <- 1
    } else{
      idx [j] <- 0
    }}

  realx <- Re(x[which(idx==1)]) # Select only the eigenvalues which are purely real


  if (numel(realx)%%2==0){
    lQ <- as.numeric(vector())
  }else{
    lQ <- realx[1]
    realx <- realx[seq_len(max(0, length(realx) - 1)) + 1]
  }


  for (h in seq_len(numel(realx)/2)){
    lQ <- rbind( lQ, 0.5*(realx[2*h-1]+realx[2*h]), (0.5*(realx[2*h-1]- realx[2*h]))^2)
  }


  imagx <- t(x[Im(x)!=0]) # Select only the eigenvalues with imaginary values
  for (h in seq_len(numel(imagx)/2)){
    lQ <- rbind(lQ, Re(imagx[2*h-1]), -abs(Im(imagx[2*h-1]))^2)
  }

  N <- numel(lQ)
  if (N==1){ K1Q <- 0  }else{ K1Q <- rbind(rep(0, N), diag(1, N - 1, N))  }
  i0 <- 0

  if (N%%2!=0){
    K1Q[1] <- lQ[1]
    lQ <- lQ[2:length(lQ)]
    i0 <- 1
  }

if (N > 1){
  for (h in seq_len(numel(lQ)/2)){
    K1Q[i0+2*h-1,i0+2*h-1] <- lQ[2*h-1]
    K1Q[i0+2*h,i0+2*h] <- lQ[2*h-1]
    K1Q[i0+2*h-1,i0+2*h] <- lQ[2*h]
  }
}

  return(K1Q)
}

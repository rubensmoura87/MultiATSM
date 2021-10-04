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
#'@importFrom stats lm splinefun
#'@importFrom powerplus Matpow
#'@importFrom pracma mldivide
#'@importFrom Jmisc  demean
#'
#'@return Risk neutral feedback matrix K1Q.
#'@examples
#'data(CM_Yields)
#'
#' Y_China <- Yields[1:6,]
#' Z_China <- Spanned_Factors(Y_China, Economies ="China", N=3)
#' mat <-c(0.25 , 0.5 , 1, 3, 5, 10)
#' dt <- 1/12
#'type <- 'Jordan'
#'Reg_K1Q(Y_China, mat, Z_China, dt, type)
#'@references
#' This function is based on the "Reg_K1Q" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@export




Reg_K1Q <- function(Y, mat, Z, dt,type){


  # Step 1: Define the boundaries of the maturities that will be interpolated
  M <- round(max(mat)/dt) # longest maturity adjusted for the time units of the model
  m <- round(min(mat)/dt) # shortest maturity adjusted for the time units of the model

  # Step 2: Interpolate yields by spline. Interpolation is made for each time unit of the model.
  b <- matrix(data=NA, ncol= ncol(Y) , nrow = M )

  for(i in 1:ncol(Y)){
    a <- splinefun(mat, t(Y[,i]), method="fmm")
    b[,i] <- t(a(seq(dt, M*dt, dt))) # each column contains one term structure in each point in time
  }

  R <- b

  # Step 3: regress by OLS interpolated yields onto pricing factors to obtain B: Y_t = A +B*Z_t
  Bn <- lm( demean(t(R))~ demean(t(Z))-1)$coefficients
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

  K1Qh <- mldivide(RHS,LHS)
  options(warn = -1)
  K1Q <- Matpow(K1Qh, numer=1, denom = h)


  #Step 7: Convert  K1Q into Jordan form
  if (type == "Jordan"){
    K1Q <- Convert2JordanForm(K1Q)
  }
  return(K1Q)

}



################################################################################################
#' Convert a generic matrix to its Jordan form
#
#'@param K1XQ  squared matrix in non-Jordan form
#'
#'@importFrom pracma numel zeros eye
#'@importFrom  wrapr seqi
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



  x <- t(eigen(K1XQ)$values) # Extract the eignevalues of matrix K1XQ.

  idx <- numeric(length(x))
  for (j in 1:length(x)){ # Exclude eignenvalues with imaginary values
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
    realx <- realx[seqi(2,length(realx))]
  }


  for (h in seqi(1,numel(realx)/2)){
    lQ <- rbind( lQ, 0.5*(realx[2*h-1]+realx[2*h]), (0.5*(realx[2*h-1]- realx[2*h]))^2)
  }


  imagx <- t(x[Im(x)!=0]) # Select only the eignenvalues with imaginary values
  for (h in seqi(1,numel(imagx)/2)){
    lQ <- rbind(lQ, Re(imagx[2*h-1]), -abs(Im(imagx[2*h-1]))^2)
  }

  N<- numel(lQ)
  K1Q <- rbind(zeros(1,N), eye(N-1,N))
  i0 <- 0

  if (N%%2!=0){
    K1Q[1] <- lQ[1]
    lQ <- lQ[2:length(lQ)]
    i0 <- 1
  }


  for (h in seqi(1,numel(lQ)/2)){
    K1Q[i0+2*h-1,i0+2*h-1] <- lQ[2*h-1]
    K1Q[i0+2*h,i0+2*h] <- lQ[2*h-1]
    K1Q[i0+2*h-1,i0+2*h] <- lQ[2*h]
  }


  return(K1Q)
}


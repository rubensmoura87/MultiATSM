# This file comprises the several inverse functions used in the computation of the Gaussian likelihood
# (see GaussianDensity function)


###############################################################################################################
#' Inverts an array of matrices so that:  inva[,,i] = inv(a[,,i])

#'@param a        matrix array (N x N x T)
#'@param whichoutput      if = 'lobabsdet' computes the log(abs(det(a))) only (text).
#'@param nargout    "nargout == 1" or "nargout == 2"(scalar)


#'@importFrom pracma numel size
#'
#'@keywords internal
#'@return "nargout == 1" returns inva -> matrix array: a^{-1} (N x N x T) \cr
#'        "nargout == 2" returns inva -> matrix array: a^{-1} (N x N x T) and  logabsdet ->  vector of log(abs(det(a)))  (1 x T)

#'@references
#' This function is modified version of the "mult__inv" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."\cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029



mult__inv <- function(a, whichoutput, nargout){


  if (!exists("whichoutput")){
    whichoutput = ""
  }

  ss <- size(a)

  if (numel(ss)<=2){
    if (nargout>1){
      inva <- solve(a)
      logabsdet <- log(abs(det(a)))
    } else {
      if (whichoutput =="logabsdet"){
        inva <- log(abs(det(a)))
      } else{
        inva <- solve(a)
      }
    }
  }

  if (numel(ss)>2){

    if (dim(a)[1]<5){
      inva <- mult_inv_small(a)$y
      dd <-  mult_inv_small(a)$dd

      if (nargout>1) {
        logabsdet <- log(abs(dd))
      } else{
        if (whichoutput == 'logabsdet'){
          inva <- log(abs(dd))
        }
      }
    } else{
      if (nargout>1){
        inva <- mult_inv_large(a)
        logabsdet <- mult_logabsdet(a)
      } else{
        if (whichoutput == 'logabsdet') {
          inva <- mult_logabsdet(a)
        } else{
          inva <- mult_inv_large(a)
        }
      }
    }
  }



  # Printing the outputs:

  if (nargout==1){
    return(inva)
  }else{
    output <- list(inva,logabsdet)
    names(output) <- c("inva", "logabsdet")
    return(output)
  }


}

#############################################################################################################
#############################################################################################################
#' Efficient computation of matrix product for arrays

#'@param a array
#'@param b array
#'@param c array

#'@keywords internal
#'
#'@details
#'#' Efficiently computes matrix product for arrays a and b and c:
#'  d[,,i] = a[,,i] b[,,i] c[,,i].
#' (efficiently = without using an inefficient loop)
#'@references
#' This function is modified version of the "mult__prod" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."\cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029



mult__prod <- function(a,b,c){


  if (!exists("c", inherits = FALSE) ||is.array(c) || is.null(c) ) { # "inherits = FALSE" ensures that we serach for variables defined with that name only.
    d <-  multiprod_2terms(a,b)
  } else {
    sa <- dim(a)
    sb <- dim(b)
    sc <- dim(c)
    if (numel(a)==1){
      d <- a*multiprod_2terms(b,c)
    } else if (numel(b)==1) {
      d <- b*multiprod_2terms(a,c)
    } else if (numel(c)==1) {
      d <- c*multiprod_2terms(a,b)
    } else if (length(sa)==2 && length(sb)==2 && length(sc)==2){
      d <- a*b*c
    } else {

      if (sa[2]!= sb[1] || sb[2]!= sc[1]) {
        stop("arrays do not commute")
      } else {
        k <- array(a, c(sa[1:2], 1, 1, sa[seqi(3,length(sa))] ))
        l <- array(b, c(1, sb[1:2], 1, sb[seqi(3,length(sb))]) )
        d <- array(0, c(dim(k),T))
        for (j in 1:T){ d[ , , , ,j] <-k*l[ , , , , j] }
        d <- colSums(d)
        m <- array(c, c(1, 1, sc))
        o <- array(0, c(dim(d)[1], 1, 1, dim(m)[length(dim(m))-1], T) )
        for (j in 1:T){
          for (p in 1:dim(d)[1]){
            for (q in 1:dim(m)[length(dim(m))-1]){
              o[ p, , , q,j] <- d[p,,,j]*m[,,,q,j] }}}
        d <- o
        sd <- dim(d)
        d <- array(d, c(sd[1], sd[seqi(4,length(sd))] ))
      }
    }
  }

  return(d)
}

#############################################################################################################
#' Inverse each 2D slice of an array (M) with arbitrary dimensions support

#'@param M multi-dimension array
#'
#'
#'@importFrom wrapr seqi
#'@keywords internal


#'@return n_D array (m x m x [p x q x  ...]), with same size as M
#'@details
#'  Inverse every 2D slice (the first two dimensions of M) for multi-dimension array M.
#'   M[,,p,q,...] * X[,,p,q,...] <- repmat(diag(m),[1,1,p,q,...])

#'@references
#' This function is based on the "mult_inv_large" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029



mult_inv_large <- function(M){


  sn <- dim(M)
  m <- sn[1]
  n <- sn[2]

  if (m!=n){
    stop("mult_inv_large: The first two dimensions of M must be m x m slices.")
  }

  p <- prod(seqi(3,length(sn)))


  RHS <- kronecker(array(1, c(p,1)), diag(m))
  X <- solve(M, RHS)


  X <- array(X, c(n, p, m))
  X <- aperm(X, c(1,3,2))
  X <- array(X, c(n,m,seqi(3,length(sn))))


  return(X)
}

#############################################################################################################
#' Inverse the (m,m,T) array of matrices for m<=4

#'@param A multi-dimensional array
#'
#'

#'@keywords internal


#'@return dd  vector of determinants (Tx1)
#'
#'@details
#' equivalent to the following for loop:
#' for (i in 1:T){ y[,,i] = inv(A[,,i]) }
#'


#'@references
#' This function is based on the "mult_inv_small" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029



mult_inv_small <- function(A){


  y <- as.numeric()
  dd <- as.numeric()
  m <- dim(A)[1]
  mm <- dim(A)


  if (m==1){
    dd <- t(t(as.vector(A)))
    y <-  1/ A
  }


  if (m==2){
    dd <- A[1,1,]*A[2,2,] - A[1,2,]*A[2,1,]
    y <- -A

    y <-  y/dd

  }

  if (m==3){
    a <- A[1,1,]
    b <- A[1,2,]
    c <- A[1,3,]
    d <- A[2,1,]
    e <- A[2,2,]
    f <- A[2,3,]
    g <- A[3,1,]
    h <- A[3,2,]
    k <- A[3,3,]
    ekfh <- e*k-f*h
    fgdk <- f*g-d*k
    dheg <- d*h-e*g
   dd <- b*fgdk + a*ekfh + c*dheg

    y <- A
    y[1,1,] <- ekfh
    y[1,2,] <- c*h-b*k
    y[1,3,] <- b*f-c*e
    y[2,1,] <- fgdk
    y[2,2,] <- a*k-c*g
    y[2,3,] <- c*d-a*f
    y[3,1,] <- dheg
    y[3,2,] <- b*g-a*h
    y[3,3,] <- a*e-b*d


    y <-  y/dd

  }


  outputs <- list(y,dd)
  names(outputs) <- c('y','dd')
  return(outputs)
}

#############################################################################################################
#' Inverse each 2D slice of an array (M) with arbitrary dimensions support

#'@param M multi-dimensional array
#'
#'@importFrom pracma lu

#'@keywords internal


#'@return X  : (n-2)_D array ([p x q x  ...])
#'
#'@details
#'  Inverse every 2D slice (the first two dimensions of M) for multi-dimension array M. In Matlab language:
#'   M[ , , p,q,...] * X[ , ,p,q,...] = repmat(diag(m),[1,1,p,q,...])
#'
#'
#'@references
#' This function is based on the "mult_logabsdet" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029



mult_logabsdet <- function(M){

  sn <- dim(M)
  m <- sn[1]
  n <- sn[2]

  if (m!=n){
    stop("mult_logabsdet: The first two dimensions of M must be m x m slices.")
  }

  p <- prod(seqi(3,length(sn)))


  u <-lu(M)$U
  X <- colSums(array(log(abs(diag(u))), c(m, p)))

  return(X)
}


#############################################################################################################
#' computes matrix product for arrays a and b:   c[,,i] = a[,,i] b[,,i]

#'@param a array (M x N x T)
#'@param b array (N x K x T)


#'@importFrom wrapr seqi
#'@return c  array (M x K x T)
#'
#'@keywords internal
#'@references
#' This function is based on the "multiprod_2terms" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029



multiprod_2terms <- function(a,b){



  if (numel(a)==1||numel(b)==1||(length(dim(a))==2 && length(dim(b))==2)){ # either a or b is scalar or matrices
    c <- a%*%b
  }else{
    sa <- dim(a)
    sb <- dim(b)

    if ((sa[1]==1 && sa[2]==1)||(sb[1]==1 && sb[2]==1)) { # either a or b is scalar
      c <- array(0, c(nrow(a), ncol(a) ,T))

      for (f in 1:T) {  c[,,f] <- a*b[,,f]   }
    } else{


      if (sa[2]!=sb[1]){
        stop('arrays do not commute')
      }else{
        d <- array(a, c(sa[1:2], 1, seqi(3,length(sa))))
        g <- array(b, c(1, sb))
        h <- array(0, c(dim(d),T))
        for (f in 1:T){    h[,,,f] <- d*g[,,,f]     } # Note: this loop does not accomodate arrays of other dimensions
        c <- colSums(h)
        sc <- dim(c)
        c <- array(c, c(sc[1:length(sc)]))
      }

    }
  }


  return(c)

}

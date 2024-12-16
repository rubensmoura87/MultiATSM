####################################################################################
#' Check whether one element is a subset of another element
#'
#' @param s1 smaller subset
#' @param s2 complete set
#'@references
#' This function is based on the "contain" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal
#'

contain <- function(s1, s2){
  y <- grepl(s1, s2)
  return(y)
}



###################################################################################
#' Transform a number bounded between a lower bound and upper bound to x by:
#'
#' @param y  Number to be transformed (scalar)
#' @param lb lower bound (scalar)
#' @param ub upper bound (scalar)
#'
#'@references
#' This function is based on the "bound2x" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal

bound2x <- function(y, lb, ub){

  x <- pos2x((y - lb)/(ub - y))
  return(x)
}

#######################################################################################################
#' Transform a positive number y to back to x by:
#'
#' @param y scalar
#'@references
#' This function is based on the "pos2x" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal

pos2x <- function(y){

  x <- y + log1p(-expm1(-y)-1)
  return(x)
}

########################################################################################################
#' Transform x to a positive number by: y = log(e^x + 1)
#'
#' @param x scalar or vector
#' @param nargout 1 or 2
#'
#'@references
#' This function is based on the "x2pos" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling)
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal

x2pos<- function(x,nargout){



  if (is.null(x) ) {
    y <- matrix(,nrow = 0, ncol=0)
  }else if (!is.null(x) & is.null(dim(x)) ) {
    y <- NA
  } else{
    y <- matrix(NA, c(dim(x)) )
  }


  if (is.null(x[x<0]) ){
    y <- NULL
  }else{
    y[x<0] <- log1p(expm1(x[x<0])+1)
  }

  if (is.null(x[x>0]) ){
    y <- NULL
  }else{
    y[x>=0] <- x[x>=0] + log1p(expm1(-x[x>=0])+1)
  }


  if (nargout>1){
    if (is.null(x) ) {
      dy <- matrix(,nrow = 0, ncol=0)
    }else if (!is.null(x) & is.null(dim(x)) ) {
      dy <- NA
    } else{
      dy <- matrix(NA, c(dim(x)) )
    }

    if (is.null(dy[x<0]) ){
      dy <- NULL
    }else{
      dy[x<0] <- (expm1(x[x<0])+1)/(expm1(x[x<0])+2)
    }


    if (is.null(dy[x>=0]) ){
      dy <- NULL
    }else{
      dy[x>=0] <- 1/(2+expm1(-x[x>=0]))
    }

  }

  if (nargout ==1){
    output <- y
  }else{    output <- list(y=y, dy=dy)  }


  return(output)
}


#########################################################################################################
#' Transform x to a number bounded btw lb and ub by:
#'
#' @param x number to be transformed (scalar)
#' @param lb lower bound (scalar)
#' @param ub upper bound (scalar)
#' @param nargout "1" or "2" (scalar)
#'
#'@references
#' This function is based on the "x2bound" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal

x2bound <-function(x, lb, ub, nargout){


  if (nargout==1){
    y <- x2pos(x,nargout = 1)
    y <- (y/(1+y))*(ub - lb) + lb
  }else{

    y  <- x2pos(x,nargout=2)$y
    dy <- x2pos(x,nargout=2)$dy

    dy <- dy*(ub-lb)/(y+1)^2
    dy <- (y/(1+y))*(ub - lb) + lb
  }


  if (nargout ==1){
    output <- y
  }else{
    output <- list(y, dy)
    names(output) <- c("y","dy")
  }

  return(output)
}


#######################################################################################################
#' Compute the square root of a matrix
#'
#'@param m squared matrix (KxK)
#'
#'@return squred matrix x (KxK) such that x%*%x = m
#'
#'@references
#' #' This function is a modified version of the "sqrtm_robust" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal



sqrtm_robust <- function(m){

  m <- as.matrix(m) # Useful for the case in which m is a scalar

  vv <- eigen(m)$vectors
  N <- nrow(vv)
  dd <- diag(N)*(sort(eigen(m)$values, decreasing = FALSE) )

  y <-pracma::sqrtm(m)$B
   if ( any(Im(y)!=0) || any(is.infinite(y))|| any(is.nan(y)) ){
    y <- mrdivide(vv%*%diag(sqrt(abs(diag(dd)))), vv ) # the y computed in this line is algebrically identical to y <- sqrtm(m).
  }


  return(y)
}


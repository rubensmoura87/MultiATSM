#' Computes numerical first order derivative of f(x)
#'
#' @param f function which contains vector (J x T) valued function handle
#' @param x parameter values
#'
#' @return
#' transformed matrix (MN x JT)
#' @references
#' This function is based on the "df__dx" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#' @export



df__dx <-function(f, x){

  h0 <- array(NaN,c(dim(x)))
  h0[] <- pmax(pmin(abs(x)*1e-3, 1e-3), 1e-6) # delta

  fx0_o <- f(x=x) # f(x)
  fis <- 0

  if (is.list(fx0_o)){
    fx0 <- f(x=x)
    fis <- 1
  }else{
    fx0 <- fx0_o
  }

  # checking if f(x+h) and f(x-h) is admissible:
  hxp <- h0; hxm <- h0

  for (i in 1:numel(x) ){
    fixed <- 0
    count <- 1
    while (!fixed && count<10){
      xp <- x
      xp[i] <- x[i]+hxp[i] # f(x+ delta)
      fxp <- f(x=xp)
      if (abs(mean((t(fxp)-t(fx0))/hxp[i]))>1e8){ # [f(x+ delta) - f(x)]/delta
        hxp[i] <- hxp[i]/2
        count <- count+1
      }else{
        fixed <- 1
      }
    }
    if (!fixed){ hxp[i] <- 0}

    fixed <- 0
    count <- 1
    while (!fixed && count<10){
      xm <- x
      xm[i] <- x[i]-hxm[i] # f(x-delta)
      fxm <- f(x=xm)
      if ( abs(mean((t(fx0) -t(fxm))/hxm[i]))>1e8){ # [f(x+ delta) - f(x)]/delta
        hxm[i] <- hxm[i]/2
        count <- count+1
      }else{
        fixed <- 1
      }
    }
    if (!fixed){ hxm[i] <-0 }
  }

  y <- matrix(NaN, nrow= numel(x), ncol= numel(f(x=x)) )


  for (i in 1:numel(x)){

    temp <- matrix(list(), 5, 5)
    for (n in 1:5){
      dxp <- matrix(0,dim(x))
      dxm <- matrix(0,dim(x))
      dxp[i] <- hxp[i]/(2^(n-1))
      dxm[i] <- hxm[i]/(2^(n-1))

      temp[[n, 1]] = (f(x=x+dxp) - f(x=x-dxm))/(dxp[i]+dxm[i])
      for (k in seqi(2,n)){
        temp[[n,k]] = ((2^(k-1))*temp[[n,k-1]] - temp[[n-1, k-1]])/(2^(k-1)-1)
      }
    }

    y[i, ] = temp[[n,n]]

  }


  if (fis==1){ # This part is needed only if fx0_o is set in a list
    fx0 <- fx0_o
    tempy <- y
    y <- NULL
    fnames <- names(fx0) # Check!
    count <- 0
    for (i in 1:numel(fnames)){
      ind <- (count+1:count+numel(fx0$fnames[i]))
      y$fnames[i] <- tempy[, ind]
    }
  }

  return(y)
}


#' Obtain the auxiliary values corresponding to each parameter, its size and its name

#'@param con If con = 'concentration' and a parameter's name contains '@@', then its auxiliary value is set to empty
#'@param varargin variable inputs used in the optimization (see "optimization" function)
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels list of necessary inputs for the estimation of JLL-based models (see \code{JLL} function)
#'@param JLLinputs     Necessary inputs for the estimation of the JLL-based models
#'
#'
#'@references
#'  This function is a modified version of the "getx" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'
#'@keywords internal



getx <- function(con, varargin, Economies, FactorLabels, JLLinputs = NULL){

  NP <- floor(sum(lengths(varargin))/4)

  # Initializing the outputs
  x <- c()
  sizex <- c()
  namex <- c()
  aIndex <- c()
  xinfofull <- list()

  countfull <- 0


  for (i in 1:NP){
    si <- strfind(varargin[[i]]$Label, ':')

    if (!isempty(si)){
      namex[i] <- killa( substr(varargin[[i]]$Label, start =1, stop =si-1) ) # remove @ from variables name, if exists
      a2t <- substr(varargin[[i]]$Label, start = si +1, stop = stringr::str_length(varargin[[i]]$Label) ) # type of the constraint of the parameter (e.g. Jordan, psd, bounded...)
      if (isempty(namex[i] )) { namex[i] <- paste('para', toString(i), sep="") }
    }else{
      namex[i] <- paste('para', toString(i), sep="")
      a2t <- varargin[[i]]$Label
    }


    temp <-  varargin[[i]]
    if ( !isempty(temp) ){

      sx <- size(temp[[1]]) # Dimensions of the variable.
      aux <- true2aux(b= varargin[[i]]$Value, ctype= a2t, lb = varargin[[i]]$LB, ub =varargin[[i]]$UB,
                      Economies, FactorLabels, JLLinputs) # Auxiliary parameters (parameters that are NOT going to be concentred out from the LLK)

      if (identical(seqi(1, numel(aux)), integer(0))  ){
        xinfofull$I[[i]] <- vector()
      }else{
        xinfofull$I[[i]] <- countfull+ seqi(1, numel(aux))
      }

      countfull <- countfull + numel(aux)

      if ( contain(con, 'concentration') && contain(temp$Label, '@') ) { # collects the parameters that will be concentrated out from the likelihood.
        a <- rbind(a, t(t((countfull+(1:numel(aux))))) )
        temp <- c()
      }else{
        if ( contain(temp$Label, 'user-defined') ){
          temp <- temp$LB(temp$Value)
        }else{
          temp <- true2aux(b = temp[[1]], ctype = temp[[2]], lb = temp[[3]], ub= temp[[4]], Economies,
                           FactorLabels, JLLinputs)
        }
      }
    }

    x <- rbind(x, t(t(as.vector(temp))) ) # collect the terms that are not concentated out from the likelihood.
    sizex <- rbind(sizex, size(temp) ) # in each row, it returns  the number of elements of parameter that are not concentrated out from the likelihood.

  }


  Outputs <- list(x, sizex, namex)
  names(Outputs) <- c("x0", "sizex", "namex")

  return(Outputs)

}


#################################################################################################################
#################################################################################################################
#################################################################################################################
#' Map constrained parameters b to unconstrained auxiliary parameters a.
#'
#' @param b             Constrained parameter
#' @param ctype         character-based vector that describes the contraints. Constraints are:
#'                      \itemize{
#'                      \item 'Jordan';
#'                      \item 'Jordan; stationary'
#'                      \item 'Jordan MultiCountry'
#'                      \item 'Jordan MultiCountry; stationary'
#'                      \item 'stationary'
#'                      \item 'psd'
#'                      \item 'BlockDiag'
#'                      \item 'bounded'
#'                      \item 'diag'
#'                      \item 'JLLstructure'
#'                        }
#' @param lb            lower bounds of b (for the bounded case). Accomodates a scalar or a matrix.
#' @param ub            upper bounds of b (for the bounded case). Accomodates a scalar or a matrix.
#' @param Economies       string-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels   string-list based which contains the labels of all the variables present in the model
#' @param JLLinputs     list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'
#' @importFrom hablar s
#' @importFrom pracma tril
#'
#'
#'@keywords internal
#'
#'@return  unconstrained auxiliary matrix.
#'
#' @references
#'  This function is a modified and extended version of the "true2aux" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029




true2aux <- function(b, ctype, lb, ub, Economies, FactorLabels, JLLinputs =NULL){



  if (contain('psd', ctype) ) {
    i <- strfind(ctype,'psd')
    i<-i[1]
    M <-as.numeric(substr(ctype, start= i+3, stop= stringr::str_length(ctype) ) )
    if (is.na(M)) { M<- 1 }
    ctype <- 'psd'
  }

  if (contain('diag', ctype) ){
    i <- strfind(ctype,'diag')
    i <- i[1]
    M <-as.numeric(substr(ctype, start= i+4, stop= stringr::str_length(ctype) ) )
    if (is.na(M)) { M<-1 }

    ctype <- 'diag'
  }

  if (contain('BlockDiag', ctype) ){
    i <- strfind(ctype,'BlockDiag')
    i <- i[1]
    M <-as.numeric(substr(ctype, start= i+9, stop= stringr::str_length(ctype) ) )
    if (is.na(M)) { M<-1 }

    ctype <- 'BlockDiag'
  }

  if (contain('JLLstructure', ctype) ){
    i <- strfind(ctype,'JLLstructure')
    i <- i[1]
    M <-as.numeric(substr(ctype, start= i+13, stop= stringr::str_length(ctype) ) )
    if (is.na(M)) { M<-1 }

    ctype <- 'JLLstructure'
  }


  # These lines of code ensure that the initial backspace from a2t is removed.
  if ( contain('bounded', ctype) ){ ctype <- 'bounded'
  } else if ( contain('Jordan MultiCountry; stationary',ctype) ) {  ctype <- 'Jordan MultiCountry; stationary'
  } else if ( contain('Jordan MultiCountry', ctype) ) {  ctype <- 'Jordan MultiCountry'
  } else if (contain('Jordan; stationary', ctype) ) { ctype <- 'Jordan; stationary'
  } else if (contain('Jordan', ctype) ) { ctype <- 'Jordan'
  } else if (contain('BlockDiag', ctype) ) { ctype <- 'BlockDiag'
  } else if (contain('JLLstructure', ctype) ) { ctype <- 'JLLstructure'
  }
  ################### CASE Jordan or Jordan stationary #########################################################
  if  (ctype ==  "Jordan" || ctype == "Jordan; stationary" ) {

    K1Q <- b
    x <- t(t(eigen(K1Q)$values))


    if ( contain('stationary', ctype) ){

      dd <- max(abs(x))
      ub <- 0.9999
      lb <- 0.0001
      maxdds <- min(max(dd, lb), ub)

      maxdd <- bound2x(maxdds, lb, ub)

      scaled <- maxdd/maxdds
      x <- x*scaled

      if (!isempty(which(is.nan(Im(x))))){
        x[which(is.nan(Im(x)))] <- Re(x[which(is.nan(Im(x)))]) + 0i # replaces the NaN in the imaginary part by 0
      }
    }


    idx <- matrix(NaN, nrow=length(x), ncol=1)
    for (j in 1:length(x)){ # Exclude eignenvalues with imaginary values
      if (Im(x[j])==0 ){
        idx [j] <- 1
      } else{
        idx [j] <- 0
      }}

    realx <- Re(x[which(idx==1)])

    # Select only the eigenvalues which are purely real
    if (numel(realx)%%2==0){  # '%%' computes the remainder from the division of two numbers
      lQ <- as.numeric(vector())
    }else{
      lQ <- realx[1]
      realx <- realx[seqi(2,length(realx))]
    }


    for (h in seqi(1,numel(realx)/2)){
      lQ <- rbind(lQ, 0.5*(realx[2*h-1]+realx[2*h]), (0.5*(realx[2*h-1]- realx[2*h]))^2)
    }

    imagx <- t(x[Im(x)!=0]) # Select only the eignenvalues with imaginary values
    for (h in seqi(1,numel(imagx)/2)){
      lQ <- rbind(lQ, Re(imagx[2*h-1]), -abs(Im(imagx[2*h-1]))^2)
    }



  a <- lQ

    ############################## Jordan MultiCountry ##################################################
  } else if  (ctype ==  "Jordan MultiCountry" || ctype == "Jordan MultiCountry; stationary" ){

    N <- length(FactorLabels$Spanned)
    C <- length(Economies)


# If the initial guess for the largest eigenvalues of all countries are greater than 1, then the numerical
# optimization of the models "JLL original" and "JLL NoDomUnit" will crash as the optimization vector will
# contain only infinities. Therefore, the code forces one of the eigenvalues to be strictly less than 1.

if (!is.null(JLLinputs)){
    IDXMaxeigen <- seq(1, N*C, by = N) # identify the indexes of the largest eigenvalue of each country
    MaxeigenCS <- diag(b)[IDXMaxeigen]

    if(all(MaxeigenCS > 0.9999)){
      BB <- which( min(MaxeigenCS) == MaxeigenCS)
      b[IDXMaxeigen[BB], IDXMaxeigen[BB]] <- 0.9998 # Replace the eigenvalue whose value is closest to 1 by 0.9998
    }
}

##

    idx0 <- 0

    for (i in 1:C){
      idx1 <-  idx0+N

      K1Q <- b[(idx0+1):idx1, (idx0+1):idx1]

      x <- t(t(eigen(K1Q)$values))


      if ( contain('stationary', ctype) ){
        dd <- max(abs(x))
        ub <- 0.9999
        lb <- 0.0001
        maxdds <- min(max(dd, lb), ub)

        maxdd <- bound2x(maxdds, lb, ub)

        scaled <- maxdd/maxdds
        x <- x*scaled

        if (!isempty(which(is.nan(Im(x))))){
          x[which(is.nan(Im(x)))] <- Re(x[which(is.nan(Im(x)))]) + 0i # replaces the NaN in the imaginary part by 0
        }

      }


      idx <- matrix(NaN, nrow=length(x), ncol=1)
      for (j in 1:length(x)){ # Exclude eignenvalues with imaginary values
        if (Im(x[j])==0){
          idx [j] <- 1
        } else{
          idx [j] <- 0
        }
      }
      realx <- Re(x[which(idx==1)])


      # Select only the eigenvalues which are purely real
      if (numel(realx)%%2==0){
        lQ <- as.numeric(vector())
      }else{
        lQ <- realx[1]
        realx <- realx[seqi(2,length(realx))]
      }


      for (h in seqi(1,numel(realx)/2)){
        lQ <- rbind(lQ, 0.5*(realx[2*h-1]+realx[2*h]), (0.5*(realx[2*h-1]- realx[2*h]))^2)

      }

      imagx <- t(x[Im(x)!=0]) # Select only the eignenvalues with imaginary values
      for (h in seqi(1,numel(imagx)/2)){
        lQ <- rbind(lQ, Re(imagx[2*h-1]), -abs(Im(imagx[2*h-1]))^2)
      }

      if (i==1) {a <- lQ
      }else{
        a <- rbind(a,lQ)
      }

      idx0 <- idx1
    }


    ############################ CASE psd ######################################################
  } else if (ctype == "psd") {

    N <- dim(b)[1]
    a <- c()

    for (i in 1:M){
      halfm <- sqrtm_robust(b[ ,(N*(i-1)+1):(N*i)]) # sqrtm (): computes a matrix square root such that Y=X^(1/2)*X^(1/2).
      MatOnes<- matrix(1 , nrow= N, ncol = N)
      a <- rbind(a, t(t(halfm[which(tril(MatOnes)==1)])) )
    }

    ############################ CASE BlockDiag ######################################################
  }  else if (ctype == "BlockDiag") {

    G <- length(FactorLabels$Global)
    K <- length(FactorLabels$Domestic)
    C <- length(Economies)

    idx0 <- 0
    a <- vector(mode="numeric")

    step <- c(G,rep(K, times=C))

    for (i in 1:(1+C)){
      idx1 <- idx0+ step[i]
      d <- as.matrix(b[seqi(idx0+1,idx1), seqi(idx0+1,idx1)])
      if (length(d)==0){ atemp <- c() } else{

      N <- dim(d)[1]
      atemp <- c()
      Mat1s <- matrix(1, nrow =N, ncol =N)
      ZeroIdx <- which(d==0) # Find the indexes of the zero elements (if any)
      for (i in 1:M){
        halfm <- sqrtm_robust(d[ ,(N*(i-1)+1):(N*i)]) # sqrtm (): computes a matrix square root such that Y=X^(1/2)*X^(1/2).
        halfm[ZeroIdx] <- 0 # if there are zero restrictions, ensure that the zeros are properly set after applying "sqrtm_robust".
        ab <- halfm[which(tril(Mat1s )==1)]
        ab <- t(t(ab[ab!=0]))
        atemp<- rbind(atemp, ab)
      }
      }
      a <- append(a,atemp)
      idx0 <- idx1
    }

    a <- t(t(a))

    ############################ JLL structure of Sigma matrix ######################################################
  }  else if (ctype == "JLLstructure") {
    N <- dim(b)[1]
    a <- c()

    C <- length(Economies)
    Nspa <- length(FactorLabels$Spanned)
    Macro <- length(FactorLabels$Domestic) - Nspa
    G <- length(FactorLabels$Global)
    K <- C*(Macro+Nspa) + G

    ZeroIdxSigmaJLL <- IDXZeroRestrictionsJLLVarCovOrtho(Macro, Nspa, G, Economies, JLLinputs$DomUnit)$VarCovOrtho
    MatOnes <- matrix(1, ncol= K, nrow = K)
    MatOnes[ZeroIdxSigmaJLL] <- 0
    IdxNONzeroSigmaJLL <- which(MatOnes!=0 & tril(MatOnes)==1)


    for (i in 1:M){
      halfm <- sqrtm_robust(b[ ,(N*(i-1)+1):(N*i)]) # sqrtm (): computes a matrix square root such that Y=X^(1/2)*X^(1/2).
      a <- rbind(a, t(t(halfm[IdxNONzeroSigmaJLL])))
    }

    ##################################### CASE 'bounded' or 'diag' ##########################################
  } else if (ctype == "bounded" || ctype== "diag"){


    if (contain('diag', ctype)) {
      J <- round(sqrt(numel(b)/M))
      bb <- matrix(NA, M*J, 1)
      for (i in 1:M){
        bb[seqi( (i-1)*J+1,i*J) ] <- diag(b[ , seqi( (i-1)*J+1, i*J) ] )
      }
      b <- bb
    }

    # make sure the constraints are not violate
    # this is important because if they are, even just by a tiny bit,
    # the outcomes may be a bunch of nan or imaginary matrices:


    if (!exists("lb") || is.null(lb) ){ lb <- -Inf }
    if (!exists("ub") || is.null(ub) ){ ub <- Inf }


    temp <- b
    temp[] <- lb[]
    lb <- temp

    temp <- b
    temp[] <- ub[]
    ub <- temp

    if (!is.null(dim(b)) || is.numeric(b)){
      d <- as.matrix(b) # d is a necessary mean to build the variables a and b with correct dimensions (see the loop below)
    }

    if (isempty(ub) || isempty(lb) || (identical(ub,Inf) & identical(lb,Inf))){
      b <- pmin(s(c(pmax(s(c(lb,b))) , ub)))
    }else{
      b <- min(s(c(max(s(c(lb,b))) , ub)))
    }

    if (is.na(b)  ){
      a <- matrix(, nrow=0, ncol=0)
    } else{
      b <- matrix(b, nrow=dim(d)[1], ncol=dim(d)[2] )
      a <- matrix(NA, nrow=dim(d)[1], ncol=dim(d)[2] )
    }

    tt <- is.infinite(lb)&  is.infinite(ub)
    a[tt] <- b[tt]

    tt <- !is.infinite(lb)  & is.infinite(ub)
    if (identical(b[tt]-lb[tt],numeric(0))){ a[tt] <- pos2x(numeric(0)) }else{ a[tt] <- pos2x(max((c(b[tt]-lb[tt],0)))) }

    tt <- is.infinite(lb) & !is.infinite(ub)
    if (identical(ub[tt]-b[tt],numeric(0))){ a[tt] <- pos2x(numeric(0)) }else{ a[tt] <- - pos2x(max((c(ub[tt]-b[tt],0)))) }

    tt<- !is.infinite(lb) & !is.infinite(ub)
    a[tt] <- bound2x(max(s(c(min(s(c(b[tt], ub[tt]))), lb[tt])) ), lb[tt], ub[tt]  )

  }


  for (j in seqi(1,numel(a)) ) {a[j] <- min(max(s(c(Re(a[j]), -1e20))), 1e20)}

  return(a)

}

##############################################################################################################
##############################################################################################################
##############################################################################################################
#' Eliminates the @
#'
#' @param s text vector containing the feature of the variable
#'
#'
#' @importFrom pracma strfind isempty
#'
#'@references
#' This function is a modified version of the "killa" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'
#'@keywords internal

killa <- function(s){

  y <- s
  ai <- strfind(s, '@')

  if (!isempty(ai) ){ y <- substr(s, start=ai+1, stop= stringr::str_length(s) ) }

  return(y)
}



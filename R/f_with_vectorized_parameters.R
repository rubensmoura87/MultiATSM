#'Use function f to generate the outputs from a ATSM
#'
#'@param x vector containing all the vectorized auxiliary parameters
#'@param sizex matrix (6x2) containing the size information of all parameters
#'@param f vector-valued objective function (function)
#'@param con if con = 'concentration', then set the value of the parameter whose name
#'          contains @@ to empty
#'@param varargin variable inputs used in the optimization (see inputs from "optimization" function)
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param JLLinputs Set of necessary inputs used in the estimation of the JLL-based models (see "JLL" function)
#'@param GVARinputs Set of necessary inputs used in the estimation of the GVAR-based models (see "GVAR" function)
#'@param nargout if nargout <- 1, returns only the values of the likelihood function.\cr
#'               If nargout <- 2, generates the entire set of outputs
#'
#'@references
#' This function is modified version of the "f_with_vectorized_parameters" function by Le and Singleton (2018).\cr
#' "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#' (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'


f_with_vectorized_parameters <- function(x, sizex, f, con, varargin, ModelType, FactorLabels, Economies,
                                         JLLinputs, GVARinputs, nargout){


  FF <- 1e9; out <- NULL; to_continue <- NULL

  if ( any(is.na(x)==0) & any(is.infinite(x) ==0) & any(Im(x)==0) ) {
    to_continue <- update_para(x, sizex, ii = NULL, con, FactorLabels, Economies, JLLinputs, GVARinputs, varargin)
  para <- getpara(to_continue) # extract only the numerical parameters of each variable that will not be concetrated out of the likelihood function.
  }


  if (nargout>1){

        out <- NULL
    tryCatch({
      if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit"){
        FF <- unlist(f(K1XQ= para[[1]],nargout = 2)$llk)
        out <- f(K1XQ= para[[1]], nargout = 2)
      }else{
        FF <- unlist(f(K1XQ= para[[1]], SSZ=para[[2]],nargout = 2)$llk)
        out <- f(K1XQ= para[[1]], SSZ=para[[2]], nargout = 2)
      }
    }, error = function(err) {
      if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit"){
        FF <- unlist(f(K1XQ= para[[1]], nargout = 1))
      }else{
        FF <- unlist(f(K1XQ= para[[1]], SSZ=para[[2]], nargout = 1))
      }
    })
  }else{

    if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit"){
      FF <- f(K1XQ= para[[1]], nargout = 1)
    }else{
      FF <- f(K1XQ= para[[1]], SSZ=para[[2]], nargout = 1)

      }
  }


  if (is.numeric(FF) & (any(is.nan(FF)) || any(is.infinite(FF)) || any(Im(FF) ==1) ) ){
    d <-numel(FF)
    FF[1:d] <- 1e9
  }

  ListOut <- list(FF, out, to_continue)
  names(ListOut) <- c("y", "out", "to_continue")

  if(nargout==1){
    return(FF)
  }else{
    return(ListOut)
  }

}

###########################################################################################################
###########################################################################################################
###########################################################################################################
#' converts the vectorized auxiliary parameter vector x to the parameters that go directly
#' into the likelihood function.


#'@param x  vector containing all the vectorized auxiliary parameters
#'@param sizex  matrix (6x2) containing the size information of all parameters
#'@param ii if empty: converts all the parameters; otherwise converts some specific parameters
#'@param con  if con = 'concentration', then set the value of the parameter whose name
#'           contains @@ to empty
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param JLLinputs Set of necessary inputs used in the estimation of the JLL-based models
#'@param GVARinputs Set of necessary inputs used in the estimation of the GVAR-based models
#'@param varargin variable inputs used in the optimization (see "Optimization" function)
#'
#'@return
#' same form as varargin, except now the parameters are updated with the values provided by the auxiliary x.
#' Importantly, by construction, all the constraints on the underlying parameters are satisfied.

#'@references
#' This function is a modified version of the "update_para" function by Le and Singleton (2018). \cr
#' "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."\cr
#' (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029
#'




update_para <-function(x, sizex, ii, con, FactorLabels, Economies, JLLinputs=NULL, GVARinputs= NULL, varargin) {


  if (!exists('ii')|| is.null(ii)){
    ud <- varargin

    NP <- floor(sum(lengths((varargin)))/4)
    counting <-0


    for (i in 1:NP){
      temp <- matrix(NaN, nrow= sizex[i,1], ncol=sizex[i,2] )
      temp[] <- x[(counting+1):(counting+numel(temp))]
      counting <- counting+numel(temp)


      if (!isempty(temp) ){
        if (contain("concentration", con) && contain("@",ud[[i]]$Label) ){
          ud[[i]]$Value <- NULL
        }else{
          ud[[i]]$Value <- aux2true(a=temp, ctype= ud[[i]]$Label, lb = ud[[i]]$LB, ub=ud[[i]]$UB,
                                    FactorLabels, Economies, JLLinputs, GVARinputs, nargout=1)
        }
      }

    }
  }else{
    pts <- t(t(c(0, cumsum(sizex[ ,1]*sizex[,2]) )))

    temp <- matrix(NaN, nrow= sizex[ii,1], ncol= sizex[ii,2]  )
    temp[1:sizex[ii,1], 1:sizex[ii,2] ] <- x[(pts[ii]+1):pts[ii+1]]


    if (!is.null(temp)){
      ud <- aux2true(a=temp, ctype= varargin[[ii]]$Label, lb = varargin[[ii]]$LB, ub=varargin[[ii]]$UB,
                     FactorLabels, Economies, JLLinputs, GVARinputs, nargout=1)
    }else{
      ud <- varargin[[ii]]$Value
    }
  }


  return(ud)
}



##############################################################################################################
#' Map auxiliary (unconstrained) parameters a to constrained parameters b

#'@param a unconstrained auxiliary parameter
#'@param ctype One of the following options:
#'             \itemize{
#'             \item 'Jordan'
#'             \item 'Jordan; stationary'
#'             \item 'Jordan MultiCountry'
#'             \item 'Jordan MultiCountry; stationary'
#'             \item 'psd';
#'             \item 'BlockDiag'
#'             \item 'bounded'
#'             \item 'diag'
#'             \item 'JLLstructure'
#'             }
#'
#'@param lb lower bounds of b (if option 'bounded' is chosen)
#'@param ub upper bounds of b (if option 'bounded' is chosen)
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param JLLinputs Inputs used in the estimation of the JLL-based models
#'@param GVARinputs Inputs used in the estimation of the GVAR-based models
#'@param nargout "nargout <- 1" returns a constrained scalar or matrix \cr
#'              "nargout <- 2" returns a list of parameters
#'
#'@importFrom pracma mod eig
#'@references
#' This function is a modified version of the "aux2true" function by Le and Singleton (2018). \cr
#' "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#' (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029
#'


aux2true <-function(a, ctype, lb, ub, FactorLabels, Economies, JLLinputs= NULL, GVARinputs= NULL, nargout){


  a <- Re(a)

  if (contain('psd', ctype) ) {
    i <- strfind(ctype,'psd')
    i<-i[1]
    M <-as.numeric(substr(ctype, start= i+3, stop= stringr::str_length(ctype) ) )
    if (is.na(M)) { M<- 1 }
    ctype <- 'psd'
  }

  if (contain('diag', ctype) ) {
    i <- strfind(ctype,'diag')
    i<-i[1]
    M <-as.numeric(substr(ctype, start= i+4, stop= stringr::str_length(ctype) ) )
    if (is.na(M)) { M<- 1 }
    ctype <- 'diag'
  }

  if (contain('BlockDiag', ctype) ) {
    i <- strfind(ctype,'BlockDiag')
    i<-i[1]
    M <-as.numeric(substr(ctype, start= i+9, stop= stringr::str_length(ctype) ) )
    if (is.na(M)) { M<- 1 }
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
  } else if ( contain('Jordan MultiCountry',ctype) ) {  ctype <- 'Jordan MultiCountry'
  } else if (contain('Jordan; stationary', ctype) ) { ctype <- 'Jordan; stationary'
  } else if (contain('Jordan', ctype) ) { ctype <- 'Jordan'
  } else if (contain('BlockDiag', ctype) ) { ctype <- 'BlockDiag'
  } else if (contain('JLLstructure', ctype) ) { ctype <- 'JLLstructure'
  }


  ################### CASE Jordan or Jordan stationary #########################################################
  if  (ctype ==  "Jordan" || ctype == "Jordan; stationary" ) {

    lQ <- a
    N <- numel(lQ)
    K1Q <- rbind(zeros(1,N), eye(N-1,N) )
    i0 <- 0
    if (mod(N,2)!=0){ K1Q[1] = lQ[1]; lQ <- t(t(lQ[2:length(lQ)])); i0<-1 }

    for (i in seqi(1,numel(lQ)/2)){
      K1Q[i0+2*i-1,i0+2*i-1] <- lQ[2*i-1]
      K1Q[i0+2*i,i0+2*i] <- lQ[2*i-1]
      K1Q[i0+2*i-1,i0+2*i] <- lQ[2*i]
    }


    if (contain('stationary', ctype)){
      dd <- eig(K1Q)
      maxdd <-  max(abs(dd))
      ub <- 0.9999
      lb <- 0.0001
      scaled <- x2bound(maxdd, lb, ub, nargout = 1)/maxdd

      lQ <- a*scaled
      if (mod(N,2)==0){
        lQ[seq(2, length(lQ), by= 2)] <- lQ[seq(2, length(lQ), by= 2)]*scaled
      }else{
        lQ[seq(3, length(lQ), by= 2)] <- lQ[seq(3, length(lQ), by= 2)]*scaled
      }

      i0 <- 0
      if (mod(N,2)!=0) { K1Q[1] <- lQ[1]; lQ <- lQ[2:length(lQ)]; i0 <- 1 }


      for (i in 1:(numel(lQ)/2) ) {
        K1Q[i0+2*i-1,i0+2*i-1] <- lQ[2*i-1]
        K1Q[i0+2*i,i0+2*i] <- lQ[2*i-1]
        K1Q[i0+2*i-1,i0+2*i] <- lQ[2*i]
      }

    }


    b <- K1Q

    ###################  Jordan MultiCountry #########################################################
  }else if (ctype ==  "Jordan MultiCountry" || ctype == "Jordan MultiCountry; stationary"){


    C <- length(Economies)
    NN <- length(FactorLabels$Spanned)

    idx0 <- 0

    for (j in 1:C){
      idx1 <- idx0 + NN

      lQ <- a[(idx0+1):idx1]
      N <- numel(lQ)
      K1Q <- rbind(zeros(1,N), eye(N-1,N) )
      i0 <- 0
      if (mod(N,2)!=0){ K1Q[1] = lQ[1]; lQ <- t(t(lQ[2:length(lQ)])); i0<-1 }

      for (i in seqi(1,numel(lQ)/2)){
        K1Q[i0+2*i-1,i0+2*i-1] <- lQ[2*i-1]
        K1Q[i0+2*i,i0+2*i] <- lQ[2*i-1]
        K1Q[i0+2*i-1,i0+2*i] <- lQ[2*i]
      }

      if (contain('stationary', ctype)){
        dd <- eig(K1Q)
        maxdd <-  max(abs(dd))
        ub <- 0.9999
        lb <- 0.0001
        scaled <- x2bound(maxdd, lb, ub, nargout = 1)/maxdd

        lQ <- a[(idx0+1):idx1]*scaled
        if (mod(N,2)==0){
          lQ[seq(2, length(lQ), by= 2)] <- lQ[seq(2, length(lQ), by= 2)]*scaled
        }else{
          lQ[seq(3, length(lQ), by= 2)] <- lQ[seq(3, length(lQ), by= 2)]*scaled
        }

        i0 <- 0
        if (mod(N,2)!=0) { K1Q[1] <- lQ[1]; lQ <- lQ[2:length(lQ)]; i0 <- 1 }


        for (i in 1:(numel(lQ)/2) ) {
          K1Q[i0+2*i-1,i0+2*i-1] <- lQ[2*i-1]
          K1Q[i0+2*i,i0+2*i] <- lQ[2*i-1]
          K1Q[i0+2*i-1,i0+2*i] <- lQ[2*i]
        }

      }


      if (j==1) {
        b <- K1Q
      }else{
        b <- magic::adiag(b, K1Q)
      }

      idx0 <- idx1
    }
    ################### CASE psd ##############################################################################
  }else if (ctype ==  "psd"){

    N <-floor(sqrt(2*numel(a)/M))
    k <- round(N*(N+1)/2)

    MatOnes <- matrix(1, nrow= N, ncol = N)
    idx <- matrix(1:(N*N), c(N, N))
    index1 <- t(t(idx[which(tril(MatOnes)==1)]))
    idx <- t(idx)
    index2 <- t(t(idx[which(tril(MatOnes)==1)]))
    convt <- zeros(N*N, N*(N+1)/2)
    convt[index1,] <- eye(N*(N+1)/2)
    convt[index2,] <- eye(N*(N+1)/2) # creates indexes to ensure that variance-covariance matrix is symmetric.

    b <- NULL
    db <- NULL

    for (i in 1:M){
      m <- matrix(convt%*%a[((i-1)*k+1):(k*i)], N, N)

      b <- cbind(b, m%*%t(m) ) # Note: the maximization is made with parameters of SSP^(1/2)
    }

    ################### CASE BlockDiag ##############################################################################
  }else if (ctype ==  "BlockDiag"){


    G <- length(FactorLabels$Global)
    K <- length(FactorLabels$Domestic)
    C <- length(Economies)


    step <- c(G*(G+1)/2,rep(K*(K+1)/2, times=C))

    # Input the zeros restrictions to the optmization vector (for the GVAR constrained models):
    if (any(GVARinputs$VARXtype == paste("constrained:", FactorLabels$Domestic))){
    # Identify the index of the constrained variable
      zz <- stringr::str_length("constrained: ")
      VarInt <- substr(GVARinputs$VARXtype, start = zz+1, stop = stringr::str_length(GVARinputs$VARXtype) )
      IdxInt <- which(FactorLabels$Domestic == VarInt)

    # Identify the indexes of the zero restrictions within the optimization vector
      MatIntVector <- matrix(NA, nrow = K, ncol=K)
      MatIntVector[IdxInt, - IdxInt] <- 0
      MatIntVector[-IdxInt, IdxInt] <- 0

      IdxZeroRest <- MatIntVector[lower.tri(MatIntVector, diag = TRUE)] # indexes of the non-zero restrictions
      IdxNonZero <- which(is.na(IdxZeroRest)) # indexes of the zero restrictions

    # Initialize the complete vector (including the zero and the non-zero restrictions)
      stepRest <- (length(a) - step[1])/C # number of parameters of the restricted vector
      ss <- rep(0, times= step[2]) # unrestricted vector of parameters per country
      tt<- rep(0, times= sum(step)) # complete unrestricted vector of parameters (per country + global)

    # Include the parameters of the marginal model
      idxI <- step[1]
       tt[seqi(1,step[1])] <- a[seqi(1,step[1])]

    # Include the parameters of the country-specific VARXs
      for (i in 1:C){
      idxF <- idxI+stepRest
      ss[IdxNonZero] <- a[(idxI+1):idxF]
      IdxFull <- sum(step[1:i])
      tt[(IdxFull+1):(IdxFull+length(ss))] <- ss
      idxI <- idxF
      }

      a <- tt
      }



    idx0 <- 0

    db <- NULL

    for (j in 1:(C+1)){

      idx1 <- idx0 + step[j]
      d <- a[seqi(idx0+1,idx1)]
if (length(d)==0 ){ btemp <- matrix(, nrow= 0, ncol= 0) } else{

      N <-floor(sqrt(2*numel(d)/M))
      k <- round(N*(N+1)/2)

      idx <- matrix(1:(N*N), c(N, N))
      MatrOne<- matrix(1, nrow = N, ncol = N)

      index1 <- t(t(idx[which(tril(MatrOne)==1)]))
      idx <- t(idx)
      index2 <- t(t(idx[which(tril(MatrOne)==1)]))
      convt <- zeros(N*N, N*(N+1)/2)
      convt[index1,] <- eye(N*(N+1)/2)
      convt[index2,] <- eye(N*(N+1)/2) # creates indexes to ensure that variance-covariance matrix is symmetric.

      btemp <- NULL
      for (i in 1:M){
        m <- matrix(convt%*%d[((i-1)*k+1):(k*i)], N, N)

        btemp <- cbind(btemp, m%*%t(m) ) # Recall that the maximization is made with parameters of SSP^(1/2)
      }
}
      if (j==1){ b <- btemp
      }else{
        b <- magic::adiag(b,btemp)
      }
      idx0 <- idx1

    }
    ################### JLL structure ##############################################################################
  }else if (ctype ==  "JLLstructure"){


    G <- length(FactorLabels$Global)
    K <- length(FactorLabels$Domestic)
    C <- length(Economies)
    Nspa <- length(FactorLabels$Spanned)
    Macro <- K - Nspa
    N <- C*K + G # Total number of factors of the model


    ZeroIdxSigmaJLL <- IDXZeroRestrictionsJLLVarCovOrtho(Macro, Nspa, G, Economies, JLLinputs$DomUnit)$VarCovOrtho
    MatOnes <- matrix(1, nrow = N, ncol = N)
    MatOnes[ZeroIdxSigmaJLL] <- 0
    IdxNONzeroSigmaJLL <- which(MatOnes!=0 & tril(MatOnes)==1)

    k <- round(N*(N+1)/2)


    # Redefine the vector of parameters to be maximized to include also the zeros restrictions
    abc <- matrix(0, N, N)
    abc[IdxNONzeroSigmaJLL] <- a # include the non-zero elements
    abc <- abc[lower.tri(abc, diag=TRUE)] # vector including the zero restrictions


    # Collect the indexes to ensure that the matrix is symmetric
    idx <- matrix(1:(N*N), c(N, N))
    Mat1s<- matrix(1, nrow = N, ncol = N)

    index1 <- t(t(idx[which(tril(Mat1s)==1)]))
    idx <- t(idx)
    index2 <- t(t(idx[which(tril( Mat1s )==1)]))
    convt <- zeros(N*N, N*(N+1)/2)
    convt[index1,] <- eye(N*(N+1)/2)
    convt[index2,] <- eye(N*(N+1)/2)


    b <- NULL
    db <- NULL

    for (i in 1:M){
      m <- matrix(convt%*%abc[((i-1)*k+1):(k*i)], N, N)


      b <- cbind(b, m%*%t(m) ) # Recall that the maximization is made with parameters of SSP^(1/2)
      # SSZ is a lower triangular matrix, hence we don't multiply m by its transpose t(m).

    }


    ################### CASE bounded ##############################################################################
  } else if (ctype ==  "bounded" || ctype == "diag"){

    if (is.null(dim(a))) {
      b <- NA
    }else{
      b <- matrix(NA, nrow=dim(a)[1], ncol=dim(a)[2] )
    }


    if (nargout>1){ db <- t(t(rep(1, times = numel(a))))  }

    if (!exists('lb') || isempty(lb) ) { lb <- -Inf    }
    if (!exists('ub') || isempty(ub) ) { ub <- Inf     }

    temp <- a; temp[]<- lb[]; lb <- temp
    temp <- a; temp[]<- ub[]; ub <- temp

    tt <- is.infinite(lb) & is.infinite(ub)
    b[tt] <- a[tt]


    tt <- (!is.infinite(lb)) & is.infinite(ub)

    if (nargout>1){
      x2p <- x2pos(a[tt], nargout=2)$y
      dx2p <- x2pos(a[tt], nargout=2)$dy

      dx2p[tt[], 1] <- dx2p

    } else{
      x2p <- x2pos(a[tt], nargout=1)
    }
    b[tt] <- x2p + lb[tt]


    tt <- is.infinite(lb) & (!is.infinite(ub))
    if (nargout>1){
      x2p  <- x2pos(-a[tt], nargout = 2)$y
      dx2p <- x2pos(-a[tt], nargout = 2)$dy
      db[tt[], 1] <- dx2p
    } else{
      x2p <- x2pos(-a[tt], nargout = 1)
    }
    b[tt] <- -x2p + ub[tt]



    tt <- (!is.infinite(lb)) & (!is.infinite(ub))
    if (nargout>1){
      b[tt] <- x2bound(a[tt], lb[tt], ub[tt], nargout=2)$y
      dx2b <- x2bound(a[tt], lb[tt], ub[tt], nargout=2)$dy

      db[tt[], 1] <- dx2b
    }else{
      b[tt] <- x2bound(a[tt], lb[tt], ub[tt], nargout = 1)
    }

  }

  if (contain('diag', ctype)){

    J <- numel(b)/M
    eyeJ <- diag(J)
    index <- which(eyeJ==1)
    bb <- matrix(NaN, nrow= J, ncol = M*J)
    if (nargout>1){ dbb <- matrix(0, nrow = J*M, ncol = M*J*J) }

    for (i in 1:M){
      bb[,((i-1)*J+1):(i*J)] <- diag(b[((i-1)*J+1):(i*J)])
      if (nargout>1){
        dbb[((i-1)*J+1):(i*J), (i-1)*J*J + index] <- diag(db[((i-1)*J+1):(i*J), 1])
      }
    }

    b <- bb
    if (nargout>1){  db <- dbb  }
  }
  else{
    if (nargout>1){
      db <- diag(db)
    }
  }


  if (nargout ==1){
    return(b)
  }else{
    return(list(b, db))
  }

}



######################################################################################################
######################################################################################################
######################################################################################################
#' Extract the parameter values from varargin
#'
#' @param varargin All parameter features
#'
#'@references
#' This function is modified version of the "getpara" function by Le and Singleton (2018). \cr
#' "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#' (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029
#'


getpara <- function(varargin){


  K <-floor(sum(lengths((varargin)))/4)
  para <- vector(mode = "list", length = K)
  for (i in 1:K){
    if (is.null(varargin[[i]]$Value)){
      break
    }
    para[[i]] <- varargin[[i]]$Value
  }

  return(para)
}

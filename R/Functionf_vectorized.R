#'Use function f to generate the outputs from a ATSM
#'
#'@param x vector containing all the vectorized auxiliary parameters
#'@param sizex matrix (6x2) containing the size information of all parameters
#'@param f vector-valued objective function (function)
#'@param con if con = 'concentration', then set the value of the parameter whose name
#'          contains @@ to empty
#'@param ListInputSet variable inputs used in the optimization (see inputs from "optimization" function)
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param JLLinputs Set of necessary inputs used in the estimation of the JLL-based models (see "JLL" function)
#'@param GVARinputs Set of necessary inputs used in the estimation of the GVAR-based models (see "GVAR" function)
#'@param WithEstimation if TRUE, returns only the values of the likelihood function, else generates the entire set of outputs
#'
#'
#'@keywords internal
#'
#'@references
#' This function is modified version of the "f_with_vectorized_parameters" function by Le and Singleton (2018).\cr
#' "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#' (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'


Functionf_vectorized <- function(x, sizex, f, con, ListInputSet, ModelType, FactorLabels, Economies, JLLinputs,
                                 GVARinputs, WithEstimation){


  FF <- 1e9; out <- NULL; ParaList_Upd <- NULL

  if ( any(is.na(x)==0) & any(is.infinite(x) ==0) & any(Im(x)==0) ) {

  ParaList_Upd <- Update_ParaList(x, sizex, con, FactorLabels, Economies, JLLinputs, GVARinputs, ListInputSet)
  Para_Upds <- getpara(ParaList_Upd) # extract only the numerical parameters of each variable that will not be concentrated out of the likelihood function.

# Gather the estimated outputs after the optimization
  if (!WithEstimation){

      if (any(ModelType == c( "JLL original", "JLL No DomUnit"))){
        FF <- unlist(f(K1XQ= Para_Upds[[1]],nargout = 2)$llk)
        out <- f(K1XQ= Para_Upds[[1]], nargout = 2)
      } else{
        FF <- unlist(f(K1XQ= Para_Upds[[1]], SSZ=Para_Upds[[2]], nargout = 2)$llk)
        out <- f(K1XQ= Para_Upds[[1]], SSZ = Para_Upds[[2]], nargout = 2)
      }

  }else{

  # Try numerical optimization.If it fails to converge, then issue an error
    tryCatch({
    if (any(ModelType == c( "JLL original", "JLL No DomUnit"))){FF <- f(K1XQ= Para_Upds[[1]], nargout = 1)
    }else{  FF <- f(K1XQ= Para_Upds[[1]], SSZ= Para_Upds[[2]], nargout = 1) }

    }, error = function(err) {
      stop("Optimization process failed to converge due to ill-defined matrices. Halting the optimization.")})


  if (is.numeric(FF) & (any(is.nan(FF)) || any(is.infinite(FF)) || any(Im(FF) ==1) ) ){
    FF[1:length(FF)] <- 1e9
  }

  }
  }

  ListOut <- list(y = FF, out = out, ParaList_Upd = ParaList_Upd)


  if(WithEstimation){ return(FF)  }else{ return(ListOut)}

}

###########################################################################################################
###########################################################################################################
#' converts the vectorized auxiliary parameter vector x to the parameters that go directly
#' into the likelihood function.


#'@param x  vector containing all the vectorized auxiliary parameters
#'@param sizex  matrix (6x2) containing the size information of all parameters
#'@param con  if con = 'concentration', then set the value of the parameter whose name
#'           contains @@ to empty
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param JLLinputs Set of necessary inputs used in the estimation of the JLL-based models
#'@param GVARinputs Set of necessary inputs used in the estimation of the GVAR-based models
#'@param ListInputSet variable inputs used in the optimization (see "Optimization" function)
#'
#'@keywords internal
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


Update_ParaList <-function(x, sizex, con, FactorLabels, Economies, JLLinputs=NULL, GVARinputs= NULL, ListInputSet) {

    COUNT <- 0

    for (i in 1:(length(ListInputSet) -1)){

      ParaTemp <- matrix(NaN, nrow= sizex[i,1], ncol=sizex[i,2] )
      ParaTemp <- x[(COUNT+1):(COUNT+length(ParaTemp))]

    if (contain("concentration", con) && contain("@",ListInputSet[[i]]$Label) ){
      ListInputSet[[i]]$Value <- NA
    }else{
      ListInputSet[[i]]$Value <- GetTruePara(ParaTemp, ListInputSet[[i]]$Label, ListInputSet[[i]]$LB,
                                             ListInputSet[[i]]$UB, FactorLabels, Economies, JLLinputs, GVARinputs)

      }
    COUNT <- COUNT + length(ParaTemp)
      }

     return(ListInputSet)
}


##############################################################################################################
#' Map auxiliary (unconstrained) parameters a to constrained parameters b
#'
#'@param ParaValue unconstrained auxiliary parameter
#'@param Const_Type_Full One of the following options:
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
#'
#'
#'@importFrom pracma mod eig
#'@references
#' This function is a modified and extended of the "aux2true" function by Le and Singleton (2018). \cr
#' "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#' (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029
#'
#'@keywords internal

GetTruePara <-function(ParaValue, Const_Type_Full, lb, ub, FactorLabels, Economies, JLLinputs= NULL, GVARinputs= NULL){

  a <- Re(ParaValue)
  Const_Type <- Adjust_Const_Type(Const_Type_Full)

  # CASE 1 : Jordan-related constraints
  if  (any(Const_Type == c("Jordan", "Jordan; stationary", "Jordan MultiCountry", "Jordan MultiCountry; stationary"))){
    b <- True_Jordan(a, Const_Type, FactorLabels, Economies)
    }

  # CASE 2: psd matrix
  else if (Const_Type ==  "psd"){  b <- True_PSD(a, Const_Type)

  # CASE 3: Block diagonal matrix
  }else if (Const_Type ==  "BlockDiag"){  b <- True_BlockDiag(a, Const_Type, FactorLabels, Economies, GVARinputs)

  # CASE 4: JLL structure of Sigma matrix
  }else if (Const_Type ==  "JLLstructure"){ b <- True_JLLstruct(a, Const_Type, FactorLabels, Economies, JLLinputs)

  # CASE 5: 'bounded' and diagonal matrix
  } else if (any(Const_Type == c( "bounded", "diag"))){ b <- True_BoundDiag(a, Const_Type, lb, ub) }


  return(b)
}

######################################################################################################
#' Extract the parameter values from the parameter list set
#'
#'@param ListInputSet All parameter features
#'
#'@references
#' This function is modified version of the "getpara" function by Le and Singleton (2018). \cr
#' "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#' (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029
#'
#'
#'@keywords internal


getpara <- function(ListInputSet){

  para <- list()
  for (i in 1:(length(ListInputSet) - 1)){para[[i]] <- ListInputSet[[i]]$Value  }

  return(para)
}

########################################################################################################
#' Transformation of the Jordan-related parameters (True form)
#'
#'@param ParaValue Constrained parameter
#'@param Const_Type Type of constraint
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal

True_Jordan <- function(ParaValue, Const_Type, FactorLabels, Economies){

  # 1) SINGLE COUNTRY SETUPS
  if  (any(Const_Type == c("Jordan", "Jordan; stationary"))) {

    lQ <- ParaValue
    N <- numel(lQ)

    if (N==1){ K1Q <- 0}else{ K1Q <- rbind(rep(0, N), diag(1, N - 1, N))}

    i0 <- 0
    if (mod(N,2)!=0){ K1Q[1] = lQ[1]; lQ <- t(t(lQ[2:length(lQ)])); i0<-1 }

    if (N > 1){
      for (i in seq_len(numel(lQ)/2)){
        K1Q[i0+2*i-1,i0+2*i-1] <- lQ[2*i-1]
        K1Q[i0+2*i,i0+2*i] <- lQ[2*i-1]
        K1Q[i0+2*i-1,i0+2*i] <- lQ[2*i]
      }
    }

    if (contain('stationary', Const_Type)){
      x <- ParaValue
      K1Q <- ImposeStat_True(x, K1Q)
}

    ##############################################
    # 2) MULTI COUNTRY SETUPS (Jordan MultiCountry)
  }else if (any(Const_Type == c("Jordan MultiCountry","Jordan MultiCountry; stationary"))){


    C <- length(Economies)
    NN <- length(FactorLabels$Spanned)

    idx0 <- 0

    for (j in 1:C){
      idx1 <- idx0 + NN

      lQ <- ParaValue[(idx0+1):idx1]
      N <- numel(lQ)
      if (N==1){ K1Q <- 0}else{ K1Q <- rbind(rep(0, N), diag(1, N - 1, N)) }
      i0 <- 0
      if (mod(N,2)!=0){ K1Q[1] = lQ[1]; lQ <- t(t(lQ[2:length(lQ)])); i0<-1 }


      if (N > 1){
        for (i in seq_len(numel(lQ)/2)){
          K1Q[i0+2*i-1,i0+2*i-1] <- lQ[2*i-1]
          K1Q[i0+2*i,i0+2*i] <- lQ[2*i-1]
          K1Q[i0+2*i-1,i0+2*i] <- lQ[2*i]
        }
      }

      if (contain('stationary', Const_Type)){
        x <- ParaValue[(idx0+1):idx1]
        K1Q <- ImposeStat_True(x, K1Q)
             }


      if (j==1) {  b <- K1Q   }else{  b <- magic::adiag(b, K1Q)   }

      idx0 <- idx1
    }

    K1Q <- b
  }

  return(K1Q)
}
###########################################################################################################
#'Makes sure that the stationary constraint under the risk-neutral measure is preserved
#'
#'@param x parameter of interest (scalar or matrix)
#'@param K1Q risk-neutral feedback matrix
#'
#'@keywords internal

ImposeStat_True <- function(x, K1Q){

  dd <- eig(K1Q)
  maxdd <-  max(abs(dd))
  ub <- 0.9999
  lb <- 0.0001
  scaled <- x2bound(maxdd, lb, ub, nargout = 1)/maxdd

  lQ <- x*scaled
  N <- numel(lQ)
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

return(K1Q)
}

##########################################################################################################
#'Transformation of a PSD matrix (true form)
#'
#'@param ParaValue Constrained parameter value
#'@param Const_Type Type of constraint
#'
#'@keywords internal

True_PSD <- function(ParaValue, Const_Type){

  # Preliminary work
  i <- strfind(Const_Type,'psd')
  i<-i[1]
  M <-as.numeric(substr(Const_Type, start= i+3, stop= stringr::str_length(Const_Type) ) )
  if (is.na(M)) { M<- 1 }

  # Rebuild true PS matrix
  N <-floor(sqrt(2*numel(ParaValue)/M))
  k <- round(N*(N+1)/2)

  MatOnes <- matrix(1, nrow= N, ncol = N)
  idx <- matrix(1:(N*N), c(N, N))
  index1 <- t(t(idx[which(tril(MatOnes)==1)]))
  idx <- t(idx)
  index2 <- t(t(idx[which(tril(MatOnes)==1)]))
  convt <- matrix(0, nrow = N*N, ncol = N*(N+1)/2)
  convt[index1,] <- diag(N*(N+1)/2)
  convt[index2,] <- diag(N*(N+1)/2) # creates indexes to ensure that variance-covariance matrix is symmetric.

  Mat_psd <- NULL
  db <- NULL

  for (i in 1:M){
    m <- matrix(convt%*%ParaValue[((i-1)*k+1):(k*i)], N, N)

    Mat_psd <- cbind(Mat_psd, m%*%t(m) ) # Note: the maximization is made with parameters of SSP^(1/2)
  }

  return(Mat_psd)
}

#########################################################################################################
#'Transformation of the block diagonal parameters (true form)
#'
#'@param ParaValue Constrained parameter
#'@param Const_Type Type of constraint
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'@param GVARinputs Inputs used in the estimation of the GVAR-based models
#'
#'@keywords internal

True_BlockDiag <- function(ParaValue, Const_Type, FactorLabels, Economies, GVARinputs){

  # Preliminary work
  i <- strfind(Const_Type,'BlockDiag')
  i<- i[1]
  M <- as.numeric(substr(Const_Type, start= i+9, stop= stringr::str_length(Const_Type) ) )
  if (is.na(M)) { M <- 1 }


  # Get the true arameter
  G <- length(FactorLabels$Global)
  K <- length(FactorLabels$Domestic)
  C <- length(Economies)


  step <- c(G*(G+1)/2,rep(K*(K+1)/2, times=C))

  # Input the zeros restrictions to the optimization vector (for the GVAR constrained models):
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
    stepRest <- (length(ParaValue) - step[1])/C # number of parameters of the restricted vector
    ss <- rep(0, times= step[2]) # unrestricted vector of parameters per country
    tt<- rep(0, times= sum(step)) # complete unrestricted vector of parameters (per country + global)

    # Include the parameters of the marginal model
    idxI <- step[1]
    tt[seq_len(step[1])] <- ParaValue[seq_len(step[1])]

    # Include the parameters of the country-specific VARXs
    for (i in 1:C){
      idxF <- idxI+stepRest
      ss[IdxNonZero] <- ParaValue[(idxI+1):idxF]
      IdxFull <- sum(step[1:i])
      tt[(IdxFull+1):(IdxFull+length(ss))] <- ss
      idxI <- idxF
    }

    ParaValue <- tt
  }



  idx0 <- 0

  db <- NULL

  for (j in 1:(C+1)){

    idx1 <- idx0 + step[j]

    if (idx0 < idx1) {seq_indices <- seq(idx0 + 1, idx1)} else {seq_indices <- integer(0) }
    d <- ParaValue[seq_indices]
    if (length(d)==0 ){ btemp <- matrix(, nrow= 0, ncol= 0) } else{

      N <-floor(sqrt(2*numel(d)/M))
      k <- round(N*(N+1)/2)

      idx <- matrix(1:(N*N), c(N, N))
      MatrOne<- matrix(1, nrow = N, ncol = N)

      index1 <- t(t(idx[which(tril(MatrOne)==1)]))
      idx <- t(idx)
      index2 <- t(t(idx[which(tril(MatrOne)==1)]))
      convt <- matrix(0, nrow = N*N, ncol = N*(N+1)/2)
      convt[index1,] <- diag(N*(N+1)/2)
      convt[index2,] <- diag(N*(N+1)/2) # creates indexes to ensure that variance-covariance matrix is symmetric.

      btemp <- NULL
      for (i in 1:M){
        m <- matrix(convt%*%d[((i-1)*k+1):(k*i)], N, N)

        btemp <- cbind(btemp, m%*%t(m) ) # Recall that the maximization is made with parameters of SSP^(1/2)
      }
    }
    if (j==1){ BD_mat <- btemp }else{  BD_mat <- magic::adiag(BD_mat,btemp)}
    idx0 <- idx1

  }

  return(BD_mat)
}

############################################################################################################
#'Transformation of the JLL-related parameters (true form)
#'
#'@param ParaValue Constrained parameter value
#'@param Const_Type Type of constraint
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param JLLinputs Inputs used in the estimation of the JLL-based models
#'
#'@keywords internal

True_JLLstruct <- function(ParaValue, Const_Type, FactorLabels, Economies, JLLinputs){

  # Preliminary work
  i <- strfind(Const_Type,'JLLstructure')
  i <- i[1]
  M <-as.numeric(substr(Const_Type, start= i+13, stop= stringr::str_length(Const_Type) ) )
  if (is.na(M)) { M<-1 }

  # Rebuild the true parameter
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
  abc[IdxNONzeroSigmaJLL] <- ParaValue # include the non-zero elements
  abc <- abc[lower.tri(abc, diag=TRUE)] # vector including the zero restrictions


  # Collect the indexes to ensure that the matrix is symmetric
  idx <- matrix(1:(N*N), c(N, N))
  Mat1s<- matrix(1, nrow = N, ncol = N)

  index1 <- t(t(idx[which(tril(Mat1s)==1)]))
  idx <- t(idx)
  index2 <- t(t(idx[which(tril( Mat1s )==1)]))
  convt <- matrix(0, nrow = N*N, ncol = N*(N+1)/2)
  convt[index1, ] <- diag(N*(N+1)/2)
  convt[index2, ] <- diag(N*(N+1)/2)


  b <- NULL
  db <- NULL

  for (i in 1:M){
    m <- matrix(convt%*%abc[((i-1)*k+1):(k*i)], N, N)
    b <- cbind(b, m%*%t(m) ) # Recall that the maximization is made with parameters of SSP^(1/2)
    # SSZ is a lower triangular matrix, hence we don't multiply m by its transpose t(m).

  }

  return(b)
}

###############################################################################################################
#'Transformation of the bounded parameters (True form)
#'
#'@param ParaValue Constrained parameter value
#'@param Const_Type Type of constraint
#'@param lb lower bound
#'@param ub upper bound
#'
#'
#'@keywords internal

True_BoundDiag <- function(ParaValue, Const_Type, lb, ub){

    # Preliminary work
    if (contain('diag', Const_Type) ) {
      i <- strfind(Const_Type,'diag')
      i<-i[1]
      M <-as.numeric(substr(Const_Type, start= i+4, stop= stringr::str_length(Const_Type) ) )
      if (is.na(M)) { M<- 1 }
    }

    # Re-build the true parameter
    if (is.null(dim(ParaValue))) { b <- NA }else{ b <- matrix(NA, nrow=dim(ParaValue)[1], ncol=dim(ParaValue)[2]) }


    if (!exists('lb') || isempty(lb) ) { lb <- -Inf    }
    if (!exists('ub') || isempty(ub) ) { ub <- Inf     }

    temp <- ParaValue; temp[]<- lb[]; lb <- temp
    temp <- ParaValue; temp[]<- ub[]; ub <- temp

    tt <- is.infinite(lb) & is.infinite(ub)
    b[tt] <- ParaValue[tt]


    tt <- (!is.infinite(lb)) & is.infinite(ub)

    x2p <- x2pos(ParaValue[tt], nargout=1)
    b[tt] <- x2p + lb[tt]


    tt <- is.infinite(lb) & (!is.infinite(ub))

    x2p <- x2pos(-ParaValue[tt], nargout = 1)
    b[tt] <- -x2p + ub[tt]


    tt <- (!is.infinite(lb)) & (!is.infinite(ub))

    b[tt] <- x2bound(ParaValue[tt], lb[tt], ub[tt], nargout = 1)


    if (contain('diag', Const_Type)){

      J <- numel(b)/M
      eyeJ <- diag(J)
      index <- which(eyeJ==1)
      bb <- matrix(NaN, nrow= J, ncol = M*J)

      for (i in 1:M){   bb[,((i-1)*J+1):(i*J)] <- diag(b[((i-1)*J+1):(i*J)])  }

      b <- bb
      }


    return(b)
  }

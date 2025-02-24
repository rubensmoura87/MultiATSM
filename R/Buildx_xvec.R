#' Obtain the auxiliary values corresponding to each parameter, its size and its name

#'@param ListInputSet variable inputs used in the optimization (see "optimization" function)
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels list of necessary inputs for the estimation of JLL-based models (see \code{JLL} function)
#'@param JLLinputs     Necessary inputs for the estimation of the JLL-based models
#'
#'
#'@references
#'  This function is a based on the "getx" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'
#'@keywords internal


Build_xvec <- function(ListInputSet, Economies, FactorLabels, JLLinputs = NULL){

  # 0) Initializing the outputs
  x <- c()
  Dim_x <- c()
  Name_x <- vector("character", length(ListInputSet) - 1)

  # 1) Build list of interest
  for (i in seq_len(length(ListInputSet) - 1)){
    si <- unlist(gregexpr(":", ListInputSet[[i]]$Label))
    Name_x[i] <- sub(".*@", "", substr(ListInputSet[[i]]$Label, start = 1, stop = si - 1)) # remove @ from variables name, if exists

    Temp_x <- GetAuxPara(ListInputSet[[i]]$Value, ListInputSet[[i]]$Label, ListInputSet[[i]]$LB,
                         ub = ListInputSet[[i]]$UB, Economies, FactorLabels, JLLinputs)

    x <- rbind(x, t(t(as.vector(Temp_x))))
    Dim_x <- rbind(Dim_x, dim(Temp_x))
  }

  Outputs <- list(x0 = x, Dim_x = Dim_x, Name_x = Name_x)
  return(Outputs)
}

#################################################################################################################
#################################################################################################################
#' Map constrained parameters b to unconstrained auxiliary parameters a.
#'
#' @param ParaValue     Constrained parameter
#' @param Const_Type_Full    character-based vector that describes the constraints. Constraints are:
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
#' @param lb            lower bounds (a scalar or a matrix) of ParaValue (for the bounded case).
#' @param ub            upper bounds (a scalar or a matrix) of ParaValue (for the bounded case).
#' @param Economies       string-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels   string-list based which contains the labels of all the variables present in the model
#' @param JLLinputs     list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
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

GetAuxPara <- function(ParaValue, Const_Type_Full, lb, ub, Economies, FactorLabels, JLLinputs = NULL){

  Const_Type <- Adjust_Const_Type(Const_Type_Full)

  # CASE 1 : Jordan-related constraints
  if (grepl("Jordan", Const_Type)){
    a <- Aux_Jordan(ParaValue, Const_Type, FactorLabels, Economies, JLLinputs)
  # CASE 2: psd matrix
  } else if (Const_Type == "psd") {
    a <- Aux_PSD(ParaValue, Const_Type)
  # CASE 3: Block diagonal matrix
  } else if (Const_Type == "BlockDiag") {
    a <- Aux_BlockDiag(ParaValue, Const_Type, FactorLabels, Economies)
  # CASE 4: JLL structure of Sigma matrix
  } else if (Const_Type == "JLLstructure") {
    a <- Aux_JLLstruct(ParaValue, Const_Type, FactorLabels, Economies, JLLinputs)
  # CASE 5: 'bounded' and diagonal matrix
  } else if (Const_Type %in% c("bounded", "diag")) {
    a <- Aux_BoundDiag(ParaValue, Const_Type, lb, ub)
  }

  for (j in seq_len(length(a)) ) {a[j] <- min(max(hablar::s(c(Re(a[j]), -1e20))), 1e20)}
  return(a)
}

#########################################################################################################
#'Adjust the constant label
#'
#'@param Const_Type_Full Type of constraint (Full Name)
#'
#'@keywords internal

Adjust_Const_Type<- function(Const_Type_Full){
  constraints <- c('Jordan MultiCountry; stationary', 'Jordan MultiCountry', 'Jordan; stationary', 'Jordan',
                   'psd', 'BlockDiag', 'JLLstructure', 'bounded', 'diag')
  for (constraint in constraints) {
    if (grepl(constraint, Const_Type_Full)) {
      return(constraint)
    }
  }
  return(Const_Type_Full)
}

##########################################################################################################
#' Transformation of the Jordan-related parameters (auxiliary form)
#'
#'@param ParaValue Constrained parameter
#'@param Const_Type Type of constraint
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'
#'@keywords internal

Aux_Jordan <- function(ParaValue, Const_Type, FactorLabels, Economies, JLLinputs){

  # 1) SINGLE COUNTRY SETUPS
  if (Const_Type %in% c("Jordan", "Jordan; stationary")){
    K1Q <- ParaValue
    x <- eigen(K1Q)$values

    # a) Impose stationary restriction (if desired)
    if (grepl('stationary', Const_Type)) {
      x <- ImposeStat_Aux(x)
    }

    # b) Compute Jordan form
    idx <- matrix(NaN, nrow=length(x), ncol=1)
    for (j in 1:length(x)){ # Exclude eigenvalues with imaginary values
      if (Im(x[j])==0 ){  idx [j] <- 1} else{  idx [j] <- 0} }

    realx <- Re(x[which(idx == 1)])

    # Select only the eigenvalues which are purely real
    if (length(realx)%%2==0){
      lQ <- c()
    }else{
      lQ <- realx[1]
      realx <- realx[seq_len(max(0, length(realx) - 1)) + 1]
    }


    for (h in seq_len(length(realx)/2)){
      lQ <- rbind(lQ, 0.5*(realx[2*h-1]+realx[2*h]), (0.5*(realx[2*h-1]- realx[2*h]))^2)
    }

    imagx <- t(x[Im(x)!=0]) # Select only the eigenvalues with imaginary values
    for (h in seq_len(length(imagx)/2)){ lQ <- rbind(lQ, Re(imagx[2*h-1]), -abs(Im(imagx[2*h-1]))^2)   }

    if(is.null(dim(lQ))){ lQ <- matrix(lQ) }

  # 2) MULTI COUNTRY SETUPS (Jordan MultiCountry)
  } else if  (any(Const_Type == c("Jordan MultiCountry", "Jordan MultiCountry; stationary"))){

    N <- length(FactorLabels$Spanned)
    C <- length(Economies)

    # Adjustment for JLL models
    if (!is.null(JLLinputs)){ ParaValue <- Jordan_JLL(ParaValue, C, N) }

    idx0 <- 0
    for (i in 1:C){
      idx1 <-  idx0 + N

      K1Q <- ParaValue[(idx0+1):idx1, (idx0+1):idx1]

      x <- t(t(eigen(K1Q)$values))

      # a) Impose stationary restriction (if desired)
      if ( grepl('stationary', Const_Type) ){ x <- ImposeStat_Aux(x) }

      # b) Compute Jordan form
      idx <- matrix(NaN, nrow=length(x), ncol=1)
      for (j in 1:length(x)){ # Exclude eigenvalues with imaginary values
        if (Im(x[j])==0){   idx [j] <- 1  } else{   idx [j] <- 0 }
      }
      realx <- Re(x[which(idx==1)])

      # Select only the eigenvalues which are purely real
      if (length(realx)%%2==0){
        lQ <- c()
      }else{
        lQ <- realx[1]
        realx <- realx[seq_len(max(0, length(realx) - 1)) + 1]
      }


      for (h in seq_len(length(realx)/2)){
        lQ <- rbind(lQ, 0.5*(realx[2*h-1]+realx[2*h]), (0.5*(realx[2*h-1]- realx[2*h]))^2)

      }

      imagx <- t(x[Im(x)!=0]) # Select only the eignenvalues with imaginary values
      for (h in seq_len(length(imagx)/2)){ lQ <- rbind(lQ, Re(imagx[2*h-1]), -abs(Im(imagx[2*h-1]))^2) }

      if (i==1) {a <- lQ }else{   a <- rbind(a,lQ) }

      idx0 <- idx1
    }

    lQ <- a
  }

  return(lQ)
}

###########################################################################################################
#'Impose stationary constraint under the risk-neutral measure
#'
#'@param yy numerical vector before imposing stationary constraint
#'
#'@keywords internal

ImposeStat_Aux <- function(yy){
  dd <- max(abs(yy))
  ub <- 0.9999
  lb <- 0.0001
  maxdds <- min(max(dd, lb), ub)
  maxdd <- bound2x(maxdds, lb, ub)

  scaled <- maxdd / maxdds
  yy <- yy * scaled

  if (any(is.nan(Im(yy)))){
    yy[is.nan(Im(yy))] <- Re(yy[is.nan(Im(yy))]) + 0i # replaces the NaN in the imaginary part by 0
  }
  return(yy)
}

##########################################################################################################
#'Check for JLL models for Jordan restrictions (auxiliary form)
#'
#'@param ParaValue Constrained parameter value
#'@param C number of countries of the economic system
#'@param N number of country-specific spanned factors
#'
#'@keywords internal

Jordan_JLL <- function(ParaValue, C, N){
  IDXMaxeigen <- seq(1, N * C, by = N) # identify the indexes of the largest eigenvalue of each country
  MaxeigenCS <- diag(ParaValue)[IDXMaxeigen]

  if (all(MaxeigenCS > 0.9999)){
    BB <- which.min(MaxeigenCS)
    ParaValue[IDXMaxeigen[BB], IDXMaxeigen[BB]] <- 0.9998 # Replace the eigenvalue whose value is closest to 1 by 0.9998
  }

  return(ParaValue)
}
##########################################################################################################
#' Transformation of a PSD matrix (auxiliary form)
#'
#'@param ParaValue Constrained parameter value
#'@param Const_Type Type of constraint
#'
#'@keywords internal

Aux_PSD <- function(ParaValue, Const_Type){

  # Preliminary work
  i <- unlist(gregexpr("psd", Const_Type))
  i <- i[1]
  M <- as.numeric(substr(Const_Type, start= i+3, stop= nchar(Const_Type) ) )
  if (is.na(M)) { M <- 1 }

  # Compute auxiliary PSD matrix
  N <- dim(ParaValue)[1]
  mat_psd <- c()

  for (i in 1:M){
    halfm <- sqrtm_robust(ParaValue[ ,(N*(i-1)+1):(N*i)]) # sqrtm (): computes a matrix square root such that Y=X^(1/2)*X^(1/2).
    MatOnes<- matrix(1 , nrow= N, ncol = N)
    mat_psd <- rbind(mat_psd, t(t(halfm[which(MatOnes == 1 & lower.tri(MatOnes, diag = TRUE))])))
  }

  return(mat_psd)
}

########################################################################################################"
#' Transformation of the JLL-related parameters (auxiliary form)
#'
#'@param ParaValue Constrained parameter value
#'@param Const_Type Type of constraint
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'
#'@keywords internal

Aux_JLLstruct <- function(ParaValue, Const_Type, FactorLabels, Economies, JLLinputs){

  # Preliminary work
  i <- unlist(gregexpr("JLLstructure", Const_Type))
  i <- i[1]
  M <-as.numeric(substr(Const_Type, start= i+13, stop= nchar(Const_Type) ) )
  if (is.na(M)) { M<-1 }

  # Transform JLL parameters
  N <- dim(ParaValue)[1]
  a <- c()

  C <- length(Economies)
  Nspa <- length(FactorLabels$Spanned)
  Macro <- length(FactorLabels$Domestic) - Nspa
  G <- length(FactorLabels$Global)
  K <- C*(Macro+Nspa) + G

  ZeroIdxSigmaJLL <- IDXZeroRestrictionsJLLVarCovOrtho(Macro, Nspa, G, Economies, JLLinputs$DomUnit)$VarCovOrtho
  MatOnes <- matrix(1, ncol= K, nrow = K)
  MatOnes[ZeroIdxSigmaJLL] <- 0
  IdxNONzeroSigmaJLL <- which(MatOnes!=0 & MatOnes == 1 & lower.tri(MatOnes, diag = TRUE))


  for (i in 1:M){
    halfm <- sqrtm_robust(ParaValue[ , (N*(i-1)+1):(N * i)]) # sqrtm (): computes a matrix square root such that Y=X^(1/2)*X^(1/2).
    a <- rbind(a, t(t(halfm[IdxNONzeroSigmaJLL])))
  }

  return(a)
}

#############################################################################################################
#' Transformation of the block diagonal parameters (auxiliary form)
#'
#'@param ParaValue Constrained parameter value
#'@param Const_Type Type of constraint
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal

Aux_BlockDiag <- function(ParaValue, Const_Type, FactorLabels, Economies) {

  # Extract M from Const_Type
  i <- regexpr("BlockDiag", Const_Type)
  M <- as.numeric(substr(Const_Type, start = i + 9, stop = nchar(Const_Type)))
  if (is.na(M)) M <- 1

  # Get dimensions
  G <- length(FactorLabels$Global)
  K <- length(FactorLabels$Domestic)
  C <- length(Economies)

  # Initialize storage
  a <- numeric(0)
  idx0 <- 0
  step <- c(G, rep(K, times = C))

  # Compute auxiliary block diagonal matrix
  for (i in seq_len(C + 1)) {
    idx1 <- idx0 + step[i]

    if (idx0 < idx1) {
      seq_indices <- seq(idx0 + 1, idx1)
      d <- as.matrix(ParaValue[seq_indices, seq_indices])

      if (length(d) > 0) {
        N <- nrow(d)
        Mat1s <- matrix(1, nrow = N, ncol = N)
        ZeroIdx <- which(d == 0)

        atemp <- do.call(rbind, lapply(seq_len(M), function(i) {
          halfm <- sqrtm_robust(d[, (N * (i - 1) + 1):(N * i)])
          halfm[ZeroIdx] <- 0
          ab <- halfm[lower.tri(Mat1s, diag = TRUE)]
          matrix(ab[ab != 0], ncol = 1)
        }))

        a <- c(a, atemp)
      }
    }

    idx0 <- idx1
  }

  return(matrix(a, ncol = 1))
}

##########################################################################################################
#'Transformation of the bounded parameters (auxiliary form)
#'
#'@param ParaValue Constrained parameter value
#'@param Const_Type Type of constraint
#'@param lb lower bound
#'@param ub upper bound
#'
#'@keywords internal

Aux_BoundDiag <- function(ParaValue, Const_Type, lb, ub){

  # Impose diagonal restriction (if desired)
  if (grepl('diag', Const_Type)) {
    i <- unlist(gregexpr("diag", Const_Type))
    i <- i[1]
    M <-as.numeric(substr(Const_Type, start= i+4, stop= nchar(Const_Type) ) )
    if (is.na(M)) { M <- 1

    J <- round(sqrt(length(ParaValue)/M))
    bb <- matrix(NA, M*J, 1)
    for (i in 1:M){
      if ((i-1)*J+1 < i*J){ IdxSeq <- ((i-1)*J+1):(i*J)}else{ IdxSeq <- numeric()}
      bb[IdxSeq ] <- diag(ParaValue[ , IdxSeq ])
    }
    ParaValue <- bb
    }
  }
  # Ensures that the constraints are preserved
  if (!exists("lb") || is.null(lb) ){ lb <- -Inf }
  if (!exists("ub") || is.null(ub) ){ ub <- Inf }


  temp <- ParaValue
  temp[] <- lb[]
  lb <- temp

  temp <- ParaValue
  temp[] <- ub[]
  ub <- temp

  if (!is.null(dim(ParaValue)) || is.numeric(ParaValue)){
    d <- as.matrix(ParaValue) # d is a necessary mean to build the variables a and b with correct dimensions (see the loop below)
  }

  if (length(ub) == 0 || length(lb) == 0  || (identical(ub,Inf) & identical(lb,Inf))){
    ParaValue <- pmin(hablar::s(c(pmax(hablar::s(c(lb,ParaValue))) , ub)))
  }else{
    ParaValue <- min(hablar::s(c(max(hablar::s(c(lb,ParaValue))) , ub)))
  }

  if (is.na(ParaValue)  ){
    a <- matrix(, nrow=0, ncol=0)
  } else{
    ParaValue <- matrix(ParaValue, nrow=dim(d)[1], ncol=dim(d)[2] )
    a <- matrix(NA, nrow=dim(d)[1], ncol=dim(d)[2] )
  }

  tt <- is.infinite(lb)&  is.infinite(ub)
  a[tt] <- ParaValue[tt]

  tt <- !is.infinite(lb)  & is.infinite(ub)
  if (identical(ParaValue[tt]-lb[tt],numeric(0))){ a[tt] <- pos2x(numeric(0)) }else{ a[tt] <- pos2x(max((c(ParaValue[tt]-lb[tt],0)))) }

  tt <- is.infinite(lb) & !is.infinite(ub)
  if (identical(ub[tt]-ParaValue[tt],numeric(0))){ a[tt] <- pos2x(numeric(0)) }else{ a[tt] <- - pos2x(max((c(ub[tt]-ParaValue[tt],0)))) }

  tt<- !is.infinite(lb) & !is.infinite(ub)
  a[tt] <- bound2x(max(hablar::s(c(min(hablar::s(c(ParaValue[tt], ub[tt]))), lb[tt])) ), lb[tt], ub[tt]  )

  return(a)
}

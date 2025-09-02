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

  if (all(!is.na(x)) & all(!is.infinite(x)) & all(Im(x) == 0)) {
    ParaList_Upd <- Update_ParaList(x, sizex, con, FactorLabels, Economies, JLLinputs, GVARinputs, ListInputSet)
    # extract only the numerical parameters of each variable that will not be concentrated out of the likelihood function.
    Para_Upds <- lapply(ParaList_Upd[-length(ParaList_Upd)], function(x) x$Value)

    # Gather the estimated outputs after the optimization
    if (!WithEstimation){
      if (any(ModelType == c("JLL original", "JLL No DomUnit"))){
        FF <- unlist(f(K1XQ = Para_Upds[[1]], nargout = 2)$llk)
        out <- f(K1XQ = Para_Upds[[1]], nargout = 2)
      } else {
        FF <- unlist(f(K1XQ = Para_Upds[[1]], SSZ = Para_Upds[[2]], nargout = 2)$llk)
        out <- f(K1XQ = Para_Upds[[1]], SSZ = Para_Upds[[2]], nargout = 2)
      }
    } else {
      # Try numerical optimization. If it fails to converge, then issue an error
      tryCatch({
        if (any(ModelType == c("JLL original", "JLL No DomUnit"))){
          FF <- f(K1XQ = Para_Upds[[1]], nargout = 1)
        } else {
          FF <- f(K1XQ = Para_Upds[[1]], SSZ = Para_Upds[[2]], nargout = 1)
        }
      }, error = function(err) {
        stop("Optimization process failed to converge due to ill-defined matrices. Halting the optimization.")
      })
    }

    if (is.numeric(FF) & (any(is.nan(FF)) || any(is.infinite(FF)) || any(Im(FF) == 1))){
      FF[1:length(FF)] <- 1e9
    }
  }

  ListOut <- list(y = FF, out = out, ParaList_Upd = ParaList_Upd)

  if (WithEstimation) {
    return(FF)
  } else {
    return(ListOut)
  }
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
#' same form as the list of parameters, except now the parameters are updated with the values provided by the auxiliary x.
#' Importantly, by construction, all the constraints on the underlying parameters are satisfied.

#'@references
#' This function is a modified version of the "update_para" function by Le and Singleton (2018). \cr
#' "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."\cr
#' (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029
#'

Update_ParaList <- function(x, sizex, con, FactorLabels, Economies, JLLinputs = NULL, GVARinputs = NULL, ListInputSet) {
  COUNT <- 0

  for (i in seq_len(length(ListInputSet) - 1)){
    ParaTemp <- matrix(NaN, nrow = sizex[i, 1], ncol = sizex[i, 2])
    ParaTemp <- x[(COUNT + 1):(COUNT + length(ParaTemp))]

    if (grepl("concentration", con) && grepl("@", ListInputSet[[i]]$Label)) {
      ListInputSet[[i]]$Value <- NA
    } else {
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
#'@references
#' This function is a modified and extended of the "aux2true" function by Le and Singleton (2018). \cr
#' "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#' (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029
#'
#'@keywords internal

GetTruePara <- function(ParaValue, Const_Type_Full, lb, ub, FactorLabels, Economies, JLLinputs = NULL, GVARinputs = NULL){

  a <- Re(ParaValue)
  Const_Type <- Adjust_Const_Type(Const_Type_Full)

  # CASE 1 : Jordan-related constraints
  if (grepl("Jordan", Const_Type)){
    b <- True_Jordan(a, Const_Type, FactorLabels, Economies)
  # CASE 2: psd matrix
  } else if (Const_Type == "psd") {
    b <- True_PSD(a, Const_Type)
  # CASE 3: Block diagonal matrix
  } else if (Const_Type == "BlockDiag") {
    b <- True_BlockDiag(a, Const_Type, FactorLabels, Economies, GVARinputs)
  # CASE 4: JLL structure of Sigma matrix
  } else if (Const_Type == "JLLstructure") {
    b <- True_JLLstruct(a, Const_Type, FactorLabels, Economies, JLLinputs)
  # CASE 5: 'bounded' and diagonal matrix
  } else if (Const_Type %in% c("bounded", "diag")) {
    b <- True_BoundDiag(a, Const_Type, lb, ub)
  }

  return(b)
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
  if (Const_Type %in% c("Jordan", "Jordan; stationary")) {
    lQ <- ParaValue
    N <- length(lQ)

    if (N == 1) {
      K1Q <- 0
    } else {
      K1Q <- rbind(rep(0, N), diag(1, N - 1, N))
    }

    i0 <- 0
    if (N %% 2 != 0) {
      K1Q[1] <- lQ[1]
      lQ <- lQ[-1]
      i0 <- 1
    }

    if (N > 1) {
      for (i in seq_len(length(lQ) / 2)) {
        K1Q[i0 + 2 * i - 1, i0 + 2 * i - 1] <- lQ[2 * i - 1]
        K1Q[i0 + 2 * i, i0 + 2 * i] <- lQ[2 * i - 1]
        K1Q[i0 + 2 * i - 1, i0 + 2 * i] <- lQ[2 * i]
      }
    }

    if (grepl('stationary', Const_Type)) {
      x <- ParaValue
      K1Q <- ImposeStat_True(x, K1Q)
    }

  # 2) MULTI COUNTRY SETUPS (Jordan MultiCountry)
  } else if (Const_Type %in% c("Jordan MultiCountry", "Jordan MultiCountry; stationary")) {
    C <- length(Economies)
    NN <- length(FactorLabels$Spanned)

    idx0 <- 0
    for (j in seq_len(C)) {
      idx1 <- idx0 + NN
      lQ <- ParaValue[(idx0 + 1):idx1]
      N <- length(lQ)

      if (N == 1) {
        K1Q <- 0
      } else {
        K1Q <- rbind(rep(0, N), diag(1, N - 1, N))
      }

      i0 <- 0
      if (N %% 2 != 0) {
        K1Q[1] <- lQ[1]
        lQ <- lQ[-1]
        i0 <- 1
      }

      if (N > 1) {
        for (i in seq_len(length(lQ) / 2)) {
          K1Q[i0 + 2 * i - 1, i0 + 2 * i - 1] <- lQ[2 * i - 1]
          K1Q[i0 + 2 * i, i0 + 2 * i] <- lQ[2 * i - 1]
          K1Q[i0 + 2 * i - 1, i0 + 2 * i] <- lQ[2 * i]
        }
      }

      if (grepl('stationary', Const_Type)) {
        x <- ParaValue[(idx0 + 1):idx1]
        K1Q <- ImposeStat_True(x, K1Q)
      }

      if (j == 1) {
        b <- K1Q
      } else {
        b <- adiag(b, K1Q)
      }

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

  dd <- eigen(K1Q)$values
  maxdd <- max(abs(dd))
  ub <- 0.9999
  lb <- 0.0001
  scaled <- x2bound(maxdd, lb, ub, nargout = 1) / maxdd

  lQ <- x * scaled
  N <- length(lQ)
  if (N %% 2 == 0) {
    lQ[seq(2, length(lQ), by = 2)] <- lQ[seq(2, length(lQ), by = 2)] * scaled
  } else {
    lQ[seq(3, length(lQ), by = 2)] <- lQ[seq(3, length(lQ), by = 2)] * scaled
  }

  i0 <- 0
  if (N %% 2 != 0) {
    K1Q[1] <- lQ[1]
    lQ <- lQ[-1]
    i0 <- 1
  }

  for (i in seq_len(length(lQ) / 2)) {
    K1Q[i0 + 2 * i - 1, i0 + 2 * i - 1] <- lQ[2 * i - 1]
    K1Q[i0 + 2 * i, i0 + 2 * i] <- lQ[2 * i - 1]
    K1Q[i0 + 2 * i - 1, i0 + 2 * i] <- lQ[2 * i]
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

  # Rebuild true PSD matrix
  N <- floor(sqrt(2 * length(ParaValue)))

  MatOnes <- matrix(1, nrow = N, ncol = N)
  idx <- matrix(1:(N * N), c(N, N))
  index1 <- t(t(idx[which(MatOnes == 1 & lower.tri(MatOnes, diag = TRUE))]))
  idx <- t(idx)
  index2 <- t(t(idx[which(MatOnes == 1 & lower.tri(MatOnes, diag = TRUE))]))
  convt <- matrix(0, nrow = N * N, ncol = N * (N + 1) / 2)
  convt[index1, ] <- diag(N * (N + 1) / 2)
  convt[index2, ] <- diag(N * (N + 1) / 2) # creates indexes to ensure that variance-covariance matrix is symmetric.

  m <- matrix(convt %*% ParaValue, N, N)
  Mat_psd <-  m %*% t(m) # Note: the maximization is made with parameters of SSP^(1/2)

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
  i <- unlist(gregexpr("BlockDiag", Const_Type))
  i <- i[1]
  M <- as.numeric(substr(Const_Type, start = i + 9, stop = nchar(Const_Type)))
  if (is.na(M)) { M <- 1 }

  # Get the true parameter
  G <- length(FactorLabels$Global)
  K <- length(FactorLabels$Domestic)
  C <- length(Economies)

  step <- c(G * (G + 1) / 2, rep(K * (K + 1) / 2, times = C))

  # Input the zeros restrictions to the optimization vector (for the GVAR constrained models):
  if (any(GVARinputs$VARXtype == paste("constrained:", FactorLabels$Domestic))){
    # Identify the index of the constrained variable
    zz <- nchar("constrained: ")
    VarInt <- substr(GVARinputs$VARXtype, start = zz + 1, stop = nchar(GVARinputs$VARXtype))
    IdxInt <- which(FactorLabels$Domestic == VarInt)

    # Identify the indexes of the zero restrictions within the optimization vector
    MatIntVector <- matrix(NA, nrow = K, ncol = K)
    MatIntVector[IdxInt, -IdxInt] <- 0
    MatIntVector[-IdxInt, IdxInt] <- 0

    IdxZeroRest <- MatIntVector[lower.tri(MatIntVector, diag = TRUE)] # indexes of the non-zero restrictions
    IdxNonZero <- which(is.na(IdxZeroRest)) # indexes of the zero restrictions

    # Initialize the complete vector (including the zero and the non-zero restrictions)
    stepRest <- (length(ParaValue) - step[1]) / C # number of parameters of the restricted vector
    ss <- rep(0, times = step[2]) # unrestricted vector of parameters per country
    tt <- rep(0, times = sum(step)) # complete unrestricted vector of parameters (per country + global)

    # Include the parameters of the marginal model
    idxI <- step[1]
    tt[seq_len(step[1])] <- ParaValue[seq_len(step[1])]

    # Include the parameters of the country-specific VARXs
    for (i in 1:C){
      idxF <- idxI + stepRest
      ss[IdxNonZero] <- ParaValue[(idxI + 1):idxF]
      IdxFull <- sum(step[1:i])
      tt[(IdxFull + 1):(IdxFull + length(ss))] <- ss
      idxI <- idxF
    }

    ParaValue <- tt
  }

  idx0 <- 0
  b <- NULL

  for (j in 1:(C + 1)){
    idx1 <- idx0 + step[j]

    if (idx0 < idx1) { seq_indices <- seq(idx0 + 1, idx1) } else { seq_indices <- integer(0) }
    d <- ParaValue[seq_indices]

    if (length(d) == 0) { btemp <- matrix(, nrow = 0, ncol = 0) } else {
      N <- floor(sqrt(2 * length(d) / M))
      k <- round(N * (N + 1) / 2)

      idx <- matrix(1:(N * N), c(N, N))
      MatrOne <- matrix(1, nrow = N, ncol = N)

      index1 <- t(t(idx[which(MatrOne == 1 & lower.tri(MatrOne, diag = TRUE))]))
      idx <- t(idx)
      index2 <- t(t(idx[which(MatrOne == 1 & lower.tri(MatrOne, diag = TRUE))]))
      convt <- matrix(0, nrow = N * N, ncol = N * (N + 1) / 2)
      convt[index1, ] <- diag(N * (N + 1) / 2)
      convt[index2, ] <- diag(N * (N + 1) / 2) # creates indexes to ensure that variance-covariance matrix is symmetric.

      btemp <- NULL
      for (i in 1:M){
        m <- matrix(convt %*% d[((i - 1) * k + 1):(k * i)], N, N)
        btemp <- cbind(btemp, m %*% t(m)) # Recall that the maximization is made with parameters of SSP^(1/2)
      }

    }
    if (j == 1) { BD_mat <- btemp } else { BD_mat <- adiag(BD_mat, btemp) }
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

  # Rebuild the true parameter
  G <- length(FactorLabels$Global)
  K <- length(FactorLabels$Domestic)
  C <- length(Economies)
  Nspa <- length(FactorLabels$Spanned)
  Macro <- K - Nspa
  R <- C * K + G # Total number of factors of the model

  ZeroIdxSigmaJLL <- IDXZeroRestrictionsJLLVarCovOrtho(Macro, Nspa, G, Economies, JLLinputs$DomUnit)$Sigma_Ye
  MatOnes <- matrix(1, nrow = R, ncol = R)
  MatOnes[ZeroIdxSigmaJLL] <- 0
  IdxNONzeroSigmaJLL <- which(MatOnes != 0 & MatOnes == 1 & lower.tri(MatOnes, diag = TRUE))

  # Redefine the vector of parameters to be maximized to include also the zeros restrictions
  abc <- matrix(0, R, R)
  abc[IdxNONzeroSigmaJLL] <- ParaValue # include the non-zero elements

  b <- abc %*% t(abc)  # NOTE: b is not identical to SSZ, but this doesn't impact the optimization

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
  if (grepl('diag', Const_Type)) {
    i <- unlist(gregexpr("diag", Const_Type))
    i <- i[1]
    M <- as.numeric(substr(Const_Type, start = i + 4, stop = nchar(Const_Type)))
    if (is.na(M)) { M <- 1 }
  }

  # Re-build the true parameter
  if (is.null(dim(ParaValue))) { b <- NA } else { b <- matrix(NA, nrow = dim(ParaValue)[1], ncol = dim(ParaValue)[2]) }

  if ( length(lb) == 0 ) { lb <- -Inf }
  if ( length(ub) == 0 ) { ub <- Inf }

  temp <- ParaValue; temp[] <- lb[]; lb <- temp
  temp <- ParaValue; temp[] <- ub[]; ub <- temp

  tt <- is.infinite(lb) & is.infinite(ub)
  b[tt] <- ParaValue[tt]

  tt <- (!is.infinite(lb)) & is.infinite(ub)
  x2p <- x2pos(ParaValue[tt], nargout = 1)
  b[tt] <- x2p + lb[tt]

  tt <- is.infinite(lb) & (!is.infinite(ub))
  x2p <- x2pos(-ParaValue[tt], nargout = 1)
  b[tt] <- -x2p + ub[tt]

  tt <- (!is.infinite(lb)) & (!is.infinite(ub))
  b[tt] <- x2bound(ParaValue[tt], lb[tt], ub[tt], nargout = 1)

  if (grepl('diag', Const_Type)){
    J <- length(b) / M
    eyeJ <- diag(J)
    index <- which(eyeJ == 1)
    bb <- matrix(NaN, nrow = J, ncol = M * J)

    for (i in 1:M) { bb[, ((i - 1) * J + 1):(i * J)] <- diag(b[((i - 1) * J + 1):(i * J)]) }

    b <- bb
  }

  return(b)
}

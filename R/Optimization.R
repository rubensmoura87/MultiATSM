#' Peform the minimization of mean(f)
#'
#' @param f vector-valued objective function (function)
#' @param ListInputSet list contain starting values and constraints:
#'                        for each input argument K (of f), we need four inputs that look like:
#'    \enumerate{
#'        \item a starting value: K0
#'        \item a variable label ('K0') followed by a ':' followed by a type of constraint. The constraint can be:
#'                  \itemize{
#'                      \item 'bounded': bounded matrix;
#'                      \item 'Jordan' or 'Jordan MultiCountry': a matrix of Jordan type;
#'                      \item 'psd': psd matrix;
#'                      \item 'stationary': largest eigenvalue of the risk-neutral feedback matrix is strictly smaller than 1;
#'                      \item 'diag' or 'BlockDiag': a diagonal or block diagonal matrix.
#'                      \item 'JLLstructure': to impose the zero-restrictions on the variance-variance matrix along
#'                              the lines of the JLL models
#'                          }
#'        \item a lower bound lb (lb <- NULL -> no lower bound)
#'        \item an upper bound ub (ub <- NULL -> no upper bound)
#'        \item Specification of the optimization settings:
#'        \itemize{
#'                  \item 'iter off': hide the printouts of the numerical optimization routines;
#'                  \item 'fminunc only': only uses fminunc for the optimization;
#'                  \item ''fminsearch only': only uses fminsearch for the optimization.
#'
#'   }
#'   }
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param ModelType A character vector indicating the model type to be estimated.
#' @param JLLinputs List. Inputs for JLL model estimation (see \code{JLL}). Default is NULL.
#' @param GVARinputs List. Inputs for GVAR model estimation (see \code{GVAR}). Default is NULL.
#' @param tol convergence tolerance (scalar). Default value is 1e-4.
#' @param TimeCount computes the required time for estimation of the model. Default is TRUE.
#'
#'@examples
#'#' # See an example of implementation in the vignette file of this package (Section 4).
#'
#'
#'
#'@references
#' This function is a conceptually based on the "LS__opt" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal

Optimization_PE <- function(f, ListInputSet, FactorLabels, Economies, ModelType, JLLinputs = NULL, GVARinputs= NULL,
                            tol= 1e-4, TimeCount = TRUE){

  # 1) Transform initial guesses of K1XQ and SSZ into auxiliary parameters that
  # will NOT be concentrated out of the log-likelihood function (llk)
  AuxVec_0 <- Build_xvec(ListInputSet, Economies, FactorLabels, JLLinputs)

  # 2) Likelihood function:
  FFvec <- function(...) Functionf_vectorized(..., sizex= AuxVec_0$Dim_x , f, con = 'concentration', ListInputSet,
                                              ModelType, FactorLabels, Economies, JLLinputs, GVARinputs,
                                              WithEstimation = TRUE)

  # 3) Optimization of the llk
  if (TimeCount) { start_time <- Sys.time() }
  AuxVec_opt <- OptimizationSetup_ATSM(AuxVec_0, FFvec, ListInputSet$OptRun, tol)
  if (TimeCount) { Optimization_Time(start_time) }

  # 4) Build the full auxiliary vector, including concentrated parameters
  Up_Temp <- Update_ParaList(AuxVec_opt$x0, sizex= AuxVec_opt$Dim_x, con= 'concentration', FactorLabels,
                         Economies, JLLinputs, GVARinputs, ListInputSet) # update the parameter set which were NOT concentrated out after the optimization.

  FF_opt <- function(...) Functionf_vectorized(..., sizex= AuxVec_opt$Dim_x, f, con = 'concentration',
                                               ListInputSet = Up_Temp, ModelType, FactorLabels, Economies,
                                               JLLinputs, GVARinputs, WithEstimation = FALSE)

  ParaLabels <- names(ListInputSet)
  Up_Temp_Full <- ParaATSM_opt_ALL(Up_Temp, FF_opt, AuxVec_opt, ParaLabels)

  # 5) Produce outputs to export:
  OutExport <- Build_xvec(Up_Temp_Full, Economies, FactorLabels, JLLinputs)
  x0_opt <-   OutExport$x0
  sizex_AllPara <- OutExport$Dim_x

  FF_export <- function(...) Functionf_vectorized(..., sizex = sizex_AllPara, f=f, con='', ListInputSet = Up_Temp_Full,
                                                  ModelType, FactorLabels, Economies, JLLinputs, GVARinputs,
                                                  WithEstimation = FALSE)

   out <- FF_export(x=x0_opt)$out

  return(out)
}


#############################################################################################################
#' Perform the optimization of the log-likelihood function of the chosen ATSM
#'
#' @param MLEinputs  A list containing the necessary inputs for building the log-likelihood function (see \code{\link{InputsForOpt}} function).
#' @param StatQ A binary variable (1 or 0) indicating whether to impose that the largest eigenvalue under Q is strictly
#'                          smaller than 1. Set to 1 to impose the restriction, or 0 otherwise.
#' @param DataFreq  A character vector specifying the data frequency. Available options: "Daily All Days", "Daily Business Days",
#'                      "Weekly", "Monthly", "Quarterly", "Annually".
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param ModelType A character vector indicating the model type to be estimated.
#' @param tol Convergence tolerance (scalar). The default is 1e-4.
#' @param TimeCount Logical. If TRUE, computes the time required for model estimation. Default is TRUE.
#' @param BS_outputs Logical. If TRUE, generates a simplified output list in the bootstrap setting. Default is FALSE.
#'
#'@examples
#' # See examples in the vignette file of this package (Section 4).
#'
#'@return
#' An object of class 'ATSMModelOutputs' containing model outputs after the optimization of the chosen ATSM specification.
#'
#' @section Available Methods:
#' - `summary(object)`
#'
#'@references
#' This function is partially adapted from the \code{LS__opt} function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@export

Optimization <- function(MLEinputs, StatQ, DataFreq, FactorLabels, Economies, ModelType, tol= 1e-4,
                            TimeCount = TRUE, BS_outputs = FALSE){

  cat("2) ATSM ESTIMATION : POINT ESTIMATE ANALYSIS \n")
  cat("2.1) Estimating ATSM ... \n")
  GVARinputs <- MLEinputs$GVARinputs
  JLLinputs <- MLEinputs$JLLinputs

# Prepare optimization
ModelParaList <- list()

if (any(ModelType ==c('JPS original', 'JPS global', "GVAR single"))){
  C <- length(Economies)
  for (i in 1:C){

    MLEinputsCS <- MLEinputs[[Economies[i]]]
    Economy <- Economies[i]
    # 1) Build the objective function
    f <- Functionf(MLEinputsCS, Economy, DataFreq, FactorLabels, ModelType, BS_outputs)
    # 2) Set the optimization settings
    VarLab <- ParaLabelsOpt(ModelType, StatQ, MLEinputsCS, BS_outputs)

    cat(paste(" ... for country:", Economies[i], ". This may take several minutes. \n"))
     ModelParaList[[ModelType]][[Economies[i]]] <- Optimization_PE(f, VarLab, FactorLabels, Economies, ModelType,
                                                                JLLinputs, GVARinputs, tol, TimeCount= TimeCount)}
} else {

  # 1) Build the objective function
  f <- Functionf(MLEinputs, Economies, DataFreq, FactorLabels, ModelType, BS_outputs)
  # 2) Set the optimization settings
  VarLab <- ParaLabelsOpt(ModelType, StatQ, MLEinputs, BS_outputs)

  cat("... This may take several minutes.\n")
  ModelParaList[[ModelType]] <- Optimization_PE(f, VarLab, FactorLabels, Economies, ModelType, JLLinputs,
                                              GVARinputs, tol, TimeCount= TimeCount)
}

# Store metadata inside the class without explicitly exporting it
attr(ModelParaList, "ModelOutInfo") <- list( Outs = ModelParaList[[ModelType]], Economies = Economies,
                                             ModelType = ModelType )

return(structure(ModelParaList, class = "ATSMModelOutputs"))
}



##############################################################################################
#' Create the variable labels used in the estimation
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param WishStationarityQ User must set "1" is she wishes to impose the largest eigenvalue under the Q to be strictly
#'                       smaller than 1. Otherwise set "0"
#'@param MLEinputs Set of inputs that are necessary to the log-likelihood function
#'@param BS_outputs Generates simplified output list in the bootstrap setting. Default is set to FALSE.
#'
#'@returns
#' list containing starting values and constraints:
#'                        for each input argument K (of f), we need four inputs that look like:
#'    \enumerate{
#'        \item a starting value: K0
#'        \item a variable label ('K0') followed by a ':' followed by a type of constraint. The constraint can be:
#'                  \itemize{
#'                      \item 'bounded': bounded matrix;
#'                      \item 'Jordan' or 'Jordan MultiCountry': a matrix of Jordan type;
#'                      \item 'psd': psd matrix;
#'                      \item 'stationary': largest eigenvalue of the risk-neutral feedback matrix is strictly smaller than 1;
#'                      \item 'diag' or 'BlockDiag': a diagonal or block diagonal matrix.
#'                      \item 'JLLstructure': to impose the zero-restrictions on the variance-voriance matrix along
#'                              the lines of the JLL models
#'                          }
#'        \item a lower bound lb (lb <- NULL -> no lower bound)
#'        \item an upper bound ub (ub <- NULL -> no upper bound)
#'        \item Specification of the optimization settings:
#'        \itemize{
#'                  \item 'iter off': hide the printouts of the numerical optimization routines;
#'                  \item 'fminunc only': only uses fminunc for the optimization;
#'                  \item ''fminsearch only': only uses fminsearch for the optimization.
#'
#'   }
#'   }
#'
#'@keywords internal

ParaLabelsOpt <- function(ModelType, WishStationarityQ, MLEinputs, BS_outputs = FALSE){

  ParaLabelsList <- list()

  ParaLabelsList[[ModelType]]$r0 <- "@r0: bounded"
  ParaLabelsList[[ModelType]]$se <- "@se: bounded"
  ParaLabelsList[[ModelType]]$K0Z <- "@K0Z: bounded"
  ParaLabelsList[[ModelType]]$K1Z <- "@K1Z: bounded"

  # K1XQ
  K1XQType <- K1XQStationary(WishStationarityQ)$SepQ

  if (ModelType %in% c("JPS original", "JPS global", 'GVAR single')) {
    ParaLabelsList[[ModelType]]$K1XQ <- K1XQStationary(WishStationarityQ)$SepQ
  } else {
    ParaLabelsList[[ModelType]]$K1XQ <- K1XQStationary(WishStationarityQ)$JointQ
  }


  # SSZ
  if (ModelType %in% c("JPS original", "JPS global", 'JPS multi')){ ParaLabelsList[[ModelType]]$SSZ <- "SSZ: psd" }
  else if (ModelType %in% c("GVAR single", 'GVAR multi')){  ParaLabelsList[[ModelType]]$SSZ <- "SSZ: BlockDiag" }
  else if (ModelType %in% c("JLL original", "JLL No DomUnit")){ ParaLabelsList[[ModelType]]$SSZ <- "@SSZ: bounded" } # Variance-covariance matrix is not estimated under Q
  else if (ModelType == "JLL joint Sigma"){ ParaLabelsList[[ModelType]]$SSZ <- "SSZ: JLLstructure" }
  # Ensures that the structure of the Variance-covariance matrix of the JLL is preserved

  # 3.2) Initial guesses for Variables that will be concentrated out of from the log-likelihood function
  K1XQ <- MLEinputs$K1XQ
  if (ModelType  %in% c("JLL original", "JLL No DomUnit")){ SSZ <- NULL} else{SSZ <- MLEinputs$SSZ}

  # Prepare list for optimization
  VarArgList <- list()
  VarArgList$K1XQ <- list(K1XQ, ParaLabelsList[[ModelType]]$K1XQ , NULL , NULL)
  VarArgList$SSZ <- list(SSZ, ParaLabelsList[[ModelType]]$SSZ, NULL, NULL)
  VarArgList$r0 <- list(NULL, ParaLabelsList[[ModelType]]$r0, NULL, NULL)
  VarArgList$se <- list(NULL, ParaLabelsList[[ModelType]]$se, 1e-6, NULL)
  VarArgList$K0Z <- list(NULL, ParaLabelsList[[ModelType]]$K0Z, NULL, NULL)
  VarArgList$K1Z <- list(NULL, ParaLabelsList[[ModelType]]$K1Z, NULL, NULL)

  if (BS_outputs){VarArgList$OptRun <- "fminunc only"} else {  VarArgList$OptRun <- "iter off"}

  LabelVar <- c('Value', 'Label', 'LB', 'UB') # Elements of each parameter
  for (d in 1:(length(VarArgList)-1)){ names(VarArgList[[d]]) <-  LabelVar}

  return(VarArgList)
}

################################################################################################################
#' Impose stationarity under the Q-measure
#'
#'@param StationaryEigenvalues Binary variable: set "1" if the user wishes the largest eigenvalue
#'                            to be strictly smaller than 1. Set "0", otherwise
#'
#'@keywords internal

K1XQStationary<- function(StationaryEigenvalues){

  K1Type <- list()

  if (StationaryEigenvalues == 1){
    K1Type$SepQ <- paste("K1XQ: ", "Jordan", "; stationary", sep="")
    K1Type$JointQ <- paste("K1XQ: ", "Jordan", " MultiCountry", "; stationary", sep="")
  }else{
    K1Type$SepQ <- paste("K1XQ: ", "Jordan", sep="")
    K1Type$JointQ <- paste("K1XQ: ", "Jordan", " MultiCountry", sep="")
  }

  return(K1Type)
}

#########################################################################################################
#'Optimization routine for the entire selected ATSM
#'
#'@param AuxVecSet List containing features for estimation of the risk-neutral parameters.
#'@param FFvec Log-likelihood function
#'@param EstType Estimation type
#'@param tol convergence tolerance (scalar). Default value is set as 1e-4.
#'
#'@keywords internal

OptimizationSetup_ATSM <- function(AuxVecSet, FFvec, EstType, tol= 1e-4){

  # 1) Optimization settings
  Max_AG_Iteration <- 1e4
  Previous_Optimal_Obj <-  -1e20
  options200 <- neldermead::optimset(MaxFunEvals = 200*length(AuxVecSet$x0), Display =  "off",
                                     MaxIter = 200, GradObj='off', TolFun= 10^-8, TolX= 10^-8)
  options1000 <- options200
  options1000$MaxIter <- 1000

  # 2) Initial checks
  converged <-(tol>1e5)
  oldF_value <- FF(AuxVecSet$x0, FFvec)

  scaling_vector <- NULL; count <- 0

  # 3) Optimization loop
  while (!converged){
    if (!grepl('fminsearch only', EstType)){

      if (!grepl('no rescaling', EstType)){
        if (length(scaling_vector) == 0){
          dFFvec <- df__dx(f=FFvec, x=AuxVecSet$x0) # first order derivative of the llK function for each point in time for the initial guess of the parameters which are NOT concentrated out of the llk.

          vv <- 1/rowMeans(abs(dFFvec))
          vv[is.infinite(vv)] <- max(vv[!is.infinite(vv)])
          vv[vv==0] <- min(vv[vv>0])
          scaling_vector <- t(t(vv))
        }

        FFtemp <- function(...) FFtemporary(..., scaling_vector = scaling_vector, FFvectorized = FFvec)
        x1 <- pracma::fminunc(x0=AuxVecSet$x0/scaling_vector, FFtemp , gr = NULL, tol = options200$TolFun,
                      maxiter = options200$MaxIter , maxfeval = options200$MaxFunEvals )

        x1 <- x1$par*scaling_vector

      } else{
        x1 <- pracma::fminunc(x0=AuxVecSet$x0, FF , gr = NULL, tol = options200$TolFun, maxiter = options200$MaxIter,
                      maxfeval = options200$MaxFunEvals)
        x1 <- x1$par
      }

      if (FF(x1, FFvec)<FF(AuxVecSet$x0, FFvec)){ AuxVecSet$x0 <- x1 }
    }
    if (!grepl('fminunc only', EstType)){
      x1<- neldermead::fminsearch(function(x) FF(x, FFvectorized = FFvec), AuxVecSet$x0, options1000)$optbase$xopt
      if (FF(x1, FFvec)<FF(AuxVecSet$x0, FFvec)){ AuxVecSet$x0 <- x1 }
    }

    newF_value <- FF(AuxVecSet$x0, FFvec)
    cat(paste("   *** Estimation round", count+1, "completed *** \n"))

    converged <-   (abs(oldF_value - newF_value)<tol) ||(count>Max_AG_Iteration && newF_value > Previous_Optimal_Obj)
    oldF_value <- newF_value

    count <- count + 1

  }
  cat('-- Done! \n')

  return(AuxVecSet)
  }
#######################################################################################################
#' Set up the vector-valued objective function (Point estimate)
#'
#' @param MLEinputs Set of inputs that are necessary to the log-likelihood function
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param DataFrequency  character-based vector: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param BS_outType Generates simplified output list in the bootstrap setting. Default is set to FALSE.
#'
#' @examples
#' # See examples in the vignette file of this package (Section 4).
#'
#' @returns
#'objective function
#'
#' @keywords internal

Functionf <- function(MLEinputs, Economies, DataFrequency, FactorLabels, ModelType, BS_outType = FALSE){

  dt <- Getdt(DataFrequency)
  mat <- MLEinputs$mat

  # 1) If one choose models in which the estimation is done country-by-country
  if (any(ModelType == c("JPS original", 'JPS global', "JPS multi", "GVAR single", "GVAR multi", "JLL joint Sigma"))){

    f <- function(...) MLEdensity(..., r0 = NULL, MLEinputs$K0Z, MLEinputs$K1Z, se= NULL, MLEinputs$Gy.0, mat,
                                  MLEinputs$Y, MLEinputs$RiskFactors, MLEinputs$SpaFact, MLEinputs$Wpca, MLEinputs$We,
                                  MLEinputs$WpcaFull, dt, Economies, FactorLabels, ModelType, MLEinputs$GVARinputs,
                                  MLEinputs$JLLinputs, BS_outType)

  } else {
    # 2) If one choose models in which the estimation is done jointly for all countries of the system with
    # the Sigma matrix estimated exclusively under the P-dynamics

    f <- function(...) MLEdensity(..., r0 = NULL, MLEinputs$SSZ, MLEinputs$K0Z, MLEinputs$K1Z, se= NULL,
                                  MLEinputs$Gy.0, mat, MLEinputs$Y, MLEinputs$RiskFactors, MLEinputs$SpaFact,
                                  MLEinputs$Wpca, MLEinputs$We, MLEinputs$WpcaFull, dt, Economies, FactorLabels,
                                  ModelType, MLEinputs$GVARinputs, MLEinputs$JLLinputs, BS_outType)
  }

  return(f)
}
#######################################################################################################
#'mean of the llk function used in the estimation of the selected ATSM
#'
#'@param x0 vector of parameters to be estimated numerically
#'@param FFvectorized log-likelihood function
#'
#'@keywords internal

FF <- function(x0, FFvectorized) {
  out <- mean(FFvectorized(x = x0))
  return(out)
  }

######################################################################################################"
#'Mean of the llk function used in the estimation of the selected ATSM
#'
#'@param xtemp temporary vector of parameters to be estimated numerically
#'@param scaling_vector scaling factor
#'@param FFvectorized log-likelihood function
#'
#'@keywords internal

FFtemporary <- function(xtemp, scaling_vector, FFvectorized){
    scaled_xtemp <- scaling_vector * xtemp

  return(FF(x0 = scaled_xtemp, FFvectorized = FFvectorized))  # Call FF with the scaled version
  }

####################################################################################################
#' Update the list of parameters
#'
#'@param Update_Temp List of model parameter features updated
#'@param FF_opt llk after optimization
#'@param AuxVecSet_opt List containing features for estimation of several parameters after optimization
#'@param ParaLabels Several variable labels
#'
#'@keywords internal

ParaATSM_opt_ALL <- function(Update_Temp, FF_opt, AuxVecSet_opt, ParaLabels){

  OptimalPara <- FF_opt(x=AuxVecSet_opt$x0)$out # computes the estimates of all parameters after the optimization.
  VarLabInt <- intersect(ParaLabels, names(OptimalPara$ests))

  Update_Temp <- stats::setNames(lapply(Update_Temp, function(elem) {

    if (is.list(elem) && grepl("@", elem$Label)) {
      namexi <- sub(".*@", "", elem$Label)  # Remove everything before '@'
      elem$Label <- namexi  # Update label

      si <- gregexpr(":", namexi)[[1]]
      if (length(si) > 0 && si[1] != -1) {
        namexi <- substr(namexi, 1, si - 1)  # Remove everything after ':'
      }

      if (namexi %in% VarLabInt) {
        elem$Value <- OptimalPara$ests[[namexi]]  # Assign correct value based on label
      }
    }

    return(elem)  # Return modified element
  }), names(Update_Temp))  # Preserve original names

  return(Update_Temp)
}

############################################################################################################
#' Computes numerical first order derivative of f(x)
#'
#' @param f function which contains vector (J x T) valued function handle
#' @param x parameter values
#'
#' @keywords internal
#'
#' @return
#' transformed matrix (MN x JT)
#' @references
#' This function is based on the "df__dx" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029

df__dx <- function(f, x) {

  h0 <- pmax(pmin(abs(x) * 1e-3, 1e-3), 1e-6) # delta

  fx0_o <- f(x = x)  # Evaluate function at original x
  fx0 <- fx0_o # Initialization

  # Checking if f(x+h) and f(x-h) is admissible
  hxp <- h0; hxm <- h0

  for (i in seq_along(x)) {
    hxp[i] <- adjust_delta(f, x, hxp[i], i, fx0, 1)
    hxm[i] <- adjust_delta(f, x, hxm[i], i, fx0, -1)
  }

  y <- matrix(NaN, nrow = length(x), ncol = length(f(x = x)))

  for (i in seq_along(x)) {
    temp <- matrix(list(), 5, 5)
    for (n in 1:5) {
      dxp <- numeric(length(x))
      dxm <- numeric(length(x))
      dxp[i] <- hxp[i] / (2^(n - 1))
      dxm[i] <- hxm[i] / (2^(n - 1))

      temp[[n, 1]] <- (f(x = x + dxp) - f(x = x - dxm)) / (dxp[i] + dxm[i])

      if (n > 1) {  # Ensure `n - 1` is always valid
        for (k in 2:n) {  # Ensure `k - 1` is valid
          temp[[n, k]] <- ((2^(k - 1)) * temp[[n, k - 1]] - temp[[n - 1, k - 1]]) / (2^(k - 1) - 1)
        }
      }
    }

    y[i, ] = temp[[5,5]]
    }

  return(y)
}

############################################################################################################
#' Adjust delta for numerical differentiation
#'
#' @param f function which contains vector (J x T) valued function handle
#' @param x parameter values
#' @param delta initial delta value
#' @param i index of the parameter being adjusted
#' @param fx0 initial function value
#' @param direction direction of adjustment (1 for positive, -1 for negative)
#'
#' @return adjusted delta value
#' @keywords internal

adjust_delta <- function(f, x, delta, i, fx0, direction) {

  fixed <- FALSE
  count <- 1
  while (!fixed && count < 10) {
    x_temp <- x
    x_temp[i] <- x[i] + direction * delta
    fx_temp <- f(x = x_temp)

    if (abs(mean((fx_temp - fx0) / delta)) > 1e8) {
      delta <- delta / 2
      count <- count + 1
    } else {
      fixed <- TRUE
    }
  }

  if (!fixed) {
    delta <- 0
  }

  return(delta)
}

############################################################################################################
#'Compute the time elapsed in the numerical optimization
#'
#'@param start_time Starting time
#'
#'@keywords internal

Optimization_Time <- function(start_time){

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
units <- c("seconds", "minutes", "hours")
scale <- c(1, 60, 3600)
idx <- max(which(elapsed >= scale))
cat(sprintf("Elapsed time: %.2f %s\n", elapsed / scale[idx], units[idx]))
}

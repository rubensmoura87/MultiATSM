#' Generates the bootstrap-related outputs
#'
#'@param ModelType   A character vector indicating the model type to be estimated.
#'@param ModelParaPE A list containing the point estimates of the model parameters. For details, refer to the outputs from the \code{\link{Optimization}} function.
#'@param NumOutPE    The point estimate derived from numerical outputs. See the outputs from the \code{\link{NumOutputs}} function for further information.
#'@param Economies   A character vector containing the names of the economies included in the system.
#'@param InputsForOutputs A list containing the necessary inputs for generating IRFs, GIRFs, FEVDs, GFEVDs and Term Premia.
#'@param FactorLabels A list of character vectors with labels for all variables in the model.
#'@param JLLlist List. Inputs for JLL model estimation (see \code{\link{JLL}} function). Default is NULL.
#'@param GVARlist List. Inputs for GVAR model estimation (see \code{\link{GVAR}} function). Default is NULL.
#'@param WishBC  Whether to estimate the physical parameter model with bias correction, based on the method by Bauer, Rudebusch and Wu (2012) (see \code{\link{Bias_Correc_VAR}} function). Default is set to 0.
#'@param BRWlist  List of necessary inputs for performing the bias-corrected estimation (see \code{\link{Bias_Correc_VAR}} function).
#'
#'
#'@examples
#' # See an example of implementation in the vignette file of this package (Section 4).
#'
#'@returns
#'list containing the following elements:
#'\itemize{
#' \item list of model parameters for one each one the draws;
#' \item list of numerical outputs (IRFs, GIRFs, FEVDs, GFEVDs and Term Premia) for each one of the draws;
#' \item Confidence bounds for the chosen level of significance.
#' }

#'@references
#' This function is a modified and extended version of the \code{VARirbound} function from "A toolbox for VAR analysis"
#' by Ambrogio Cesa-Bianchi (https://github.com/ambropo/VAR-Toolbox)
#' @export


Bootstrap <- function(ModelType, ModelParaPE, NumOutPE, Economies, InputsForOutputs, FactorLabels,
                      JLLlist = NULL, GVARlist = NULL, WishBC = 0, BRWlist = NULL){

  cat("3) BOOTSTRAP ANALYSIS \n")
  WishBoot<- InputsForOutputs[[ModelType]]$Bootstrap$WishBoot

  if (WishBoot ==0){ cat("No Bootstrap analysis was generated \n\n")
    }
  else{
    cat("3.1) Estimating bootstrap setup. This may take several hours. \n")

    StatQ <- InputsForOutputs$StationaryQ
    UMatY <- InputsForOutputs$UnitMatYields

    # 1) Pre-allocation of list of outputs
    ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
    DataFreq <- InputsForOutputs$DataFreq
    N <- length(FactorLabels$Spanned)

    ModelBootstrap <- list()
    dt <- Getdt(DataFreq)
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
      mat <- ModelParaPE[[ModelType]][[Economies[1]]]$inputs$mat }else{mat <- ModelParaPE[[ModelType]]$inputs$mat}

    ModelBootstrap <- list(GeneralInputs = list(mat= mat, dt = dt , N= N))

    # 2) Obtain the residuals from the original model
    # a) P-dynamics residuals
    residPdynOriginal <- PdynResid_BS(ModelType, Economies, ModelParaPE)

    # b) Yield residuals
      BFull_Original <- Get_BFull(ModelParaPE, FactorLabels, mat, Economies, ModelType)
      residYieOriginal <- residY_original(ModelParaPE, BFull_Original, ModelType, Economies)

      # 3) Bootstrap samples
     ## Loop over the number of draws
      tt <- 1 # numbers of accepted draws
      ww <- 1 # index for printing on screen
      start_time <- Sys.time()
      while (tt<=ndraws){ # Display number of loops
          if (tt==10*ww){
          cat(paste('Loop ', tt, ' / ', ndraws, ' draws \n'))
          ww <- ww +1
        }

        # a) generate the artificial data
    invisible(utils::capture.output(Series_artificial <- Gen_Artificial_Series(ModelParaPE, residPdynOriginal,
                                                        residYieOriginal, ModelType, BFull_Original,
                                                        InputsForOutputs, Economies, FactorLabels,
                                                        GVARlist, JLLlist, WishBC, BRWlist)))

      # b) Prepare the inputs the model estimation
      Y_BS <- t(Series_artificial$Y_BS)
      Global_BS <- t(Series_artificial$GlobalMacro_BS)
      Dom_BS <- t(Series_artificial$DomesticMacro_BS)
      t0_BS <- colnames(Global_BS)[1]
      tF_BS <- utils::tail(colnames(Global_BS), 1)

      invisible(utils::capture.output(ATSMInputs_BS <- InputsForOpt(t0_BS, tF_BS, ModelType, Y_BS, Global_BS, Dom_BS, FactorLabels,
                                           Economies, DataFreq, GVARlist, JLLlist, WishBC,
                                           BRWlist, UMatY, CheckInputs= FALSE, BS_Adj= TRUE)))

      # c) Run the optimization
      invisible(utils::capture.output(Draw_Opt <- Optimization(ATSMInputs_BS, StatQ, DataFreq, FactorLabels,
                                                         Economies, ModelType, tol = 1e-1, TimeCount=F,
                                                         BS_outputs = TRUE)))

      ModelBootstrap <- AdjustOptm_BS(ModelType, ModelBootstrap, Draw_Opt, Economies, tt)

      # if the optimization crashes after some draws, we can keep the outputs of draws before
      saveRDS(ModelBootstrap, paste(tempdir(),"/Bootstrap_", InputsForOutputs$'Label Outputs','.rds',sep=""))
      tt<-tt+1
      }

    cat('-- Done! \n')
    Optimization_Time(start_time)

    # 4) Compute the numerical outputs from the bootstrap samples
    cat("3.2) Computing numerical outputs. \n")
    ModelBootstrap$NumOutDraws <- NumOutputs_Bootstrap(ModelType, ModelBootstrap, InputsForOutputs, FactorLabels,
                                                       Economies)

    # 5) Compute confidence bounds
    cat("3.3) Computing confidence bounds and producing graphical outputs. \n")
    ModelBootstrap$ConfBounds <- BootstrapBoundsSet(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies)
    ## To save space, clean the repeated outputs from the JLL outputs
    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){ModelBootstrap <- CleanOrthoJLL_Boot(ModelBootstrap, ndraws, ModelType)}

    saveRDS(ModelBootstrap, paste(tempdir(),"/Bootstrap_", InputsForOutputs$'Label Outputs','.rds',sep=""))

    return(ModelBootstrap)
  }
}
################################################################################################################
#' Compute some key parameters from the P-dynamics (Bootstrap set)
#'
#'@param ModelType    string-vector containing the label of the model to be estimated
#'@param Economies    A character vector containing the names of the economies included in the system.
#'@param ModelPara_PE point estimate from the model parameters
#'
#'@keywords internal

PdynResid_BS <- function(ModelType, Economies, ModelPara_PE){

  # SepQ models
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))){

    eZ <- list()

    for (i in 1:length(Economies)){
      ZZ <- ModelPara_PE[[ModelType]][[Economies[i]]]$inputs$AllFactors # K x T
      K0Z <- ModelPara_PE[[ModelType]][[Economies[i]]]$ests$K0Z
      K1Z <- ModelPara_PE[[ModelType]][[Economies[i]]]$ests$K1Z
      T <- ncol(ZZ)
      eZ_row <- ZZ[ ,2:T] - matrix(K0Z, nrow= nrow(K0Z),  ncol= T-1) - K1Z%*%ZZ[ , 1:(T-1)] # T x K
      eZ[[Economies[i]]] <- t(eZ_row)
    }
    # JointQ models
  }else{
    ZZ <- ModelPara_PE[[ModelType]]$inputs$AllFactors
    K0Z <- ModelPara_PE[[ModelType]]$ests$K0Z
    K1Z <- ModelPara_PE[[ModelType]]$ests$K1Z
    T <- ncol(ZZ)
    eZ <- ZZ[ ,2:T] - matrix(K0Z, nrow= nrow(K0Z),  ncol= T-1) - K1Z%*%ZZ[ , 1:(T-1)] # T x K
    eZ <- t(eZ)
  }

  return(eZ)
}


#################################################################################################################
#' Compute the residuals from the original model
#'
#'@param residPdynOriginal Time-series of the residuals from the P-dynamics equation (T x F)
#'@param residYieOriginal Time-series of the residuals from the observational equation (T x J or T x CJ)
#'@param InputsForOutputs list containing the desired inputs for the construction of the numerical outputs.
#'@param ModelType A character vector indicating the model type to be estimated
#'@param nlag Number of lags in the P-dynamics. Default is set to 1.
#'
#'@keywords internal

ResampleResiduals_BS <- function(residPdynOriginal, residYieOriginal, InputsForOutputs, ModelType, nlag= 1){

  methodBS <- InputsForOutputs[[ModelType]]$Bootstrap$methodBS
  BlockLength <- InputsForOutputs[[ModelType]]$Bootstrap$BlockLength

  T <- nrow(residYieOriginal)
  K <- ncol(residPdynOriginal)

  if (methodBS == "bs"){
    # Use the residuals to bootstrap: generate a random number bounded
    # between 0 and the number of residuals, then use the ceiling function to select
    # that row of the residuals (this is equivalent to sampling with replacement)
    Rand <- matrix(stats::runif(T, min = 0, max = 1), ncol = 1)
    rr <- ceiling((T-nlag)*Rand) # T x 1
    uPdyn <- residPdynOriginal[rr[1:(T-nlag)], ] # (T-1) x K
    uYiel <- residYieOriginal[rr, ] # TxJ  or T x CJ
  } else if ( methodBS == "wild"){

    # Wild bootstrap based on simple distribution (~Rademacher)
    Rand <- matrix(stats::runif(T, min = 0, max = 1), ncol = 1)
    rr <- 1 - 2*(Rand > 0.5)
    uPdyn<- residPdynOriginal*(rr[1:(T-nlag)]%*%matrix(1,nrow=1, ncol=K))
    uYiel <- residYieOriginal*(rr%*%matrix(1,nrow=1, ncol= ncol(residYieOriginal)) )

  } else if (methodBS == "block"){

    # Blocks overlap and are drawn with replacement
    FullBlocksSet <- dim(residPdynOriginal)[1] - BlockLength + 1 # all possible blocks that can be drawn
    SampleBlock <- ceiling((T-nlag)/BlockLength) #

    Rand <- matrix(stats::runif(SampleBlock, min = 0, max = 1), ncol = 1)
    bb <- ceiling(SampleBlock*Rand)
    IdxBlocks <- matrix(NA, nrow= BlockLength, ncol= FullBlocksSet)
    for (mm in 1:FullBlocksSet){
      IdxBlocks[ , mm] <- mm:(mm + BlockLength-1)
    }
    rr <- as.vector(IdxBlocks[,bb])[1:T]
    uPdyn <- residPdynOriginal[rr[1:(T-nlag)],]
    uYiel <- residYieOriginal[rr, ]
  }else{
    stop(paste('The method ', methodBS, ' is not available'))
  }

  return(list(residFact = uPdyn, residYields = uYiel))
}

##############################################################################################################
#' Compute the residuals from the observational equation
#'
#'@param ModelParaPE list of point estimates of the model parameter
#'@param BFull matrix B of loadings (CJ x F or J x F)
#'@param ModelType A character vector indicating the model type to be estimated
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal

residY_original <- function(ModelParaPE, BFull, ModelType, Economies){


  # For models estimated on a country-by-country basis
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))){

    residYie <- list()

    for (i in 1:length(Economies)){
    RiskFactorLabels <- rownames(ModelParaPE[[ModelType]][[Economies[i]]]$inputs$AllFactors)
    YY <- ModelParaPE[[ModelType]][[Economies[i]]]$inputs$Y
    ZZ <- ModelParaPE[[ModelType]][[Economies[i]]]$inputs$AllFactors
    A <- ModelParaPE[[ModelType]][[Economies[i]]]$rot$P$A
    B__Full <- BFull[[Economies[i]]]
    YYhat <- matrix(A, nrow= nrow(YY), ncol=ncol(YY)) + B__Full%*%ZZ # Model-implied yields

    residYie[[Economies[i]]] <- t(YY - YYhat)
    }
    } else{
    # For models estimated jointly
    YY <- ModelParaPE[[ModelType]]$inputs$Y
    ZZ <- ModelParaPE[[ModelType]]$inputs$AllFactors # K x T
    A <- ModelParaPE[[ModelType]]$rot$P$A
    YYhat <- matrix(A, nrow= nrow(YY), ncol=ncol(YY)) + BFull%*%ZZ # Model-implied yields
    residYie <- t(YY - YYhat)
  }


  return(residYie)
}
#########################################################################################################""
#'Compute the B matrix of loadings
#'
#'@param ModelParaPE list of point estimates of the model parameter
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param mat vector of bond yield maturities
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param ModelType A character vector indicating the model type to be estimated
#'
#'
#'@keywords internal

Get_BFull <- function(ModelParaPE, FactorLabels, mat, Economies, ModelType){

  J <- length(mat)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)

  # For models estimated separately
if(any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
  K <- nrow(ModelParaPE[[ModelType]][[Economies[1]]]$inputs$AllFactors)

  BFull <- list()
  for (i in 1:length(Economies)){
  B_CS <- matrix(0, nrow = J, ncol= K)
  LabelSpannedCS <- c(FactorLabels$Tables[[Economies[i]]][-(1:M)])
  RiskFactorLabels <- rownames(ModelParaPE[[ModelType]][[Economies[i]]]$inputs$AllFactors)
  idxSpanned <- match(LabelSpannedCS, RiskFactorLabels)
  B <- ModelParaPE[[ModelType]][[Economies[i]]]$rot$P$B
  B_CS[ , idxSpanned] <- B
  BFull[[Economies[i]]] <- B_CS
  }
  # For models estimated jointly
}else{
  K <- nrow(ModelParaPE[[ModelType]]$inputs$AllFactors)
  B <- ModelParaPE[[ModelType]]$rot$P$B
BFull <- BUnspannedAdapJoint(G,M,N,C, J, B)
}

return(BFull)
}
###################################################################################################################
#' Build the time-series of the risk factors in each bootstrap draw
#'
#'@param ModelParaPE list of point estimates of the model parameter
#'@param residPdynOriginal Time-series of the residuals from the P-dynamics equation (T x F)
#'@param residYieOriginal Time-series of the residuals from the observational equation (T x J or T x CJ)
#'@param InputsForOutputs list containing the desired inputs for the construction
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param ModelType Desired model to be estimated
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param GVARlist list of necessary inputs for the estimation of GVAR-based models
#'@param JLLlist list of necessary inputs for the estimation of JLL-based models
#'@param WishBRW Whether the user wishes to estimate the physical parameter model with the Bias correction model from BRW (2012) (see "Bias_Correc_VAR" function).\cr
#'              Default is set to 0.
#'@param BRWlist list of necessary inputs for performing the bias-corrected estimation (see "Bias_Correc_VAR" function)
#'@param nlag Number of lags in the P-dynamics. Default is set to 1.
#'
#'@keywords internal

BuildRiskFactors_BS <- function(ModelParaPE, residPdynOriginal, residYieOriginal, InputsForOutputs, Economies,
                                ModelType, FactorLabels, GVARlist, JLLlist, WishBRW, BRWlist, nlag = 1){

  # 1) Initialization
  if(any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
    RiskFact_Temp <- list()
    ZZ_list <- list()
    resid_list <- list()
    }


  for (i in 1:length(Economies)){
    if ((any(ModelType ==c("GVAR jointQ", "VAR jointQ","JLL original", "JLL NoDomUnit","JLL jointSigma")))
        & i >1 ){break}

    MaxEigen <- 1.1 # Initialization
    while (MaxEigen > 1){ #Test whether the VAR is stationary (if not, drop the draw)

      # Separately estimated models
    if(any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
      RiskFactors <- ModelParaPE[[ModelType]][[Economies[i]]]$inputs$AllFactors
      T <- ncol(RiskFactors)
      K <- nrow(RiskFactors)
      ZZ_row <- t(RiskFactors) # T x K
      D0Z <- ModelParaPE[[ModelType]][[Economies[i]]]$ests$K0Z
      D1Z <- ModelParaPE[[ModelType]][[Economies[i]]]$ests$K1Z
      resids_BS <- ResampleResiduals_BS(residPdynOriginal[[Economies[i]]], residYieOriginal[[Economies[i]]],
                                        InputsForOutputs, ModelType)
    # Jointly estimated models
    }else{
      RiskFactors <- ModelParaPE[[ModelType]]$inputs$AllFactors
      T <- ncol(RiskFactors)
      K <- nrow(RiskFactors)
      ZZ_row <- t(RiskFactors) # T x K
      D0Z <- ModelParaPE[[ModelType]]$ests$K0Z
      D1Z <- ModelParaPE[[ModelType]]$ests$K1Z
      resids_BS <- ResampleResiduals_BS(residPdynOriginal, residYieOriginal, InputsForOutputs, ModelType)
    }

    ZZ_Boot <- matrix(NA, nrow= T - 1 + nlag, ncol = K)
    dimnames(ZZ_Boot) <- list(rownames(ZZ_row), colnames(ZZ_row))

    # 2)  Compute artificial time-series
    # 2.1) initial values for the artificial data
    LAG <- c()
    for (jj in 1:nlag){
      ZZ_Boot[jj,] <- ZZ_row[jj,]
      LAG <- rbind(ZZ_Boot[jj,], LAG)
    }
    # Initialize the artificial series and the LAGplus vector
    LAGplus <- LAG
    LAGplus <- cbind(1, LAG)

    # 2.2) generate artificial series
    Ft <- rbind(t(D0Z), t(D1Z))
    # From observation nlag+1 to nobs, compute the artificial data
    for (jj in (nlag+1):(T-1+nlag)){
      for (mm in 1:K){
        # Compute the value for time=jj
        ZZ_Boot[jj,mm] = LAGplus %*% as.matrix(Ft[,mm]) + resids_BS$residFact[jj-nlag,mm]
      }
      # now update the LAG matrix
      if (jj<T-1+nlag){
        LAG <- rbind(ZZ_Boot[jj, ], LAG[1,seq_len((nlag-1)*K) ] )
        LAGplus <- cbind(1, LAG)
      }
    }

    # 3) Test whether the VAR is stationary (if not, drop the draw)
    if (any(ModelType ==c("JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
      K1Z_BS <- FeedbackMat_BS(ModelType, t(ZZ_Boot), FactorLabels, Economies, GVARlist, JLLlist, WishBRW, BRWlist)
    }else{
      RiskFact_Temp[[Economies[i]]] <- t(ZZ_Boot)
      if (any(ModelType ==c("JPS original", "JPS global"))){Economies_temp <- Economies[i]
      } else{ Economies_temp  <- Economies}
      K1Z_BS <- FeedbackMat_BS(ModelType, RiskFact_Temp, FactorLabels, Economies_temp, GVARlist, JLLlist,
                               WishBRW, BRWlist)
    }

    MaxEigen <- max(abs(eigen(K1Z_BS)$value))

    }

    # 4) Store outputs to export for sepQ models
    if(any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
    ZZ_list[[Economies[i]]] <- ZZ_Boot
    resid_list[[Economies[i]]] <- resids_BS
    }

  }

  if(any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
    Out <- list(ZZ_BS= ZZ_list, resids_BS = resid_list)
  }else{
  Out <- list(ZZ_BS= ZZ_Boot, resids_BS = resids_BS)}

  return(Out)
}

###################################################################################################################
#'Generate artificial time-series in the bootstrap setup
#'
#'@param ModelParaPE list of point estimates of the model parameter
#'@param residPdynOriginal Time-series of the residuals from the P-dynamics equation (T x F)
#'@param residYieOriginal Time-series of the residuals from the observational equation (T x J or T x CJ)
#'@param ModelType Desired model to be estimated
#'@param BFull matrix B of loadings (CJ x F or J x F)
#'@param InputsForOutputs list containing the desired inputs for the construction
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param GVARlist list of necessary inputs for the estimation of GVAR-based models
#'@param JLLlist list of necessary inputs for the estimation of JLL-based models
#'@param WishBRW Whether the user wishes to estimate the physical parameter model with the Bias correction model from BRW (2012) (see "Bias_Correc_VAR" function).\cr
#'              Default is set to 0.
#'@param BRWlist list of necessary inputs for performing the bias-corrected estimation (see "Bias_Correc_VAR" function)
#'@param nlag Number of lags in the P-dynamics. Default is set to 1.
#'
#'
#'@keywords internal

Gen_Artificial_Series <- function(ModelParaPE, residPdynOriginal, residYieOriginal, ModelType, BFull,
                                  InputsForOutputs, Economies, FactorLabels, GVARlist, JLLlist, WishBRW, BRWlist,
                                  nlag = 1){

  # 1) Artificial time-series of the risk factors
  BS_Set <- BuildRiskFactors_BS(ModelParaPE, residPdynOriginal, residYieOriginal, InputsForOutputs, Economies,
                               ModelType, FactorLabels, GVARlist, JLLlist, WishBRW, BRWlist, nlag)

  # Extract Unspanned factors from the full risk factor set
  ZZ_BS <- BS_Set$ZZ_BS

  if(any(ModelType == c("JPS original", "JPS global", "GVAR single"))){

    G <- length(FactorLabels$Global)
    N <- length(FactorLabels$Spanned)

    UnspannedFactors_CS_BS <- function(Economy, G, N) {
      AllFactors <- BS_Set$ZZ_BS[[Economy]]
      AllFactors[ , (G + 1):(ncol(AllFactors) - N)]
    }

    DomesticMacro_BS <- do.call(cbind, lapply(Economies, UnspannedFactors_CS_BS, G, N))
    GlobalMacro_BS <- do.call(cbind, lapply(Economies, function(country) ZZ_BS[[country]][ , seq_len(G), drop = FALSE]))


  }else{

  Idxs <- Idx_UnspanFact(t(ZZ_BS), FactorLabels, Economies) # indexes of the variable of interest
  GlobalMacro_BS <- ZZ_BS[ , Idxs$IdxGlobal, drop = FALSE]
  DomesticMacro_BS <- ZZ_BS[ , Idxs$IdxunSpa, drop = FALSE]

}
  # 2) Artificial time-series of bond yields
  Y_BS <- BuildYields_BS(ModelParaPE, ModelType, ZZ_BS, BFull, BS_Set, Economies)


  return(list(ZZ_BS = ZZ_BS, GlobalMacro_BS = GlobalMacro_BS, DomesticMacro_BS = DomesticMacro_BS,
              Y_BS= Y_BS))
}
#################################################################################################################
#'Clean unnecessary outputs of JLL models in the bootstrap setup
#'
#'@param ModelBootstrap List of outputs to store bootstrap draws
#'@param ndraws Total number of bootstrap draws
#'@param ModelType A character vector indicating the model type to be estimated
#'
#'@keywords internal


CleanOrthoJLL_Boot <- function(ModelBootstrap, ndraws, ModelType){
  for (tt in 1:ndraws){
    ModelBootstrap$NumOutDraws$IRF[[ModelType]][[tt]]$Yields$Ortho <- NULL
    ModelBootstrap$NumOutDraws$FEVD[[ModelType]][[tt]]$Yields$Ortho <- NULL
  }
  ModelBootstrap$ConfBounds$IRF[[ModelType]]$Yields$Ortho <- NULL
  ModelBootstrap$ConfBounds$FEVD[[ModelType]]$Yields$Ortho <- NULL
  return(ModelBootstrap)
}

##############################################################################################################
#'Compute the Feedback matrix of each bootstrap draw
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param RiskFactors_TS Time-series of risk factors of the bootstrap (F x T)
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param GVARlist list of necessary inputs for the estimation of GVAR-based models
#'@param JLLlist list of necessary inputs for the estimation of JLL-based models
#'@param WishBRW Whether the user wishes to estimate the physical parameter model with the Bias correction model from BRW (2012) (see "Bias_Correc_VAR" function).\cr
#'              Default is set to 0.
#'@param BRWlist list of necessary inputs for performing the bias-corrected estimation (see \code{\link{Bias_Correc_VAR}} function)
#'
#'@keywords internal

FeedbackMat_BS <- function(ModelType, RiskFactors_TS, FactorLabels, Economies, GVARlist, JLLlist,
                            WishBRW, BRWlist){

    # Model-specific inputs
    SpeInputs <- SpecificMLEInputs(ModelType, Economies, RiskFactors_TS, FactorLabels, GVARlist, JLLlist,
                                   WishBRW, BRWlist)


    if(WishBRW) SpeInputs$BRWinputs[c("checkBRW", "checkSigma")] <- 0

    # 1) Two special cases
    # a) JLL models without bias correction (avoid the unnecessary numerical optimization from the Sigma matrix)
    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) & !WishBRW){
      SpeInputs$JLLinputs$WishSigmas <- 0
      N <- length(FactorLabels$Spanned)
      PdynPara <- JLL(RiskFactors_TS, N, SpeInputs$JLLinputs)
      K1Z_BS <- PdynPara$k1

    # b) GVAR single model and  bias correction
      }else if(ModelType ==  "GVAR single" & WishBRW){
      N <- length(FactorLabels$Spanned)
      PdynPara <- Bias_Correc_VAR(ModelType, SpeInputs$BRWinputs, t(RiskFactors_TS[[1]]), N, Economies,
                                  FactorLabels, SpeInputs$GVARinputs)

      K1Z_BS <- PdynPara$Phi_tilde

    }else{
      # 2) All other specifications
      PdynPara <- GetPdynPara(RiskFactors_TS, FactorLabels, Economies, ModelType, SpeInputs$BRWinputs,
                            SpeInputs$GVARinputs, SpeInputs$JLLinputs)

    if(any(ModelType %in% c("JPS original", "JPS global", "GVAR single"))){
      K1Z_BS <- PdynPara[[Economies[1]]]$K1Z
    }else{K1Z_BS <- PdynPara$K1Z}
    }

  return(K1Z_BS)
}

################################################################################################################
#'Build the time-series of bond yields for each bootstrap draw
#'
#'@param ModelParaPE list of point estimates of the model parameter
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param RiskFactors_BS Time-series of the risk factors (F x T)
#'@param BFull B matrix of loadings
#'@param BS_Set Set of bootstrap inputs
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal

BuildYields_BS <- function(ModelParaPE, ModelType, RiskFactors_BS, BFull, BS_Set, Economies){


  # Models estimated jointly
  if(any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
    T <- nrow(RiskFactors_BS[[Economies[1]]])
    TS_Labels <- rownames(RiskFactors_BS[[Economies[1]]])
    Y_listBS <- list()

    for (i in 1:length(Economies)){
      YieldsLabels  <- rownames(ModelParaPE[[ModelType]][[Economies[i]]]$inputs$Y)

      Y_CS <- matrix(NA, nrow= T  , ncol = length(YieldsLabels))
      dimnames(Y_CS) <- list(TS_Labels, YieldsLabels)

      A <- ModelParaPE[[ModelType]][[Economies[i]]]$rot$P$A
      ZZ_BS <- RiskFactors_BS[[Economies[i]]]
      Y_CS <-  BS_Set$resids_BS[[Economies[i]]]$residYields + matrix(A, nrow= T, ncol=length(YieldsLabels), byrow=T) + ZZ_BS%*%t(BFull[[Economies[i]]])
      rownames(Y_CS) <-  TS_Labels

      Y_listBS[[Economies[i]]] <- Y_CS
    }

    Y_BS <- do.call(cbind, lapply(1:length(Economies), function(i) {Y_listBS[[Economies[i]]]}))

    # Models estimated separetely
  }else{

    T <- nrow(RiskFactors_BS)
    TS_Labels <- rownames(RiskFactors_BS)
    YieldsLabels  <- rownames(ModelParaPE[[ModelType]]$inputs$Y)

    Y_BS <- matrix(NA, nrow= T  , ncol = length(YieldsLabels))
    dimnames(Y_BS) <- list(TS_Labels, YieldsLabels)

    A <- ModelParaPE[[ModelType]]$rot$P$A
    Y_BS <-  BS_Set$resids_BS$residYields + matrix(A, nrow= T, ncol=length(YieldsLabels), byrow=T) + RiskFactors_BS%*%t(BFull)
    rownames(Y_BS) <-  TS_Labels
  }

  return(Y_BS)
}

##############################################################################################################
#'Gathers the estimate of the bootstrap draws
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelBootstrap List to store the bootstrap set
#'@param Draw_Opt List of model estimated parameters
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param tt number of the bootstrap draw
#'
#'@keywords internal

AdjustOptm_BS <- function(ModelType, ModelBootstrap, Draw_Opt, Economies, tt){
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))){

    for (i in 1:length(Economies)){
      ModelBootstrap$ParaDraws[[ModelType]][[Economies[i]]][[tt]] <- Draw_Opt[[ModelType]][[Economies[i]]]
    }

  }else{  ModelBootstrap$ParaDraws[[ModelType]][[tt]] <- Draw_Opt[[ModelType]]}

  return(ModelBootstrap)
}

#############################################################################################################
####################################################################################################
#' Prepare the factor set for GVAR models (Bootstrap version)
#'
#' @param ModelType A character vector containing the label of the model to be estimated.
#' @param RiskFactors A matrix of the complete set of risk factors (K x T).
#' @param Wgvar A transition matrix from GVAR models (C x C).
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#'
#' @return A list containing the factor set for GVAR models.
#' @keywords internal
DataSet_BS <- function(ModelType, RiskFactors, Wgvar, Economies, FactorLabels) {

  if (!any(ModelType %in% c("GVAR single", "GVAR multi"))) {
    return(NULL)
  }

  # 1) Pre-allocate list of factors
  T <- ncol(RiskFactors) # length of model's time dimension
  C <- length(Economies) # number of economies in the economic system
  N <- length(FactorLabels$Spanned) # number of country-specific spanned factors
  M <- length(FactorLabels$Domestic) - N # number of country-specific macro variables
  G <- length(FactorLabels$Global) # number of global variables

  ListFactors <- vector(mode = 'list', length = length(Economies) + 1) # length = all countries + global factors
  names(ListFactors) <- c(Economies, 'Global')

  # Initialize country-specific factors (CSF) and star factors (SF)
  CSF <- vector(mode = 'list', length = length(FactorLabels$Domestic))
  names(CSF) <- FactorLabels$Domestic
  SF <- vector(mode = 'list', length = length(FactorLabels$Star))
  names(SF) <- FactorLabels$Star

  for (i in 1:C) {
    ListFactors[[Economies[i]]] <- list(Factors = c(CSF, SF))
  }

  # Initialize global factors (GF)
  GF <- vector(mode = 'list', length = length(FactorLabels$Global))
  names(GF) <- FactorLabels$Global
  ListFactors[[length(Economies) + 1]] <- GF

  # 2) Fill in list with the corresponding factors
  # A) Country-specific variables (economy-related variables)
  for (i in 1:C) {
    for (j in 1:M) {
      ListFactors[[Economies[i]]]$Factors[[j]] <- as.matrix(RiskFactors[(c(FactorLabels$Tables[[Economies[i]]][j])), ])
    }
  }

  # B) Country-specific variables (pricing-related variables)
  idx0 <- M
  for (i in 1:C) {
    for (j in 1:N) {
      ListFactors[[Economies[i]]]$Factors[[idx0 + j]] <- as.matrix(RiskFactors[(c(FactorLabels$Tables[[Economies[i]]][idx0 + j])), ])
    }
  }

  # C) Foreign country-specific variables (economy and pricing-related)
  idx1 <- M + N
  Z <- lapply(1:(M + N), function(j) {
    X <- matrix(NA, nrow = C, ncol = T)
    for (i in 1:C) {
      X[i, ] <- ListFactors[[Economies[i]]]$Factors[[j]]
    }
    X
  })
  names(Z) <- FactorLabels$Domestic

  # C.1) If star variables are computed with time-varying weights
  if (is.list(Wgvar)) {
    # Use only the transition matrices that are included in the sample span
    t_First <- as.Date(colnames(RiskFactors)[1], format = "%d-%m-%Y")
    t_Last <- as.Date(colnames(RiskFactors)[T], format = "%d-%m-%Y")

    t_First_Wgvar <- format(t_First, "%Y")
    t_Last_Wgvar <- format(t_Last, "%Y")

    Wgvar_subset <- Wgvar[names(Wgvar) >= t_First_Wgvar & names(Wgvar) <= t_Last_Wgvar]

    # Add common column label (i.e. the year of the observation) to all variables
    Dates <- as.Date(colnames(RiskFactors), format = "%d-%m-%Y")
    YearLabels <- substr(Dates, 1, 4)
    Z <- lapply(Z, function(x) {
      colnames(x) <- paste0(YearLabels, seq_along(colnames(x)))
      x
    })

    # Compute the star variables with time-varying weights
    for (i in 1:C) {
      for (j in 1:(M + N)) {
        StarTimeVarTemp <- matrix(NA, nrow = T, ncol = 1)

        for (k in 1:length(Wgvar_subset)) {
          YearRef <- names(Wgvar_subset)[k] # year of reference
          IdxYear <- grep(YearRef, colnames(Z[[j]]))
          WgvarYear <- Wgvar_subset[[k]]
          StarTimeVarTemp[IdxYear] <- t(WgvarYear[i, ] %*% Z[[j]][, IdxYear])
        }
        # If the last year of the transition matrix happens earlier than the year of the last observation from the sample,
        # then use the last transition matrix for the remaining observations
        if (anyNA(StarTimeVarTemp)) {
          LenlastYear <- length(IdxYear)
          IdxLastObs <- IdxYear[LenlastYear] + 1
          StarTimeVarTemp[IdxLastObs:T] <- t(WgvarYear[i, ] %*% Z[[j]][, (IdxLastObs):T])
        }

        ListFactors[[Economies[i]]]$Factors[[idx1 + j]] <- StarTimeVarTemp
      }
    }

    # C.2) If star variables are computed with time fixed weights
  } else {
    if (any(ModelType == c("GVAR single", "GVAR multi"))) {
      for (i in 1:C) {
        for (j in 1:(M + N)) {
          ListFactors[[Economies[i]]]$Factors[[idx1 + j]] <- t(Wgvar[i, ] %*% Z[[j]])
        }
      }
    }
  }

  # D) Global Factors
  for (i in seq_len(G)) {
    ListFactors[[length(Economies) + 1]][[i]] <- as.matrix(RiskFactors[(c(FactorLabels$Global[i])), ])
  }

  return(ListFactors)
}

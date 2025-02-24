#' Numerical outputs (IRFs, GIRFs, FEVD, and GFEVD) for bootstrap
#'
#'@param ModelType A character vector indicating the model type to be estimated.
#'@param ModelParaBoot A list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param InputsForOutputs A list containing the necessary inputs for generating IRFs, GIRFs, FEVDs, GFEVDs and Term Premia.
#'@param FactorLabels  A list of character vectors with labels for all variables in the model.
#'@param Economies  A character vector containing the names of the economies included in the system.
#'
#'@keywords internal


NumOutputs_Bootstrap <- function(ModelType, ModelParaBoot, InputsForOutputs, FactorLabels, Economies){

  # IRF and GIRF
  IRFout <- IRFandGIRF_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)

  # FEVD and GFEVD
  FEVDout <- FEVDandGFEVD_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels, Economies)

  # Output summary
  NumericalOutputs <- list(IRF = IRFout$IRFs, FEVD = FEVDout$FEVDs, GIRF = IRFout$GIRFs, GFEVD = FEVDout$GFEVDs)

  return(NumericalOutputs)
}


######################################################################################################
########################################### 1) IRFs and GIRFs ########################################
######################################################################################################
#' IRFs and GIRFs after bootstrap for all models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param IRFhoriz single numerical vector containing the desired horizon of analysis for the IRFs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal


IRFandGIRF_BS <- function(ModelType, ModelParaBoot, IRFhoriz, FactorLabels, Economies){

  C <- length(Economies)

  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N


  # Pre-allocation
  GIRFoutputs <- list()
  IRFoutputs <- list()

  # 1) SINGLE COUNTRY MODELS
  if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){

    for (i in 1:C){
      ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]][[1]])
      K <- nrow(ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[1]]$ests$K1Z)
      J <- length(ModelParaBoot$GeneralInputs$mat)

      YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[1]]$inputs$Y) # Yield labels
      for (tt in 1:ndraws){

        # Summarize inputs for the IRFs
        SIGMA <- ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[tt]]$ests$SSZ # KxK (variance-covariance matrix)
        K1Z <- ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[tt]]$ests$K1Z # KxK (feedback matrix)
        B  <- BUnspannedAdapSep_BS(G, M, ModelParaBoot, Economies, Economies[i], ModelType,tt)

        # a) Compute IRFs
        IRFs <- ComputeIRFs(SIGMA, K1Z, B, FactorLabels, K, J, IRFhoriz, YieldsLabel, ModelType, Economies[i])
        IRFoutputs[[ModelType]][[Economies[i]]][[tt]] <- IRFs # Store Country specific IRFs

        # b) Compute GIRFs
        G0.y <- ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[tt]]$ests$Gy.0
        GIRFs <- ComputeGIRFs(SIGMA, K1Z, B, G0.y, FactorLabels, K, J, IRFhoriz, YieldsLabel, ModelType, Economies[i])
        GIRFoutputs[[ModelType]][[Economies[i]]][[tt]] <- GIRFs # Store Country specific GIRFs
        }

    }
  } else{

    # 2) JOINT COUNTRY MODELS
    ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])
    for (tt in 1:ndraws){
    J <- length(ModelParaBoot$GeneralInputs$mat)
    K <- nrow(ModelParaBoot$ParaDraws[[ModelType]][[1]]$ests$K1Z)
    YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelType]][[1]]$inputs$Y) # Yield labels

    # a) Summarize inputs for the IRFs
    if ( any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
      SIGMA <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$JLLoutcomes$Sigmas$Sigma_Y # For JLL models, we selected the cholesky factor, which won't be compute inside "the"ComputeIRFs"
    }else{  SIGMA <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$SSZ} # KxK (variance-covariance matrix)

    K1Z <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$K1Z # KxK (feedback matrix)
    BSpanned <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$rot$P$B
    B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

    # b) Compute IRFs
    IRFoutputs[[ModelType]][[tt]] <- ComputeIRFs(SIGMA, K1Z, B, FactorLabels, K, C*J, IRFhoriz,
                                                 YieldsLabel, ModelType)

    # c) Compute GIRFs
    G0.y <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$Gy.0
    GIRFs <- ComputeGIRFs(SIGMA, K1Z, B, G0.y, FactorLabels, K, C*J, IRFhoriz, YieldsLabel, ModelType)
    GIRFoutputs[[ModelType]][[tt]] <- GIRFs # Store Country specific GIRFs

    # 3) JLL-BASED MODELS (orthogonalized outputs)
    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){

      # Summarize inputs for the IRFs
      K1Ze <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$JLLoutcomes$k1_e # KxK (feedback matrix)
      PI <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$JLLoutcomes$PI
      Se <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$JLLoutcomes$Sigmas$Sigma_Ye

      # a) Compute IRFs orthogonalized
      IRFOrtho <- list()
      IRFOrtho[[ModelType]][[tt]] <- ComputeIRFs(Se, K1Ze, B, FactorLabels, K, C*J, IRFhoriz, YieldsLabel,
                                                ModelType, PI = PI, Mode= "Ortho")

      # Gather Outputs
      IRFoutputs[[ModelType]][[tt]]$Factors <- list(NonOrtho = IRFoutputs[[ModelType]][[tt]]$Factors,
                                              Ortho = IRFOrtho[[ModelType]][[tt]]$Factors)
      IRFoutputs[[ModelType]][[tt]]$Yields <- list(NonOrtho = IRFoutputs[[ModelType]][[tt]]$Yields,
                                             Ortho = IRFOrtho[[ModelType]][[tt]]$Yields)

      # b) Compute GIRFs orthogonalized
      GIRFsOrtho <- list()
      GIRFsOrtho[[ModelType]][[tt]] <- ComputeGIRFs(SIGMA, K1Z, B, G0.y, FactorLabels, K, C*J, IRFhoriz, YieldsLabel,
                                              ModelType, PI = PI, Mode = "Ortho")

      # Gather Outputs
      GIRFoutputs[[ModelType]][[tt]]$Factors <- list(NonOrtho = GIRFoutputs[[ModelType]][[tt]]$Factors,
                                               Ortho = GIRFsOrtho[[ModelType]][[tt]]$Factors)
      GIRFoutputs[[ModelType]][[tt]]$Yields <- list(NonOrtho = GIRFoutputs[[ModelType]][[tt]]$Yields,
                                              Ortho = GIRFsOrtho[[ModelType]][[tt]]$Yields)

    }
    }
}
  Out <- list(IRFs = IRFoutputs, GIRFs = GIRFoutputs)
  return(Out)
}

#########################################################################################################
################################### 2) FEVD and GFEVDs ##################################################
#########################################################################################################
#' FEVDs and GFEVDs after bootstrap for all models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param FEVDhoriz single numerical vector conataining the desired horizon of analysis for the FEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@keywords internal


FEVDandGFEVD_BS <- function(ModelType, ModelParaBoot, FEVDhoriz, FactorLabels, Economies){

  C <- length(Economies)

  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N

  # Pre-allocation
  GFEVDoutputs <- list()
  FEVDoutputs <- list()

  # 1) SINGLE COUNTRY MODELS
  if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){

    for (i in 1:C){
      ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]][[1]])
      K <- nrow(ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[1]]$ests$K1Z)
      J <- length(ModelParaBoot$GeneralInputs$mat)

      YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[1]]$inputs$Y) # Yield labels
      for (tt in 1:ndraws){

        # Summarize inputs for the IRFs
        SIGMA <- ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[tt]]$ests$SSZ # KxK (variance-covariance matrix)
        K1Z <- ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[tt]]$ests$K1Z # KxK (feedback matrix)
        G0 <-  ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[tt]]$ests$Gy.0
        B  <- BUnspannedAdapSep_BS(G, M, ModelParaBoot, Economies, Economies[i], ModelType,tt)

        # a) Compute FEVDs
        FEVDs <- ComputeFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, J, FEVDhoriz, YieldsLabel, ModelType, Economies[i])
        FEVDoutputs[[ModelType]][[Economies[i]]][[tt]] <- FEVDs # Store Country specific IRFs

        # b) Compute GFEVDs
        GFEVDs <- ComputeGFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, J, FEVDhoriz, YieldsLabel, ModelType, Economies[i])
        GFEVDoutputs[[ModelType]][[Economies[i]]][[tt]] <- GFEVDs # Store Country specific GIRFs
      }

    }
  }else{

    # 2) JOINT COUNTRY MODELS
    ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])
    for (tt in 1:ndraws){
      J <- length(ModelParaBoot$GeneralInputs$mat)
      K <- nrow(ModelParaBoot$ParaDraws[[ModelType]][[1]]$ests$K1Z)
      YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelType]][[1]]$inputs$Y) # Yield labels

      # a) Summarize inputs for the IRFs
      if ( any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
        Chol_Fac_JLL <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$JLLoutcomes$Sigmas$Sigma_Y
      }

      K1Z <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$K1Z
      SIGMA <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$SSZ
      BSpanned <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$rot$P$B
      G0 <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$Gy.0
      B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

      # b) Compute FEVDs
      FEVDoutputs[[ModelType]][[tt]] <- ComputeFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, C*J, FEVDhoriz,
                                                     YieldsLabel, ModelType, CholFac_JLL = Chol_Fac_JLL)

      # c) Compute GFEVDs
      GFEVDs <- ComputeGFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, C*J, FEVDhoriz, YieldsLabel, ModelType)
      GFEVDoutputs[[ModelType]][[tt]] <- GFEVDs

      # 3) JLL-BASED MODELS (orthogonalized outputs)
      if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){

        # Summarize inputs for the IRFs
        K1Ze <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$JLLoutcomes$k1_e # KxK (feedback matrix)
        PI <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$JLLoutcomes$PI
        SIGMA_Ortho <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$JLLoutcomes$Sigmas$VarCov_Ortho
        Chol_Fac_JLL_Ortho <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ests$JLLoutcomes$Sigmas$Sigma_Ye

        # a) Compute FEVDs orthogonalized
        FEVDOrtho <- list()
        FEVDOrtho[[ModelType]][[tt]] <- ComputeFEVDs(SIGMA_Ortho, K1Ze, G0, B, FactorLabels, K, C*J, FEVDhoriz,
                                                     YieldsLabel, ModelType, CholFac_JLL = Chol_Fac_JLL_Ortho, PI = PI,
                                                     Mode= "Ortho")

        # Gather Outputs
        FEVDoutputs[[ModelType]][[tt]]$Factors <- list(NonOrtho = FEVDoutputs[[ModelType]][[tt]]$Factors,
                                                 Ortho = FEVDOrtho[[ModelType]][[tt]]$Factors)
        FEVDoutputs[[ModelType]][[tt]]$Yields <- list(NonOrtho = FEVDoutputs[[ModelType]][[tt]]$Yields,
                                                Ortho = FEVDOrtho[[ModelType]][[tt]]$Yields)

        # b) Compute GFEVDs orthogonalized
        GFEVDsOrtho <- list()
        GFEVDsOrtho[[ModelType]][[tt]] <- ComputeGFEVDs(SIGMA_Ortho, K1Ze, G0, B, FactorLabels, K, C*J, FEVDhoriz,
                                                  YieldsLabel, ModelType, PI = PI, Mode = "Ortho")

        # Gather Outputs
        GFEVDoutputs[[ModelType]][[tt]]$Factors <- list(NonOrtho = GFEVDoutputs[[ModelType]][[tt]]$Factors,
                                                  Ortho = GFEVDsOrtho[[ModelType]][[tt]]$Factors)
        GFEVDoutputs[[ModelType]][[tt]]$Yields <- list(NonOrtho = GFEVDoutputs[[ModelType]][[tt]]$Yields,
                                                 Ortho = GFEVDsOrtho[[ModelType]][[tt]]$Yields)
      }
}
    }
  Out <- list(FEVDs = FEVDoutputs, GFEVDs = GFEVDoutputs)
  return(Out)
}

######################################################################################################
######################################## AUXILIARY FUNCTIONS #########################################
######################################################################################################
#' Obtain the full form of B unspanned for "sep Q" models within the bootstrap setting
#'
#'@param G number of global unspanned factors
#'@param M number of country-specific domestic unspanned factors
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param Economy  string-vector containing the names of the economy under study
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param tt number of the bootstrap draw
#'
#'
#'@keywords internal


BUnspannedAdapSep_BS <- function(G, M, ModelParaBoot, Economies, Economy, ModelType, tt){

  C <- length(Economies)
  J <- length(ModelParaBoot$GeneralInputs$mat)
  i <- match(Economy, Economies)
  N <- ModelParaBoot$GeneralInputs$N


  if( ModelType== "JPS original"){
    K <- nrow(ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[tt]]$ests$K1Z)
    BUnspanned <- matrix(0, nrow=J, ncol= K)
    BSpanned <- ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[tt]]$rot$P$B
    BUnspanned[ , (K-N+1):K] <-  BSpanned
  }


  else if( any(ModelType == c("JPS global", "GVAR single"))){
    K <- nrow(ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[tt]]$ests$K1Z)
    BUnspanned <- matrix(0, nrow=J, ncol= K)
    BSpanned <- ModelParaBoot$ParaDraws[[ModelType]][[Economies[i]]][[tt]]$rot$P$B

    IDX <- list()
    idx0 <- G+M
    for (h in 1:C){
      idx1 <- idx0 + N
      IDX[[h]] <- (idx0+1):idx1
      idx0 <- idx1 + M
    }
    BUnspanned[ , IDX[[i]]] <-  BSpanned
  }

  return(BUnspanned)
}

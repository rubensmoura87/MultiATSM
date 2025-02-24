#' Constructs the model numerical outputs (model fit, IRFs, GIRFs, FEVDs, GFEVDs, and risk premia decomposition)
#'
#'@param ModelType A character vector indicating the model type to be estimated.
#'@param ModelPara A list containing the point estimates of the model parameters. For details, refer to the outputs from the \code{\link{Optimization}} function.
#'@param InputsForOutputs A list containing the necessary inputs for generating IRFs, GIRFs, FEVDs, GFEVDs and Term Premia.
#'@param FactorLabels  A list of character vectors with labels for all variables in the model.
#'@param Economies A character vector containing the names of the economies included in the system.
#'
#'
#'@examples
#' # See an example of implementation in the vignette file of this package (Section 4).
#'
#'
#'@returns
#'List of the model numerical outputs, namely
#'\enumerate{
#'\item Model fit of bond yields
#'\item IRFs
#'\item FEVDs
#'\item GIRFs
#'\item GFEVDs
#'\item Bond yield decomposition
#'}
#'
#'@details
#'Both IRFs and FEVDs are computed using the Cholesky decomposition method. The risk factors are ordered as follows: (i) global unspanned factors, and (ii) domestic unspanned and spanned factors for each country. The order of countries follows the sequence defined in the \code{Economies} vector.
#'
#'@references
#' Pesaran, H. Hashem, and Shin, Yongcheol. "Generalized impulse response analysis in linear multivariate models." Economics letters 58.1 (1998): 17-29.
#'@export


NumOutputs <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies){

  cat("2.2) Computing numerical outputs \n")

  AllNumOutputs <- list()
  AllNumOutputs <- OutputConstruction(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies)

  # Save relevant numerical outputs
  PEoutputs<- list(ModelPara, AllNumOutputs)
  names(PEoutputs) <- c("Model Parameters", "Numerical Outputs")
  saveRDS(PEoutputs, paste(tempdir(),"/PEoutputs_", InputsForOutputs$'Label Outputs','.rds',sep=""))

  # Generate graphs, if previously selected
  GraphicalOutputs(ModelType, ModelPara, AllNumOutputs, InputsForOutputs, Economies, FactorLabels)
  return(AllNumOutputs)
}

######################################################################################################
######################################################################################################
#' Numerical outputs (variance explained, model fit, IRFs, GIRFs, FEVDs, GFEVDs, and risk premia decomposition)
#' for all models
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (See the "Optimization" function)
#'@param InputsForOutputs list conataining the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs,
#'                        and risk premia decomposition
#'@param FactorLabels string-list based which contains all the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal

OutputConstruction <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies){

  # Output summary
  # Total Variance Explained and Model fit
  Total_Var_exp <- VarianceExplained(ModelType, ModelPara, FactorLabels, Economies)
  ModFit <- YieldsFit(ModelType, ModelPara, FactorLabels, Economies)
  cat(" ** Model Fit \n")

  # IRF and GIRF
  IRFout <- IRFandGIRF(ModelType, ModelPara, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)
  cat(" ** IRFs and GIRFs  \n")

  # FEVD and GFEVD
  FEVDout <- FEVDandGFEVD(ModelType, ModelPara, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels, Economies)
  cat(" ** FEVDs and GFEVDs  \n")

  # Risk Premia Decomposition
  TermPremia <- TermPremiaDecomp(ModelPara, FactorLabels, ModelType, InputsForOutputs, Economies)
  cat(" ** Term Premia \n")

  NumericalOutputs <- list(Total_Var_exp, ModFit, IRFout$IRFs, FEVDout$FEVDs, IRFout$GIRFs, FEVDout$GFEVDs, TermPremia)
  names(NumericalOutputs) <- c("PC var explained", "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremiaDecomp")

  return(NumericalOutputs)
}

######################################################################################################
####################### 1) Total Variance Explained #########################################
######################################################################################################
#' Percentage explained by the spanned factors of the variations in the set of observed yields for all models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara List of model parameter estimates (see the "Optimization" function)
#'@param FactorLabels   string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal
#'

VarianceExplained <- function(ModelType, ModelPara, FactorLabels, Economies){

C <- length(Economies)
N <- length(FactorLabels$Spanned)

Total_Var_exp <- list()


# Models estimated individually
if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
  for (i in 1:C){
    H <- eigen(stats::cov(t(ModelPara[[ModelType]][[Economies[i]]]$inputs$Y)))$values
    percentages_explained <- cumsum(H)/sum(H)
    Total_Var_exp[[i]] <-percentages_explained[1:N]
  }

}else{
  # Models estimated jointly
  J <- length(ModelPara[[ModelType]]$inputs$mat)
  idx0 <- 0
  for (i in 1:C){
    idx1 <- idx0 + J
    H <- eigen(stats::cov(t(ModelPara[[ModelType]]$inputs$Y[(idx0+1):idx1,])))$values
    percentages_explained <- cumsum(H)/sum(H)
    Total_Var_exp[[i]] <-percentages_explained[1:N]
    idx0 <- idx1
  }
}
  names(Total_Var_exp) <- Economies

  return(Total_Var_exp)
}


######################################################################################################
########################################### 2) Model Fit #############################################
######################################################################################################
#' Computes two measures of model fit for bond yields (all models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param ModelPara List of model parameter estimates (See the "Optimization" function)
#'@param FactorLabels a string-list based which contains the labels of all the variables present in the model
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'@details
#' "Model-implied yields" is the measure of fit based exclusively on the risk-neutral parameters, whereas the
#' "Model-Fit" takes into account both the risk-neutral and the physical paameters.
#'
#' @references
#' See, for instance, Jotiskhatira, Le and Lundblad (2015). "Why do interest rates in different currencies co-move?" (Journal of Financial Economics)
#'
#'@keywords internal


YieldsFit <- function(ModelType, ModelPara, FactorLabels, Economies){

  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)

  Output <- list()

  # I) Models estimated individually
  if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){

    mat <- ModelPara[[ModelType]][[Economies[1]]]$inputs$mat
    J <- length(mat)

    for (i in 1:C){

      YieldData <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y
      Z <- ModelPara[[ModelType]][[Economies[i]]]$inputs$AllFactors

      T <- ncol(Z)

      Afull <- ModelPara[[ModelType]][[Economies[i]]]$rot$P$A
      Bspanned <- ModelPara[[ModelType]][[Economies[i]]]$rot$P$B
      K0Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K0Z
      K1Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K1Z

      # 2) MODEL FIT MEASURES
      # a) Model Fit (Yields)
      # Extract spanned factors from the list of unspanned factors
      if (ModelType == "JPS original"){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]]) }
      else { AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

      LabelSpannedCS <- c(FactorLabels$Tables[[Economies[i]]][-(1:M)])
      IdxSpanned <- match(LabelSpannedCS, AllLabels)

      P <- Z[IdxSpanned,] # Set of spanned factors

      # Compute model fit
      Yieldfit <- Y_Fit(Afull, Bspanned, P, J, T, dimnames(YieldData))

      # b) Model-implied Yields
      Bfull <- BUnspannedAdapSep(G, M, ModelPara, Economies, Economies[i], ModelType)
      dimnames(Bfull) <- list(rownames(YieldData), rownames(Z))

      YieldModelImplied <- Y_ModImp(Afull, Bfull, K0Z, K1Z, Z, J, T, dimnames(YieldData))

      # 3) Prepare outputs
      fits<- list(Yieldfit, YieldModelImplied)
      names(fits) <- c("Yield Fit","Yield Model Implied")

      Output[[ModelType]][[Economies[i]]] <-  fits
    }

    # II) Models estimated jointly
    }else{

    Z <- ModelPara[[ModelType]]$inputs$AllFactors
    YieldData <- ModelPara[[ModelType]]$inputs$Y

    mat <- ModelPara[[ModelType]]$inputs$mat
    Afull <- ModelPara[[ModelType]]$rot$P$A
    Bspanned <- ModelPara[[ModelType]]$rot$P$B
    K0Z <- ModelPara[[ModelType]]$ests$K0Z
    K1Z <- ModelPara[[ModelType]]$ests$K1Z

    J <- length(mat)
    T <- ncol(Z)

    # 2) MODEL FIT MEASURES
    # a) Model Fit (Yields)
    # Extract spanned factors from the list of unspanned factors
    IdxSpanned <- c()
    idxSpa0 <- G + M
    for (j in 1:C){
      idxSpa1 <- idxSpa0 + N

      if (j ==1){ IdxSpanned <- (idxSpa0+1):idxSpa1
      }  else{      IdxSpanned <- c(IdxSpanned, (idxSpa0+1):idxSpa1) }

      idxSpa0 <- idxSpa1 + M
    }

    P <- Z[IdxSpanned, ] # Set of spanned factors

    # Compute model fit
    Yieldfit <- Y_Fit(Afull, Bspanned, P, C*J, T, dimnames(YieldData))

    # b) Model-implied Yields
    Bfull <- BUnspannedAdapJoint(G,M,N,C, J, Bspanned)
    dimnames(Bfull) <- list(rownames(YieldData), rownames(Z))

    YieldModelImplied <- Y_ModImp(Afull, Bfull, K0Z, K1Z, Z, C*J, T, dimnames(YieldData))

    Output<- list(Yieldfit, YieldModelImplied)
    names(Output) <- c("Yield Fit","Yield Model Implied")

  }

  return(Output)
}

######################################################################################################
########################################### 3) IRFs and GIRFs ########################################
######################################################################################################
#' IRFs and GIRFs for all models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (See the "Optimization" function)
#'@param IRFhoriz single numerical vector containing the desired horizon of analysis for the IRFs
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@details
#' The Structural shocks from the IRFs are identified via Cholesky decomposition
#'
#'
#'@keywords internal


IRFandGIRF <- function(ModelType, ModelPara, IRFhoriz, FactorLabels, Economies){

  # Pre-allocation
  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N

  IRFoutputs <- list()
  GIRFoutputs <- list()

  # 1) SINGLE COUNTRY MODELS
  if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){

    J <- length(ModelPara[[ModelType]][[Economies[1]]]$inputs$mat)
    K <- nrow(ModelPara[[ModelType]][[Economies[1]]]$inputs$AllFactors)

    for (i in 1:C){

      YieldsLabel<- rownames(ModelPara[[ModelType]][[Economies[i]]]$inputs$Y) # Yield labels
      # a) Summarize inputs for the IRFs
      SIGMA <- ModelPara[[ModelType]][[Economies[i]]]$ests$SSZ # KxK (variance-covariance matrix)
      K1Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K1Z # KxK (feedback matrix)
      B <- BUnspannedAdapSep(G,M, ModelPara, Economies, Economies[i], ModelType)

      # b) Compute IRFs
      IRFs <- ComputeIRFs(SIGMA, K1Z, B, FactorLabels, K, J, IRFhoriz, YieldsLabel, ModelType, Economies[i])
      IRFoutputs[[ModelType]][[Economies[i]]] <- IRFs # Store Country specific IRFs

      # c) Compute GIRFs
      G0.y <- ModelPara[[ModelType]][[Economies[i]]]$ests$Gy.0
      GIRFs <- ComputeGIRFs(SIGMA, K1Z, B, G0.y, FactorLabels, K, J, IRFhoriz, YieldsLabel, ModelType, Economies[i])
      GIRFoutputs[[ModelType]][[Economies[i]]] <- GIRFs # Store Country specific GIRFs
      }

  } else{

  # 2) JOINT COUNTRY MODELS
    J <- length(ModelPara[[ModelType]]$inputs$mat)
    K <- nrow(ModelPara[[ModelType]]$inputs$AllFactors)
    YieldsLabel<- rownames(ModelPara[[ModelType]]$inputs$Y) # Yield labels

    # a) Summarize inputs for the IRFs
    if ( any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
      SIGMA <- ModelPara[[ModelType]]$ests$JLLoutcomes$Sigmas$Sigma_Y # For JLL models, we selected the cholesky factor, which won't be compute inside "the"ComputeIRFs"
    }else{  SIGMA <- ModelPara[[ModelType]]$ests$SSZ} # KxK (variance-covariance matrix)

    K1Z <- ModelPara[[ModelType]]$ests$K1Z # KxK (feedback matrix)
    BSpanned <- ModelPara[[ModelType]]$rot$P$B
    B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

    # b) Compute IRFs
    IRFoutputs[[ModelType]] <- ComputeIRFs(SIGMA, K1Z, B, FactorLabels, K, C*J, IRFhoriz, YieldsLabel, ModelType)

    # c) Compute GIRFs
    G0.y <- ModelPara[[ModelType]]$ests$Gy.0
    GIRFs <- ComputeGIRFs(SIGMA, K1Z, B, G0.y, FactorLabels, K, C*J, IRFhoriz, YieldsLabel, ModelType)
    GIRFoutputs[[ModelType]] <- GIRFs # Store Country specific GIRFs


    # 3) JLL-BASED MODELS (orthogonalized outputs)
    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){

      # Summarize inputs for the IRFs
      K1Ze <- ModelPara[[ModelType]]$ests$JLLoutcomes$k1_e # KxK (feedback matrix)
      PI <- ModelPara[[ModelType]]$ests$JLLoutcomes$PI
      Se <- ModelPara[[ModelType]]$ests$JLLoutcomes$Sigmas$Sigma_Ye

      # a) Compute IRFs orthogonalized
      IRFOrtho <- list()
      IRFOrtho[[ModelType]] <- ComputeIRFs(Se, K1Ze, B, FactorLabels, K, C*J, IRFhoriz, YieldsLabel,
                                           ModelType, PI = PI, Mode= "Ortho")

      # Gather Outputs
      IRFoutputs[[ModelType]]$Factors <- list(NonOrtho = IRFoutputs[[ModelType]]$Factors,
                                              Ortho = IRFOrtho[[ModelType]]$Factors)
      IRFoutputs[[ModelType]]$Yields <- list(NonOrtho = IRFoutputs[[ModelType]]$Yields,
                                             Ortho = IRFOrtho[[ModelType]]$Yields)

      # b) Compute GIRFs orthogonalized
      GIRFsOrtho <- list()
      GIRFsOrtho[[ModelType]] <- ComputeGIRFs(SIGMA, K1Z, B, G0.y, FactorLabels, K, C*J, IRFhoriz, YieldsLabel,
                                              ModelType, PI = PI, Mode = "Ortho")

      # Gather Outputs
      GIRFoutputs[[ModelType]]$Factors <- list(NonOrtho = GIRFoutputs[[ModelType]]$Factors,
                                               Ortho = GIRFsOrtho[[ModelType]]$Factors)
      GIRFoutputs[[ModelType]]$Yields <- list(NonOrtho = GIRFoutputs[[ModelType]]$Yields,
                                              Ortho = GIRFsOrtho[[ModelType]]$Yields)

      }
    }

  Out <- list(IRFs = IRFoutputs, GIRFs = GIRFoutputs)
  return(Out)
}


#########################################################################################################
################################### 4) FEVD and GFEVD ###################################################
#########################################################################################################
#' FEVDs and GFEVDs for all models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param FEVDhoriz single numerical vector containing the desired horizon of analysis for the FEVDs and GFEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@details
#' Structural shocks are identified via Cholesky decomposition
#'
#'@keywords internal

FEVDandGFEVD <- function(ModelType, ModelPara, FEVDhoriz, FactorLabels, Economies){

  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N

  FEVDoutputs <- list()
  GFEVDoutputs <- list()

  # 1) SINGLE COUNTRY MODELS
  if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){

    K <- nrow(ModelPara[[ModelType]][[Economies[1]]]$inputs$AllFactors)
    J <- length(ModelPara[[ModelType]][[Economies[1]]]$inputs$mat)

    for (i in 1:C){

    YieldsLabel <- rownames(ModelPara[[ModelType]][[Economies[i]]]$inputs$Y) # Yield labels

    # a) Summarize inputs for the FEVDs
    SIGMA <- ModelPara[[ModelType]][[Economies[i]]]$ests$SSZ
    K1Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K1Z
    G0 <-  ModelPara[[ModelType]][[Economies[i]]]$ests$Gy.0
    B <- BUnspannedAdapSep(G, M, ModelPara, Economies, Economies[i], ModelType)

    # b) Compute FEVDs
    FEVDs <- ComputeFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, J, FEVDhoriz, YieldsLabel, ModelType, Economies[i])
    FEVDoutputs[[ModelType]][[Economies[i]]] <- FEVDs # Store Country specific FEVDs

    # c) Compute GFEVDs
    GFEVDs <- ComputeGFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, J, FEVDhoriz, YieldsLabel, ModelType, Economies[i])
    GFEVDoutputs[[ModelType]][[Economies[i]]] <- GFEVDs # Store Country specific GFEVDs
  }
  } else{

    # 2) JOINT COUNTRY MODELS
    J <- length(ModelPara[[ModelType]]$inputs$mat)
    K <- nrow(ModelPara[[ModelType]]$inputs$AllFactors)
    YieldsLabel <- rownames(ModelPara[[ModelType]]$inputs$Y) # Yield labels

    # a) Summarize inputs for the FEVD
    if ( any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
      Chol_Fac_JLL <- ModelPara[[ModelType]]$ests$JLLoutcomes$Sigmas$Sigma_Y
      }

    K1Z <- ModelPara[[ModelType]]$ests$K1Z # KxK (feedback matrix)
    SIGMA <- ModelPara[[ModelType]]$ests$SSZ # KxK (variance-covariance matrix)
    BSpanned <- ModelPara[[ModelType]]$rot$P$B
    B <- BUnspannedAdapJoint(G, M, N, C, J, BSpanned)
    G0 <- ModelPara[[ModelType]]$ests$Gy.0

    # b) Compute FEVD
    FEVDoutputs[[ModelType]] <- ComputeFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, C*J, FEVDhoriz, YieldsLabel,
                                             ModelType, CholFac_JLL = Chol_Fac_JLL)

    # c) Compute GFEVD
    GFEVDoutputs[[ModelType]] <- ComputeGFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, C*J, FEVDhoriz, YieldsLabel,
                                               ModelType)

    # 3) JLL-BASED MODELS (orthogonalized outputs)
    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){

      # Summarize inputs for the IRFs
      K1Ze <- ModelPara[[ModelType]]$ests$JLLoutcomes$k1_e # KxK (feedback matrix)
      PI <- ModelPara[[ModelType]]$ests$JLLoutcomes$PI
      SIGMA_Ortho <- ModelPara[[ModelType]]$ests$JLLoutcomes$Sigmas$VarCov_Ortho
      Chol_Fac_JLL_Ortho <- ModelPara[[ModelType]]$ests$JLLoutcomes$Sigmas$Sigma_Ye

      # a) Compute FEVDs orthogonalized
      FEVDOrtho <- list()
      FEVDOrtho[[ModelType]] <- ComputeFEVDs(SIGMA_Ortho, K1Ze, G0, B, FactorLabels, K, C*J, FEVDhoriz, YieldsLabel,
                                             ModelType, CholFac_JLL = Chol_Fac_JLL_Ortho, PI = PI, Mode= "Ortho")

      # Gather Outputs
      FEVDoutputs[[ModelType]]$Factors <- list(NonOrtho = FEVDoutputs[[ModelType]]$Factors,
                                              Ortho = FEVDOrtho[[ModelType]]$Factors)
      FEVDoutputs[[ModelType]]$Yields <- list(NonOrtho = FEVDoutputs[[ModelType]]$Yields,
                                             Ortho = FEVDOrtho[[ModelType]]$Yields)

      # b) Compute GFEVDs orthogonalized
      GFEVDsOrtho <- list()
      GFEVDsOrtho[[ModelType]] <- ComputeGFEVDs(SIGMA_Ortho, K1Ze, G0, B, FactorLabels, K, C*J, FEVDhoriz,
                                                YieldsLabel, ModelType, PI = PI, Mode = "Ortho")

      # Gather Outputs
      GFEVDoutputs[[ModelType]]$Factors <- list(NonOrtho = GFEVDoutputs[[ModelType]]$Factors,
                                                Ortho = GFEVDsOrtho[[ModelType]]$Factors)
      GFEVDoutputs[[ModelType]]$Yields <- list(NonOrtho = GFEVDoutputs[[ModelType]]$Yields,
                                               Ortho = GFEVDsOrtho[[ModelType]]$Yields)

  }
}

  Out <- list(FEVDs = FEVDoutputs, GFEVDs = GFEVDoutputs)
  return(Out)
}

######################################################################################################
####################### 5) Risk Premia Decomposition #################################################
######################################################################################################
#' Decomposition of yields into the average of expected future short-term interest rate and risk premia for all models
#'
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param InputsForOutputs list containing the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs,
#'                        and risk premia decomposition
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'
#'@keywords internal

TermPremiaDecomp <- function(ModelPara, FactorLabels, ModelType, InputsForOutputs, Economies){

  WishFP <- InputsForOutputs$ForwardPremia

  # 1) Compute the country-specific expected components
  EP_list <- ExpectedComponent(ModelPara, InputsForOutputs, ModelType, Economies, FactorLabels, WishFP)

  # 2) Compute Term Premium
  OutputTP <- TermPremia(ModelPara, EP_list$avexp, ModelType, Economies)

  # 3) Forward Premia (if desired)
  if( WishFP){
    OutputFP <- ForwardPremia(ModelPara, EP_list$avexpFP, ModelType, FactorLabels, InputsForOutputs, Economies)
    }else {  OutputFP <- NA  }

  Output <- list(RiskPremia = OutputTP, ForwardPremia = OutputFP)


  return(Output)
}

######################################################################################################
######################################## AUXILIARY FUNCTIONS #########################################
######################################################################################################
#' Model-implied yields (cross-section)
#'
#'@param ALoad  A loadings
#'@param BLoad B loadings
#'@param Spa_TS time series of spanned factors
#'@param MatLength length of the vector of maturities
#'@param TDim Time-series dimension
#'@param YieldLab Label of yields
#'
#'@keywords internal

Y_Fit <- function(ALoad, BLoad, Spa_TS , MatLength, TDim, YieldLab){

  Yieldfit<- matrix(NA, nrow = MatLength, ncol = TDim)
  if(is.null(dim(Spa_TS))){
    for (h in 1:TDim){   Yieldfit[,h] <- ALoad + BLoad%*%Spa_TS[h] }
  }else{
  for (h in 1:TDim){   Yieldfit[,h] <- ALoad + BLoad%*%Spa_TS[,h] }
  }
    dimnames(Yieldfit) <- YieldLab

  return(Yieldfit)
}

##############################################################
#' Model-implied yields (P-dynamics)
#'
#'@param ALoad  A loadings
#'@param BLoad B loadings
#'@param K0Z intercept from the P-dynamics
#'@param K1Z feedback matrix from the P-dynamics
#'@param PdynFact time series of the risk-factors spanned factors
#'@param MatLength length of the vector of maturities
#'@param TDim Time-series dimension
#'@param YieldLab Label of yields
#'
#'
#'@keywords internal

Y_ModImp <- function(ALoad, BLoad, K0Z, K1Z, PdynFact, MatLength, TDim, YieldLab){

  YieldModelImplied <- matrix(NA, nrow= MatLength, ncol = TDim)
  for (h in 2:TDim){ #  first observation is discarded
    YieldModelImplied[,h] <- ALoad + BLoad%*%(K0Z + K1Z%*% PdynFact[,h-1])
  }
  dimnames(YieldModelImplied) <- YieldLab


  return(YieldModelImplied)
}


#####################################################################################################
#' Transform B_spanned into B_unspanned for jointQ models
#'
#'
#'@param G number of global unspanned factors
#'@param M number of domestic unspanned factors
#'@param N number of domestic spanned factors
#'@param C number of economies of the economic system
#'@param J number of country-specific observed bond yields
#'@param BSpanned B that accomodates only the map to the spanned factors only
#'
#'
#'@keywords internal

BUnspannedAdapJoint <- function(G,M,N,C, J, BSpanned){

  K <- C*(N+M) +G
  CJ <- C*J

  BUnspanned <- matrix(0, nrow=CJ, ncol= K)

  idxA <- 0
  idxB <- G + M
  idxC <- 0

  for (i in 1:C){
    idxAA <- idxA + J
    idxBB <- idxB + N
    idxCC <- idxC + N
    BUnspanned[(idxA+1):idxAA, (idxB+1):idxBB] <- BSpanned[(idxA+1):idxAA, (idxC+1):idxCC]
    idxA <- idxAA
    idxB <- idxBB +M
    idxC <- idxCC
  }


  return(BUnspanned)
}

#################################################################################################
#' Transform B_spanned into B_unspanned for sepQ models
#'
#' @param G number of global unspanned factors
#' @param M number of domestic unspanned factors per country
#' @param ModelPara list of model parameter estimates (See the "Optimization" function)
#' @param Economies complet set of economies of the economic system
#' @param Economy  specific economy under study
#' @param ModelType a string-vector containing the label of the model to be estimated
#'
#'@keywords internal

BUnspannedAdapSep <- function(G,M, ModelPara, Economies, Economy, ModelType){

  i <- match(Economy, Economies)
  C <- length(Economies)
  N <- ModelPara[[ModelType]][[Economy]]$inputs$N
  J <- length(ModelPara[[ModelType]][[Economy]]$inputs$mat)


  if( ModelType == "JPS original"){
    K <- N+M + G
    BUnspanned <- matrix(0, nrow=J, ncol= K)
    BSpanned <- ModelPara[[ModelType]][[Economies[i]]]$rot$P$B
    BUnspanned[ , (K-N+1):K] <-  BSpanned
  }



  else if( any(ModelType == c("JPS global","GVAR single" ))){
    K <- C*(N+M) + G
    BUnspanned <- matrix(0, nrow=J, ncol= K)
    BSpanned <- ModelPara[[ModelType]][[Economies[i]]]$rot$P$B

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
#####################################################################################################
#'Compute IRFs of all models
#'
#'@param SIGMA Variance-covariance matrix
#'@param K1Z Loading As
#'@param BLoad Loading Bs
#'@param FactorLabels List containing the label of factors
#'@param FacDim Dimension of the P-dynamics
#'@param MatLength Length of the maturity vector
#'@param IRFhoriz Horizon of the analysis
#'@param YieldsLabel Label of bond yields
#'@param ModelType Desired model type
#'@param Economy specific economy under study
#'@param PI matrix PI for JLL-based models
#'@param Mode allows for the orthogonalized version in the case of JLL-based models
#'
#'@keywords internal

ComputeIRFs <- function(SIGMA, K1Z, BLoad, FactorLabels, FacDim, MatLength, IRFhoriz, YieldsLabel, ModelType,
                        Economy = NULL, PI = NULL, Mode = FALSE){


  # 1) Initialization of IRFs of interest
  tempFactors <- array(0, c(FacDim, FacDim, IRFhoriz))
  tempYields  <- array(0, c(MatLength, FacDim,IRFhoriz))


  # 2) Compute the IRFs
  if (Mode == "Ortho" & any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {AdjTerm <- PI
  } else{ AdjTerm <- diag(FacDim) }

  # Choleski term
  if ( any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
     S <- SIGMA  } else{ S <- t(chol(SIGMA))}

  # Shock at t=0:
  tempFactors[ ,  , 1] <- S
  tempYields[ , , 1]  <- BLoad %*%AdjTerm%*%S
  # Shock at t=1:
  for (r in 2:IRFhoriz){
    if (r == 2) { A1h <- K1Z} else { A1h <- A1h %*% K1Z}
    tempFactors[ , , r] <- A1h%*%S # IRF (t+h) = A1^h*S
    tempYields[ , , r]  <- BLoad%*%AdjTerm%*%A1h%*%S
  }
  IRFRiskFactors <- aperm(tempFactors, c(3,1,2))
  IRFYields <- aperm(tempYields, c(3,1,2))

  Horiz <- t(t(0:(IRFhoriz-1))) #Add a column for horizon of interest


  # 3) Adjust the variable labels
  # Factor Labels
  if ( ModelType == "JPS original") {AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]])}
  else if(any(ModelType == c("JPS global", "GVAR single"))) {
    AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  }else{ AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

  dimnames(IRFRiskFactors) <- list(Horiz, AllFactorsLabels, AllFactorsLabels)
  dimnames(IRFYields) <- list(Horiz, YieldsLabel, AllFactorsLabels)

  Out <- list(Factors = IRFRiskFactors, Yields = IRFYields)

  return(Out)
}
##########################################################################################################
##########################################################################################################
#'Compute GIRFs for all models
#'
#'@param Sigma.y Variance-covariance matrix
#'@param F1 Feedback matrix
#'@param BLoad Loading Bs
#'@param G0.y matrix of contemporaneous terms
#'@param FactorLabels List containing the labels of the factors
#'@param FacDim Dimension of the P-dynamics
#'@param MatLength Length of the maturity vector
#'@param GIRFhoriz Horizon of the analysis
#'@param YieldsLabel Label o yields
#'@param ModelType desired Model type
#'@param Economy Economy under study
#'@param PI matrix PI for JLL-based models
#'@param Mode allows for the orthogonalized version in the case of JLL-based models
#'
#'
#'#'#' @references
#' \itemize{
#' \item This function is a partially based on the version of the "irf" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#'@keywords internal

ComputeGIRFs <- function(Sigma.y, F1, BLoad, G0.y, FactorLabels,  FacDim, MatLength, GIRFhoriz, YieldsLabel,
                         ModelType, Economy = NULL, PI = NULL, Mode = FALSE){

  # 1) Dynamic multiplier:
  Ry.h <- array(NA, c(FacDim,FacDim, GIRFhoriz))
  Ry.h[, ,1] <- diag(FacDim) # dynamic multiplier at t=0

  for (w in 2:GIRFhoriz) {  Ry.h[, ,w] <- F1%*%Ry.h[, ,w-1] }

  # 2) Build the vector containing the one unit-shock for each variable of the system
  ey.j <- diag(FacDim)

  # 3) GIRFs:
  if (Mode == "Ortho" & any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {AdjTerm <- PI
  } else{ AdjTerm <- diag(FacDim) }

  # 3.1) Factors
  AllResponsesToAllShocksFactors <- array(NA, c(FacDim,GIRFhoriz,FacDim))
  AllResponsesToShockOfOneVariableFactors <- matrix(NA, ncol= GIRFhoriz , nrow = FacDim)
  for (g in 1:FacDim){
    for (w in 1:GIRFhoriz){
      numFactors <- AdjTerm%*%(Ry.h[,,w]%*% solve(G0.y)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
      demFactors <- 1/sqrt((t(ey.j[,g])%*%Sigma.y%*%ey.j[,g])) # denominator from equation at the bottom of the page 22 (PS, 1998)
      AllResponsesToShockOfOneVariableFactors[,w] <- numFactors*drop(demFactors)
    }
    AllResponsesToAllShocksFactors[,,g] <- AllResponsesToShockOfOneVariableFactors
  }

  GIRFFactors <- aperm(AllResponsesToAllShocksFactors, c(2,1,3))

  # 3.2) Yields
  AllResponsesToAllShocksYields <- array(NA, c(MatLength, GIRFhoriz, FacDim))
  AllResponsesToShockOfOneVariableYields <- matrix(NA, ncol= GIRFhoriz , nrow = MatLength)
  for (g in 1:FacDim){
    for (w in 1:GIRFhoriz){
      numYields <- BLoad%*%AdjTerm%*%(Ry.h[,,w]%*% solve(G0.y)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
      demYields <- 1/sqrt((t(ey.j[,g])%*%Sigma.y%*%ey.j[,g])) # denominator from equation at the bottom of the page 22 (PS, 1998)
      AllResponsesToShockOfOneVariableYields[,w] <- numYields*drop(demYields)
    }
    AllResponsesToAllShocksYields[,,g] <- AllResponsesToShockOfOneVariableYields
  }

  GIRFYields <- aperm(AllResponsesToAllShocksYields, c(2,1,3))

  # 4) Prepare labels for the output
  if(ModelType == "JPS original" ){  labelsGIRF <- c(FactorLabels$Global,FactorLabels$Tables[[Economy]]) }
  else if(any(ModelType == c("JPS global", "GVAR single"))) {
    labelsGIRF <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  }else{ labelsGIRF <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

  # 4.1) Add columns containig the horizons
  Horiz <- t(t(0:(GIRFhoriz-1))) #Add a column for horizon of interest

  # 4.2) Labels
  dimnames(GIRFFactors) <- list(Horiz, labelsGIRF, labelsGIRF)
  dimnames(GIRFYields) <- list(Horiz, YieldsLabel, labelsGIRF)


  GIRFoutputs <- list(Factors = GIRFFactors, Yields = GIRFYields)


  return(GIRFoutputs)
}

#######################################################################################################
#'Compute FEVDs for all models
#'
#'@param SIGMA Variance-covariance matrix
#'@param K1Z Loading As
#'@param G0 contemporaneous terms
#'@param BLoad Loading Bs
#'@param FactorLabels List containing the label of factors
#'@param FacDim Dimension of the P-dynamics
#'@param MatLength Length of the maturity vector
#'@param FEVDhoriz Horizon of the analysis
#'@param YieldsLabel Label of bond yields
#'@param ModelType Desired model type
#'@param Economy specific economy under study
#'@param CholFac_JLL Cholesky factorization term from JLL models
#'@param PI matrix PI for JLL-based models
#'@param Mode allows for the orthogonalized version in the case of JLL-based models
#'
#'@references
#' \itemize{
#' \item This function is a modified and extended version of the "fevd" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#'@keywords internal

ComputeFEVDs <- function(SIGMA, K1Z, G0, BLoad, FactorLabels, FacDim, MatLength, FEVDhoriz, YieldsLabel,
                         ModelType, Economy = NULL, CholFac_JLL = NULL, PI = NULL, Mode = FALSE){

  # 1) Dynamic multipliers
  Ry.h <- array(NA, c(FacDim, FacDim, FEVDhoriz))
  Ry.h[, , 1] <- diag(FacDim) # dynamic multiplier at t=0

  for (l in 2:FEVDhoriz) { Ry.h[, ,l] <- K1Z%*%Ry.h[, ,l-1]   }

  # 2) Initialization
  vslct <- diag(FacDim)
  eslct <- diag(FacDim)

  # 2.1) Minor preliminary work
  invG <- diag(nrow(G0))/G0
  invG[!is.finite(invG)] <- 0
  invGSigmau <- solve(G0)%*%SIGMA

  # Choleski term
  if ( any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
    P <- CholFac_JLL }else{ P <- t(chol(invGSigmau))}

  # Adjustment term
  if (Mode == "Ortho" & any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {AdjTerm <- PI
  } else{ AdjTerm <- diag(FacDim) }

  scale <- 1

  # 3) FEVD
  # 3.1) Factors
  FEVDresFactors <- array(NA, c(nrow = FacDim, ncol = FacDim, FEVDhoriz))
  num <- matrix(0, nrow = FacDim, ncol = FacDim)
  den <- rep(0, times = FacDim)

  for (l in 1:FEVDhoriz){
    acc1 <- (eslct%*%Ry.h[ , , l]%*%P%*%vslct)^2
    num <- num + acc1
    acc2 <- diag(eslct%*%Ry.h[ , , l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%eslct)
    den<- den + acc2
    FEVDresFactors[ , , l] <- scale*num/den
  }

  FEVDFactors <- aperm(FEVDresFactors, c(3,2,1))

  # 3.2) Yields
  eslctCJ <- diag(MatLength)
  vslctCJ <- diag(FacDim)

  FEVDresYields <- array(NA, c(nrow = MatLength, ncol= FacDim, FEVDhoriz))
  num <- matrix(0, nrow = MatLength, ncol= FacDim)
  den <- matrix(0, nrow = MatLength, ncol= FacDim)

  for (l in 1:FEVDhoriz){
    acc1 <- (eslctCJ%*%BLoad%*%AdjTerm%*%Ry.h[,,l]%*%P%*%vslctCJ)^2
    num <- num + acc1
    acc2 <- diag(eslctCJ%*%BLoad%*%AdjTerm%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%t(AdjTerm)%*%t(BLoad)%*%eslctCJ)
    den<- den + acc2
    FEVDresYields[ ,,l] <- scale*num/den
  }

  FEVDYields <- aperm(FEVDresYields, c(3, 2, 1))

  # 4) Prepare labels
  if ( ModelType == "JPS original") {AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]])}
  else if(any(ModelType == c("JPS global", "GVAR single"))) {
    AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  }else{ AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

  # 5) Export outputs
  Horiz <- 1:(FEVDhoriz) # # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error

  dimnames(FEVDFactors) <- list(Horiz, AllFactorsLabels, AllFactorsLabels)
  dimnames(FEVDYields) <- list(Horiz, AllFactorsLabels, YieldsLabel)

  Out <- list(Factors = FEVDFactors, Yields = FEVDYields)

  return(Out)
}

########################################################################################################
#'Compute GFEVDs for all models
#'
#'@param SIGMA Variance-covariance matrix
#'@param K1Z Loading As
#'@param G0 contemporaneous terms
#'@param BLoad Loading Bs
#'@param FactorLabels List containing the label of factors
#'@param FacDim Dimension of the P-dynamics
#'@param MatLength Length of the maturity vector
#'@param GFEVDhoriz Horizon of the analysis
#'@param YieldsLabel Label of bond yields
#'@param ModelType Desired model type
#'@param Economy specific economy under study
#'@param PI matrix PI for JLL-based models
#'@param Mode allows for the orthogonalized version in the case of JLL-based models
#'
#'@references
#' \itemize{
#' \item This function is a modified and extended version of the "fevd" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#'@keywords internal

ComputeGFEVDs <- function(SIGMA, K1Z, G0, BLoad, FactorLabels, FacDim, MatLength, GFEVDhoriz, YieldsLabel,
                          ModelType, Economy, PI = NULL, Mode = FALSE){

  # 1) Dynamic multipliers
  Ry.h <- array(NA, c(FacDim, FacDim, GFEVDhoriz))
  Ry.h[, , 1] <- diag(FacDim) # dynamic multiplier at t=0

  for (l in 2:GFEVDhoriz) { Ry.h[, ,l] <- K1Z%*%Ry.h[, ,l-1]  }

  # 2) Initialization/ Minor preliminary work
  GFEVDresFac <- array(NA, c(nrow = FacDim, ncol= GFEVDhoriz, FacDim))
  vslct <- diag(FacDim)
  eslct <- diag(FacDim)

  invG <- diag(nrow(G0))/G0
  invG[!is.finite(invG)] <- 0
  invGSigmau <- solve(G0)%*%SIGMA

  scale <- 1/diag(SIGMA)

  # Adjustment term
  if (Mode == "Ortho" & any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {AdjTerm <- PI
  } else{ AdjTerm <- diag(FacDim) }


  # 3) GFEVD
  # 3.1) Factors
  for(h in 1:FacDim){
    n <- 1
    num <- matrix(0, nrow = FacDim, ncol= GFEVDhoriz)
    den <- matrix(0, nrow = FacDim, ncol= GFEVDhoriz)
    while (n <= GFEVDhoriz){
      for (l in 1:n){
        acc1 <- t((eslct[,h]%*%AdjTerm%*%Ry.h[,,l]%*%invGSigmau%*%vslct)^2) # Contribution of all j variables to explain i
        num[,n] <- num[ ,n] + acc1
        acc2 <- eslct[,h]%*%AdjTerm%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%eslct[,h]
        den[,n]<- den[,n] + matrix(1, nrow = FacDim)%*%acc2
      }
      GFEVDresFac[ , n, h] <- t(t(scale*num[,n]))/den[,n]
      n <- n+1
    }
  }

  GFEVDFactors <- aperm(GFEVDresFac, c(2,1,3)) # Non-normalized GFEVD (i.e. rows need not sum up to 1)

  #  Normalization of the GFEVD for the factors
  # (Make sure that the sum of the errors equal to one in each period)
  DEM <- array(NA, c(nrow = GFEVDhoriz, ncol=1, FacDim))
  GFEVDFactorsNormalized <- array(NA, c(nrow = GFEVDhoriz, ncol= FacDim, FacDim))

  for (h in 1:FacDim){
    for (n in 1:GFEVDhoriz){
      DEM[n, 1, h] <- sum(GFEVDFactors[n,,h])
      GFEVDFactorsNormalized[n,,h] <- GFEVDFactors[n,,h]/DEM[n,,h]
    }
  }

  # 3.2) Yields
  # Initialization
  GFEVDresYie <- array(NA, c(nrow = MatLength, ncol= FacDim, GFEVDhoriz))
  vslctYie <- diag(FacDim)
  eslctYie <- diag(MatLength)

  num <- matrix(0, nrow = MatLength, ncol=FacDim)
  den <- matrix(0, nrow = MatLength, ncol=FacDim)

  for (l in 1:GFEVDhoriz){
    acc1 <- (eslctYie%*%BLoad%*%AdjTerm%*%Ry.h[,,l]%*%invGSigmau%*%vslctYie)^2
    num <- num + acc1
    acc2 <- diag(eslctYie%*%BLoad%*%AdjTerm%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%t(AdjTerm)%*%t(BLoad)%*%eslctYie)
    den <- den + acc2
    for (q in 1:FacDim){
      GFEVDresYie[ ,q,l] <- scale[q]*(num/den)[,q] # note: unlike the GFEVD of the factors, note that the "scale" variable is now at the acc1
    }
  }

  GFEVDYields <- aperm(GFEVDresYie, c(3,2,1)) # Non-normalized GFEVD (i.e. rows need not sum up to 1)

  #  Normalization of the GFEVD for the factors
  # (Make sure that the sum of the errors equal to one in each period)
  DEM <- array(NA, c(nrow = GFEVDhoriz, ncol=1, MatLength))
  GFEVDYieldsNormalized <- array(NA, c(nrow = GFEVDhoriz, ncol= FacDim, MatLength))

  for (h in 1:MatLength){
    for (n in 1:GFEVDhoriz){
      DEM[n, 1, h] <- sum(GFEVDYields[n,,h])
      GFEVDYieldsNormalized[n,,h] <- GFEVDYields[n,,h]/DEM[n,,h]
    }
  }

  # 4) Prepare labels
  if ( ModelType == "JPS original") {AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]])}
  else if(any(ModelType == c("JPS global", "GVAR single"))) {
    AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  }else{ AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

  # 5) Outputs to export
  Horiz <- 1:(GFEVDhoriz) # # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error

  dimnames(GFEVDFactorsNormalized) <- list(Horiz, AllFactorsLabels, AllFactorsLabels)
  dimnames(GFEVDYieldsNormalized) <- list(Horiz, AllFactorsLabels, YieldsLabel)

  Out <- list(Factors = GFEVDFactorsNormalized, Yields = GFEVDYieldsNormalized)

  return(Out)
}
########################################################################################################
########################################################################################################
#' Fit yields for all maturities of interest
#'
#'@param MatInt numerical vector containing the fit maturities of interest
#'@param ModelPara List of model parameter estimates (See the "Optimization" function)
#'@param FactorLabels a string-list based which contains all the labels of all the variables present in the model
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'@param YLab Label of yields ("Months" or "Yields")
#'
#'@keywords internal

YieldsFitAllSep <- function(MatInt, ModelPara, FactorLabels, ModelType, Economies, YLab){


  dt <- ModelPara[[ModelType]][[Economies[1]]]$inputs$dt
  mat <- ModelPara[[ModelType]][[Economies[1]]]$inputs$mat

  C <- length(Economies)
  N <- ModelPara[[ModelType]][[Economies[1]]]$inputs$N
  T <- ncol(ModelPara[[ModelType]][[Economies[1]]]$inputs$AllFactors)
  J <- length(mat)
  M <- length(FactorLabels$Domestic) - N
  G <- length(FactorLabels$Global)

# Initialization of inputs to store outputs
  FittedYieldsPerMat <- list()

  FitLat <- matrix(NA, nrow= length(MatInt), ncol = T)
  colnames(FitLat) <- colnames(ModelPara[[ModelType]][[Economies[1]]]$inputs$AllFactors)
  rownames(FitLat) <- paste(MatInt,YLab, sep="")


  for (i in 1:C){
  BnX <- ModelPara[[ModelType]][[Economies[i]]]$rot$X$B
  AnX <- ModelPara[[ModelType]][[Economies[i]]]$rot$X$A
  K1XQ <- ModelPara[[ModelType]][[Economies[i]]]$ests$K1XQ
  SSX <- ModelPara[[ModelType]][[Economies[i]]]$rot$X$Q$SS
  r0 <- ModelPara[[ModelType]][[Economies[i]]]$ests$r0
  Wpca <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Wpca
  ZZ <- ModelPara[[ModelType]][[Economies[i]]]$inputs$AllFactors


if (ModelType == "JPS original"){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]]) }
else if (any(ModelType == c("JPS global", 'GVAR single'))){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

  LabelSpannedCS <- c(FactorLabels$Tables[[Economies[i]]][-(1:M)])
  b <- match(LabelSpannedCS, AllLabels)

  PP <- ZZ[b , ]

  X <- matrix(NA, nrow=nrow(PP), ncol = ncol(PP))

  for (t in 1:T){    X[,t] <- solve(Wpca%*%BnX, tol = 1e-50)%*%(PP[,t] - Wpca%*%AnX)  }

  # Loadings of the other maturities
  MatTU <- mat/dt
  MatAll <- 1:max(MatTU)

  LoadingsLat <- A0N__BnAn(MatAll, K1XQ, ModelType, dX= NULL, r0, SSX, Economies)
  AnXAll <- (LoadingsLat$AnX)/dt
  BnXAll <- (LoadingsLat$BnX)/dt

  for(h in 1:length(MatInt)){

    IdxMatInt <- seq(MatInt[h], length(AnXAll), by = max(MatTU))

    AnXInt <- AnXAll[IdxMatInt]
    BnXInt <- BnXAll[IdxMatInt,]

    for (t in 1:T){      FitLat[h,t] <- AnXInt + BnXInt%*%X[,t]    }
  }

  FittedYieldsPerMat[[i]] <- FitLat
  }

  names(FittedYieldsPerMat) <- Economies

  return(FittedYieldsPerMat)
}


########################################################################################################
########################################################################################################
#' Fit yields for all maturities of interest
#'
#'@param MatInt numerical vector containing the fit maturities of interest
#'@param ModelPara List of model parameter estimates (See the "Optimization" function)
#'@param FactorLabels a string-list based which contains all the labels of all the variables present in the model
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'@param YLab Label of yields ("Months" or "Yields")
#'
#'@keywords internal

YieldsFitAllJoint <- function(MatInt, ModelPara, FactorLabels, ModelType, Economies, YLab){

  BnX <- ModelPara[[ModelType]]$rot$X$B
  AnX <- ModelPara[[ModelType]]$rot$X$A
  K1XQ <- ModelPara[[ModelType]]$ests$K1XQ
  SSX <- ModelPara[[ModelType]]$rot$X$Q$SS
  r0 <- ModelPara[[ModelType]]$ests$r0
  Wpca <- ModelPara[[ModelType]]$inputs$Wpca
  ZZ <- ModelPara[[ModelType]]$inputs$AllFactors

  dt <- ModelPara[[ModelType]]$inputs$dt
  mat <- ModelPara[[ModelType]]$inputs$mat


  C <- length(Economies)
  N <- ModelPara[[ModelType]]$inputs$N
  T <- ncol(ZZ)
  J <- length(mat)
  M <- length(FactorLabels$Domestic) - N
  G <- length(FactorLabels$Global)

  b <- IdxSpanned(G,M,N,C)
  PP <- ZZ[b , ]

  X <- matrix(NA, nrow=nrow(PP), ncol = ncol(PP))

  for (t in 1:T){   X[,t] <- solve(Wpca%*%BnX, tol = 1e-50)%*%(PP[,t] - Wpca%*%AnX)  }

  # Loadings of the other maturities
  MatTU <- mat/dt
  MatAll <- 1:max(MatTU)

  LoadingsLat <- A0N__BnAn(MatAll, K1XQ, ModelType, dX= NULL, r0, SSX, Economies)
  AnXAll <- (LoadingsLat$AnX)/dt
  BnXAll <- (LoadingsLat$BnX)/dt

  FittedYieldsPerMat <- list()

  m <- length(MatInt)
  for(h in 1:m){

    IdxMatInt <- seq(MatInt[h], length(AnXAll), by = max(MatTU))

    AnXInt <- AnXAll[IdxMatInt]
    BnXInt <- BnXAll[IdxMatInt,]

    FitLat <- matrix(NA, nrow=C, ncol = ncol(PP))
    colnames(FitLat) <- colnames(ZZ)
    rownames(FitLat) <- Economies


    for (t in 1:T){      FitLat[,t] <- AnXInt + BnXInt%*%X[, t]    }


    FittedYieldsPerMat[[h]] <- FitLat
  }


  names(FittedYieldsPerMat) <- paste(MatInt,YLab, sep="")


  # Reorganize fitted yields per country
  FittedYieldsCS <- list()

  for (i in 1:C){
    FitCS <- do.call("rbind", lapply(FittedYieldsPerMat, "[", i, ))
    FittedYieldsCS[[Economies[i]]] <- FitCS
  }

  return(FittedYieldsCS)
}
#################################################################################################
#'Adjust vector of maturities
#'
#'@param mat vector of maturities (J x 1)
#'@param UnitYields Available options: "Month" and "Year"
#'
#'@keywords internal

MatAdjusted <- function(mat, UnitYields){

  if (UnitYields== "Month"){  k <- 12
  } else  if (UnitYields== "Year"){  k <- 1}
  matAdjUnit <- mat*k

  return(matAdjUnit)
}
########################################################################################################
#'Get the expected component of all models
#'
#'@param ModelPara list of model parameter estimates
#'@param InputsForOutputs list containing the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs,
#'                        and risk premia decomposition
#'@param ModelType desired model type
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels string-list based which contains all the labels of all the variables present in the model
#'@param WishFP If users wants to compute the forward premia. Default is FALSE.
#'
#'@keywords internal

ExpectedComponent <- function(ModelPara, InputsForOutputs, ModelType, Economies, FactorLabels, WishFP = FALSE){

  # 0) Preliminary work
  N <- length(FactorLabels$Spanned)
  UnitYields <- InputsForOutputs$UnitMatYields

  if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
    mat <- ModelPara[[ModelType]][[1]]$inputs$mat} else{ mat <- ModelPara[[ModelType]]$inputs$mat}

  matAdjUnit <- MatAdjusted(mat, UnitYields)

  if( WishFP){
    matMIN <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[1]
    matMAX <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[2]
  }

  #1) Compute risk-neutral parameters
  rhos <- rhoParas(ModelPara, N, ModelType, Economies)

  # 2) Compute the expectation component
  # Pure expected component
  avexp <- Compute_EP(ModelPara, ModelType, UnitYields, matAdjUnit, N, rhos, Economies, FactorLabels)

  # expected component related to the forward premia (if desired)
  if( WishFP){
    avexpFP <- Compute_EP(ModelPara, ModelType, UnitYields, matAdjUnit, N, rhos, Economies, FactorLabels,
                          WishFP, matMIN, matMAX)
  } else{ avexpFP <- list() }

  return(list(avexp = avexp, avexpFP = avexpFP))
}

#################################################################################################
#'Compute  risk-neutral intercept and slope
#'
#'@param ModelPara list of model parameter estimates
#'@param N number of country-specific spanned factors
#'@param ModelType desired model type
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal


rhoParas <- function(ModelPara, N, ModelType, Economies){

  # Compute the intercept and slope coefficients of the short rate expressed as a function of the spanned factors
  # By definition: r_t = r0 + rho1_X* X_t
  # But X_t = (W*Bx)^(-1)(P- WAx)
  # so r_t = rho0_PP + rho1_PP*P_t
  # where (i) rho0_PP = r0 - rho1_X*(W*BX)^(-1)W*AX and (ii) rho1_PP = rho1_X (W*BX)^(-1)

  # 1) Models estimated for countries SEPARETLY
  if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){

  rho0_PP <- vector(mode = "list", length = length(Economies))
  rho1_PP <- vector(mode = "list", length = length(Economies))
  names(rho0_PP) <- Economies
  names(rho1_PP) <- Economies

  dt <- ModelPara[[ModelType]][[Economies[1]]]$inputs$dt

  for (i in 1:length(Economies)){
    BnX <- ModelPara[[ModelType]][[Economies[i]]]$rot$X$B
    AnX <- ModelPara[[ModelType]][[Economies[i]]]$rot$X$A
    Wpca <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Wpca
    r0 <- ModelPara[[ModelType]][[Economies[i]]]$ests$r0

    # Compute rhos
    rho1_X <- rep(1,N)
    rho0_PP[[i]] <- as.numeric((r0 - rho1_X%*%solve(Wpca%*%BnX,tol = 1e-50)%*%Wpca%*%AnX)/dt)
    rho1_PP[[i]] <- (rho1_X%*%solve(Wpca%*%BnX, tol = 1e-50))/dt
  }

  }else{

    # 2) Models estimated for countries JOINTLY
    dt <- ModelPara[[ModelType]]$inputs$dt

    BnX <- ModelPara[[ModelType]]$rot$X$B
    AnX <- ModelPara[[ModelType]]$rot$X$A
    Wpca <- ModelPara[[ModelType]]$inputs$Wpca
    r0 <- ModelPara[[ModelType]]$ests$r0

    # Compute rhos
    rho1_X_CS <- rep(1,N)
    rho1_X <- matrix(0, nrow=length(Economies), ncol=N*length(Economies))

    idx0 <-0
    for (j in 1:length(Economies)){
      idx1<- idx0+N
      rho1_X[j,(idx0+1):idx1] <- rho1_X_CS
      idx0<- idx1
    }

    rho0_PP <- (r0 - rho1_X%*%solve(Wpca%*%BnX, tol = 1e-50)%*%Wpca%*%AnX)/dt
    rho1_PP <- (rho1_X%*%solve(Wpca%*%BnX, tol = 1e-50))/dt
  }
  return(list(rho0_PP = rho0_PP, rho1_PP = rho1_PP))
  }

###################################################################################################
#'Compute the expected component for all models
#'
#'@param ModelPara list of model parameter estimates
#'@param ModelType Desired model type
#'@param UnitYields Available options: "Month" and "Year"
#'@param matAdjUnit Adjusted vector of matutities
#'@param N number of country-specific spanned factors
#'@param rhoList List of risk-neutral parameters
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels List of factor labels
#'@param WishFP If users wants to compute the forward premia. Default is FALSE.
#'@param matMIN For the forward premia, the shortest maturity of the remium of interest
#'@param matMAX For the forward premia, the longest maturity of the remium of interest
#'
#'@keywords internal

Compute_EP <- function(ModelPara, ModelType, UnitYields, matAdjUnit, N, rhoList, Economies, FactorLabels,
                       WishFP = FALSE, matMIN = FALSE, matMAX = FALSE){

  # 0) Preliminary work
  SepQ_Lab <- c("JPS original", "JPS global", "GVAR single")
  M <- length(FactorLabels$Domestic) - N

  if (ModelType %in% SepQ_Lab){
    dt <- ModelPara[[ModelType]][[Economies[1]]]$inputs$dt
  }else{dt <- ModelPara[[ModelType]]$inputs$dt}

  if (UnitYields == "Month"){ YLab <- "M"; k <- 12} else { YLab <- "Y"; k <- 1} #  second case is for UnitYields = "Year"

  if( WishFP){ # Forward premia features
    EP_Lab <- "FP_" ;  ExpecCompLength <- 1
    matMINAdj <- round((matMIN/k)/dt);    matMAXAdj <- round((matMAX/k)/dt)
    }else{
    EP_Lab <- "RP_";  ExpecCompLength <-  round((matAdjUnit/k)/dt)
  }

  # 1) Compute the expected component
  avexp <- list()
  for (i in 1:length(Economies)){

    # a) Pre-allocation
    if (ModelType %in% SepQ_Lab){ # SepQ models
    K0Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K0Z
    K1Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K1Z
    ZZ <- ModelPara[[ModelType]][[Economies[i]]]$inputs$AllFactors
    }else{ # jointQ models
    K0Z <- ModelPara[[ModelType]]$ests$K0Z
    K1Z <- ModelPara[[ModelType]]$ests$K1Z
    ZZ <- ModelPara[[ModelType]]$inputs$AllFactors
    }

    # b) Extract spanned factors from the list of unspanned factors
    if (ModelType %in% SepQ_Lab){
    if (ModelType == "JPS original"){
      AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]])
    }else{  AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries) }
    LabelSpannedCS <- c(FactorLabels$Tables[[Economies[i]]][-(1:M)])
    IdxSpanned <- match(LabelSpannedCS, AllLabels)
    }else{
    IdxSpanned <- IdxAllSpanned(ModelType, FactorLabels, Economies)
    }

    # c) Get the expected component
    avexpCS <- matrix(0, ncol(ZZ), length(ExpecCompLength))
    dimnames(avexpCS) <- list(colnames(ZZ), paste(EP_Lab, ExpecCompLength, YLab, sep=""))

    for (h in 1:length(ExpecCompLength)){ # per bond maturity
      for (t in 1:ncol(ZZ)){ # Per point in time

        # Initialiation
        if( WishFP){g <- matrix(NA, nrow(K0Z), matMAXAdj) }else{g <- matrix(NA, nrow(K0Z), ExpecCompLength[h])}
        rownames(g) <- rownames(ZZ)

        # Fitted P-dynamics
        g[ ,1] <- ZZ[, t]
        if( WishFP){loopLim <- matMAXAdj}else{loopLim <- ExpecCompLength[h]}
        for (j in 2:loopLim){g[ ,j] <- K0Z + K1Z%*%g[ ,j-1]}

        # Adjust `g_lim` dynamically
        g_lim <- if (WishFP) {round((matMIN / k) / dt):round((matMAX / k) / dt)
          } else {  seq_len(ncol(g)) }

        g <- g[IdxSpanned, g_lim] # extract relevant variables

        # Adjustment term
        if (ModelType %in% SepQ_Lab){ MaxExpec <- pmax(rhoList$rho0_PP[[i]] + (rhoList$rho1_PP[[i]]%*%g),0)
        }else{ MaxExpec <- pmax(rhoList$rho0_PP[i] + (rhoList$rho1_PP[i, ]%*%g),0)}

        avexpCS[t,h] <- mean(MaxExpec)
      }
    }

    avexp[[Economies[i]]] <- avexpCS*100
    }

  return(avexp)
}
###################################################################################################
#'Compute the term premia
#'
#'@param ModelPara list of model parameter estimates
#'@param avexp list containing the country-specific pure expected component
#'@param ModelType desired model type
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal

TermPremia <- function(ModelPara, avexp, ModelType, Economies){

  # 0) Preliminary work
  YieldData <- list()
  TermPremium <- list()

  SepQ_Lab <- c("JPS original", "JPS global", "GVAR single")

  if(!ModelType %in% SepQ_Lab){
    mat <- ModelPara[[ModelType]]$inputs$mat
    J <- length(mat)
  }

 # 1) Compute the term premia
  for (i in 1:length(Economies)){

    if(ModelType %in% SepQ_Lab){ # SepQ models
      Y <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y
      YieldData[[Economies[i]]] <- t(Y)*100
    } else{ # jointQ models
      IdxRP <- (1:J) +J*(i-1)
      Y <- ModelPara[[ModelType]]$inputs$Y
      YieldData[[Economies[i]]] <- t(Y[IdxRP, ]*100)
    }
    TermPremium[[Economies[i]]] <- YieldData[[Economies[i]]] - avexp[[Economies[i]]]
  }

  Output <- list(TermPremium, avexp)
  names(Output) <- c("Term Premia","Expected Component")

  return(Output)
}
###############################################################################################
#' Compute the forward premia for all models
#'
#'@param ModelPara list of model parameter estimates
#'@param avexpFP list containing the country-specific expected component of the forward period
#'@param ModelType desired model type
#'@param FactorLabels List of factor labels
#'@param InputsForOutputs list containing the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs,
#'                        and risk premia decomposition
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal

ForwardPremia <- function(ModelPara, avexpFP, ModelType, FactorLabels, InputsForOutputs, Economies){

  # 0) Preliminary work: redefine necessary inputs
  SepQ_Lab <- c("JPS original", "JPS global", "GVAR single")

  if (ModelType %in% SepQ_Lab){    mat <- ModelPara[[ModelType]][[1]]$inputs$mat
  }else{ mat <- ModelPara[[ModelType]]$inputs$mat}

  J <- length(mat)

  # Yield labels
  UnitYields <- InputsForOutputs$UnitMatYields
  if (UnitYields== "Month"){  YLab <- "M"} else{   YLab <- "Y"  }

  # Inputs from the forward premia specification
  matAdjUnit <- MatAdjusted(mat, UnitYields)

  matMIN <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[1]
  matMAX <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[2]
  IDXMatMIN <- match(matMIN,matAdjUnit)
  IDXMatMAX <- match(matMAX,matAdjUnit)
  IDXMatBoth <- c(IDXMatMIN,IDXMatMAX)

  FR <- list()
  ForwardPremium <- list()

  # CASE 1: If any of the data of the specific maturity were not used in the estimation, then compute the fitted value
  if( anyNA(IDXMatBoth)){

    MatBoth <- c(matMIN, matMAX)
    IdxNA <- which(is.na(IDXMatBoth))

    # a) Missing maturity (model-implied)
    MissingMat <-  MatBoth[IdxNA] # Maturity not available
    if ( ModelType %in% SepQ_Lab){
      YieldMissing <- YieldsFitAllSep(MissingMat, ModelPara, FactorLabels, ModelType, Economies, YLab)
    }else{ YieldMissing <- YieldsFitAllJoint(MissingMat, ModelPara, FactorLabels, ModelType, Economies, YLab) }


    for (i in 1:length(Economies)){

      if ( ModelType %in% SepQ_Lab){ Y <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y
      } else{ Y <- ModelPara[[ModelType]]$inputs$Y }

      # b) If both maturities are missing
      if (length(MissingMat) == 2){
        YieldMIN <- t(t(YieldMissing[[Economies[i]]][1,]))*100
        YieldMAX <- t(t(YieldMissing[[Economies[i]]][2,]))*100
      } else{
        # c) If only one maturity is missing
        # Available maturity
        IdxNotMissing0 <- IDXMatBoth[!is.na(IDXMatBoth)]
        if ( ModelType %in% SepQ_Lab){ IdxNotMissingCS <- IdxNotMissing0
        }else{IdxNotMissingCS <- IdxNotMissing0 + J*(i-1) }

        YieldNotMissing <- t(t(Y[IdxNotMissingCS, ]))

        if (MissingMat ==1){
          YieldMIN <- YieldNotMissing*100
          YieldMAX <- t(YieldMissing[[Economies[i]]])*100
        }else{
          YieldMIN <- t(YieldMissing[[Economies[i]]])*100
          YieldMAX <- YieldNotMissing*100
        }
      }

      FR[[Economies[i]]] <- (matMAX*YieldMAX -  matMIN*YieldMIN)/(matMAX -  matMIN) # Fitted forward rate
      ForwardPremium[[Economies[i]]] <- FR[[Economies[i]]] - avexpFP[[Economies[i]]] # Forward Premia
      colnames(ForwardPremium[[Economies[i]]]) <- paste("FP_", matMIN, "-", matMAX, YLab, sep="")
      colnames(FR[[Economies[i]]]) <- paste("Mat", matMIN, "-", matMAX, YLab, sep="")
    }

  }else{

    # CASE 2: when all data is available
    YieldData <- list()

    for (i in 1:length(Economies)){

      # SepQ models
      if (ModelType %in% SepQ_Lab){
      Y <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y
      YieldData[[Economies[i]]] <- t(Y)*100

      YieldMIN <- t(t(YieldData[[Economies[i]]][ , IDXMatMIN]))
      YieldMAX <- t(t(YieldData[[Economies[i]]][ , IDXMatMAX]))

      }else{
      # jointQ models
      Y <- ModelPara[[ModelType]]$inputs$Y

      IdxMinCS <- IDXMatMIN + J*(i-1)
      IdxMaxCS <- IDXMatMAX + J*(i-1)
      YieldMIN <- t(t(Y[IdxMinCS, ]*100))
      YieldMAX <- t(t(Y[IdxMaxCS, ]*100))
      }

      FR[[Economies[i]]] <- (matMAX*YieldMAX -  matMIN*YieldMIN)/(matMAX -  matMIN) # Fitted forward rate
      ForwardPremium[[Economies[i]]] <- FR[[Economies[i]]] - avexpFP[[Economies[i]]] # Forward Premia
      colnames(ForwardPremium[[Economies[i]]]) <- paste("FP_", matMIN, "-", matMAX, YLab, sep="")
      colnames(FR[[Economies[i]]]) <- paste("Mat", matMIN, "-", matMAX, YLab, sep="")
    }
  }

  OutputFP <- list(ForwardPremium, avexpFP, FR)
  names(OutputFP) <- c("Forward Premia","Expected Component", "Forward Rate")

  return(OutputFP)
}

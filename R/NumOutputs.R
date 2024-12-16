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

  # If one chooseS models in which the estimation is done country-by-country
  if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
    NumOutSep <- OutputConstructionSep(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies)

  }
  # If one chooseS models in which the estimation is done jointly for all countries
  else { NumOutJoint <- OutputConstructionJoint(ModelType, ModelPara,InputsForOutputs, FactorLabels, Economies) }



  # Prepare final list of outputs
  if (!exists("NumOutJoint")){ AllNumOutputs <- NumOutSep}
  if (!exists("NumOutSep")){ AllNumOutputs <- NumOutJoint}
  if ( exists("NumOutSep") & exists("NumOutJoint")){

    for (i in 1:length(NumOutSep)){AllNumOutputs[[i]] <- append(NumOutSep[[i]], NumOutJoint[[i]]) }

    names(AllNumOutputs) <- names(NumOutSep)
  }

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
####################### OUTPUTS FOR MODELS IN WHICH THE ESTIMATION ###################################
########################       IS DONE COUNTRY-BY-COUNTRY      #######################################
######################################################################################################
######################################################################################################

#' Numerical outputs (variance explained, model fit, IRFs, GIRFs, FEVDs, GFEVDs, and risk premia decomposition)
#' for "sep Q" models
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (See the "Optimization" function)
#'@param InputsForOutputs list conataining the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs,
#'                        and risk premia decomposition
#'@param FactorLabels string-list based which contains all the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal



OutputConstructionSep <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies){


  # Output summary
  # Total Variance Explained and Model fit
  Total_Var_exp <- VarianceExplained(ModelType, ModelPara, FactorLabels, Economies)
  ModFit <- YieldsFit(ModelType, ModelPara, FactorLabels, Economies)
  cat(" ** Model Fit \n")

  # IRF and GIRF
  IRFout <- IRFandGIRF(ModelType, ModelPara, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)
  cat(" ** IRFs and GIRFs  \n")

  FEVDoutputs <- FEVDsep(ModelType, ModelPara, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels, Economies)
  cat(" ** FEVD \n")

  # GIRF and GFEVD
  GFEVDoutputs <- GFEVDsep(ModelType, ModelPara, InputsForOutputs[[ModelType]]$GFEVD$horiz, FactorLabels, Economies)
  cat(" ** GFEVD \n")

  # Risk Premia Decomposition
  TermPremia <- TermPremiaDecompSep(ModelPara, FactorLabels, ModelType, InputsForOutputs, Economies)
  cat(" ** Term Premia \n")

  NumericalOutputs <- list(Total_Var_exp, ModFit, IRFout$IRFs, FEVDoutputs, IRFout$GIRFs, GFEVDoutputs, TermPremia)
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
  J <- numel(ModelPara[[ModelType]]$inputs$mat)
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
########################################### 3) IRFs ##################################################
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
################################### 3) FEVD #############################################################
#########################################################################################################
#' FEVDs for "sep Q" models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param FEVDhoriz single numerical vector conataining the desired horizon of analysis for the FEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@details
#' Structural shocks are identified via Cholesky decomposition
#'
#'@keywords internal


FEVDsep <- function(ModelType, ModelPara, FEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N

  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR single")] # Exclude all models in which the estimation is made jointly.

  FEVDoutputs <- list()

  for (j in  idxIndividual){
    for (i in 1:C){
      N <- length(FactorLabels$Spanned)
      M <- length(FactorLabels$Domestic) - N
      K <- nrow(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$AllFactors)
      J <- length(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$mat)

      G0 <-  ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$Gy.0
      Sigma_y <- ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$SSZ
      F1 <- ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$K1Z


      # 1) Dynamic multipliers
      Ry.h <- array(NA, c(K,K,FEVDhoriz))
      Ry.h[, ,1] <- diag(K) # dynamic multiplier at t=0

      for (l in 2:FEVDhoriz) {
        Ry.h[, ,l] <- F1%*%Ry.h[, ,l-1]
      }

      # 2) Initialization
      vslct <- diag(K)
      eslct <- diag(K)

      # 2.1) Minor preliminary work
      invG <- diag(nrow(G0))/G0
      invG[!is.finite(invG)] <- 0
      invGSigmau <- solve(G0)%*%Sigma_y

      P <- t(chol(invGSigmau))
      scale <- 1

      # 2.2) Factor loadings preparation
      BallFac <- BUnspannedAdapSep(G,M, ModelPara, Economies, Economies[i], ModelTypeSet[j])


      # 3) FEVD
      # 3.1) Factors
      FEVDresFactors <- array(NA, c(nrow = K, ncol=K, FEVDhoriz))
      num <- matrix(0, nrow =K, ncol=K)
      den <- rep(0, times = K)

      for (l in 1:FEVDhoriz){
        acc1 <- (eslct%*%Ry.h[,,l]%*%P%*%vslct)^2
        num <- num + acc1
        acc2 <- diag(eslct%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%eslct)
        den<- den + acc2
        FEVDresFactors[ ,,l] <- scale*num/den
      }

      FEVDFactors <- aperm(FEVDresFactors, c(3,2,1))

      # 3.2) Yields
      eslctCJ <- diag(J)
      vslctCJ <- diag(K)


      FEVDresYields <- array(NA, c(nrow = J, ncol=K, FEVDhoriz))
      num <- matrix(0, nrow =J, ncol=K)
      den <- matrix(0, nrow =J, ncol=K)

      for (l in 1:FEVDhoriz){
        acc1 <- (eslctCJ%*%BallFac%*%Ry.h[,,l]%*%P%*%vslctCJ)^2
        num <- num + acc1
        acc2 <- diag(eslctCJ%*%BallFac%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%t(BallFac)%*%eslctCJ)
        den<- den + acc2
        FEVDresYields[ ,,l] <- scale*num/den
      }


      FEVDYields <- aperm(FEVDresYields, c(3, 2, 1))

      # 4) Prepare labels
      if ( ModelTypeSet[j] == "JPS original") {AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]])}
      else {AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries) }

      YieldsLabel<- rownames(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$Y)

      dimnames(FEVDFactors)[[1]] <- 1:(FEVDhoriz) # # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error
      dimnames(FEVDFactors)[[2]] <- AllFactorsLabels
      dimnames(FEVDFactors)[[3]] <- AllFactorsLabels


      dimnames(FEVDYields)[[1]] <- 1:(FEVDhoriz) # # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error
      dimnames(FEVDYields)[[2]] <- AllFactorsLabels
      dimnames(FEVDYields)[[3]] <- YieldsLabel


      FEVDoutputs[[ModelTypeSet[j]]][[Economies[i]]] <- list(FEVDFactors,FEVDYields)
      names(FEVDoutputs[[ModelTypeSet[j]]][[Economies[i]]]) <- c("Factors","Yields")


    }
  }
  # Clean empty lists
  FEVDoutputs  <- FEVDoutputs[unlist(lapply(FEVDoutputs, length) != 0)]

  names(FEVDoutputs) <- ModelTypeSet[idxIndividual]


  return(FEVDoutputs)

}

#########################################################################################################
################################### 5) GFEVD #############################################################
#########################################################################################################
#' GFEVDs for "sep Q" models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param GFEVDhoriz single numerical vector conataining the desired horizon of analysis for the GFEVDs
#'@param FactorLabels string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#' @references
#' \itemize{
#' \item This function is a modified and extended version of the "fevd" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#'@keywords internal




GFEVDsep <- function(ModelType, ModelPara, GFEVDhoriz, FactorLabels, Economies){


  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N

  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR single")] # Exclude all models in which the estimation is made jointly.

  GFEVDoutputs <- list()

  for (j in  idxIndividual){
    for (i in 1:C){
      N <- length(FactorLabels$Spanned)
      M <- length(FactorLabels$Domestic) - N
      K <- nrow(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$AllFactors)
      J <- length(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$mat)

      G0 <-  ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$Gy.0
      Sigma_y <- ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$SSZ
      F1 <- ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$K1Z

      # 1) Dynamic multipliers
      Ry.h <- array(NA, c(K,K,GFEVDhoriz))
      Ry.h[, ,1] <- diag(K) # dynamic multiplier at t=0

      for (l in 2:GFEVDhoriz) {
        Ry.h[, ,l] <- F1%*%Ry.h[, ,l-1]
      }

      # 2) Initialization/ Minor preliminary work
      GFEVDresFac <- array(NA, c(nrow = K, ncol=GFEVDhoriz,K))
      vslct <- diag(K)
      eslct <- diag(K)

      invG <- diag(nrow(G0))/G0
      invG[!is.finite(invG)] <- 0
      invGSigmau <- solve(G0)%*%Sigma_y

      scale <- 1/diag(Sigma_y)

      # 3) GFEVD
      # 3.1) Factors
      for(h in 1:K){
        n<-1
        num <- matrix(0, nrow =K, ncol=GFEVDhoriz)
        den <- matrix(0, nrow =K, ncol=GFEVDhoriz)
        while (n <= GFEVDhoriz){
          for (l in 1:n){
            acc1 <- t((eslct[,h]%*%Ry.h[,,l]%*%invGSigmau%*%vslct)^2) # Contribution of all j variables to explain i
            num[,n] <- num[ ,n] + acc1
            acc2 <- eslct[,h]%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%eslct[,h]
            den[,n]<- den[,n] + matrix(1, nrow=K)%*%acc2
          }
          GFEVDresFac[ ,n,h] <- t(t(scale*num[,n]))/den[,n]
          n <- n+1
        }
      }

      GFEVDFactors <- aperm(GFEVDresFac, c(2,1,3)) # Non-normalized GFEVD (i.e. rows need not sum up to 1)

      #  Normalization of the GFEVD for the factors
      # (Make sure that the sum of the errors equal to one in each period)
      DEM <- array(NA, c(nrow = GFEVDhoriz, ncol=1, K))
      GFEVDFactorsNormalized <- array(NA, c(nrow = GFEVDhoriz, ncol=K, K))

      for (h in 1:K){
        for (n in 1:GFEVDhoriz){
          DEM[n, 1, h] <- sum(GFEVDFactors[n,,h])
          GFEVDFactorsNormalized[n,,h] <- GFEVDFactors[n,,h]/DEM[n,,h]
        }
      }

      # 3.2) Yields
      # Get the full B
      Bfull <- BUnspannedAdapSep(G,M, ModelPara, Economies, Economies[i], ModelTypeSet[j])


      # Initialization
      GFEVDresYie <- array(NA, c(nrow = J, ncol=K, GFEVDhoriz))
      vslctYie <- diag(K)
      eslctYie <- diag(J)


      num <- matrix(0, nrow =J, ncol=K)
      den <- matrix(0, nrow =J, ncol=K)

      for (l in 1:GFEVDhoriz){
        acc1 <- (eslctYie%*%Bfull%*%Ry.h[,,l]%*%invGSigmau%*%vslctYie)^2
        num <- num + acc1
        acc2 <- diag(eslctYie%*%Bfull%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%t(Bfull)%*%eslctYie)
        den <- den + acc2
        for (q in 1:K){
          GFEVDresYie[ ,q,l] <- scale[q]*(num/den)[,q] # note: unlike the GFEVD of the factors, note that the "scale" variable is now at the acc1
        }
      }


      GFEVDYields <- aperm(GFEVDresYie, c(3,2,1)) # Non-normalized GFEVD (i.e. rows need not sum up to 1)

      #  Normalization of the GFEVD for the factors
      # (Make sure that the sum of the errors equal to one in each period)
      DEM <- array(NA, c(nrow = GFEVDhoriz, ncol=1, J))
      GFEVDYieldsNormalized <- array(NA, c(nrow = GFEVDhoriz, ncol=K, J))

      for (h in 1:J){
        for (n in 1:GFEVDhoriz){
          DEM[n, 1, h] <- sum(GFEVDYields[n,,h])
          GFEVDYieldsNormalized[n,,h] <- GFEVDYields[n,,h]/DEM[n,,h]
        }
      }


      # 4) Adjust labels:
      if ( ModelTypeSet[j] == "JPS original") {AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]])}
      if ( ModelTypeSet[j] != "JPS original") {AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries) }

      YieldsLabel<- rownames(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$Y)

      dimnames(GFEVDFactorsNormalized)[[1]] <- 1:(GFEVDhoriz) # We subtract 1, because the first element is the contemporaneous one
      dimnames(GFEVDFactorsNormalized)[[2]] <- AllFactorsLabels
      dimnames(GFEVDFactorsNormalized)[[3]] <- AllFactorsLabels

      dimnames(GFEVDYieldsNormalized)[[1]] <- 1:(GFEVDhoriz) # We subtract 1, because the first element is the contemporaneous one
      dimnames(GFEVDYieldsNormalized)[[2]] <- AllFactorsLabels
      dimnames(GFEVDYieldsNormalized)[[3]] <- YieldsLabel

      GFEVDoutputs[[ModelTypeSet[j]]][[Economies[i]]] <- list(GFEVDFactorsNormalized,GFEVDYieldsNormalized)
      names(GFEVDoutputs[[ModelTypeSet[j]]][[Economies[i]]]) <- c("Factors","Yields")


    }
  }

  # Clean empty lists
  GFEVDoutputs  <- GFEVDoutputs[unlist(lapply(GFEVDoutputs, length) != 0)]

  names(GFEVDoutputs) <- ModelTypeSet[idxIndividual]

  return(GFEVDoutputs)

}


######################################################################################################
####################### 7) Risk Premia Decomposition #################################################
######################################################################################################
#' Decomposition of yields into the average of expected future short-term interest rate and risk premia for "joint Q" models
#'
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param InputsForOutputs list conataining the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs,
#'                        and risk premia decomposition
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'
#'@keywords internal

TermPremiaDecompSep <- function(ModelPara, FactorLabels, ModelType, InputsForOutputs, Economies){



  # 1) Redefine inputs of interest
  # a) General Inputs
  T <- ncol(ModelPara[[ModelType]][[1]]$inputs$AllFactors)
  mat <- ModelPara[[ModelType]][[1]]$inputs$mat
  J <- length(mat)

  C <- length(Economies)
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N

  UnitYields <- InputsForOutputs$UnitMatYields
  if (UnitYields== "Month"){
    k <- 12
    YLab <- "M"
  } else  if (UnitYields== "Year"){
    k <- 1
    YLab <- "Y"
  }
  matAdjUnit <- mat*k

  WishFP <- InputsForOutputs$ForwardPremia

  matMIN <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[1]
  matMAX <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[2]


  #b) Specific Inputs
  avexp <- list()

  rho0_PP <- vector(mode = "list", length = C)
  rho1_PP <- vector(mode = "list", length = C)
  names(rho0_PP) <- Economies
  names(rho1_PP) <- Economies

  for (i in 1:C){
    ZZ <- ModelPara[[ModelType]][[Economies[i]]]$inputs$AllFactors
    BnX <- ModelPara[[ModelType]][[Economies[i]]]$rot$X$B
    AnX <- ModelPara[[ModelType]][[Economies[i]]]$rot$X$A
    Wpca <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Wpca
    K0Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K0Z
    K1Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K1Z
    SSZ <- ModelPara[[ModelType]][[Economies[i]]]$ests$SSZ

    r0 <- ModelPara[[ModelType]][[Economies[i]]]$ests$r0

    mat <- ModelPara[[ModelType]][[Economies[i]]]$inputs$mat
    dt <- ModelPara[[ModelType]][[Economies[i]]]$inputs$dt

    K <- nrow(K0Z)

    # 2) Expectations component
    # a) Extract spanned factors from the list of unspanned factors

    if (ModelType == "JPS original"){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]]) }
    if (ModelType == "JPS global" || ModelType == 'GVAR single'){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

    rownames(SSZ) <- AllLabels
    colnames(SSZ) <- AllLabels
    LabelSpannedCS <- c(FactorLabels$Tables[[Economies[i]]][-(1:M)])
    IdxSpanned <- match(LabelSpannedCS, AllLabels)

    # b) Compute the intercept and slope coefficients of the short rate expressed as a function of the spanned factors
    # By definition: r_t = r0 + rho1_X* X_t
    # But X_t = (W*Bx)^(-1)(P- WAx)
    # so r_t = rho0_PP + rho1_PP*P_t
    # where (i) rho0_PP = r0 - rho1_X*(W*BX)^(-1)W*AX and (ii) rho1_PP = rho1_X (W*BX)^(-1)
    rho1_X <- rep(1,N)

    rho0_PP[[i]] <- as.numeric((r0 - rho1_X%*%solve(Wpca%*%BnX,tol = 1e-50)%*%Wpca%*%AnX)/dt)
    rho1_PP[[i]] <- (rho1_X%*%solve(Wpca%*%BnX, tol = 1e-50))/dt


    # c) Compute expectations component
    ExpecCompLength <-  round((matAdjUnit/k)/dt)


    # c.1) Pure expected component
    # Per country
    avexpCS <- matrix(0, T, length(ExpecCompLength))
    rownames(avexpCS) <- colnames(ZZ)
    colnames(avexpCS) <- paste("RP_",ExpecCompLength, YLab, sep="")

    for (h in 1:length(ExpecCompLength)){ # per bond maturity

      for (t in 1:T){ # Per point in time

        g <- matrix(NA, K, ExpecCompLength[h])
        rownames(g) <- rownames(ZZ)

        g[ ,1] <- ZZ[, t]
        for (j in 2:ExpecCompLength[h]){g[ ,j] <- K0Z + K1Z%*%g[ ,j-1]} # Fitted P-dynamics

        g <- g[IdxSpanned, ] # extract relevant variables

        MaxExpec <- pmax(rho0_PP[[i]] + (rho1_PP[[i]]%*%g),0)
        avexpCS[t,h] <- mean(MaxExpec)
      }
    }

    avexp[[Economies[i]]] <- avexpCS*100
  }

  ############################################## FORWARD PREMIA PART ####################################
  if( WishFP == 1){

    avexpFP <- list()
    # c.1) Expected forward rate component
    for(i in 1:C){ # Per country
      ZZ <- ModelPara[[ModelType]][[Economies[i]]]$inputs$AllFactors
      K0Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K0Z
      K1Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K1Z

      avexpCS <- matrix(0, T, 1)
      rownames(avexpCS) <- colnames(ZZ)
      colnames(avexpCS) <- paste("FP_", matMIN, "-", matMAX, YLab, sep="")

      matMINAdj <- round((matMIN/k)/dt)
      matMAXAdj <- round((matMAX/k)/dt)
      for (t in 1:T){ # Per point in time

        g <- matrix(NA, K, matMAXAdj)
        rownames(g) <- rownames(ZZ)

        g[ ,1] <- ZZ[, t]
        for (j in 2:matMAXAdj){g[ ,j] <- K0Z + K1Z%*%g[ ,j-1]} # Fitted P-dynamics

        g <- g[IdxSpanned, matMINAdj:matMAXAdj] # extract relevant variables

        MaxExpec <- pmax(rho0_PP[[i]] + (rho1_PP[[i]]%*%g),0)
        avexpCS[t] <- mean(MaxExpec)
      }

      avexpFP[[Economies[i]]] <- avexpCS*100
    }
  }
  ##########################################################################################################

  # 3) Compute Term Premium
  YieldData <- list()
  TermPremium <- list()

  # a) Pure term Premia expected component
  for (i in 1:C){
    Y <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y
    YieldData[[Economies[i]]] <- t(Y)*100
    TermPremium[[Economies[i]]] <- YieldData[[Economies[i]]] - avexp[[Economies[i]]]
  }

  Output <- list(TermPremium, avexp)
  names(Output) <- c("Term Premia","Expected Component")

  ################################################## FORWARD PREMIA ####################################
  if( WishFP == 1){
    # d) Forward Premia

    IDXMatMIN <- match(matMIN,matAdjUnit)
    IDXMatMAX <- match(matMAX,matAdjUnit)
    IDXMatBoth <- c(IDXMatMIN,IDXMatMAX)

    FR <- list()
    ForwardPremium <- list()

    # If any of the data of the specific maturity were not used in the estimation, then compute the fitted value
    if( anyNA(IDXMatBoth)){

      MatBoth <- c(matMIN, matMAX)
      IdxNA <- which(is.na(IDXMatBoth))

      # Missing maturity (model-implied)
      MissingMat <-  MatBoth[IdxNA] # Maturity not available
      YieldMissing <- YieldsFitAllSep(MissingMat, ModelPara, FactorLabels, ModelType, Economies, YLab)

      for (i in 1:C){
        # If both matrutities are missing
        if (length(MissingMat) == 2){
          YieldMIN <- t(t(YieldMissing[[Economies[i]]][1,]*100))
          YieldMAX <- t(t(YieldMissing[[Economies[i]]][2,]*100))
        } else{
          # If only one maturity is missing
          # Available maturity
          IdxNotMissing0 <- IDXMatBoth[!is.na(IDXMatBoth)]
          IdxNotMissingCS <- IdxNotMissing0 + J*(i-1)
          YieldNotMissing <- t(t(Y[IdxNotMissingCS, ]))

          if (MissingMat ==1){
            YieldMIN <- YieldNotMissing*100
            YieldMAX <- t(YieldMissing[[Economies[i]]]*100)
          }else{
            YieldMIN <- t(YieldMissing[[Economies[i]]]*100)
            YieldMAX <- YieldNotMissing*100
          }
        }

        FR[[Economies[i]]] <- (matMAX*YieldMAX -  matMIN*YieldMIN)/(matMAX -  matMIN) # Fitted forward rate
        ForwardPremium[[Economies[i]]] <- FR[[Economies[i]]] - avexpFP[[Economies[i]]] # Forward Premia
        colnames(ForwardPremium[[Economies[i]]]) <- paste("FP_", matMIN, "-", matMAX, YLab, sep="")
        colnames(FR[[Economies[i]]]) <- paste("Mat", matMIN, "-", matMAX, YLab, sep="")
      }


      # Otherwise just use available data
    }else{

      for (i in 1:C){

        YieldMIN <- t(t(YieldData[[Economies[i]]][ , IDXMatMIN]))
        YieldMAX <- t(t(YieldData[[Economies[i]]][ , IDXMatMAX]))

        FR[[Economies[i]]] <- (matMAX*YieldMAX -  matMIN*YieldMIN)/(matMAX -  matMIN) # Fitted forward rate
        ForwardPremium[[Economies[i]]] <- FR[[Economies[i]]] - avexpFP[[Economies[i]]] # Forward Premia
        colnames(ForwardPremium[[Economies[i]]]) <- paste("FP_", matMIN, "-", matMAX, YLab, sep="")
        colnames(FR[[Economies[i]]]) <- paste("Mat", matMIN, "-", matMAX, YLab, sep="")
      }
    }

    OutputFP <- list(ForwardPremium, avexpFP, FR)
    names(OutputFP) <- c("Forward Premia","Expected Component", "Forward Rate")

    Output <- list(Output,OutputFP)
    names(Output) <- c("RiskPremia", "ForwardPremia")

    ###############################################################################################################
  }else {
    OutputFP <- NA
    Output <- list(Output,OutputFP)
    names(Output) <- c("RiskPremia", "ForwardPremia")
  }

  return(Output)
}




######################################################################################################
######################################################################################################
####################### OUTPUTS FOR MODELS IN WHICH THE ESTIMATION ###################################
########################       IS DONE FOR ALL COUNTRIES JOINTLY      ################################
######################################################################################################
######################################################################################################
#' Numerical outputs (variance explained, model fit, IRFs, GIRFs, FEVDs, GFEVDs and risk premia decomposition)
#' for "joint Q" models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param InputsForOutputs list conataining the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs
#'                        and risk premia decomposition
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@keywords internal


OutputConstructionJoint <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels,
                                    Economies){


  # Output summary
  # Total Variance Explained and Model Fit
  Total_Var_exp <- VarianceExplained(ModelType, ModelPara, FactorLabels, Economies)
  ModFit <- YieldsFit(ModelType, ModelPara, FactorLabels, Economies)

  # IRF and FEVD
  IRFout <- IRFandGIRF(ModelType, ModelPara, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)
  FEVDoutputs <- FEVDjoint(ModelType, ModelPara, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels, Economies)

  # GIRFS and GFEVD
  GFEVDoutputs <- GFEVDjoint(ModelType, ModelPara, InputsForOutputs[[ModelType]]$GFEVD$horiz, FactorLabels, Economies)

  # Risk Premia decomposition
  TermPremia <- TermPremiaDecompJoint(ModelPara, FactorLabels, ModelType, InputsForOutputs, Economies)


  # FOR JLL models: Non-orthogonalized IRFs
  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){

    # Generate the outputs with orthogonalized factros:
    FEVDOrtho <- FEVDjointOrthogoJLL(ModelType, ModelPara, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels,Economies)
    GFEVDOrtho <- GFEVDjointOrthoJLL(ModelType, ModelPara, InputsForOutputs[[ModelType]]$GFEVD$horiz, FactorLabels, Economies)


      # Merge the lists of orthogonalized and non-orthogonalized factors
       # FEVD
      FEVDoutputs[[ModelType]]$Factors <- list(FEVDoutputs[[ModelType]]$Factors, FEVDOrtho[[ModelType]]$Factors)
      FEVDoutputs[[ModelType]]$Yields <- list(FEVDoutputs[[ModelType]]$Yields, FEVDOrtho[[ModelType]]$Yields)
      names(FEVDoutputs[[ModelType]]$Factors) <- c("NonOrtho", "Ortho")
      names(FEVDoutputs[[ModelType]]$Yields) <- c("NonOrtho", "Ortho")
      # GFEVD
      GFEVDoutputs[[ModelType]]$Factors <- list(GFEVDoutputs[[ModelType]]$Factors, GFEVDOrtho[[ModelType]]$Factors)
      GFEVDoutputs[[ModelType]]$Yields <- list(GFEVDoutputs[[ModelType]]$Yields, GFEVDOrtho[[ModelType]]$Yields)
      names(GFEVDoutputs[[ModelType]]$Factors) <- c("NonOrtho", "Ortho")
      names(GFEVDoutputs[[ModelType]]$Yields) <- c("NonOrtho", "Ortho")

    }

  NumericalOutputs <- list(Total_Var_exp, ModFit, IRFout$IRFs, FEVDoutputs, IRFout$GIRFs, GFEVDoutputs, TermPremia)
  names(NumericalOutputs) <- c("PC var explained","Fit","IRF", 'FEVD', "GIRF", "GFEVD", "TermPremiaDecomp")

  return(NumericalOutputs)

}


#########################################################################################################
################################### 4) FEVD #############################################################
#########################################################################################################
#' FEVDs for "joint Q" models
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param FEVDhoriz single numerical vector conataining the desired horizon of analysis for the FEVDs
#'@param FactorLabels   string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#
#'
#'@details
#' Structural shocks are identified via Cholesky decomposition
#'
#'@keywords internal


FEVDjoint <- function(ModelType, ModelPara, FEVDhoriz, FactorLabels, Economies){


  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic)-N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G


  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR single") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  FEVDoutputs <- list()

  for (j in  jointModelIdx){
    J <- length(ModelPara[[ModelTypeSet[[j]]]]$inputs$mat)

    G0 <- ModelPara[[ModelTypeSet[[j]]]]$ests$Gy.0
    Sigma_y <-  ModelPara[[ModelTypeSet[[j]]]]$ests$SSZ
    F1 <- ModelPara[[ModelTypeSet[[j]]]]$ests$K1Z


    # 1) Dynamic multipliers
    Ry.h <- array(NA, c(K,K,FEVDhoriz))
    Ry.h[, ,1] <- diag(K) # dynamic multiplier at t=0

    for (i in 2:FEVDhoriz) {
      Ry.h[, ,i] <- F1%*%Ry.h[, ,i-1]
    }

    # 2) Initialization
    vslct <- diag(K)
    eslct <- diag(K)

    # 2.1) Minor preliminary work
    invG <- diag(nrow(G0))/G0
    invG[!is.finite(invG)] <- 0
    invGSigmau <- solve(G0)%*%Sigma_y

    # Choleski term
    if ( any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
      P <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$Sigmas$Sigma_Y
    }else{ P <- t(chol(Sigma_y))}

    scale <- 1

    # 2.2) Factor loadings preparation
    BSpanned <- ModelPara[[ModelTypeSet[[j]]]]$rot$P$B
    B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

    # 3) FEVD
    # 3.1) Factors
    FEVDresFactors <- array(NA, c(nrow = K, ncol=K, FEVDhoriz))
    num <- matrix(0, nrow =K, ncol=K)
    den <- rep(0, times = K)

    for (l in 1:FEVDhoriz){
      acc1 <- (eslct%*%Ry.h[,,l]%*%P%*%vslct)^2
      num <- num + acc1
      acc2 <- diag(eslct%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%eslct)
      den<- den + acc2
      FEVDresFactors[ ,,l] <- scale*num/den
    }

    FEVDFactors <- aperm(FEVDresFactors, c(3,2,1))


    # 3.2) Yields
    eslctCJ <- diag(C*J)
    vslctCJ <- diag(K)


    FEVDresYields <- array(NA, c(nrow = C*J, ncol=K, FEVDhoriz))
    num <- matrix(0, nrow =C*J, ncol=K)
    den <- matrix(0, nrow =C*J, ncol=K)

    for (l in 1:FEVDhoriz){
      acc1 <- (eslctCJ%*%B%*%Ry.h[,,l]%*%P%*%vslctCJ)^2
      num <- num + acc1
      acc2 <- diag(eslctCJ%*%B%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%t(B)%*%eslctCJ)
      den<- den + acc2
      FEVDresYields[ ,,l] <- scale*num/den
    }


    FEVDYields <- aperm(FEVDresYields, c(3, 2, 1))

    # 4) Prepare labels
    labelsFEVDFactors <- c(FactorLabels$Global,FactorLabels$Tables$AllCountries)
    YieldsLabel<- rownames(ModelPara[[ModelTypeSet[[j]]]]$inputs$Y)

    dimnames(FEVDFactors)[[1]] <- 1:FEVDhoriz # We subtract 1, because the first element is the contemporaneous one
    dimnames(FEVDFactors)[[2]] <- labelsFEVDFactors
    dimnames(FEVDFactors)[[3]] <- labelsFEVDFactors


    dimnames(FEVDYields)[[1]] <- 1:FEVDhoriz # We subtract 1, because the first element is the contemporaneous one
    dimnames(FEVDYields)[[2]] <- labelsFEVDFactors
    dimnames(FEVDYields)[[3]] <- YieldsLabel


    FEVDoutputsAllCountries <- list(FEVDFactors,FEVDYields)
    names(FEVDoutputsAllCountries) <- c("Factors","Yields")

    FEVDoutputs[[j-L]] <- FEVDoutputsAllCountries


  }
  # Clean empty lists
  FEVDoutputs  <- FEVDoutputs[unlist(lapply(FEVDoutputs, length) != 0)]
  names(FEVDoutputs) <- ModelTypeSet[jointModelIdx]


  return(FEVDoutputs)

}

#########################################################################################################
################################### 6) GFEVD ############################################################
#########################################################################################################
#' GFEVDs for "joint Q" models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param GFEVDhoriz single numerical vector conataining the desired horizon of analysis for the GFEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@references
#' \itemize{
#' \item This function is a modified and extended version of the "fevd" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#'
#'@keywords internal


GFEVDjoint <- function(ModelType, ModelPara, GFEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic)-N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G


  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR single") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  GFEVDoutputs <- list()

  for (j in  jointModelIdx){

    J <- numel(ModelPara[[ModelTypeSet[[j]]]]$inputs$mat)
    Y <- ModelPara[[ModelTypeSet[[j]]]]$inputs$Y
    ZZ <- ModelPara[[ModelTypeSet[[j]]]]$inputs$AllFactors

    G0 <- ModelPara[[ModelTypeSet[[j]]]]$ests$Gy.0
    Sigma_y <- ModelPara[[ModelTypeSet[[j]]]]$ests$SSZ
    F1 <- ModelPara[[ModelTypeSet[[j]]]]$ests$K1Z


    # 1) Dynamic multipliers
    Ry.h <- array(NA, c(K,K,GFEVDhoriz))
    Ry.h[, ,1] <- diag(K) # dynamic multiplier at t=0

    for (i in 2:GFEVDhoriz) {
      Ry.h[, ,i] <- F1%*%Ry.h[, ,i-1]
    }

    # 2) Initialization/ Minor preliminary work
    GFEVDresFac <- array(NA, c(nrow = K, ncol=GFEVDhoriz,K))
    vslct <- diag(K)
    eslct <- diag(K)

    invG <- diag(nrow(G0))/G0
    invG[!is.finite(invG)] <- 0
    invGSigmau <- solve(G0)%*%Sigma_y

    scale <- 1/diag(Sigma_y)

    # 3) GFEVD
    # 3.1) Factors
    for(i in 1:K){
      n<-1
      num <- matrix(0, nrow =K, ncol=GFEVDhoriz)
      den <- matrix(0, nrow =K, ncol=GFEVDhoriz)
      while (n <= GFEVDhoriz){
        for (l in 1:n){
          acc1 <- t((eslct[,i]%*%Ry.h[,,l]%*%invGSigmau%*%vslct)^2)
          num[,n] <- num[ ,n] + acc1
          acc2 <- eslct[,i]%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%eslct[,i]
          den[,n]<- den[,n] + matrix(1, nrow=K)%*%acc2
        }
        GFEVDresFac[ ,n,i] <- t(t(scale*num[,n]))/den[,n]
        n <- n+1
      }
    }


    GFEVDFactors <- aperm(GFEVDresFac, c(2,1,3)) # Non-normalized GFEVD (i.e. rows need not sum up to 1)

    #  Normalization of the GFEVD for the factors
    # (Make sure that the sum of the errors equal to one in each period)
    DEM <- array(NA, c(nrow = GFEVDhoriz, ncol=1, K))
    GFEVDFactorsNormalized <- array(NA, c(nrow = GFEVDhoriz, ncol=K, K))

    for (h in 1:K){
      for (n in 1:GFEVDhoriz){
        DEM[n, 1, h] <- sum(GFEVDFactors[n,,h])
        GFEVDFactorsNormalized[n,,h] <- GFEVDFactors[n,,h]/DEM[n,,h]
      }
    }

    # 3.2) Yields
    # Get the full B
    BSpanned <- ModelPara[[ModelTypeSet[[j]]]]$rot$P$B
    B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)


    # Initialization
    GFEVDresYie <- array(NA, c(nrow = C*J, ncol=K, GFEVDhoriz))
    vslctYie <- diag(K)
    eslctYie <- diag(C*J)


    num <- matrix(0, nrow =C*J, ncol=K)
    den <- matrix(0, nrow =C*J, ncol=K)

    for (l in 1:GFEVDhoriz){
      acc1 <- (eslctYie%*%B%*%Ry.h[,,l]%*%invGSigmau%*%vslctYie)^2
      num <- num + acc1
      acc2 <- diag(eslctYie%*%B%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%t(B)%*%eslctYie)
      den <- den + acc2
      for (q in 1:K){
        GFEVDresYie[ ,q,l] <- scale[q]*(num/den)[,q] # note: unlike the GFEVD of the factors, note that the "scale" variable is now at the acc1
      }
    }


    GFEVDYields <- aperm(GFEVDresYie, c(3,2,1))

    #  Normalization of the GFEVD for the factors
    # (Make sure that the sum of the errors equal to one in each period)
    DEM <- array(NA, c(nrow = GFEVDhoriz, ncol=1, C*J))
    GFEVDYieldsNormalized <- array(NA, c(nrow = GFEVDhoriz, ncol=K, C*J))

    for (h in 1:(C*J)){
      for (n in 1:GFEVDhoriz){
        DEM[n, 1, h] <- sum(GFEVDYields[n,,h])
        GFEVDYieldsNormalized[n,,h] <- GFEVDYields[n,,h]/DEM[n,,h]
      }
    }

    # 4) Adjust labels:
    labelsGFEVDFactors <- c(FactorLabels$Global,FactorLabels$Tables$AllCountries)
    YieldsLabel <- rownames(ModelPara[[ModelTypeSet[[j]]]]$inputs$Y)


    dimnames(GFEVDFactorsNormalized)[[1]] <- 1:(GFEVDhoriz) # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error
    dimnames(GFEVDFactorsNormalized)[[2]] <- labelsGFEVDFactors
    dimnames(GFEVDFactorsNormalized)[[3]] <- labelsGFEVDFactors

    dimnames(GFEVDYieldsNormalized)[[1]] <- 1:(GFEVDhoriz) # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error
    dimnames(GFEVDYieldsNormalized)[[2]] <- labelsGFEVDFactors
    dimnames(GFEVDYieldsNormalized)[[3]] <- YieldsLabel

    GFEVDoutputsAllCountries <- list(GFEVDFactorsNormalized,GFEVDYieldsNormalized)
    names(GFEVDoutputsAllCountries) <- c("Factors","Yields")

    GFEVDoutputs[[j-L]] <- GFEVDoutputsAllCountries

  }
  # Clean empty lists
  GFEVDoutputs  <- GFEVDoutputs[unlist(lapply(GFEVDoutputs, length) != 0)]

  names(GFEVDoutputs) <- ModelTypeSet[jointModelIdx]


  return(GFEVDoutputs)

}

######################################################################################################
####################### 7) Risk Premia Decomposition #################################################
######################################################################################################
#' Decomposition of yields into the average of expected future short-term interest rate and risk premia for "joint Q" models
#'
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param InputsForOutputs list conataining the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs and risk premia decomposition
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@keywords internal

TermPremiaDecompJoint <- function(ModelPara, FactorLabels, ModelType, InputsForOutputs, Economies){


  # 1) Redefine inputs of interest
  ZZ <- ModelPara[[ModelType]]$inputs$AllFactors
  Y <- ModelPara[[ModelType]]$inputs$Y


  BnX <- ModelPara[[ModelType]]$rot$X$B
  AnX <- ModelPara[[ModelType]]$rot$X$A
  Wpca <- ModelPara[[ModelType]]$inputs$Wpca
  K0Z <- ModelPara[[ModelType]]$ests$K0Z
  K1Z <- ModelPara[[ModelType]]$ests$K1Z
  r0 <- ModelPara[[ModelType]]$ests$r0

  mat <- ModelPara[[ModelType]]$inputs$mat
  dt <- ModelPara[[ModelType]]$inputs$dt

  K <- nrow(ZZ)
  T <- ncol(ZZ)
  J <-length(mat)
  C <- length(Economies)
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- (K-G)/C - N

  UnitYields <- InputsForOutputs$UnitMatYields
  if (UnitYields== "Month"){
    k <- 12
    YLab <- "M"
  }  else if (UnitYields== "Year"){
    k <- 1
    YLab <- "Y"
  }
  matAdjUnit <- mat*k

  WishFP <- InputsForOutputs$ForwardPremia

  matMIN <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[1]
  matMAX <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[2]

  # 2) Expectations component
  # a) Extract spanned factors from the list of unspanned factors
  IdxSpanned <- c()
  idxSpa0 <- G + M
  for (j in 1:C){
    idxSpa1 <- idxSpa0 + N

    if (j ==1){ IdxSpanned <- (idxSpa0+1):idxSpa1
    }  else{
      IdxSpanned <- c(IdxSpanned, (idxSpa0+1):idxSpa1)
    }
    idxSpa0 <- idxSpa1 + M
  }


  # b) Compute the intercept and slope coefficients of the short rate expressed as a function of the spanned factors
  # By definition: r_t = r0 + rho1_X* X_t
  # But X_t = (W*Bx)^(-1)(P- WAx)
  # so r_t = rho0_PP + rho1_PP*P_t
  # where (i) rho0_PP = r0 - rho1_X*(W*BX)^(-1)W*AX and (ii) rho1_PP = rho1_X (W*BX)^(-1)
  rho1_X_CS <- rep(1,N)
  rho1_X <- matrix(0, nrow=C, ncol=N*C)
  idx0 <-0

  for (j in 1:C){
    idx1<- idx0+N
    rho1_X[j,(idx0+1):idx1] <- rho1_X_CS
    idx0<- idx1
  }


  rho0_PP <- (r0 - rho1_X%*%solve(Wpca%*%BnX, tol = 1e-50)%*%Wpca%*%AnX)/dt
  rho1_PP <- (rho1_X%*%solve(Wpca%*%BnX, tol = 1e-50))/dt

  # c) Compute expectations component
  ExpecCompLength <-  round((matAdjUnit/k)/dt)


  avexp <- list()
  # c.1) Pure term Premia expected component
  for(i in 1:C){ # Per country
    avexpCS <- matrix(0, T, length(ExpecCompLength))
    rownames(avexpCS) <- colnames(ZZ)
    colnames(avexpCS) <- paste("RP_",ExpecCompLength, YLab, sep="")

    for (h in 1:length(ExpecCompLength)){ # per bond maturity

      for (t in 1:T){ # Per point in time

        g <- matrix(NA, K, ExpecCompLength[h])
        rownames(g) <- rownames(ZZ)

        g[ ,1] <- ZZ[, t]
        for (j in 2:ExpecCompLength[h]){g[ ,j] <- K0Z + K1Z%*%g[ ,j-1]} # Fitted P-dynamics

        g <- g[IdxSpanned, ] # extract relevant variables

        MaxExpec <- pmax(rho0_PP[i] + (rho1_PP[i, ]%*%g),0)
        avexpCS[t,h] <- mean(MaxExpec)
      }
    }

    avexp[[Economies[i]]] <- avexpCS*100
  }

  ############################################## FORWARD PREMIA ####################################
  if( WishFP == 1){

    avexpFP <- list()
    # c.1) Forward Premia expected component
    for(i in 1:C){ # Per country
      avexpCS <- matrix(0, T, 1)
      rownames(avexpCS) <- colnames(ZZ)
      colnames(avexpCS) <- paste("FP_", matMIN, "-", matMAX, YLab, sep="")

      matMINAdj <- round((matMIN/k)/dt)
      matMAXAdj <- round((matMAX/k)/dt)
      for (t in 1:T){ # Per point in time

        g <- matrix(NA, K, matMAXAdj)
        rownames(g) <- rownames(ZZ)

        g[ ,1] <- ZZ[, t]
        for (j in 2:matMAXAdj){g[ ,j] <- K0Z + K1Z%*%g[ ,j-1]} # Fitted P-dynamics

        g <- g[IdxSpanned, matMINAdj:matMAXAdj] # extract relevant variables

        MaxExpec <- pmax(rho0_PP[i] + (rho1_PP[i, ]%*%g),0)
        avexpCS[t] <- mean(MaxExpec)
      }

      avexpFP[[Economies[i]]] <- avexpCS*100
    }
  }
  ##########################################################################################################

  # 3) Compute Term Premium
  YieldData <- list()
  TermPremium <- list()

  # a) Pure term Premia expected component

  for (i in 1:C){
    IdxRP <- (1:J) +J*(i-1)
    YieldData[[Economies[i]]] <- t(Y[IdxRP, ]*100)
    TermPremium[[Economies[i]]] <- YieldData[[Economies[i]]] - avexp[[Economies[i]]]
  }

  Output <- list(TermPremium, avexp)
  names(Output) <- c("Term Premia","Expected Component")

  ################################################## FORWARD PREMIA ####################################
  if( WishFP == 1){
    # d) Forward Premia

    IDXMatMIN <- match(matMIN,matAdjUnit)
    IDXMatMAX <- match(matMAX,matAdjUnit)
    IDXMatBoth <- c(IDXMatMIN,IDXMatMAX)

    FR <- list()

    # If any of the data of the specific maturity were not used in the estimation, then compute the fitted value
    if( anyNA(IDXMatBoth)){

      MatBoth <- c(matMIN, matMAX)
      IdxNA <- which(is.na(IDXMatBoth))

      # Missing maturity (model-implied)
      MissingMat <-  MatBoth[IdxNA] # Maturity not avaiable
      YieldMissing <- YieldsFitAllJoint(MissingMat, ModelPara, FactorLabels, ModelType, Economies, YLab)

      for (i in 1:C){
        # If both matrutities are missing
        if (length(MissingMat) == 2){
          YieldMIN <- t(t(YieldMissing[[Economies[i]]][1,]*100))
          YieldMAX <- t(t(YieldMissing[[Economies[i]]][2,]*100))
        } else{
          # If only one maturity is missing
          # Available maturity
          IdxNotMissing0 <- IDXMatBoth[!is.na(IDXMatBoth)]
          IdxNotMissingCS <- IdxNotMissing0 + J*(i-1)
          YieldNotMissing <- t(t(Y[IdxNotMissingCS, ]))

          if (MissingMat ==1){
            YieldMIN <- YieldNotMissing*100
            YieldMAX <- t(YieldMissing[[Economies[i]]]*100)
          }else{
            YieldMIN <- t(YieldMissing[[Economies[i]]]*100)
            YieldMAX <- YieldNotMissing*100
          }
        }

        FR[[Economies[i]]] <- (matMAX*YieldMAX -  matMIN*YieldMIN)/(matMAX -  matMIN) # Fitted forward rate
        TermPremium[[Economies[i]]] <- FR[[Economies[i]]] - avexpFP[[Economies[i]]] # Forward Premia
        colnames(TermPremium[[Economies[i]]]) <- paste("FP_", matMIN, "-", matMAX, YLab, sep="")
        colnames(FR[[Economies[i]]]) <- paste("Mat", matMIN, "-", matMAX, YLab, sep="")
      }


      # Otherwise just use available data
    }else{

      for (i in 1:C){
        IdxMinCS <- IDXMatMIN + J*(i-1)
        IdxMaxCS <- IDXMatMAX + J*(i-1)

        YieldMIN <- t(t(Y[IdxMinCS, ]*100))
        YieldMAX <- t(t(Y[IdxMaxCS, ]*100))

        FR[[Economies[i]]] <- (matMAX*YieldMAX -  matMIN*YieldMIN)/(matMAX -  matMIN) # Fitted forward rate
        TermPremium[[Economies[i]]] <- FR[[Economies[i]]] - avexpFP[[Economies[i]]] # Forward Premia
        colnames(TermPremium[[Economies[i]]]) <- paste("FP_", matMIN, "-", matMAX, YLab, sep="")
        colnames(FR[[Economies[i]]]) <- paste("Mat", matMIN, "-", matMAX, YLab, sep="")
      }
    }

    OutputFP <- list(TermPremium, avexpFP, FR)
    names(OutputFP) <- c("Forward Premia","Expected Component", "Forward Rate")

    Output <- list(Output,OutputFP)
    names(Output) <- c("RiskPremia", "ForwardPremia")

    ###############################################################################################################
  }else {
    OutputFP <- NA
    Output <- list(Output,OutputFP)
    names(Output) <- c("RiskPremia", "ForwardPremia")
  }

  return(Output)
}


#########################################################################################################
################################### 9) FEVD with orthogonalized factors #################################
#########################################################################################################
#' Orthogonalized FEVDs for JLL models
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param FEVDhoriz single numerical vector conataining the desired horizon of analysis for the FEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'@details
#' Structural shocks are identified via Cholesky decomposition
#'
#'@keywords internal


FEVDjointOrthogoJLL <- function(ModelType, ModelPara, FEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic)-N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G

  jointModelIdxJLL <- idxWishModels[idxWishModels >= which(ModelTypeSet == "JLL original") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  FEVDoutputs <- list()

  for (j in  jointModelIdxJLL){
    J <- length(ModelPara[[ModelTypeSet[[j]]]]$inputs$mat)

    G0 <- ModelPara[[ModelTypeSet[[j]]]]$ests$Gy.0
    Sigma_y <-  ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$Sigmas$VarCov_Ortho
    F1e <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$k1_e
    PI <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$PI

    # 1) Dynamic multipliers
    Ry.h <- array(NA, c(K,K,FEVDhoriz))
    Ry.h[, ,1] <- diag(K) # dynamic multiplier at t=0

    for (i in 2:FEVDhoriz) {
      Ry.h[, ,i] <- F1e%*%Ry.h[, ,i-1]
    }

    # 2) Initialization
    vslct <- diag(K)
    eslct <- diag(K)

    # 2.1) Minor preliminary work
    invG <- diag(nrow(G0))/G0
    invG[!is.finite(invG)] <- 0
    invGSigmau <- solve(G0)%*%Sigma_y

    # Choleski term
    P <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$Sigmas$Sigma_Ye

    scale <- 1

    # 2.2) Factor loadings preparation
    BSpanned <- ModelPara[[ModelTypeSet[[j]]]]$rot$P$B
    B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

    # 3) FEVD
    # 3.1) Factors
    FEVDresFactors <- array(NA, c(nrow = K, ncol=K, FEVDhoriz))
    num <- matrix(0, nrow =K, ncol=K)
    den <- rep(0, times = K)

    for (l in 1:FEVDhoriz){
      acc1 <- (eslct%*%Ry.h[,,l]%*%P%*%vslct)^2
      num <- num + acc1
      acc2 <- diag(eslct%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%eslct)
      den<- den + acc2
      FEVDresFactors[ ,,l] <- scale*num/den
    }

    FEVDFactors <- aperm(FEVDresFactors, c(3,2,1))

    # 3.2) Yields
    eslctCJ <- diag(C*J)
    vslctCJ <- diag(K)


    FEVDresYields <- array(NA, c(nrow = C*J, ncol=K, FEVDhoriz))
    num <- matrix(0, nrow =C*J, ncol=K)
    den <- matrix(0, nrow =C*J, ncol=K)

    for (l in 1:FEVDhoriz){
      acc1 <- (eslctCJ%*%B%*%PI%*%Ry.h[,,l]%*%P%*%vslctCJ)^2
      num <- num + acc1
      acc2 <- diag(eslctCJ%*%B%*%PI%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%t(PI)%*%t(B)%*%eslctCJ)
      den<- den + acc2
      FEVDresYields[ ,,l] <- scale*num/den
    }


    FEVDYields <- aperm(FEVDresYields, c(3, 2, 1))

    # 4) Prepare labels
    labelsFEVDFactors <- c(FactorLabels$Global,FactorLabels$Tables$AllCountries)
    YieldsLabel<- rownames(ModelPara[[ModelTypeSet[[j]]]]$inputs$Y)

    dimnames(FEVDFactors)[[1]] <- 1:(FEVDhoriz) # # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error
    dimnames(FEVDFactors)[[2]] <- labelsFEVDFactors
    dimnames(FEVDFactors)[[3]] <- labelsFEVDFactors


    dimnames(FEVDYields)[[1]] <- 1:(FEVDhoriz) # # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error
    dimnames(FEVDYields)[[2]] <- labelsFEVDFactors
    dimnames(FEVDYields)[[3]] <- YieldsLabel


    FEVDoutputsAllCountries <- list(FEVDFactors,FEVDYields)
    names(FEVDoutputsAllCountries) <- c("Factors","Yields")

    FEVDoutputs[[j-L]] <- FEVDoutputsAllCountries


  }
  # Clean empty lists
  FEVDoutputs  <- FEVDoutputs[unlist(lapply(FEVDoutputs, length) != 0)]

  names(FEVDoutputs) <- ModelTypeSet[jointModelIdxJLL]


  return(FEVDoutputs)

}

#########################################################################################################
################################### 11) GFEVD With orthogonalized factor s################################
#########################################################################################################
#' Orthogonalized GFEVDs for JLL models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param ModelPara List of model parameter estimates (See the "Optimization" function)
#'@param GFEVDhoriz single numerical vector conataining the desired horizon of analysis for the GFEVDs
#'@param FactorLabels  a string-list based which contains all the labels of all the variables present in the model
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@references
#' \itemize{
#' \item This function is a modified and extended version of the "fevd" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#'
#'@keywords internal


GFEVDjointOrthoJLL <- function(ModelType, ModelPara, GFEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic)-N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G


  jointModelIdxJLL <- idxWishModels[idxWishModels >= which(ModelTypeSet == "JLL original") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  GFEVDoutputs <- list()

  for (j in  jointModelIdxJLL){

    J <- numel(ModelPara[[ModelTypeSet[[j]]]]$inputs$mat)
    Y <- ModelPara[[ModelTypeSet[[j]]]]$inputs$Y
    ZZ <- ModelPara[[ModelTypeSet[[j]]]]$inputs$AllFactors

    G0 <- ModelPara[[ModelTypeSet[[j]]]]$ests$Gy.0
    Sigma_y <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$Sigmas$VarCov_Ortho
    F1e <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$k1_e
    PI <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$PI


    # 1) Dynamic multipliers
    Ry.h <- array(NA, c(K,K,GFEVDhoriz))
    Ry.h[, ,1] <- diag(K) # dynamic multiplier at t=0

    for (i in 2:GFEVDhoriz) {
      Ry.h[, ,i] <- F1e%*%Ry.h[, ,i-1]
    }

    # 2) Initialization/ Minor preliminary work
    GFEVDresFac <- array(NA, c(nrow = K, ncol=GFEVDhoriz,K))
    vslct <- diag(K)
    eslct <- diag(K)

    invG <- diag(nrow(G0))/G0
    invG[!is.finite(invG)] <- 0
    invGSigmau <- solve(G0)%*%Sigma_y

    scale <- 1/diag(Sigma_y)

    # 3) GFEVD
    # 3.1) Factors
    for(i in 1:K){
      n<-1
      num <- matrix(0, nrow =K, ncol=GFEVDhoriz)
      den <- matrix(0, nrow =K, ncol=GFEVDhoriz)
      while (n <= GFEVDhoriz){
        for (l in 1:n){
          acc1 <- t((eslct[,i]%*%PI%*%Ry.h[,,l]%*%invGSigmau%*%vslct)^2)
          num[,n] <- num[ ,n] + acc1
          acc2 <- eslct[,i]%*%PI%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%eslct[,i]
          den[,n]<- den[,n] + matrix(1, nrow=K)%*%acc2
        }
        GFEVDresFac[ ,n,i] <- t(t(scale*num[,n]))/den[,n]
        n <- n+1
      }
    }


    GFEVDFactors <- aperm(GFEVDresFac, c(2,1,3)) # Non-normalized GFEVD (i.e. rows need not sum up to 1)

    #  Normalization of the GFEVD for the factors
    # (Make sure that the sum of the errors equal to one in each period)
    DEM <- array(NA, c(nrow = GFEVDhoriz, ncol=1, K))
    GFEVDFactorsNormalized <- array(NA, c(nrow = GFEVDhoriz, ncol=K, K))

    for (h in 1:K){
      for (n in 1:GFEVDhoriz){
        DEM[n, 1, h] <- sum(GFEVDFactors[n,,h])
        GFEVDFactorsNormalized[n,,h] <- GFEVDFactors[n,,h]/DEM[n,,h]
      }
    }

    # 3.2) Yields
    # Get the full B
    BSpanned <- ModelPara[[ModelTypeSet[[j]]]]$rot$P$B
    B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)




    # Initialization
    GFEVDresYie <- array(NA, c(nrow = C*J, ncol=K, GFEVDhoriz))
    vslctYie <- diag(K)
    eslctYie <- diag(C*J)

    num <- matrix(0, nrow =C*J, ncol=K)
    den <- matrix(0, nrow =C*J, ncol=K)

    for (l in 1:GFEVDhoriz){
      acc1 <- (eslctYie%*%B%*%PI%*%Ry.h[,,l]%*%invGSigmau%*%vslctYie)^2
      num <- num + acc1
      acc2 <- diag(eslctYie%*%B%*%PI%*%Ry.h[,,l]%*%invGSigmau%*%t(invG)%*%t(Ry.h[,,l])%*%t(PI)%*%t(B)%*%eslctYie)
      den <- den + acc2
      for (q in 1:K){
        GFEVDresYie[ ,q,l] <- scale[q]*(num/den)[,q] # note: unlike the GFEVD of the factors, note that the "scale" variable is now at the acc1
      }
    }


    GFEVDYields <- aperm(GFEVDresYie, c(3,2,1))

    #  Normalization of the GFEVD for the factors
    # (Make sure that the sum of the errors equal to one in each period)
    DEM <- array(NA, c(nrow = GFEVDhoriz, ncol=1, C*J))
    GFEVDYieldsNormalized <- array(NA, c(nrow = GFEVDhoriz, ncol=K, C*J))

    for (h in 1:(C*J)){
      for (n in 1:GFEVDhoriz){
        DEM[n, 1, h] <- sum(GFEVDYields[n,,h])
        GFEVDYieldsNormalized[n,,h] <- GFEVDYields[n,,h]/DEM[n,,h]
      }
    }

    # 4) Adjust labels:
    labelsGFEVDFactors <- c(FactorLabels$Global,FactorLabels$Tables$AllCountries)
    YieldsLabel <- rownames(ModelPara[[ModelTypeSet[[j]]]]$inputs$Y)


    dimnames(GFEVDFactorsNormalized)[[1]] <- 1:(GFEVDhoriz) # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error
    dimnames(GFEVDFactorsNormalized)[[2]] <- labelsGFEVDFactors
    dimnames(GFEVDFactorsNormalized)[[3]] <- labelsGFEVDFactors

    dimnames(GFEVDYieldsNormalized)[[1]] <- 1:(GFEVDhoriz) # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error
    dimnames(GFEVDYieldsNormalized)[[2]] <- labelsGFEVDFactors
    dimnames(GFEVDYieldsNormalized)[[3]] <- YieldsLabel

    GFEVDoutputsAllCountries <- list(GFEVDFactorsNormalized,GFEVDYieldsNormalized)
    names(GFEVDoutputsAllCountries) <- c("Factors","Yields")

    GFEVDoutputs[[j-L]] <- GFEVDoutputsAllCountries

  }
  # Clean empty lists
  GFEVDoutputs  <- GFEVDoutputs[unlist(lapply(GFEVDoutputs, length) != 0)]

  names(GFEVDoutputs) <- ModelTypeSet[jointModelIdxJLL]


  return(GFEVDoutputs)

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
  for (h in 1:TDim){   Yieldfit[,h] <- ALoad + BLoad%*%Spa_TS[,h] }
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
    tempYields[ , , r]      <- BLoad%*%AdjTerm%*%A1h%*%S
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
#'@param G0.y Contemporaneous matrix
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
#' \item This function is a ppartially based on the version of the "irf" function from
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

  #3.2) Yields
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
if (any(ModelType == c("JPS global", 'GVAR single'))){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

  LabelSpannedCS <- c(FactorLabels$Tables[[Economies[i]]][-(1:M)])
  b <- match(LabelSpannedCS, AllLabels)

  PP <- ZZ[b , ]

  X <- matrix(NA, nrow=nrow(PP), ncol = ncol(PP))

  for (t in 1:T){    X[,t] <- solve(Wpca%*%BnX, tol = 1e-50)%*%(PP[,t] - Wpca%*%AnX)  }

  # Loadings of the other maturities
  MatTU <- mat/dt
  MatAll <- 1:max(MatTU)

  LoadingsLat <- A0N__BnAn(MatAll, K1XQ, ModelType, dX= NULL, r0, SSX, Economies)
  AnXAll <- LoadingsLat[[2]]/dt
  BnXAll <- LoadingsLat[[1]]/dt

  for(h in 1:length(MatInt)){

    IdxMatInt <- seq(MatInt[h], length(AnXAll), by = max(MatTU))

    AnXInt <- AnXAll[IdxMatInt]
    BnXInt <- BnXAll[IdxMatInt,]

    for (t in 1:T){
      FitLat[h,t] <- AnXInt + BnXInt%*%X[,t]
    }
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

  for (t in 1:T){
    X[,t] <- solve(Wpca%*%BnX, tol = 1e-50)%*%(PP[,t] - Wpca%*%AnX)
  }


  # Loadings of the other maturities
  MatTU <- mat/dt
  MatAll <- 1:max(MatTU)

  LoadingsLat <- A0N__BnAn(MatAll, K1XQ, ModelType, dX= NULL, r0, SSX, Economies)
  AnXAll <- LoadingsLat[[2]]/dt
  BnXAll <- LoadingsLat[[1]]/dt



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

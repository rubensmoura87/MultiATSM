#' Construct the model numerical outputs (model fit, IRFs, GIRFs, FEVDs, and GFEVDs)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param ModelPara List of model parameter estimates (See the "Optimization" function)
#'@param InputsForOutputs list conataining the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, and GFEVDs
#'@param FactorLabels  a string-list based which contains all the labels of all the variables present in the model
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@examples
#' # See examples in the vignette file of this package (Section 4).
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
#'}
#'
#'
#'@export



NumOutputs <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies){


  AllNumOutputs <- list()

  # If one chooseS models in which the estimation is done country-by-country
  if ( "JPS" %in% ModelType || "JPS jointP" %in% ModelType ||  "GVAR sepQ" %in% ModelType){
    NumOutSep <- OutputConstructionSep(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies)

  }
  # If one chooseS models in which the estimation is done jointly for all countries
  if ( "GVAR jointQ" %in% ModelType || "VAR jointQ" %in% ModelType || "JLL original" %in% ModelType
       || "JLL NoDomUnit" %in% ModelType || "JLL jointSigma" %in% ModelType){

    NumOutJoint <- OutputConstructionJoint(ModelType, ModelPara,InputsForOutputs, FactorLabels, Economies)


  }



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

#' Numerical outputs (variance explained, model fit, IRFs, GIRFs, FEVDs, and GFEVDs) for "sep Q" models
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (See the "Optimization" function)
#'@param InputsForOutputs list conataining the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, and GFEVDs
#'@param FactorLabels string-list based which contains all the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system



OutputConstructionSep <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies){


  # Output summary
  # Total Variance Explained and Model fit
  Total_Var_exp <- VarianceExplainedSep(ModelType, ModelPara, FactorLabels, Economies)
  ModFit <- YieldsFitsep(ModelType, ModelPara, FactorLabels, Economies)

  # IRF and FEVD
  IRFoutputs <- IRFsep(ModelType, ModelPara, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)
  FEVDoutputs <- FEVDsep(ModelType, ModelPara, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels, Economies)

  # GIRF and GFEVD
  GIRFoutputs <- GIRFSep(ModelType, ModelPara, InputsForOutputs[[ModelType]]$GIRF$horiz, FactorLabels, Economies)
  GFEVDoutputs <- GFEVDsep(ModelType, ModelPara, InputsForOutputs[[ModelType]]$GFEVD$horiz, FactorLabels, Economies)

  NumericalOutputs <- list(Total_Var_exp, ModFit, IRFoutputs, FEVDoutputs, GIRFoutputs, GFEVDoutputs)
  names(NumericalOutputs) <- c("PC var explained", "Fit", "IRF", "FEVD", "GIRF", "GFEVD")

  return(NumericalOutputs)

}

######################################################################################################
####################### 1) Total Variance Explained #########################################
######################################################################################################
#' Percentage explained by the spanned factors of the variations in the set of observed yields for "sep Q" models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara List of model parameter estimates (see the "Optimization" function)
#'@param FactorLabels   string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system

#'


VarianceExplainedSep<-function(ModelType, ModelPara, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)


  C <- length(Economies)
  N <- length(FactorLabels$Spanned)

  Total_Var_exp_per_model <- list()
  Total_Var_exp_All <- list()
  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR sepQ")] # Exclude all models in which the estimation is made jointly

  for(j in idxIndividual){
    for (i in 1:C){
      H <- eigen(stats::cov(t(ModelPara[[ModelTypeSet[j]]][[Economies[i]]]$inputs$Y)))$values
      percentages_explained <- cumsum(H)/sum(H)
      Total_Var_exp_per_model[[i]] <-percentages_explained[1:N]
    }
    names(Total_Var_exp_per_model) <- Economies
    Total_Var_exp_All[[j]] <- Total_Var_exp_per_model
  }
  Total_Var_exp_All  <- Total_Var_exp_All[unlist(lapply(Total_Var_exp_All, length) != 0)]

  names(Total_Var_exp_All) <- ModelTypeSet[idxIndividual]

  return(Total_Var_exp_All)
}


######################################################################################################
########################################### 2) Model Fit #############################################
######################################################################################################
#' Computes two measures of model fit for bond yields
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


YieldsFitsep <- function(ModelType, ModelPara, FactorLabels, Economies){

  # 1) Redefining variables

  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)

  Output <- list()

  for (i in 1:C){
    mat <- ModelPara[[ModelType]][[Economies[i]]]$inputs$mat
    Z <- ModelPara[[ModelType]][[Economies[i]]]$inputs$AllFactors
    YieldData <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y

    J <- length(mat)
    T <- ncol(Z)

    Afull <- ModelPara[[ModelType]][[Economies[i]]]$rot$P$A
    Bspanned <- ModelPara[[ModelType]][[Economies[i]]]$rot$P$B
    K0Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K0Z
    K1Z <- ModelPara[[ModelType]][[Economies[i]]]$ests$K1Z

    # 2) MODEL FIT MEASURES
    # a) Model Fit (Yields)
    # Extract spanned factors from the list of unspanned factors
    if (ModelType == "JPS"){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]]) }
    if (ModelType == "JPS jointP" || ModelType == 'GVAR sepQ'){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

    LabelSpannedCS <- c(FactorLabels$Tables[[Economies[i]]][-(1:M)])
    IdxSpanned <- match(LabelSpannedCS, AllLabels)

    P <- Z[IdxSpanned,] # Set of spanned factors

    # Compute model fit
    Yieldfit<- matrix(NA, nrow= J, ncol =T)
    for (h in 1:T){
      Yieldfit[,h] <- Afull + Bspanned%*%P[,h]
    }
    rownames(Yieldfit) <- rownames(YieldData)
    colnames(Yieldfit) <- colnames(YieldData)


    # b) Model-implied Yields
    Bfull <- BUnspannedAdapSep(G,M, ModelPara, Economies, Economies[i], ModelType)

    row.names(Bfull) <- rownames(YieldData)
    colnames(Bfull) <- rownames(Z)

    YieldModelImplied<- matrix(NA, nrow= J, ncol =T)
    for (h in 2:T){ #  first observation is discarded
      YieldModelImplied[,h] <- Afull + Bfull%*%(K0Z + K1Z%*% Z[,h-1])
    }
    rownames(YieldModelImplied) <- rownames(YieldData)
    colnames(YieldModelImplied) <- colnames(YieldData)

    # 3) Prepare outputs
    fits<- list(Yieldfit, YieldModelImplied)
    names(fits) <- c("Yield Fit","Yield Model Implied")

    Output[[ModelType]][[Economies[i]]] <-  fits
  }


  return(Output)
}

######################################################################################################
########################################### 2) IRFs ##################################################
######################################################################################################
#' IRFs for "sep Q" models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (See the "Optimization" function)
#'@param IRFhoriz single numerical vector conataining the desired horizon of analysis for the IRFs
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@details
#' Structural shocks are identified via Cholesky decomposition


IRFsep <- function(ModelType, ModelPara, IRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  # Pre-allocation
  IRFoutputs <- list()
  IRFCS <- list()

  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR sepQ")] # Exclude all models in which the estimation is made jointly

  for(j in idxIndividual){
    for (i in 1:C){
      J <- numel(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$mat)
      K <- nrow(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$AllFactors)


      # Generate factor labels depending on the models to be estimated
      if ( ModelTypeSet[j] == "JPS") {AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]])}
      if ( ModelTypeSet[j] != "JPS") {AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries) }

      # Summarize inputs for the IRFs
      SIGMA <- ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$SSZ # KxK (variance-covariance matrix)
      A0 <- ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$K0Z # Kx1 (matrix of intercepts)
      A1 <- ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$K1Z # KxK (feedback matrix)


      B <- BUnspannedAdapSep(G,M, ModelPara, Economies, Economies[i], ModelTypeSet[j])


      # Initialization of IRFs of interest
      tempFactors <- array(0, c(K, K, IRFhoriz))
      tempYields  <- array(0, c(J,K,IRFhoriz))

      # Compute the IRFs
      S <- t(chol(SIGMA)) # Choleski term
      # Shock at t=0:
      tempFactors[ ,  , 1] <- S
      tempYields[ , , 1]  <- B %*% S
      # Shock at t=1:
      for (r in 2:IRFhoriz){
        tempFactors[ , , r] <- powerplus::Matpow(A1,numer=r-1, denom=1)%*%S # IRF (t+h) = A1^h*S
        tempYields[ , , r]      <- B%*%powerplus::Matpow(A1,numer=r-1, denom=1)%*%S
      }
      IRFRiskFactors <- aperm(tempFactors, c(3,1,2))
      IRFYields <- aperm(tempYields, c(3,1,2))

      Horiz <- t(t(0:(IRFhoriz-1))) #Add a column for horizon of interest


      # Adjust the variable labels
      dimnames(IRFRiskFactors)[[1]] <- Horiz
      dimnames(IRFRiskFactors)[[2]] <- AllFactorsLabels
      dimnames(IRFRiskFactors)[[3]] <- AllFactorsLabels

      YieldsLabel<- rownames(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$Y)
      dimnames(IRFYields)[[1]] <- Horiz
      dimnames(IRFYields)[[2]] <- YieldsLabel
      dimnames(IRFYields)[[3]] <- AllFactorsLabels

      # Store Country specific IRFs
      IRFCS[[i]] <- list(IRFRiskFactors, IRFYields)
      names(IRFCS[[i]]) <- c("Factors", "Yields")


    }

    # Prepare output per country
    names(IRFCS) <- Economies

    IRFoutputs[[j]] <- IRFCS

  }
  IRFoutputs  <- IRFoutputs[unlist(lapply(IRFoutputs, length) != 0)]

  names(IRFoutputs) <- ModelTypeSet[idxIndividual]

  return(IRFoutputs)

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


FEVDsep <- function(ModelType, ModelPara, FEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N

  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR sepQ")] # Exclude all models in which the estimation is made jointly.

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
      if ( ModelTypeSet[j] == "JPS") {AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]])}
      if ( ModelTypeSet[j] != "JPS") {AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries) }

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
################################### 4) GIRF #############################################################
#########################################################################################################
#' GIRFs for "sep Q" models
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param GIRFhoriz single numerical vector conataining the desired horizon of analysis for the GIRFs
#'@param FactorLabels   string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#' @references
#' \itemize{
#' \item This function is a modified and extended version of the "irf" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }




GIRFSep <- function(ModelType, ModelPara, GIRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)


  N <- length(FactorLabels$Spanned)
  C <- length(Economies) # Number of economies in the system
  M <- length(FactorLabels$Domestic) - N # Number of country-specific domestic variables
  G <- length(FactorLabels$Global) # Number of global variables

  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR sepQ")] # Exclude all models in which the estimation is made jointly

  GIRFoutputs <- list()

  for (j in  idxIndividual){
    for (i in 1:C){
      J <- numel(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$mat)
      K <- nrow(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$AllFactors)

      G0.y <- ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$Gy.0
      Sigma.y <- ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$SSZ
      F1 <- ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$ests$K1Z

      B <- BUnspannedAdapSep(G,M, ModelPara, Economies, Economies[i], ModelTypeSet[j])

      # 1) Dynamic multiplier:
      Ry.h <- array(NA, c(K,K,GIRFhoriz))
      Ry.h[, ,1] <- diag(K) # dynamic multiplier at t=0

      for (w in 2:GIRFhoriz) {
        Ry.h[, ,w] <- F1%*%Ry.h[, ,w-1]
      }

      # 2) Build the vector containing the one unit-shock for each variable of the system
      ey.j <- diag(K)

      # 3) GIRFs:
      # 3.1) Factors
      AllResponsesToAllShocksFactors <- array(NA, c(K,GIRFhoriz,K))
      AllResponsesToShockOfOneVariableFactors <- matrix(NA, ncol= GIRFhoriz , nrow = K)
      for (g in 1:K){
        for (w in 1:GIRFhoriz){
          numFactors <- (Ry.h[,,w]%*% solve(G0.y)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
          demFactors <- 1/sqrt((t(ey.j[,g])%*%Sigma.y%*%ey.j[,g])) # denominator from equation at the bottom of the page 22 (PS, 1998)
          AllResponsesToShockOfOneVariableFactors[,w] <- numFactors*drop(demFactors)
        }
        AllResponsesToAllShocksFactors[,,g] <- AllResponsesToShockOfOneVariableFactors
      }

      GIRFFactors <- aperm(AllResponsesToAllShocksFactors, c(2,1,3))

      #3.2) Yields
      AllResponsesToAllShocksYields <- array(NA, c(J,GIRFhoriz,K))
      AllResponsesToShockOfOneVariableYields <- matrix(NA, ncol= GIRFhoriz , nrow = J)
      for (g in 1:K){
        for (w in 1:GIRFhoriz){
          numYields <- B%*%(Ry.h[,,w]%*% solve(G0.y)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
          demYields <- 1/sqrt((t(ey.j[,g])%*%Sigma.y%*%ey.j[,g])) # denominator from equation at the bottom of the page 22 (PS, 1998)
          AllResponsesToShockOfOneVariableYields[,w] <- numYields*drop(demYields)
        }
        AllResponsesToAllShocksYields[,,g] <- AllResponsesToShockOfOneVariableYields
      }

      GIRFYields <- aperm(AllResponsesToAllShocksYields, c(2,1,3))

      # 4) Prepare labels for the output
      if(ModelTypeSet[j] == "JPS" ){  labelsGIRF <- c(FactorLabels$Global,FactorLabels$Tables[[Economies[i]]]) }
      if(ModelTypeSet[j] == "JPS jointP" || ModelTypeSet[j] == "GVAR sepQ" ){  labelsGIRF <- c(FactorLabels$Global,FactorLabels$Tables$AllCountries) }

      # 4.1) Add columns containig the horizons
      Horiz <- t(t(0:(GIRFhoriz-1))) #Add a column for horizon of interest

      # 4.2) Labels
      dimnames(GIRFFactors)[[1]] <- Horiz # We subtract 1, because the first element is the contemporaneous one
      dimnames(GIRFFactors)[[2]] <- labelsGIRF
      dimnames(GIRFFactors)[[3]] <- labelsGIRF

      YieldsLabel<- rownames(ModelPara[[ModelTypeSet[[j]]]][[Economies[i]]]$inputs$Y)
      dimnames(GIRFYields)[[1]] <- Horiz # We subtract 1, because the first element is the contemporaneous one
      dimnames(GIRFYields)[[2]] <- YieldsLabel
      dimnames(GIRFYields)[[3]] <- labelsGIRF


      GIRFoutputs[[ModelTypeSet[j]]][[Economies[i]]] <- list(GIRFFactors,GIRFYields)
      names(GIRFoutputs[[ModelTypeSet[j]]][[Economies[i]]]) <- c("Factors","Yields")


    }
  }

  GIRFoutputs  <- GIRFoutputs[unlist(lapply(GIRFoutputs, length) != 0)]

  names(GIRFoutputs) <- ModelTypeSet[idxIndividual]

  return(GIRFoutputs)
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




GFEVDsep <- function(ModelType, ModelPara, GFEVDhoriz, FactorLabels, Economies){


  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N

  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR sepQ")] # Exclude all models in which the estimation is made jointly.

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
      if ( ModelTypeSet[j] == "JPS") {AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]])}
      if ( ModelTypeSet[j] != "JPS") {AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries) }

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
######################################################################################################
####################### OUTPUTS FOR MODELS IN WHICH THE ESTIMATION ###################################
########################       IS DONE FOR ALL COUNTRIES JOINTLY      ################################
######################################################################################################
######################################################################################################
#' Numerical outputs (variance explained, model fit, IRFs, GIRFs, FEVDs, and GFEVDs) for "joint Q" models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param InputsForOutputs list conataining the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, and GFEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'


OutputConstructionJoint <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels,
                                    Economies){

ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
idxWishModels <- which(ModelTypeSet == ModelType)


  # Output summary
    # Total Variance Explained and Model Fit
  Total_Var_exp <- VarianceExplainedJoint(ModelType, ModelPara, FactorLabels, Economies)
  ModFit <- YieldsFitJoint(ModelType, ModelPara, FactorLabels, Economies)

  # IRF and FEVD
  IRFoutputs <- IRFjoint(ModelType, ModelPara, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)
  FEVDoutputs <- FEVDjoint(ModelType, ModelPara, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels, Economies)

  # GIRFS and GFEVD
  GIRFoutputs <- GIRFjoint(ModelType, ModelPara, InputsForOutputs[[ModelType]]$GIRF$horiz, FactorLabels, Economies)
  GFEVDoutputs <- GFEVDjoint(ModelType, ModelPara, InputsForOutputs[[ModelType]]$GFEVD$horiz, FactorLabels, Economies)


  # FOR JLL models: Non-orthogonalized IRFs
  if ("JLL original" %in% ModelType || "JLL NoDomUnit" %in% ModelType ||"JLL jointSigma" %in% ModelType){

    # Generate the outputs with orthogonalized factros:
    IRFOrtho<- IRFjointOrthoJLL(ModelType, ModelPara, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)
    GIRFOrtho <- GIRFjointOrthoJLL(ModelType, ModelPara, InputsForOutputs[[ModelType]]$GIRF$horiz, FactorLabels, Economies)
    FEVDOrtho <- FEVDjointOrthogoJLL(ModelType, ModelPara, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels,Economies)
    GFEVDOrtho <- GFEVDjointOrthoJLL(ModelType, ModelPara, InputsForOutputs[[ModelType]]$GFEVD$horiz, FactorLabels, Economies)


    jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR jointQ") ]
    for (j in jointModelIdx){
      # Merge the lists of orthogonalized and non-orthogonalized factors
      # IRF
      IRFoutputs[[ModelTypeSet[j]]]$Factors <- list(IRFoutputs[[ModelTypeSet[j]]]$Factors, IRFOrtho[[ModelTypeSet[j]]]$Factors)
      IRFoutputs[[ModelTypeSet[j]]]$Yields <- list(IRFoutputs[[ModelTypeSet[j]]]$Yields, IRFOrtho[[ModelTypeSet[j]]]$Yields)
      names(IRFoutputs[[ModelTypeSet[j]]]$Factors) <- c("NonOrtho", "Ortho")
      names(IRFoutputs[[ModelTypeSet[j]]]$Yields) <- c("NonOrtho", "Ortho")
      # GIRF
      GIRFoutputs[[ModelTypeSet[j]]]$Factors <- list(GIRFoutputs[[ModelTypeSet[j]]]$Factors, GIRFOrtho[[ModelTypeSet[j]]]$Factors)
      GIRFoutputs[[ModelTypeSet[j]]]$Yields <- list(GIRFoutputs[[ModelTypeSet[j]]]$Yields, GIRFOrtho[[ModelTypeSet[j]]]$Yields)
      names(GIRFoutputs[[ModelTypeSet[j]]]$Factors) <- c("NonOrtho", "Ortho")
      names(GIRFoutputs[[ModelTypeSet[j]]]$Yields) <- c("NonOrtho", "Ortho")
      # FEVD
      FEVDoutputs[[ModelTypeSet[j]]]$Factors <- list(FEVDoutputs[[ModelTypeSet[j]]]$Factors, FEVDOrtho[[ModelTypeSet[j]]]$Factors)
      FEVDoutputs[[ModelTypeSet[j]]]$Yields <- list(FEVDoutputs[[ModelTypeSet[j]]]$Yields, FEVDOrtho[[ModelTypeSet[j]]]$Yields)
      names(FEVDoutputs[[ModelTypeSet[j]]]$Factors) <- c("NonOrtho", "Ortho")
      names(FEVDoutputs[[ModelTypeSet[j]]]$Yields) <- c("NonOrtho", "Ortho")
      # GFEVD
      GFEVDoutputs[[ModelTypeSet[j]]]$Factors <- list(GFEVDoutputs[[ModelTypeSet[j]]]$Factors, GFEVDOrtho[[ModelTypeSet[j]]]$Factors)
      GFEVDoutputs[[ModelTypeSet[j]]]$Yields <- list(GFEVDoutputs[[ModelTypeSet[j]]]$Yields, GFEVDOrtho[[ModelTypeSet[j]]]$Yields)
      names(GFEVDoutputs[[ModelTypeSet[j]]]$Factors) <- c("NonOrtho", "Ortho")
      names(GFEVDoutputs[[ModelTypeSet[j]]]$Yields) <- c("NonOrtho", "Ortho")

    }
  }



  NumericalOutputs <- list(Total_Var_exp, ModFit, IRFoutputs, FEVDoutputs, GIRFoutputs, GFEVDoutputs)
  names(NumericalOutputs) <- c("PC var explained","Fit","IRF", 'FEVD', "GIRF", "GFEVD")

  return(NumericalOutputs)

}


######################################################################################################
####################### 1) Total Variance Explained #########################################
######################################################################################################
#' Percentage explained by the spanned factors of the variations in the set of observed yields for "joint Q" models
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'



VarianceExplainedJoint<-function(ModelType, ModelPara, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)


  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  J <- numel(ModelPara[[ModelType]]$inputs$mat)
  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR sepQ") ]
  L <- which(ModelTypeSet == "GVAR sepQ") # index of the model right before the first desired joint model.

  Total_Var_exp_per_country <- list()
  Total_Var_exp_All <- list()

  for(j in jointModelIdx){
    idx0 <- 0
    for (i in 1:C){
      idx1 <- idx0 + J
      H <- eigen(stats::cov(t(ModelPara[[ModelTypeSet[j]]]$inputs$Y[(idx0+1):idx1,])))$values
      percentages_explained <- cumsum(H)/sum(H)
      Total_Var_exp_per_country[[i]] <-percentages_explained[1:N]
      idx0 <- idx1
    }
    names(Total_Var_exp_per_country) <- Economies
    Total_Var_exp_All[[j-L]] <- Total_Var_exp_per_country
  }

  # Clean empty lists
  Total_Var_exp_All  <- Total_Var_exp_All[unlist(lapply(Total_Var_exp_All, length) != 0)]

  names(Total_Var_exp_All) <- ModelTypeSet[jointModelIdx]


  return(Total_Var_exp_All)
}

######################################################################################################
########################################### 2) Model Fit #############################################
######################################################################################################
#' Computes two measures of model fit for bond yields
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@details
#' "Model-implied yields" is the measure of fit based exclusively on the risk-neutral parameters, whereas the
#' "Model-Fit" takes into account both the risk-neutral and the physical parameters.
#'
#' @references
#' See, for instance, Jotiskhatira, Le and Lundblad (2015). "Why do interest rates in different currencies co-move?" (Journal of Financial Economics)


YieldsFitJoint <- function(ModelType, ModelPara, FactorLabels, Economies){

  # 1) Redefining variables
  Z <- ModelPara[[ModelType]]$inputs$AllFactors
  YieldData <- ModelPara[[ModelType]]$inputs$Y

  mat <- ModelPara[[ModelType]]$inputs$mat
  Afull <- ModelPara[[ModelType]]$rot$P$A
  Bspanned <- ModelPara[[ModelType]]$rot$P$B
  K0Z <- ModelPara[[ModelType]]$ests$K0Z
  K1Z <- ModelPara[[ModelType]]$ests$K1Z

  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)
  J <- length(mat)
  CJ <- C*J
  T <- ncol(Z)



  # 2) MODEL FIT MEASURES
  # a) Model Fit (Yields)
  # Extract spanned factors from the list of unspanned factors
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

  P <- Z[IdxSpanned,] # Set of spanned factors

  # Compute model fit
  Yieldfit<- matrix(NA, nrow= CJ, ncol =T)
  for (i in 1:T){
    Yieldfit[,i] <- Afull + Bspanned%*%P[,i]
  }
  rownames(Yieldfit) <- rownames(YieldData)
  colnames(Yieldfit) <- colnames(YieldData)


  # b) Model-implied Yields
  Bfull <- BUnspannedAdapJoint(G,M,N,C, J, Bspanned)
  row.names(Bfull) <- rownames(YieldData)
  colnames(Bfull) <- rownames(Z)

  YieldModelImplied<- matrix(NA, nrow= CJ, ncol =T)
  for (i in 2:T){ #  first observation is discarded
    YieldModelImplied[,i] <- Afull + Bfull%*%(K0Z + K1Z%*% Z[,i-1])
  }
  rownames(YieldModelImplied) <- rownames(YieldData)
  colnames(YieldModelImplied) <- colnames(YieldData)


  Output<- list(Yieldfit, YieldModelImplied)
  names(Output) <- c("Yield Fit","Yield Model Implied")

  return(Output)
}

######################################################################################################
########################################### 3) IRFs ##################################################
######################################################################################################
#' IRFs for "joint Q" models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param IRFhoriz single numerical vector conataining the desired horizon of analysis for the IRFs
#'@param FactorLabels string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@details
#' Structural shocks are identified via Cholesky decomposition

#'
#'
IRFjoint <- function(ModelType, ModelPara, IRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies)
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  # Pre-allocation
  IRFoutputs <- list()

  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR sepQ") ]
  L <- which(ModelTypeSet == "GVAR sepQ") # index of the model right before the first desired joint model.

  for(j in jointModelIdx){

    J <- numel(ModelPara[[ModelType]]$inputs$mat)
    K <- nrow(ModelPara[[ModelType]]$inputs$AllFactors)
    M <- (K-G)/C - N
    CJ <- C*J

    # Generate factor labels depending on the models to be estimated
    AllCountryLabels <- c()
    AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries)

    # Summarize inputs for the IRFs
    SIGMA <- ModelPara[[ModelTypeSet[[j]]]]$ests$SSZ # KxK (variance-covariance matrix)
    A0 <- ModelPara[[ModelTypeSet[[j]]]]$ests$K0Z # Kx1 (matrix of intercepts)
    A1 <- ModelPara[[ModelTypeSet[[j]]]]$ests$K1Z # KxK (feedback matrix)

    BSpanned <- ModelPara[[ModelTypeSet[[j]]]]$rot$P$B
    B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

    # Initialization of IRFs of interest
    tempFactors <- array(0, c(K, K, IRFhoriz))
    tempYields  <- array(0, c(CJ,K,IRFhoriz))

    # Compute the IRFs
   if ( ModelType == "JLL original" || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma" ){
     S <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$Sigmas$Sigma_Y
   }else{ S <- t(chol(SIGMA))} # Choleski term
    # Shock at t=0:
    tempFactors[ ,  , 1] <- S
    tempYields[ , , 1]  <- B %*% S
    # Shock at t=1:
    for (r in 2:IRFhoriz){
      tempFactors[ , , r] <- powerplus::Matpow(A1,numer=r-1, denom=1)%*%S # IRF (t+h) = A1^h*S
      tempYields[ , , r]      <- B%*%powerplus::Matpow(A1,numer=r-1, denom=1)%*%S
    }
    IRFRiskFactors <- aperm(tempFactors, c(3,1,2))
    IRFYields <- aperm(tempYields, c(3,1,2))

    Horiz <- t(t(0:(IRFhoriz-1))) #Add a column for horizon of interest


    # Adjust the variable labels
    dimnames(IRFRiskFactors)[[1]] <- Horiz
    dimnames(IRFRiskFactors)[[2]] <- AllFactorsLabels
    dimnames(IRFRiskFactors)[[3]] <- AllFactorsLabels

    YieldsLabel<- rownames(ModelPara[[ModelTypeSet[[j]]]]$inputs$Y)
    dimnames(IRFYields)[[1]] <- Horiz
    dimnames(IRFYields)[[2]] <- YieldsLabel
    dimnames(IRFYields)[[3]] <- AllFactorsLabels


    IRFoutputsAllCountries <- list(IRFRiskFactors,IRFYields)
    names(IRFoutputsAllCountries) <- c("Factors","Yields")

    # Prepare output per country
    IRFoutputs[[j-L]] <- IRFoutputsAllCountries

  }
  # Clean empty lists
  IRFoutputs  <- IRFoutputs[unlist(lapply(IRFoutputs, length) != 0)]

  names(IRFoutputs) <- ModelTypeSet[jointModelIdx]


  return(IRFoutputs)

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


FEVDjoint <- function(ModelType, ModelPara, FEVDhoriz, FactorLabels, Economies){


  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic)-N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G


  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR sepQ") ]
  L <- which(ModelTypeSet == "GVAR sepQ") # index of the model right before the first desired joint model.

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
    if ( ModelType == "JLL original" || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma" ){
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
################################### 5) GIRF #############################################################
#########################################################################################################
#' GIRFs for "joint Q" models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param ModelPara List of model parameter estimates (See the "Optimization" function)
#'@param GIRFhoriz single numerical vector conataining the desired horizon of analysis for the GIRFs
#'@param FactorLabels  a string-list based which contains all the labels of all the variables present in the model
#'@param Economies a string-vector containing the names of the economies which are part of the economic system

#'
#'@references
#' \itemize{
#' \item This function is a modified and extended version of the "irf" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }



GIRFjoint <- function(ModelType, ModelPara, GIRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies) # Number of economies in the system
  M <- length(FactorLabels$Domestic) - N # Number of country-specific domestic variables
  G <- length(FactorLabels$Global) # Number of global variables
  K <- C*(M+N)+G # All factors of the system


  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR sepQ") ]
  L <- which(ModelTypeSet == "GVAR sepQ") # index of the model right before the first desired joint model.

  GIRFoutputs <- list()

  for (j in  jointModelIdx){
    J <- numel(ModelPara[[ModelTypeSet[[j]]]]$inputs$mat)
    CJ <- C*J

    G0.y <- ModelPara[[ModelTypeSet[[j]]]]$ests$Gy.0
    Sigma.y <- ModelPara[[ModelTypeSet[[j]]]]$ests$SSZ
    F1 <- ModelPara[[ModelTypeSet[[j]]]]$ests$K1Z

    BSpanned <- ModelPara[[ModelTypeSet[[j]]]]$rot$P$B
    B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

    # 1) Dynamic multiplier:
    Ry.h <- array(NA, c(K,K,GIRFhoriz))
    Ry.h[, ,1] <- diag(K) # dynamic multiplier at t=0

    for (i in 2:GIRFhoriz) {
      Ry.h[, ,i] <- F1%*%Ry.h[, ,i-1]
    }

    # 2) Build the vector containing the one unit-shock for each variable of the system
    ey.j <- diag(K)

    # 3) GIRFs:
    # 3.1) Factors
    AllResponsesToAllShocksFactors <- array(NA, c(K,GIRFhoriz,K))
    AllResponsesToShockOfOneVariableFactors <- matrix(NA, ncol= GIRFhoriz , nrow = K)
    for (g in 1:K){
      for (i in 1:GIRFhoriz){
        numFactors <- (Ry.h[,,i]%*% solve(G0.y)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
        demFactors <- 1/sqrt((t(ey.j[,g])%*%Sigma.y%*%ey.j[,g])) # denominator from equation at the bottom of the page 22 (PS, 1998)
        AllResponsesToShockOfOneVariableFactors[,i] <- numFactors*drop(demFactors)
      }
      AllResponsesToAllShocksFactors[,,g] <- AllResponsesToShockOfOneVariableFactors
    }


    GIRFFactors <- aperm(AllResponsesToAllShocksFactors, c(2,1,3))


    #3.2) Yields
    AllResponsesToAllShocksYields <- array(NA, c(CJ,GIRFhoriz,K))
    AllResponsesToShockOfOneVariableYields <- matrix(NA, ncol= GIRFhoriz , nrow = CJ)
    for (g in 1:K){
      for (i in 1:GIRFhoriz){
        numYields <- B%*%(Ry.h[,,i]%*% solve(G0.y)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
        demYields <- 1/sqrt((t(ey.j[,g])%*%Sigma.y%*%ey.j[,g])) # denominator from equation at the bottom of the page 22 (PS, 1998)
        AllResponsesToShockOfOneVariableYields[,i] <- numYields*drop(demYields)
      }
      AllResponsesToAllShocksYields[,,g] <- AllResponsesToShockOfOneVariableYields
    }

    GIRFYields <- aperm(AllResponsesToAllShocksYields, c(2,1,3))

    # 4) Prepare labels for the output
    # 4.1) Add columns containig the horizons
    Horiz <- t(t(0:(GIRFhoriz-1))) #Add a column for horizon of interest

    # 4.2) Labels
    labelsGIRF <- c(FactorLabels$Global,FactorLabels$Tables$AllCountries)

    dimnames(GIRFFactors)[[1]] <- 0:(GIRFhoriz-1) # We subtract 1, because the first element is the contemporaneous one
    dimnames(GIRFFactors)[[2]] <- labelsGIRF
    dimnames(GIRFFactors)[[3]] <- labelsGIRF

    YieldsLabel<- rownames(ModelPara[[ModelTypeSet[[j]]]]$inputs$Y)
    dimnames(GIRFYields)[[1]] <- 0:(GIRFhoriz-1) # We subtract 1, because the first element is the contemporaneous one
    dimnames(GIRFYields)[[2]] <- YieldsLabel
    dimnames(GIRFYields)[[3]] <- labelsGIRF


    GIRFoutputsAllCountries <- list(GIRFFactors,GIRFYields)
    names(GIRFoutputsAllCountries) <- c("Factors","Yields")

    GIRFoutputs[[j-L]] <- GIRFoutputsAllCountries

  }
  # Clean empty lists
  GIRFoutputs  <- GIRFoutputs[unlist(lapply(GIRFoutputs, length) != 0)]

  names(GIRFoutputs) <- ModelTypeSet[jointModelIdx]

  return(GIRFoutputs)
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


GFEVDjoint <- function(ModelType, ModelPara, GFEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic)-N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G


  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR sepQ") ]
  L <- which(ModelTypeSet == "GVAR sepQ") # index of the model right before the first desired joint model.

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
############################ 7) IRFs with orthogonalized factors #####################################
######################################################################################################
#' Orthogonalized IRFs for JLL models
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara list of model parameter estimates (see the "Optimization" function)
#'@param IRFhoriz single numerical vector conataining the desired horizon of analysis for the IRFs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@details
#' Structural shocks are identified via Cholesky decomposition




IRFjointOrthoJLL <- function(ModelType, ModelPara, IRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  # Pre-allocation
  IRFoutputs <- list()

  jointModelIdxJLL <- idxWishModels[idxWishModels >= which(ModelTypeSet == "JLL original") ]
  L <- which(ModelTypeSet == "GVAR sepQ") # index of the model right before the first desired joint model.

  for(j in jointModelIdxJLL){

    J <- numel(ModelPara[[ModelTypeSet[[j]]]]$inputs$mat)
    K <- nrow(ModelPara[[ModelTypeSet[[j]]]]$inputs$AllFactors)
    M <- (K-G)/C - N
    CJ <- C*J

    # Generate factor labels depending on the models to be estimated
    AllCountryLabels <- c()
    AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountriesJLL)

    # Summarize inputs for the IRFs
    SIGMAe <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$Sigmas$VarCov_Ortho # KxK (variance-covariance matrix)
    A0e <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$k0_e # Kx1 (matrix of intercepts)
    A1e <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$k1_e # KxK (feedback matrix)
    PI <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$PI

    BSpanned <- ModelPara[[ModelTypeSet[[j]]]]$rot$P$B
    B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

    # Initialization of IRFs of interest
    tempFactors <- array(0, c(K, K, IRFhoriz))
    tempYields  <- array(0, c(CJ,K,IRFhoriz))

    # Compute the IRFs
    Se <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$Sigmas$Sigma_Ye

    # Shock at t=0:
    tempFactors[ ,  , 1] <- Se
    tempYields[ , , 1]  <- B%*%PI%*% Se
    # Shock at t=1:
    for (r in 2:IRFhoriz){
      tempFactors[ , , r] <- powerplus::Matpow(A1e,numer=r-1, denom=1)%*%Se # IRF (t+h) = A1^h*S
      tempYields[ , , r]      <- B%*%PI%*%powerplus::Matpow(A1e,numer=r-1, denom=1)%*%Se
    }
    IRFRiskFactors <- aperm(tempFactors, c(3,1,2))
    IRFYields <- aperm(tempYields, c(3,1,2))

    Horiz <- t(t(0:(IRFhoriz-1))) #Add a column for horizon of interest

    # Adjust the variable labels
    dimnames(IRFRiskFactors)[[1]] <- Horiz
    dimnames(IRFRiskFactors)[[2]] <- AllFactorsLabels
    dimnames(IRFRiskFactors)[[3]] <- AllFactorsLabels

    YieldsLabel<- rownames(ModelPara[[ModelTypeSet[[j]]]]$inputs$Y)
    dimnames(IRFYields)[[1]] <- Horiz
    dimnames(IRFYields)[[2]] <- YieldsLabel
    dimnames(IRFYields)[[3]] <- AllFactorsLabels


    IRFoutputsAllCountries <- list(IRFRiskFactors,IRFYields)
    names(IRFoutputsAllCountries) <- c("Factors","Yields")

    # Prepare output per country
    IRFoutputs[[j-L]] <- IRFoutputsAllCountries

  }
  # Clean empty lists
  IRFoutputs  <- IRFoutputs[unlist(lapply(IRFoutputs, length) != 0)]

  names(IRFoutputs) <- ModelTypeSet[jointModelIdxJLL]


  return(IRFoutputs)

}

#########################################################################################################
################################### 8) FEVD with orthogonalized factors #################################
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


FEVDjointOrthogoJLL <- function(ModelType, ModelPara, FEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic)-N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G

  jointModelIdxJLL <- idxWishModels[idxWishModels >= which(ModelTypeSet == "JLL original") ]
  L <- which(ModelTypeSet == "GVAR sepQ") # index of the model right before the first desired joint model.

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
################################### 9) GIRF With orthogonalized factors #################################
#########################################################################################################
#' Orthogonalized GIRFs for JLL models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param ModelPara List of model parameter estimates (See the "Optimization" function)
#'@param GIRFhoriz single numerical vector conataining the desired horizon of analysis for the GIRFs
#'@param FactorLabels  a string-list based which contains all the labels of all the variables present in the model
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@references
#' \itemize{
#' \item This function is a modified and extended version of the "irf" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }



GIRFjointOrthoJLL <- function(ModelType, ModelPara, GIRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies) # Number of economies in the system
  M <- length(FactorLabels$Domestic) - N # Number of country-specific domestic variables
  G <- length(FactorLabels$Global) # Number of global variables
  K <- C*(M+N)+G # All factors of the system


  jointModelIdxJLL <- idxWishModels[idxWishModels >= which(ModelTypeSet == "JLL original") ]
  L <- which(ModelTypeSet == "GVAR sepQ") # index of the model right before the first desired joint model.

  GIRFoutputs <- list()

  for (j in  jointModelIdxJLL){
    J <- numel(ModelPara[[ModelTypeSet[[j]]]]$inputs$mat)
    CJ <- C*J

    G0.y <- ModelPara[[ModelTypeSet[[j]]]]$ests$Gy.0
    Sigma.y <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$Sigmas$VarCov_Ortho
    F1e <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$k1_e
    PI <- ModelPara[[ModelTypeSet[[j]]]]$ests$JLLoutcomes$PI

    BSpanned <- ModelPara[[ModelTypeSet[[j]]]]$rot$P$B
    B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

    # 1) Dynamic multiplier:
    Ry.h <- array(NA, c(K,K,GIRFhoriz))
    Ry.h[, ,1] <- diag(K) # dynamic multiplier at t=0

    for (i in 2:GIRFhoriz) {
      Ry.h[, ,i] <- F1e%*%Ry.h[, ,i-1]
    }

    # 2) Build the vector containing the one unit-shock for each variable of the system
    ey.j <- diag(K)

    # 3) GIRFs:
    # 3.1) Factors
    AllResponsesToAllShocksFactors <- array(NA, c(K,GIRFhoriz,K))
    AllResponsesToShockOfOneVariableFactors <- matrix(NA, ncol= GIRFhoriz , nrow = K)
    for (g in 1:K){
      for (i in 1:GIRFhoriz){
        numFactors <- (PI%*%Ry.h[,,i]%*% solve(G0.y)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22
        demFactors <- 1/sqrt((t(ey.j[,g])%*%Sigma.y%*%ey.j[,g])) # denominator from equation at the bottom of the page 22
        AllResponsesToShockOfOneVariableFactors[,i] <- numFactors*drop(demFactors)
      }
      AllResponsesToAllShocksFactors[,,g] <- AllResponsesToShockOfOneVariableFactors
    }

    GIRFFactors <- aperm(AllResponsesToAllShocksFactors, c(2,1,3))

    #3.2) Yields
    AllResponsesToAllShocksYields <- array(NA, c(CJ,GIRFhoriz,K))
    AllResponsesToShockOfOneVariableYields <- matrix(NA, ncol= GIRFhoriz , nrow = CJ)
    for (g in 1:K){
      for (i in 1:GIRFhoriz){
        numYields <- B%*%(PI%*%Ry.h[,,i]%*% solve(G0.y)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22
        demYields <- 1/sqrt((t(ey.j[,g])%*%Sigma.y%*%ey.j[,g])) # denominator from equation at the bottom of the page 22
        AllResponsesToShockOfOneVariableYields[,i] <- numYields*drop(demYields)
      }
      AllResponsesToAllShocksYields[,,g] <- AllResponsesToShockOfOneVariableYields
    }


    GIRFYields <- aperm(AllResponsesToAllShocksYields, c(2,1,3))


    # 4) Prepare labels for the output
    # 4.1) Add columns containig the horizons
    Horiz <- t(t(0:(GIRFhoriz-1))) #Add a column for horizon of interest

    # 4.2) Labels
    labelsGIRF <- c(FactorLabels$Global,FactorLabels$Tables$AllCountries)

    dimnames(GIRFFactors)[[1]] <- Horiz # We subtract 1, because the first element is the contemporaneous one
    dimnames(GIRFFactors)[[2]] <- labelsGIRF
    dimnames(GIRFFactors)[[3]] <- labelsGIRF

    YieldsLabel<- rownames(ModelPara[[ModelTypeSet[[j]]]]$inputs$Y)
    dimnames(GIRFYields)[[1]] <- Horiz # We subtract 1, because the first element is the contemporaneous one
    dimnames(GIRFYields)[[2]] <- YieldsLabel
    dimnames(GIRFYields)[[3]] <- labelsGIRF


    GIRFoutputsAllCountries <- list(GIRFFactors,GIRFYields)
    names(GIRFoutputsAllCountries) <- c("Factors","Yields")

    GIRFoutputs[[j-L]] <- GIRFoutputsAllCountries

  }
  # Clean empty lists
  GIRFoutputs  <- GIRFoutputs[unlist(lapply(GIRFoutputs, length) != 0)]

  names(GIRFoutputs) <- ModelTypeSet[jointModelIdxJLL]

  return(GIRFoutputs)
}

#########################################################################################################
################################### 10) GFEVD With orthogonalized factor s################################
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


GFEVDjointOrthoJLL <- function(ModelType, ModelPara, GFEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS", "JPS jointP", "GVAR sepQ", "VAR jointQ" , "GVAR jointQ",
                    "JLL original", "JLL NoDomUnit", "JLL jointSigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic)-N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G


  jointModelIdxJLL <- idxWishModels[idxWishModels >= which(ModelTypeSet == "JLL original") ]
  L <- which(ModelTypeSet == "GVAR sepQ") # index of the model right before the first desired joint model.

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



#####################################################################################################
#' Transform B_spanned into B_unspanned for jointQ models


#'@param G number of global unspanned factors
#'@param M number of domestic unspanned factors
#'@param N number of domestic spanned factors
#'@param C number of economies of the economic system
#'@param J number of country-specific observed bond yields
#'@param BSpanned B that accomodates only the map to the spanned factors only

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

BUnspannedAdapSep <- function(G,M, ModelPara, Economies, Economy, ModelType){

  i <- match(Economy, Economies)
  C <- length(Economies)
  N <- ModelPara[[ModelType]][[Economy]]$inputs$N
  J <- length(ModelPara[[ModelType]][[Economy]]$inputs$mat)


  if( ModelType == "JPS"){
    K <- N+M + G
    BUnspanned <- matrix(0, nrow=J, ncol= K)
    BSpanned <- ModelPara[[ModelType]][[Economies[i]]]$rot$P$B
    BUnspanned[ , (K-N+1):K] <-  BSpanned
  }



  if( ModelType == "JPS jointP" || ModelType == "GVAR sepQ" ){
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

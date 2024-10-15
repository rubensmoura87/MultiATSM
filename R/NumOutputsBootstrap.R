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


  AllNumOutputs <- list()


  # If one chooseS models in which the estimation is done country-by-country
  if ( any(ModelType == c("JPS original" , "JPS global", "GVAR single"))){
    NumOutSep <- OutputConstructionSep_BS(ModelType, ModelParaBoot, InputsForOutputs, FactorLabels, Economies)
  } else{
  # If one chooseS models in which the estimation is done jointly for all countries
    NumOutJoint <- OutputConstructionJoint_BS(ModelType, ModelParaBoot, InputsForOutputs, FactorLabels, Economies)
  }


  # Prepare final list of outputs
  if (!exists("NumOutJoint")){ AllNumOutputs <- NumOutSep}
  if (!exists("NumOutSep")){ AllNumOutputs <- NumOutJoint}
  if ( exists("NumOutSep") & exists("NumOutJoint")){

    for (i in 1:length(NumOutSep)){AllNumOutputs[[i]] <- append(NumOutSep[[i]], NumOutJoint[[i]]) }
    names(AllNumOutputs) <- names(NumOutSep)
  }

  return(AllNumOutputs)
}


######################################################################################################
######################################################################################################
####################### OUTPUTS FOR MODELS IN WHICH THE ESTIMATION ###################################
########################       IS DONE COUNTRY-BY-COUNTRY      #######################################
######################################################################################################
######################################################################################################
#' Gathers all the model numerical ouputs after bootstrap for "sep Q" models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param InputsForOutputs list conataining the desired inputs for the cunstruction of the model fit, IRFs, GIRFs, FEVDs, and GFEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@keywords internal


OutputConstructionSep_BS <- function(ModelType, ModelParaBoot, InputsForOutputs, FactorLabels, Economies){


  # Output summary
  # IRF and FEVD
  IRFoutputs <- IRFsep_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)
  FEVDoutputs <- FEVDsep_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels, Economies)

  # GIRF and GFEVD
  GIRFoutputs <- GIRFSep_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$GIRF$horiz, FactorLabels, Economies)
  GFEVDoutputs <- GFEVDsep_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$GFEVD$horiz, FactorLabels, Economies)

  NumericalOutputs <- list(IRF = IRFoutputs, FEVD = FEVDoutputs, GIRF = GIRFoutputs, GFEVD = GFEVDoutputs)

  return(NumericalOutputs)

}



######################################################################################################
########################################### 2) IRFs ##################################################
######################################################################################################
#' IRFs after bootstrap for "sep Q" models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param IRFhoriz single numerical vector conataining the desired horizon of analysis for the IRFs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal



IRFsep_BS <- function(ModelType, ModelParaBoot, IRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies)
  J <- numel(ModelParaBoot$GeneralInputs$mat)
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N

  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]][[1]])

  # Pre-allocation
  IRFCS <- list()
  IRFoutputsCS <- list()
  IRFoutputsAllCountries <- list()
  IRFoutputs <- list()


  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR single")] # Exclude all models in which the estimation is made jointly

  for(j in idxIndividual){
    for (i in 1:C){
      for (tt in 1:ndraws){

        K <- nrow(ModelParaBoot$ParaDraws[[ModelTypeSet[j]]][[Economies[i]]][[tt]]$ests$K1Z)


        # Generate factor labels depending on the models to be estimated
        if ( ModelTypeSet[j] == "JPS original") {AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]])}
        if ( ModelTypeSet[j] != "JPS original") {AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries) }

        # Summarize inputs for the IRFs
        SIGMA <- ModelParaBoot$ParaDraws[[ModelTypeSet[j]]][[Economies[i]]][[tt]]$ests$SSZ # KxK (variance-covariance matrix)
        A1 <- ModelParaBoot$ParaDraws[[ModelTypeSet[j]]][[Economies[i]]][[tt]]$ests$K1Z # KxK (feedback matrix)

        B  <- BUnspannedAdapSep_BS(G, M, ModelParaBoot, Economies, Economies[i], ModelTypeSet[j],tt)

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
          if (r == 2) { A1h <- A1} else { A1h <- A1h %*% A1}

          tempFactors[ , , r] <- A1h%*%S # IRF (t+h) = A1^h*S
          tempYields[ , , r]      <- B%*%A1h%*%S
        }
        IRFRiskFactors <- aperm(tempFactors, c(3,1,2))
        IRFYields <- aperm(tempYields, c(3,1,2))

        Horiz <- t(t(0:(IRFhoriz-1))) # Full horizon of interest


        # Adjust the variable labels
        dimnames(IRFRiskFactors)[[1]] <- Horiz
        dimnames(IRFRiskFactors)[[2]] <- AllFactorsLabels
        dimnames(IRFRiskFactors)[[3]] <- AllFactorsLabels


        YieldsLabel <- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$inputs$Y)
        dimnames(IRFYields)[[1]] <- Horiz
        dimnames(IRFYields)[[2]] <- YieldsLabel
        dimnames(IRFYields)[[3]] <- AllFactorsLabels

        # Store Country specific IRFs
        IRFCS <- list(IRFRiskFactors, IRFYields)
        names(IRFCS) <- c("Factors", "Yields")


        # Prepare output per country per draw
        IRFoutputsCS[[tt]] <- IRFCS

      }
      # All draws of all counntries for ONE model
      IRFoutputsAllCountries[[Economies[i]]] <- IRFoutputsCS

    }
    # All draws of all counntries for ALL models
    IRFoutputs[[ModelTypeSet[j]]] <- IRFoutputsAllCountries
  }

  IRFoutputs  <- IRFoutputs[unlist(lapply(IRFoutputs, length) != 0)] # clean empty lists
  #names(IRFoutputs) <- ModelTypeSet[idxIndividual]

  return(IRFoutputs)

}

#########################################################################################################
################################### 3) FEVD #############################################################
#########################################################################################################
#' FEVDs after bootstrap for "sep Q" models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param FEVDhoriz single numerical vector conataining the desired horizon of analysis for the FEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@keywords internal


FEVDsep_BS <- function(ModelType, ModelParaBoot, FEVDhoriz, FactorLabels, Economies){


  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  J <- numel(ModelParaBoot$GeneralInputs$mat)
  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]][[1]])

  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR single")] # Exclude all models in which the estimation is made jointly.

  FEVDoutputsCS <- list()
  FEVDoutputsAllCountries <- list()
  FEVDoutputs <- list()


  for (j in  idxIndividual){
    for (i in 1:C){
      for (tt in 1:ndraws){

        K <- nrow(ModelParaBoot$ParaDraws[[ModelTypeSet[j]]][[Economies[i]]][[tt]]$ests$K1Z)

        G0 <-  ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$ests$Gy.0
        Sigma_y <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$ests$SSZ
        F1 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$ests$K1Z


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
        BallFac <- BUnspannedAdapSep_BS(G, M, ModelParaBoot, Economies, Economies[i], ModelTypeSet[j], tt)


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
        if ( ModelTypeSet[j] != "JPS original") {AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries) }

        YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$inputs$Y)

        dimnames(FEVDFactors)[[1]] <- 1:(FEVDhoriz)
        dimnames(FEVDFactors)[[2]] <- AllFactorsLabels
        dimnames(FEVDFactors)[[3]] <- AllFactorsLabels


        dimnames(FEVDYields)[[1]] <- 1:(FEVDhoriz)
        dimnames(FEVDYields)[[2]] <- AllFactorsLabels
        dimnames(FEVDYields)[[3]] <- YieldsLabel


        FEVDoutputsCS[[tt]] <- list(FEVDFactors,FEVDYields)
        names(FEVDoutputsCS[[tt]]) <- c("Factors","Yields")

      }

      FEVDoutputsAllCountries[[Economies[i]]] <- FEVDoutputsCS

    }
    FEVDoutputs[[ModelTypeSet[j]]] <- FEVDoutputsAllCountries
  }
  # Clean empty lists
  FEVDoutputs  <- FEVDoutputs[unlist(lapply(FEVDoutputs, length) != 0)]


  return(FEVDoutputs)

}


#########################################################################################################
################################### 4) GIRF #############################################################
#########################################################################################################
#' GIRFs after bootstrap for "sep Q" models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param GIRFhoriz single numerical vector conataining the desired horizon of analysis for the GIRFs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@references
#' \itemize{
#' \item This function is a modified and extended version of the "irf" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#'
#'@keywords internal



GIRFSep_BS <- function(ModelType, ModelParaBoot, GIRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  N <- length(FactorLabels$Spanned)
  C <- length(Economies) # Number of economies in the system
  M <- length(FactorLabels$Domestic) - N   # Number of country-specific domestic variables
  G <- length(FactorLabels$Global) # Number of global variables
  J <- numel(ModelParaBoot$GeneralInputs$mat)
  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]][[1]])

  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR single")] # Exclude all models in which the estimation is made jointly

  GIRFoutputsCS <- list()
  GIRFoutputsAllCountries <- list()
  GIRFoutputs <- list()

  for (j in  idxIndividual){
    for (i in 1:C){
      for (tt in 1:ndraws){


        K <- nrow(ModelParaBoot$ParaDraws[[ModelTypeSet[j]]][[Economies[i]]][[tt]]$ests$K1Z)
        Gy.0 <- ModelParaBoot$ParaDraws[[ModelTypeSet[j]]][[Economies[i]]][[tt]]$ests$Gy.0
        Sigma.y <- ModelParaBoot$ParaDraws[[ModelTypeSet[j]]][[Economies[i]]][[tt]]$ests$SSZ
        F1 <- ModelParaBoot$ParaDraws[[ModelTypeSet[j]]][[Economies[i]]][[tt]]$ests$K1Z

        B <- BUnspannedAdapSep_BS(G, M, ModelParaBoot, Economies, Economies[i], ModelTypeSet[j], tt)

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
            numFactors <- (Ry.h[,,w]%*% solve(Gy.0)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
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
            numYields <- B%*%(Ry.h[,,w]%*% solve(Gy.0)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
            demYields <- 1/sqrt((t(ey.j[,g])%*%Sigma.y%*%ey.j[,g])) # denominator from equation at the bottom of the page 22 (PS, 1998)
            AllResponsesToShockOfOneVariableYields[,w] <- numYields*drop(demYields)
          }
          AllResponsesToAllShocksYields[,,g] <- AllResponsesToShockOfOneVariableYields
        }

        GIRFYields<- aperm(AllResponsesToAllShocksYields, c(2,1,3))

        # 4) Prepare labels for the output
        if(ModelTypeSet[j] == "JPS original" ){  labelsGIRF <- c(FactorLabels$Global,FactorLabels$Tables[[Economies[i]]]) }
        if(ModelTypeSet[j] == "JPS global" || ModelTypeSet[j] == "GVAR single" ){  labelsGIRF <- c(FactorLabels$Global,FactorLabels$Tables$AllCountries) }

        # 4.1) Labels
        Horiz <- t(t(0:(GIRFhoriz-1))) # horizon of interest

        dimnames(GIRFFactors)[[1]] <- Horiz # We subtract 1, because the first element is the contemporaneous one
        dimnames(GIRFFactors)[[2]] <- labelsGIRF
        dimnames(GIRFFactors)[[3]] <- labelsGIRF

        YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$inputs$Y)
        dimnames(GIRFYields)[[1]] <- Horiz # We subtract 1, because the first element is the contemporaneous one
        dimnames(GIRFYields)[[2]] <- YieldsLabel
        dimnames(GIRFYields)[[3]] <- labelsGIRF


        GIRFoutputsCS[[tt]] <- list(GIRFFactors,GIRFYields)
        names(GIRFoutputsCS[[tt]]) <- c("Factors","Yields")

      }

      GIRFoutputsAllCountries[[Economies[i]]] <- GIRFoutputsCS
    }

    GIRFoutputs[[ModelTypeSet[j]]] <- GIRFoutputsAllCountries
  }
  GIRFoutputs  <- GIRFoutputs[unlist(lapply(GIRFoutputs, length) != 0)]


  return(GIRFoutputs)
}

#########################################################################################################
################################### 5) GFEVD #############################################################
#########################################################################################################
#' GFEVDs after bootstrap for "sep Q" models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param GFEVDhoriz single numerical vector conataining the desired horizon of analysis for the GFEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
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
#'@keywords internal




GFEVDsep_BS <- function(ModelType, ModelParaBoot, GFEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  J <- length(ModelParaBoot$GeneralInputs$mat)
  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]][[1]])

  idxIndividual <- idxWishModels[ idxWishModels <= which(ModelTypeSet== "GVAR single")] # Exclude all models in which the estimation is made jointly.

  GFEVDoutputsCS <- list()
  GFEVDoutputsAllCountires <- list()
  GFEVDoutputs <- list()

  for (j in  idxIndividual){
    for (i in 1:C){
      for (tt in 1:ndraws){

        K <- nrow(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$ests$K1Z)

        G0 <-  ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$ests$Gy.0
        Sigma_y <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$ests$SSZ
        F1 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$ests$K1Z


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
              acc1 <- t((eslct[,h]%*%Ry.h[,,l]%*%invGSigmau%*%vslct)^2)
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
        Bfull <- BUnspannedAdapSep_BS(G, M, ModelParaBoot, Economies, Economies[i], ModelTypeSet[j], tt)



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

        YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[Economies[i]]][[tt]]$inputs$Y)

        dimnames(GFEVDFactorsNormalized)[[1]] <- 1:(GFEVDhoriz)
        dimnames(GFEVDFactorsNormalized)[[2]] <- AllFactorsLabels
        dimnames(GFEVDFactorsNormalized)[[3]] <- AllFactorsLabels

        dimnames(GFEVDYieldsNormalized)[[1]] <- 1:(GFEVDhoriz)
        dimnames(GFEVDYieldsNormalized)[[2]] <- AllFactorsLabels
        dimnames(GFEVDYieldsNormalized)[[3]] <- YieldsLabel

        GFEVDoutputsCS[[tt]] <- list(GFEVDFactorsNormalized,GFEVDYieldsNormalized)
        names(GFEVDoutputsCS[[tt]]) <- c("Factors","Yields")

      }

      GFEVDoutputsAllCountires[[Economies[i]]] <- GFEVDoutputsCS

    }
    GFEVDoutputs[[ModelTypeSet[j]]]<- GFEVDoutputsAllCountires

  }
  # Clean empty lists
  GFEVDoutputs  <- GFEVDoutputs[unlist(lapply(GFEVDoutputs, length) != 0)]


  return(GFEVDoutputs)

}


######################################################################################################
######################################################################################################
####################### OUTPUTS FOR MODELS IN WHICH THE ESTIMATION ###################################
########################       IS DONE FOR ALL COUNTRIES JOINTLY      ################################
######################################################################################################
######################################################################################################
#'  Gathers all the model numerical ouputs after bootstrap for "joint Q" models
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param InputsForOutputs list conataining the desired inputs for the cunstruction of IRFs, GIRFs, FEVDs, and GFEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@keywords internal


OutputConstructionJoint_BS <- function(ModelType, ModelParaBoot, InputsForOutputs, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)
  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])

  # Output summary
  # IRF and FEVD
  IRFoutputs <- IRFjoint_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$IRF$horiz,
                            FactorLabels, Economies)
  FEVDoutputs <- FEVDjoint_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$FEVD$horiz,
                              FactorLabels, Economies)

  # GIRFS and GFEVD
  GIRFoutputs <- GIRFjoint_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$GIRF$horiz,
                              FactorLabels, Economies)
  GFEVDoutputs <- GFEVDjoint_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$GFEVD$horiz,
                                FactorLabels, Economies)

  # FOR JLL models: Non-orthogonalized IRFs
  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){

    # Generate the outputs with orthogonalized factros:
    IRFOrtho<- IRFjointOrthoJLL_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$IRF$horiz,
                                   FactorLabels, Economies)
    GIRFOrtho <- GIRFjointOrthoJLL_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$FEVD$horiz,
                                      FactorLabels, Economies)
    FEVDOrtho <- FEVDjointOrthogoJLL_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$GIRF$horiz,
                                        FactorLabels, Economies)
    GFEVDOrtho <- GFEVDjointOrthoJLL_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$GFEVD$horiz,
                                        FactorLabels, Economies)


    jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR multi") ]

    for (j in jointModelIdx){
      # Merge the lists of orthogonalized and non-orthogonalized factors
      # IRF
      for (tt in 1:ndraws){
        IRFoutputs[[ModelTypeSet[j]]][[tt]]$Factors <- list(IRFoutputs[[ModelTypeSet[j]]][[tt]]$Factors, IRFOrtho[[ModelTypeSet[j]]][[tt]]$Factors)
        IRFoutputs[[ModelTypeSet[j]]][[tt]]$Yields <- list(IRFoutputs[[ModelTypeSet[j]]][[tt]]$Yields, IRFOrtho[[ModelTypeSet[j]]][[tt]]$Yields)
        names(IRFoutputs[[ModelTypeSet[j]]][[tt]]$Factors) <- c("NonOrtho", "Ortho")
        names(IRFoutputs[[ModelTypeSet[j]]][[tt]]$Yields) <- c("NonOrtho", "Ortho")
        # GIRF
        GIRFoutputs[[ModelTypeSet[j]]][[tt]]$Factors <- list(GIRFoutputs[[ModelTypeSet[j]]][[tt]]$Factors, GIRFOrtho[[ModelTypeSet[j]]][[tt]]$Factors)
        GIRFoutputs[[ModelTypeSet[j]]][[tt]]$Yields <- list(GIRFoutputs[[ModelTypeSet[j]]][[tt]]$Yields, GIRFOrtho[[ModelTypeSet[j]]][[tt]]$Yields)
        names(GIRFoutputs[[ModelTypeSet[j]]][[tt]]$Factors) <- c("NonOrtho", "Ortho")
        names(GIRFoutputs[[ModelTypeSet[j]]][[tt]]$Yields) <- c("NonOrtho", "Ortho")
        # FEVD
        FEVDoutputs[[ModelTypeSet[j]]][[tt]]$Factors <- list(FEVDoutputs[[ModelTypeSet[j]]][[tt]]$Factors, FEVDOrtho[[ModelTypeSet[j]]][[tt]]$Factors)
        FEVDoutputs[[ModelTypeSet[j]]][[tt]]$Yields <- list(FEVDoutputs[[ModelTypeSet[j]]][[tt]]$Yields, FEVDOrtho[[ModelTypeSet[j]]][[tt]]$Yields)
        names(FEVDoutputs[[ModelTypeSet[j]]][[tt]]$Factors) <- c("NonOrtho", "Ortho")
        names(FEVDoutputs[[ModelTypeSet[j]]][[tt]]$Yields) <- c("NonOrtho", "Ortho")
        # GFEVD
        GFEVDoutputs[[ModelTypeSet[j]]][[tt]]$Factors <- list(GFEVDoutputs[[ModelTypeSet[j]]][[tt]]$Factors, GFEVDOrtho[[ModelTypeSet[j]]][[tt]]$Factors)
        GFEVDoutputs[[ModelTypeSet[j]]][[tt]]$Yields <- list(GFEVDoutputs[[ModelTypeSet[j]]][[tt]]$Yields, GFEVDOrtho[[ModelTypeSet[j]]][[tt]]$Yields)
        names(GFEVDoutputs[[ModelTypeSet[j]]][[tt]]$Factors) <- c("NonOrtho", "Ortho")
        names(GFEVDoutputs[[ModelTypeSet[j]]][[tt]]$Yields) <- c("NonOrtho", "Ortho")
      }

    }
  }



  NumericalOutputs <- list(IRFoutputs, FEVDoutputs, GIRFoutputs, GFEVDoutputs)
  names(NumericalOutputs) <- c("IRF", 'FEVD', "GIRF", "GFEVD")

  return(NumericalOutputs)

}



######################################################################################################
########################################### 2) IRFs ##################################################
######################################################################################################
#' IRFs after bootstrap for "joint Q" models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param IRFhoriz single numerical vector conataining the desired horizon of analysis for the IRFs
#'@param FactorLabels a string-list based which contains all the labels of all the variables present in the model
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'
#'
#'@keywords internal



IRFjoint_BS <- function(ModelType, ModelParaBoot, IRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  # Pre-allocation
  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  J <- length(ModelParaBoot$GeneralInputs$mat)
  CJ <- C*J
  K <- (M+N)*C + G

  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])

  IRFoutputs <- list()
  IRFoutputsAllDraws <- list()

  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR single") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  for(j in jointModelIdx){
    for (tt in 1:ndraws){

      # Generate factor labels depending on the models to be estimated
      AllCountryLabels <- c()
      AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountries)

      # Summarize inputs for the IRFs
      SIGMA <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$SSZ # KxK (variance-covariance matrix)
      A0 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$K0Z # Kx1 (matrix of intercepts)
      A1 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$K1Z # KxK (feedback matrix)

      BSpanned <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$rot$P$B
      B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

      # Initialization of IRFs of interest
      tempFactors <- array(0, c(K, K, IRFhoriz))
      tempYields  <- array(0, c(CJ,K,IRFhoriz))

      # Compute the IRFs
      if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){  # Choleski term
        S <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$Sigmas$Sigma_Y }
      else{   S <- t(chol(SIGMA)) }
      # Shock at t=0:
      tempFactors[ ,  , 1] <- S
      tempYields[ , , 1]  <- B %*% S
      # Shock at t=1:
      for (r in 2:IRFhoriz){
        if (r == 2) { A1h <- A1} else { A1h <- A1h %*% A1}

        tempFactors[ , , r] <- A1h%*%S # IRF (t+h) = A1^h*S
        tempYields[ , , r]      <- B%*%A1h%*%S
      }
      IRFRiskFactors <- aperm(tempFactors, c(3,1,2))
      IRFYields <- aperm(tempYields, c(3,1,2))

      Horiz <- t(t(0:(IRFhoriz-1))) #Add a column for horizon of interest


      # Adjust the variable labels
      dimnames(IRFRiskFactors)[[1]] <- Horiz
      dimnames(IRFRiskFactors)[[2]] <- AllFactorsLabels
      dimnames(IRFRiskFactors)[[3]] <- AllFactorsLabels


      YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$inputs$Y)
      dimnames(IRFYields)[[1]] <- Horiz
      dimnames(IRFYields)[[2]] <-  YieldsLabel
      dimnames(IRFYields)[[3]] <- AllFactorsLabels


      IRFoutputsAllCountries <- list(IRFRiskFactors,IRFYields)
      names(IRFoutputsAllCountries) <- c("Factors","Yields")

      # Prepare output per country
      IRFoutputsAllDraws[[tt]]  <- IRFoutputsAllCountries
      #


    }

    # All draws of all counntries for ALL models
    IRFoutputs[[j-L]] <- IRFoutputsAllDraws

  }



  # Clean empty lists
  IRFoutputs  <- IRFoutputs[unlist(lapply(IRFoutputs, length) != 0)]

  names(IRFoutputs) <- ModelTypeSet[jointModelIdx]


  return(IRFoutputs)

}


#########################################################################################################
################################### 3) FEVD #############################################################
#########################################################################################################
#' FEVDs after bootstrap for "joint Q" models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param FEVDhoriz single numerical vector conataining the desired horizon of analysis for the FEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@keywords internal


FEVDjoint_BS <- function(ModelType, ModelParaBoot, FEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic) - N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G
  J <- length(ModelParaBoot$GeneralInputs$mat)

  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR single") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  FEVDoutputs <- list()
  FEVDoutputsAllDraws <- list()

  for (j in  jointModelIdx){
    for (tt in 1: ndraws){

      G0 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$Gy.0
      Sigma_y <-  ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$SSZ
      F1 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$K1Z

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
      if ( any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma" ))){
        P <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$Sigmas$Sigma_Y
      }else{ P <- t(chol(invGSigmau)) }

      scale <- 1

      # 2.2) Factor loadings preparation
      BSpanned <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$rot$P$B
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
      YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$inputs$Y)

      dimnames(FEVDFactors)[[1]] <- 1:(FEVDhoriz)
      dimnames(FEVDFactors)[[2]] <- labelsFEVDFactors
      dimnames(FEVDFactors)[[3]] <- labelsFEVDFactors


      dimnames(FEVDYields)[[1]] <- 1:(FEVDhoriz)
      dimnames(FEVDYields)[[2]] <- labelsFEVDFactors
      dimnames(FEVDYields)[[3]] <- YieldsLabel


      FEVDoutputsAllCountries <- list(FEVDFactors,FEVDYields)
      names(FEVDoutputsAllCountries) <- c("Factors","Yields")

      FEVDoutputsAllDraws[[tt]] <- FEVDoutputsAllCountries


    }

    FEVDoutputs[[j-L]] <- FEVDoutputsAllDraws

  }
  # Clean empty lists
  FEVDoutputs  <- FEVDoutputs[unlist(lapply(FEVDoutputs, length) != 0)]

  names(FEVDoutputs) <- ModelTypeSet[jointModelIdx]


  return(FEVDoutputs)

}


#########################################################################################################
################################### 4) GIRF #############################################################
#########################################################################################################
#' GIRFs after bootstrap for "joint Q" models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param GIRFhoriz single numerical vector conataining the desired horizon of analysis for the GIRFs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@references
#' \itemize{
#' \item This function is a modified and extended version of the "irf" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#'
#'@keywords internal




GIRFjoint_BS <- function(ModelType, ModelParaBoot, GIRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])

  N <- length(FactorLabels$Spanned)
  C <- length(Economies) # Number of economies in the system
  M <- length(FactorLabels$Domestic) - N # Number of country-specific domestic variables
  G <- length(FactorLabels$Global) # Number of global variables
  K <- C*(M+N)+G # All factors of the system
  J <- length(ModelParaBoot$GeneralInputs$mat)
  CJ <- C*J

  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR single") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  GIRFoutputs <- list()
  GIRFoutputsAllDraws <- list()

  for (j in  jointModelIdx){
    for (tt in 1: ndraws){


      Gy.0 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$Gy.0
      Sigma.y <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$SSZ
      F1 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$K1Z

      BSpanned <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$rot$P$B
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
          numFactors <- (Ry.h[,,i]%*% solve(Gy.0)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
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
          numYields <- B%*%(Ry.h[,,i]%*% solve(Gy.0)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
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

      dimnames(GIRFFactors)[[1]] <- Horiz # We subtract 1, because the first element is the contemporaneous one
      dimnames(GIRFFactors)[[2]] <- labelsGIRF
      dimnames(GIRFFactors)[[3]] <- labelsGIRF

      YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$inputs$Y)
      dimnames(GIRFYields)[[1]] <- Horiz # We subtract 1, because the first element is the contemporaneous one
      dimnames(GIRFYields)[[2]] <- YieldsLabel
      dimnames(GIRFYields)[[3]] <- labelsGIRF


      GIRFoutputsAllCountries <- list(GIRFFactors,GIRFYields)
      names(GIRFoutputsAllCountries) <- c("Factors","Yields")

      GIRFoutputsAllDraws[[tt]] <- GIRFoutputsAllCountries

    }

    GIRFoutputs[[j-L]] <- GIRFoutputsAllDraws
  }
  # Clean empty lists
  GIRFoutputs  <- GIRFoutputs[unlist(lapply(GIRFoutputs, length) != 0)]

  names(GIRFoutputs) <- ModelTypeSet[jointModelIdx]

  return(GIRFoutputs)
}

#########################################################################################################
################################### 5) GFEVD #############################################################
#########################################################################################################
#' GFEVDs after bootstrap for "joint Q" models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot List of model parameter estimates (See the "Optimization" function) after a bootstrap draw
#'@param GFEVDhoriz single numerical vector conataining the desired horizon of analysis for the GFEVDs
#'@param FactorLabels string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
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


GFEVDjoint_BS <- function(ModelType, ModelParaBoot, GFEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic) - N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G
  J <- length(ModelParaBoot$GeneralInputs$mat)

  jointModelIdx <- idxWishModels[idxWishModels > which(ModelTypeSet == "GVAR single") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  GFEVDoutputs <- list()
  GFEVDoutputsAllDraws <- list()

  for (j in  jointModelIdx){
    for (tt in 1: ndraws){

      G0 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$Gy.0
      Sigma_y <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$SSZ
      F1 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$K1Z

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
      BSpanned <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$rot$P$B
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
      YieldsLabel <- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$inputs$Y)


      dimnames(GFEVDFactorsNormalized)[[1]] <- 1:(GFEVDhoriz)
      dimnames(GFEVDFactorsNormalized)[[2]] <- labelsGFEVDFactors
      dimnames(GFEVDFactorsNormalized)[[3]] <- labelsGFEVDFactors

      dimnames(GFEVDYieldsNormalized)[[1]] <- 1:(GFEVDhoriz)
      dimnames(GFEVDYieldsNormalized)[[2]] <- labelsGFEVDFactors
      dimnames(GFEVDYieldsNormalized)[[3]] <- YieldsLabel

      GFEVDoutputsAllCountries <- list(GFEVDFactorsNormalized,GFEVDYieldsNormalized)
      names(GFEVDoutputsAllCountries) <- c("Factors","Yields")

      GFEVDoutputsAllDraws[[tt]] <- GFEVDoutputsAllCountries

    }

    GFEVDoutputs[[j-L]] <- GFEVDoutputsAllDraws
  }


  # Clean empty lists
  GFEVDoutputs  <- GFEVDoutputs[unlist(lapply(GFEVDoutputs, length) != 0)]

  names(GFEVDoutputs) <- ModelTypeSet[jointModelIdx]


  return(GFEVDoutputs)

}


######################################################################################################
############################ 6) IRFs with orthogonalized factors #####################################
######################################################################################################
#' IRFs after bootstrap for JLL-based models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param IRFhoriz single numerical vector conataining the desired horizon of analysis for the IRFs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'
#'@keywords internal


IRFjointOrthoJLL_BS <- function(ModelType, ModelParaBoot, IRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])

  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  K <- C*(M+N) + G
  J <- length(ModelParaBoot$GeneralInputs$mat)
  CJ <- C*J


  # Pre-allocation
  IRFoutputs <- list()


  jointModelIdxJLL <- idxWishModels[idxWishModels >= which(ModelTypeSet == "JLL original") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  for(j in jointModelIdxJLL){
    IRFAllDraws <- list()
    for (tt in 1:ndraws){

      # Generate factor labels depending on the models to be estimated
      AllCountryLabels <- c()
      AllFactorsLabels <-  c(FactorLabels$Global, FactorLabels$Tables$AllCountriesJLL)

      # Summarize inputs for the IRFs
      SIGMAe <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$Sigmas$VarCov_Ortho # KxK (variance-covariance matrix)
      A0e <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$k0_e # Kx1 (matrix of intercepts)
      A1e <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$k1_e # KxK (feedback matrix)
      PI <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$PI

      BSpanned <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$rot$P$B
      B <- BUnspannedAdapJoint(G,M,N,C, J, BSpanned)

      # Initialization of IRFs of interest
      tempFactors <- array(0, c(K, K, IRFhoriz))
      tempYields  <- array(0, c(CJ,K,IRFhoriz))

      # Compute the IRFs
      Se <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$Sigmas$Sigma_Ye # Choleski term
      # Shock at t=0:
      tempFactors[ ,  , 1] <- Se
      tempYields[ , , 1]  <- B%*%PI%*% Se
      # Shock at t=1:
      for (r in 2:IRFhoriz){
        if (r == 2) { A1h <- A1e} else { A1h <- A1h %*% A1e}

        tempFactors[ , , r] <- A1h%*%Se # IRF (t+h) = A1^h*S
        tempYields[ , , r]      <- B%*%PI%*%A1h%*%Se
      }
      IRFRiskFactors <- aperm(tempFactors, c(3,1,2))
      IRFYields <- aperm(tempYields, c(3,1,2))

      Horiz <- t(t(0:(IRFhoriz-1))) #Add a column for horizon of interest


      # Adjust the variable labels
      dimnames(IRFRiskFactors)[[1]] <- Horiz
      dimnames(IRFRiskFactors)[[2]] <- AllFactorsLabels
      dimnames(IRFRiskFactors)[[3]] <- AllFactorsLabels

      YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$inputs$Y)
      dimnames(IRFYields)[[1]] <- Horiz
      dimnames(IRFYields)[[2]] <- YieldsLabel
      dimnames(IRFYields)[[3]] <- AllFactorsLabels


      IRFoutputsAllCountries <- list(IRFRiskFactors,IRFYields)
      names(IRFoutputsAllCountries) <- c("Factors","Yields")

      # Prepare output per country
      IRFAllDraws[[tt]] <- IRFoutputsAllCountries

    }

    IRFoutputs[[j-L]] <- IRFAllDraws
  }
  # Clean empty lists
  IRFoutputs  <- IRFoutputs[unlist(lapply(IRFoutputs, length) != 0)]

  names(IRFoutputs) <- ModelTypeSet[jointModelIdxJLL]


  return(IRFoutputs)

}



#########################################################################################################
################################### 7) FEVD with orthogonalized factors #################################
#########################################################################################################
#' FEVDs after bootstrap for JLL-based models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param FEVDhoriz single numerical vector conataining the desired horizon of analysis for the FEVDs
#'@param FactorLabels string-list based which contains all the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@keywords internal


FEVDjointOrthogoJLL_BS <- function(ModelType, ModelParaBoot, FEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])


  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic) - N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G
  J <- length(ModelParaBoot$GeneralInputs$mat)


  jointModelIdxJLL <- idxWishModels[idxWishModels >= which(ModelTypeSet == "JLL original") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  FEVDoutputs <- list()
  FEVDAllDraws <- list()

  for (j in  jointModelIdxJLL){
    for (tt in 1:ndraws){

      G0 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$Gy.0
      Sigma_y <-  ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$Sigmas$VarCov_Ortho
      F1e <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$k1_e
      PI <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$PI

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

      P <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$Sigmas$Sigma_Ye
      scale <- 1

      # 2.2) Factor loadings preparation
      BSpanned <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$rot$P$B
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
      YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$inputs$Y)

      dimnames(FEVDFactors)[[1]] <- 1:(FEVDhoriz) # We subtract 1, because the first element is the contemporaneous one
      dimnames(FEVDFactors)[[2]] <- labelsFEVDFactors
      dimnames(FEVDFactors)[[3]] <- labelsFEVDFactors


      dimnames(FEVDYields)[[1]] <- 1:(FEVDhoriz) # We subtract 1, because the first element is the contemporaneous one
      dimnames(FEVDYields)[[2]] <- labelsFEVDFactors
      dimnames(FEVDYields)[[3]] <- YieldsLabel


      FEVDoutputsAllCountries <- list(FEVDFactors,FEVDYields)
      names(FEVDoutputsAllCountries) <- c("Factors","Yields")

      FEVDAllDraws[[tt]] <- FEVDoutputsAllCountries


    }

    FEVDoutputs[[j-L]] <- FEVDAllDraws
  }
  # Clean empty lists
  FEVDoutputs  <- FEVDoutputs[unlist(lapply(FEVDoutputs, length) != 0)]

  names(FEVDoutputs) <- ModelTypeSet[jointModelIdxJLL]


  return(FEVDoutputs)

}

#########################################################################################################
################################### 8) GIRF With orthogonalized factors #################################
#########################################################################################################
#' GIRFs after bootstrap for JLL-based models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param GIRFhoriz single numerical vector conataining the desired horizon of analysis for the GIRFs
#'@param FactorLabels  string-list based which contains the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@references
#' \itemize{
#' \item This function is a modified and extended version of the "irf" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#'@keywords internal


GIRFjointOrthoJLL_BS <- function(ModelType, ModelParaBoot, GIRFhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  C <- length(Economies) # Number of economies in the system
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N # Number of country-specific domestic variables
  G <- length(FactorLabels$Global) # Number of global variables
  K <- C*(M+N)+G # All factors of the system
  J <- length(ModelParaBoot$GeneralInputs$mat)
  CJ <- C*J

  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])

  jointModelIdxJLL <- idxWishModels[idxWishModels >= which(ModelTypeSet == "JLL original") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  GIRFoutputs <- list()

  for (j in  jointModelIdxJLL){

    GIRFallDraws <- list()

    for (tt in 1:ndraws){

      Gy.0 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$Gy.0
      Sigma.y <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$Sigmas$VarCov_Ortho
      F1e <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$k1_e
      PI <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$PI

      BSpanned <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$rot$P$B
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
          numFactors <- (PI%*%Ry.h[,,i]%*% solve(Gy.0)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
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
          numYields <- B%*%(PI%*%Ry.h[,,i]%*% solve(Gy.0)%*%Sigma.y%*%ey.j[,g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
          demYields <- 1/sqrt((t(ey.j[,g])%*%Sigma.y%*%ey.j[,g])) # denominator from equation at the bottom of the page 22 (PS, 1998)
          AllResponsesToShockOfOneVariableYields[,i] <- numYields*drop(demYields)
        }
        AllResponsesToAllShocksYields[,,g] <- AllResponsesToShockOfOneVariableYields
      }

      GIRFYields <- aperm(AllResponsesToAllShocksYields, c(2,1,3))


      # 4.2) Labels
      labelsGIRF <- c(FactorLabels$Global,FactorLabels$Tables$AllCountries)

      dimnames(GIRFFactors)[[1]] <- 0:(GIRFhoriz-1) # We subtract 1, because the first element is the contemporaneous one
      dimnames(GIRFFactors)[[2]] <- labelsGIRF
      dimnames(GIRFFactors)[[3]] <- labelsGIRF

      YieldsLabel<- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$inputs$Y)
      dimnames(GIRFYields)[[1]] <- 0:(GIRFhoriz-1) # We subtract 1, because the first element is the contemporaneous one
      dimnames(GIRFYields)[[2]] <- YieldsLabel
      dimnames(GIRFYields)[[3]] <- labelsGIRF


      GIRFoutputsAllCountries <- list(GIRFFactors,GIRFYields)
      names(GIRFoutputsAllCountries) <- c("Factors","Yields")

      GIRFallDraws[[tt]] <- GIRFoutputsAllCountries

    }

    GIRFoutputs[[j-L]] <- GIRFallDraws
  }
  # Clean empty lists
  GIRFoutputs  <- GIRFoutputs[unlist(lapply(GIRFoutputs, length) != 0)]

  names(GIRFoutputs) <- ModelTypeSet[jointModelIdxJLL]

  return(GIRFoutputs)
}

#########################################################################################################
################################### 9) GFEVD With orthogonalized factor s################################
#########################################################################################################
#' GFEVDs after bootstrap for JLL-based models

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#'@param GFEVDhoriz single numerical vector conataining the desired horizon of analysis for the GFEVDs
#'@param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
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


GFEVDjointOrthoJLL_BS <- function(ModelType, ModelParaBoot, GFEVDhoriz, FactorLabels, Economies){

  ModelTypeSet <- c("JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original",
                    "JLL No DomUnit", "JLL joint Sigma")
  idxWishModels <- which(ModelTypeSet == ModelType)

  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)
  M <- length(FactorLabels$Domestic)-N
  G <- length(FactorLabels$Global)
  K <- C*(M+N) + G
  J <- length(ModelParaBoot$GeneralInputs$mat)
  CJ <- C*J

  jointModelIdxJLL <- idxWishModels[idxWishModels >= which(ModelTypeSet == "JLL original") ]
  L <- which(ModelTypeSet == "GVAR single") # index of the model right before the first desired joint model.

  GFEVDoutputs <- list()

  for (j in  jointModelIdxJLL){

    GFEVDAllDraws <- list()

    for (tt in 1:ndraws){

      G0 <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$Gy.0
      Sigma_y <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$Sigmas$VarCov_Ortho
      F1e <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$k1_e
      PI <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$ests$JLLoutcomes$PI

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
      BSpanned <- ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$rot$P$B
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
      YieldsLabel <- rownames(ModelParaBoot$ParaDraws[[ModelTypeSet[[j]]]][[tt]]$inputs$Y)


      dimnames(GFEVDFactorsNormalized)[[1]] <- 1:(GFEVDhoriz) # We subtract 1, because the first element is the contemporaneous one
      dimnames(GFEVDFactorsNormalized)[[2]] <- labelsGFEVDFactors
      dimnames(GFEVDFactorsNormalized)[[3]] <- labelsGFEVDFactors

      dimnames(GFEVDYieldsNormalized)[[1]] <- 1:(GFEVDhoriz) # We subtract 1, because the first element is the contemporaneous one
      dimnames(GFEVDYieldsNormalized)[[2]] <- labelsGFEVDFactors
      dimnames(GFEVDYieldsNormalized)[[3]] <- YieldsLabel

      GFEVDoutputsAllCountries <- list(GFEVDFactorsNormalized,GFEVDYieldsNormalized)
      names(GFEVDoutputsAllCountries) <- c("Factors","Yields")

      GFEVDAllDraws[[tt]] <- GFEVDoutputsAllCountries

    }

    GFEVDoutputs[[j-L]] <- GFEVDAllDraws
  }
  # Clean empty lists
  GFEVDoutputs  <- GFEVDoutputs[unlist(lapply(GFEVDoutputs, length) != 0)]

  names(GFEVDoutputs) <- ModelTypeSet[jointModelIdxJLL]


  return(GFEVDoutputs)

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

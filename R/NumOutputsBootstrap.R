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
  IRFout <- IRFandGIRF_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)
  FEVDoutputs <- FEVDsep_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels, Economies)

  # GIRF and GFEVD
  GFEVDoutputs <- GFEVDsep_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$GFEVD$horiz, FactorLabels, Economies)

  NumericalOutputs <- list(IRF = IRFout$IRFs, FEVD = FEVDoutputs, GIRF = IRFout$GIRFs, GFEVD = GFEVDoutputs)

  return(NumericalOutputs)

}



######################################################################################################
########################################### 2) IRFs ##################################################
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
      J <- numel(ModelParaBoot$GeneralInputs$mat)

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
    J <- numel(ModelParaBoot$GeneralInputs$mat)
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


  ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])

  # Output summary
  # IRF and FEVD
  IRFout <- IRFandGIRF_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$IRF$horiz,
                            FactorLabels, Economies)
  FEVDoutputs <- FEVDjoint_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$FEVD$horiz,
                              FactorLabels, Economies)

  # GIRFS and GFEVD
  GFEVDoutputs <- GFEVDjoint_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$GFEVD$horiz,
                                FactorLabels, Economies)

  # FOR JLL models: Non-orthogonalized IRFs
  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){

    FEVDOrtho <- FEVDjointOrthogoJLL_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$GIRF$horiz,
                                        FactorLabels, Economies)
    GFEVDOrtho <- GFEVDjointOrthoJLL_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$GFEVD$horiz,
                                        FactorLabels, Economies)


      for (tt in 1:ndraws){
        # FEVD
        FEVDoutputs[[ModelType]][[tt]]$Factors <- list(FEVDoutputs[[ModelType]][[tt]]$Factors, FEVDOrtho[[ModelType]][[tt]]$Factors)
        FEVDoutputs[[ModelType]][[tt]]$Yields <- list(FEVDoutputs[[ModelType]][[tt]]$Yields, FEVDOrtho[[ModelType]][[tt]]$Yields)
        names(FEVDoutputs[[ModelType]][[tt]]$Factors) <- c("NonOrtho", "Ortho")
        names(FEVDoutputs[[ModelType]][[tt]]$Yields) <- c("NonOrtho", "Ortho")
        # GFEVD
        GFEVDoutputs[[ModelType]][[tt]]$Factors <- list(GFEVDoutputs[[ModelType]][[tt]]$Factors, GFEVDOrtho[[ModelType]][[tt]]$Factors)
        GFEVDoutputs[[ModelType]][[tt]]$Yields <- list(GFEVDoutputs[[ModelType]][[tt]]$Yields, GFEVDOrtho[[ModelType]][[tt]]$Yields)
        names(GFEVDoutputs[[ModelType]][[tt]]$Factors) <- c("NonOrtho", "Ortho")
        names(GFEVDoutputs[[ModelType]][[tt]]$Yields) <- c("NonOrtho", "Ortho")
      }


  }


  NumericalOutputs <- list(IRF = IRFout$IRFs, FEVD = FEVDoutputs, GIRF = IRFout$GIRFs, GFEVD = GFEVDoutputs)

  return(NumericalOutputs)

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

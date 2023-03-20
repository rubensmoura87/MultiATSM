#' Compute some key parameters from the P-dynamics (Bootstrap set)
#'
#'@param ModelType           string-vector containing the label of the model to be estimated
#'@param AllFactorsUnderP   complete set of factors that may be used in the estimation of P (KxT)
#'@param FactorLabels       string-list based which contains the labels of all the variables present in the model (see "LabFac" function)
#'@param Economies          string-vector containing the names of the economies which are part of the economic system
#'@param JLLinputs         List containing the necessary inputs for the estimation of the JLL-based models (see "JLL" function). Default is set to NULL.
#'@param GVARinputs        List containing the necessary inputs for the estimation of the GVAR-based models (see "GVAR" function). Default is set to NULL.
#'
#'


PdynamicsSet_BS<- function(ModelType, AllFactorsUnderP, FactorLabels, Economies, JLLinputs = NULL, GVARinputs= NULL){

   T <- ncol(AllFactorsUnderP)
   N <- length(FactorLabels$Spanned)

  # JPS-related models
  if (ModelType == 'JPS' || ModelType == 'JPS jointP' || ModelType == 'VAR jointQ'){
    K <- nrow(AllFactorsUnderP)
    VAR <-stats::lm( t(AllFactorsUnderP[,2:T])~ t(AllFactorsUnderP[,1:(T-1)])) # VAR(1) under the P.
    K0Z <- t(t(VAR$coefficients[1,]))
    K1Z <- t(VAR$coefficients[2:(K+1),])
    eZ<- VAR$residuals # (T-1) x K
  }


  # GVAR-related models
  if (ModelType ==  'GVAR sepQ'|| ModelType == 'GVAR jointQ'){
    GVARpara <- GVAR(GVARinputs, N)
    K0Z <- GVARpara$F0
    K1Z <- GVARpara$F1

    eZ <- AllFactorsUnderP[,2:T] - matrix(K0Z, nrow= nrow(K0Z),  ncol= T-1) - K1Z%*%AllFactorsUnderP[,1:(T-1)] # T x K
    eZ <- t(eZ)
  }

  # JLL-related models
  if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma"){
    JLLinputs$WishSigmas <- 0
    JLLPara <- JLL(AllFactorsUnderP, N, JLLinputs)
    K0Z <-JLLPara$k0
    K1Z <-JLLPara$k1
    eZ <- AllFactorsUnderP[,2:T] - matrix(K0Z, nrow= nrow(K0Z),  ncol= T-1) - K1Z%*%AllFactorsUnderP[,1:(T-1)] # T x K
    eZ <- t(eZ)
  }


  Outputs <- list(K0Z, K1Z, eZ)
  names(Outputs) <- c("K0Z", "K1Z", "eZ")

  return(Outputs)
}

####################################################################################################
#' Prepare the factor set for GVAR models  (Bootstrap version)
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param RiskFactors Complete set of risk factors (KxT)
#' @param Wgvar  transition matrix from GVAR models (CxC)
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#'


DataSet_BS <- function(ModelType, RiskFactors, Wgvar, Economies, FactorLabels){


  if (ModelType == "GVAR sepQ" || ModelType == "GVAR jointQ"){

  # 1) Pre-allocate list of factors
  T <- ncol(RiskFactors) # length of model's time dimension
  C <- length(Economies) # number of economies in of the economic system
  N <- length(FactorLabels$Spanned) # number of countrey-specific spanned factors
  M <- length(FactorLabels$Domestic) - N # Number of country-specific macro variables
  M.star <- length(FactorLabels$Star) - N # Number of foreign-country-specific macro variables
  G <- length(FactorLabels$Global) # Number of global variables

  ListFactors <- vector(mode='list', length = length(Economies)+1) # length = all countries + global factors
  names(ListFactors) <- c(Economies, 'Global')

  # Country-specific factors (CSF)
  CSF <- vector(mode='list', length = length(FactorLabels$Domestic))
  names(CSF) <- FactorLabels$Domestic
  for (i in 1:C){  ListFactors[[Economies[i]]] <- CSF }

  #  Star factors (SF)
  SF <- vector(mode='list', length = length(FactorLabels$Star))
  names(SF) <- FactorLabels$Star
  for (i in 1:length(Economies)){
    ListFactors[[Economies[i]]] <- list(append(CSF,SF))
    names(ListFactors[[Economies[i]]]) <- 'Factors'
  }

  # Global Factors (GF)
  GF <- vector(mode='list', length = length(FactorLabels$Global))
  names(GF) <- FactorLabels$Global
  ListFactors[[ length(Economies)+1 ]] <- GF

  # Yields
  YieldsSeries <- vector(mode='list', length = C)
  Wpca <- vector(mode='list', length = C)
  names(Wpca) <- rep("Wpca", times=C)
  YieldsList <- vector(mode='list', length = C)



  # 2) Fill in list with the corresponding factors
  # A) Country-specific variables (economy-related variables)

  for (i in 1:C) {
    for (j in 1:M){
      ListFactors[[Economies[i]]]$Factors[[j]]<- as.matrix(RiskFactors[(c(FactorLabels$Tables[[Economies[i]]][j])),])
    }
  }

  # B) Country-specific variables (pricing-related variables)
  idx0 <- M
  for (i in 1:C) {
    for (j in 1:N){
      ListFactors[[Economies[i]]]$Factors[[idx0+j]] <- as.matrix(RiskFactors[(c(FactorLabels$Tables[[Economies[i]]][idx0+j])),])
    }
  }

  # C) Foreign country-specific variables (economy and pricing-related)
  idx1 <- M+N
  Z <- list()

  for (j in 1:(M+N)){
    X <- matrix(NA, nrow= C, ncol=T)
    for (i in 1:C){
      X[i,] <- ListFactors[[Economies[i]]]$Factors[[j]]
      Z[[j]] <- X # Each element of the list contains the same country-specific variable of all countries
    }
  }
  names(Z) <- FactorLabels$Domestic

  if ("GVAR sepQ" %in% ModelType || "GVAR jointQ" %in% ModelType ){
    for (i in 1:C){
      for (j in 1:(M+N)){
        ListFactors[[Economies[i]]]$Factors[[idx1+j]] <- t(Wgvar[i,]%*%Z[[j]])
      }
    }
  }

  # D) Global Factors
  for (i in seqi(1,G)){
    ListFactors[[length(Economies)+1]][[i]] <-as.matrix( RiskFactors[(c(FactorLabels$Global[i])),])
  }

  } else{
    ListFactors <- NULL
    }


  return(ListFactors)
}

#' Generates several  inputs that are necessary to build the likelihood function
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param Yields time series of yields (JxT or CJ x T)
#'@param PdynamicsFactors time series of the risk factors (K x T)
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param mat  vector of maturities (in years) used in the estimation
#'@param Economies  string-vector containing the names of the economies of the system. \cr
#'                 If the ModelType selected is "JPS", "JPS jointP", "GVAR sepQ", then only one economy can be selected. \cr
#'                  For the other models, more than one economy must be selected.
#'@param DataFrequency  character-based-vector. Avaialable options are: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"
#'@param JLLinputs    list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'@param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#'
#'
#'@importFrom pracma null
#'
#'@return List of necessary inputs  for constructing the model's log-likelihood function
#'@examples
#'\donttest{
#' # Example 1:
#' data(CM_Factors)
#' data(CM_Yields)
#'
#' ModelType <- "JPS"
#' Economies <- "Mexico"
#' Factors <- RiskFactors
#' N <- 3
#' GlobalVar <- c("GBC", "CPI_OECD") # Global Variables
#' DomVar <- c("Eco_Act", "Inflation") # Domestic Variables
#' FactorLabels <- LabFac(N, DomVar,GlobalVar, Economies, ModelType)
#'
#' mat <- c(0.25, 0.5, 1, 3, 5, 10)
#' DataFrequency <- "Monthly"
#'
#' i <- length(Economies)
#' ATSMInputs <- InputsForMLEdensity(ModelType, Yields, Factors, FactorLabels, mat,
#'                                  Economies, DataFrequency)


#'# Example 2:
#' data(CM_Factors)
#' data(CM_Yields)
#' data(CM_Factors_GVAR)
#'
#' ModelType <- "GVAR jointQ"

#' Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#' mat <- c(0.25, 0.5, 1, 3, 5, 10)
#' DataFrequency <- "Monthly"
#' Factors  <- RiskFactors
#' N <- 3
#' GlobalVar <- c("GBC", "CPI_OECD") # Global Variables
#' DomVar <- c("Eco_Act", "Inflation") # Domestic Variables
#' FactorLabels <- LabFac(N, DomVar,GlobalVar, Economies, ModelType)

#'
#' GVARinputs <- list()
#' GVARinputs$Economies <- Economies
#' GVARinputs$GVARFactors <- FactorsGVAR
#' GVARinputs$VARXtype <- "unconstrained"
#' GVARinputs$Wgvar <- matrix( c(0, 0.83, 0.86, 0.38,
#'                0.65, 0, 0.13, 0.55,
#'                0.32, 0.12, 0, 0.07,
#'                0.03, 0.05, 0.01, 0), nrow = 4, ncol = 4)

#'ATSMInputs <- InputsForMLEdensity(ModelType, Yields, Factors, FactorLabels, mat, Economies,
#'                                  DataFrequency, JLLinputs= NULL , GVARinputs)
#'
#' # Example 3:
#' if (requireNamespace('neldermead', quietly = TRUE)) {
#'
#' data(CM_Factors)
#' data(CM_Yields)
#' ModelType <- "JLL jointSigma"
#' GlobalVar <- c("GBC", "CPI_OECD") # Global Variables
#' DomVar <- c("Eco_Act", "Inflation") # Domestic Variables
#' N <- 3
#' Economies <- c( "China", "Brazil", "Mexico", "Uruguay")
#' FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#'
#' Factors <- RiskFactors
#' mat <- c(0.25, 0.5, 1, 3, 5, 10)
#' DataFrequency <- "Monthly"

#' JLLinputs <- list()
#' JLLinputs$Economies <- Economies
#' JLLinputs$DomUnit <- "China"
#' JLLinputs$WishSigmas <- 1
#' JLLinputs$SigmaNonOrtho <- NULL
#' JLLinputs$JLLModelType <- ModelType
#'
#' ATSMInputs <- InputsForMLEdensity(ModelType, Yields, Factors, FactorLabels, mat, Economies,
#'                                  DataFrequency, JLLinputs)
#'} else {
#'  message("skipping functionality due to missing Suggested dependency")
#'}
#'
#'



#'}
#' @details
#' To ensure that the risk factors matrix is correctly built for the model "JPS", the global factors should be allocated on the first G rows of this matrix.
#'@export




InputsForMLEdensity <- function(ModelType, Yields, PdynamicsFactors, FactorLabels,  mat, Economies, DataFrequency,
                                JLLinputs = NULL, GVARinputs  = NULL){

# Check whether the model choice is compatible with the number of countries selected
if (length(Economies) == 1 &
    ( ModelType == "GVAR jointQ" || ModelType == "VAR jointQ" || ModelType == "JLL original"
      || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma" )){

  stop("The models 'GVAR jointQ', 'VAR jointQ', 'JLL original', 'JLL NoDomUnit', and 'JLL jointSigma'
       require the estimation of several countries")
}


  # for the cases in which the estimation is done on a country-by-country basis
  if (ModelType == 'JPS' || ModelType == 'JPS jointP' || ModelType == "GVAR sepQ"){
    i <- get("i", globalenv())
    Economies <- Economies[i]
  }


  # Pre-allocation
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  C <- length(Economies)
  J <- nrow(Yields)/C

  if (C==1){ # For the cases in which the model estimation is done on a country-by-country basis
  IdxCountry <- which(grepl(Economies, rownames(Yields)))
  J <- length(IdxCountry)
  }

if (DataFrequency == "Daily All Days"){ dt <- 1/365}
if (DataFrequency == "Daily Business Days"){ dt <- 1/252}
if (DataFrequency == "Weekly"){ dt <- 1/52}
if (DataFrequency == "Monthly"){ dt <- 1/12}
if (DataFrequency == "Quartely"){ dt <- 1/4}
if (DataFrequency == "Annually"){ dt <- 1}

    idxJ0 <- 0
    idxN0 <- 0

    if (C == 1){idxJ0 <- IdxCountry[1]-1 }

    # Compute the inputs
    for (i in 1:C){ # Country-specific inputs
      idxJ1 <- idxJ0 + J
      idxN1 <- idxN0 + N

       # Spanned Factors (N x T)
      YCS <- Yields[(idxJ0+1):idxJ1,] # Yields (JxT)
      WpcaCS <- 100*pca_weights_one_country(YCS, Economy= Economies[i])[1:N,] # matrix of weigths for the portfolio without errors (N x J)
      WeCS <- t(null(WpcaCS)) # matrix of weigths for the yield portfolios priced with errors
      WpcaFullCS <- rbind(WpcaCS, WeCS)
      PPCS <- Spanned_Factors(YCS, Economies = Economies[i], N)
      K1XQCS <-  Reg_K1Q(YCS, mat, PPCS, dt, type="Jordan")

      if (i==1){
        K1XQ <- K1XQCS
        Wpca <- WpcaCS
        WpcaFull <- WpcaFullCS
        We <- WeCS
        PP <- PPCS
        Y <- YCS
      }else{
        K1XQ <- magic::adiag(K1XQ,K1XQCS)
        Wpca <- magic::adiag(Wpca,WpcaCS)
        We <- magic::adiag(We,WeCS)
        WpcaFull <- magic::adiag(WpcaFull, WpcaFullCS)
        PP <- rbind(PP, PPCS)
        Y <- rbind(Y, YCS)
        }
      idxJ0<- idxJ1
      idxN0 <- idxN1
    }

# Time series of the risk factors
  if (ModelType == 'JPS'){
    idxCountryFactors <- which(grepl(Economies, rownames(PdynamicsFactors))) # index of the rows of the country-specific factors
    idxFactors <- c(seqi(1,G), idxCountryFactors) # index global + country-specific factors
    ZZ <- PdynamicsFactors[idxFactors,]
  }else{
    ZZ <- PdynamicsFactors
  }

  # Time series labels
    colnames(Y) <- colnames(ZZ)

# Matrix of contemporaneous terms
  K <- nrow(ZZ)
  Gy.0 <- diag(K)

# Parameters of the of the P-dynamics
  # JPS-related models
  if (ModelType == 'JPS' || ModelType == 'JPS jointP' || ModelType == 'VAR jointQ'){
  VARpara <- VAR(ZZ, VARtype= 'unconstrained', Bcon = NULL)
  K0Z <- VARpara$K0Z
  K1Z <- VARpara$K1Z
  SSZ <- VARpara$SSZ
  }

  # GVAR-related models
  if (ModelType == 'GVAR sepQ'|| ModelType == 'GVAR jointQ'){
    GVARpara <- GVAR(GVARinputs, N)
    K0Z <- GVARpara$F0
    K1Z <- GVARpara$F1
    SSZ <- GVARpara$Sigma_y
  }

  # JLL-related models
  if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma"){
    JLLinputs$WishSigmas <- 1
    JLLPara <- JLL(ZZ, N, JLLinputs)
    K0Z <-JLLPara$k0
    K1Z <-JLLPara$k1
    JLLinputs$WishSigmas <- 0 # Ensures that the variance-covariance matrix will no longer be estimated within the JLL function

    if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit"){  SSZ <- JLLPara$Sigmas$VarCov_NonOrtho }
    if (ModelType == "JLL jointSigma"){ SSZ <- JLLPara$Sigmas$VarCov_Ortho
    # NOTE: the required vectorization that preceeds the numerical optimization of the variance-covariance matrix
    # is done on the orthogonalized variance-covariance matrix because we want to preserve the restrictions imposed in this matrix.
    # However, to compute the loadings A and B from Y= A + B*P we have to use the non-orthogonalized variance-covariance
    # matrix. We make this adjustment in the function 'A0N_MLEdensity_WOE...' by redefining SSZ before computing A and B.
    # (this avoids overcomplicating the overall structure of the code")
    }

    }


  # Outputs
  Output <- list(Wpca, We, WpcaFull, Y, PP, Gy.0, K1XQ, ZZ, SSZ, K0Z, K1Z, JLLinputs, GVARinputs)
  names(Output) <- c("Wpca","We", "WpcaFull", "Y","PP", "Gy.0","K1XQ", "ZZ", "SSZ", "K0Z", "K1Z", "JLLinputs", "GVARinputs")

  return(Output)
}


#####################################################################################################################
#' Generates several  inputs that are necessary to build the likelihood function - Bootstrap version
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param Y_artificial time series of yields (CJ x T or JxT)
#'@param Z_artificial time series of the risk factors (F x T)
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param mat  vector of maturities (in years) used in the estimation
#'@param Economies    string-vector containing the names of the economies of the system. \cr
#'                 If the ModelType selected is "JPS", "JPS jointP", "GVAR sepQ", then only one economy can be selected. \cr
#'                  For the other models, more than one economy must be selected.
#'@param DataFrequency  character-based-vector. Avaialable options are: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"
#'@param JLLinputs    list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'@param GVARinputs   list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#'
#'@importFrom pracma null
#'




InputsForMLEdensity_BS <- function(ModelType, Y_artificial, Z_artificial, FactorLabels, mat,
                                   Economies, DataFrequency, JLLinputs = NULL, GVARinputs= NULL){

  # for the cases in which the estimation is done on a country-by-country basis
  if (ModelType == 'JPS' || ModelType == 'JPS jointP' || ModelType == "GVAR sepQ"){
    i <- get("i", globalenv())
    Economies <- Economies[i]
  }

  J <- length(mat)
  C <- length(Economies)
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N

  if (DataFrequency == "Daily All Days"){ dt <- 1/365}
  if (DataFrequency == "Daily Business Days"){ dt <- 1/252}
  if (DataFrequency == "Weekly"){ dt <- 1/52}
  if (DataFrequency == "Monthly"){ dt <- 1/12}
  if (DataFrequency == "Quartely"){ dt <- 1/4}
  if (DataFrequency == "Annually"){ dt <- 1}


  # Set/Initialize artificial variables

  ZZ_artificial <- t(Z_artificial)
  YY_artificial <- t(Y_artificial)

  # Set artificial parameter values
  idxJ0 <- 0
  idxN0 <- 0

  for ( i in 1:C){
    idxJ1 <- idxJ0 + J
    idxN1 <- idxN0 + N

    # Country-specific variables
    P_CS <- ZZ_artificial[(c(FactorLabels$Tables[[Economies[i]]] [(M+1):(M+N)])),]
    Y_CS <-  YY_artificial[(idxJ0+1):idxJ1,]


    K1XQ__CS <-  Reg_K1Q(Y_CS, mat, P_CS, dt, type="Jordan")

    Wpca_CS <- mrdivide(P_CS,Y_CS) # To ensure that we will have reasonable estimates for Wpca, we use the real data for PP

    ######################################################
    # Ensures that the weghts can be properly interpreted
    if (all(Wpca_CS[1,] < 0)){   Wpca_CS[1,] <- Wpca_CS[1,]*(-1)   }
    if (N ==2){   if (Wpca_CS[2,1] > Wpca_CS[2,J]){   Wpca_CS[2,] <- Wpca_CS[2,]*(-1)   }}
    if (N > 3){   if (Wpca_CS[2,1] > Wpca_CS[2,J]){   Wpca_CS[2,] <- Wpca_CS[2,]*(-1)}
      Med <- round(stats::median(1:J))
      if (Wpca_CS[3,1] > Wpca_CS[3,Med] & Wpca_CS[3,J] > Wpca_CS[3,Med]){ Wpca_CS[3,] <- Wpca_CS[3,]*(-1) } }
    ####################################################
    We_CS  <- t(null(Wpca_CS))
    WpcaFull_CS <- rbind(Wpca_CS, We_CS)

    # Build the initial guess from the artificial parameters
    if (i == 1 ){
      P_artificial <- P_CS
      Y_artificial <- Y_CS
      Wpca_artificial <- Wpca_CS
      We_artificial <- We_CS
      WpcaFull_artificial <- WpcaFull_CS
      K1XQ_artificial <- K1XQ__CS
    }else{
      P_artificial <-rbind(P_artificial, P_CS)
      Y_artificial <- rbind(Y_artificial, Y_CS)
      Wpca_artificial <- magic::adiag(Wpca_artificial, Wpca_CS)
      We_artificial <- magic::adiag(We_artificial, We_CS)
      WpcaFull_artificial <- magic::adiag(WpcaFull_artificial, WpcaFull_CS)
      K1XQ_artificial <- magic::adiag(K1XQ_artificial, K1XQ__CS)
    }
    idxJ0<- idxJ1
    idxN0 <- idxN1
  }

  K <- nrow(ZZ_artificial)
  Gy.0 <- diag(K)

  # SSZ, K0Z and K1Z
  if (ModelType == 'JPS' || ModelType == 'JPS jointP' || ModelType == 'VAR jointQ'){
    VARpara <- VAR(ZZ_artificial, VARtype= 'unconstrained', Bcon = NULL)
    K0Z_artificial <- VARpara$K0Z
    K1Z_artificial <- VARpara$K1Z
    SSZ_artificial <- VARpara$SSZ
  }

  if (ModelType == 'GVAR sepQ'|| ModelType == 'GVAR jointQ'){
    GVARpara <- GVAR(GVARinputs, N)
    K0Z_artificial <- GVARpara$F0
    K1Z_artificial <- GVARpara$F1
    SSZ_artificial <- GVARpara$Sigma_y
  }

  if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma"){
    JLLinputs$WishSigmas <- 1
    JLLPara <- JLL(ZZ_artificial, N, JLLinputs)
    K0Z_artificial <-JLLPara$k0
    K1Z_artificial <-JLLPara$k1
    JLLinputs$WishSigmas <- 0 # Ensures that the variance-covariance matrix will no longer be estimated within the JLL function

    if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit"){  SSZ_artificial <- JLLPara$Sigmas$VarCov_NonOrtho }
    if (ModelType == "JLL jointSigma"){ SSZ_artificial <- JLLPara$Sigmas$VarCov_Ortho
    # NOTE: the required vectorization that preceeds the numerical optimization of the variance-covariance matrix
    # is done on the orthogonalized variance-covariance matrix because we want to preserve the restrictions imposed in this matrix.
    # However, to compute the loadings A and B from Y= A + B*P we have to use the non-orthogonalized variance-covariance
    # matrix. We make this adjustment in the function 'A0N_MLEdensity_WOE...' by redefining SSZ before computing A and B.
    # (this avoids overcomplicating the overall structure of the code")
    }

  }
  # Prepare the list of outputs:
  ListOutputs <- list(Wpca_artificial, We_artificial, WpcaFull_artificial, Y_artificial, P_artificial, ZZ_artificial,
                      Gy.0, SSZ_artificial, K1XQ_artificial, K0Z_artificial, K1Z_artificial, JLLinputs, GVARinputs)

  names(ListOutputs) <- c("Wpca","We", "WpcaFull", "Y", "PP", "ZZ", "Gy.0", "SSZ", "K1XQ", "K0Z", "K1Z",
                          "JLLinputs", "GVARinputs")


  return(ListOutputs)
}

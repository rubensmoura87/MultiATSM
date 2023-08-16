#' Concatenate the model-specific inputs in a list
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param Data dataset generated from the "DataForEstimation" function
#'@param Economies string-vector containing the names of the economies of the system
#'@param VARXtype  string-vector containing the VARX feature (see "GVAR" function) (GVAR-based models)
#'@param t_First_Wgvar Sample starting date (year) (GVAR-based models)
#'@param t_Last_Wgvar  Sample last date (year) (GVAR-based models)
#'@param W_type  Criterion used in the computation of the star variables (see "Transition_Matrix" function)
#'               (GVAR-based models)
#'@param DomUnit name of the economy which is assigned as the dominant unit (JLL-based models)
#'@param WishSigmas equal to "1" if one wishes the variance-covariance matrices and the Cholesky factorizations (JLL-based models)
#'@param SigmaNonOrtho NULL or some F x F matrix from the non-orthogonalized dynamics (JLL-based models)
#'@param BiasCorrection binary variable. it takes value equal to 1 if the user whishes the estimates to be bias-corrected
#'                      and 0, otherwise. (BRW model)
#'@param flag_mean flag whether mean- (TRUE) or median- (FALSE) unbiased estimation is desired
#'@param gamma adjustment parameter (BRW model)
#'@param N_iter number of iterations (BRW model)
#'@param N_burn number of burn-in iterations (BRW model)
#'@param B  number of bootstrap samples (BRW model)
#'@param checkBRW flag whether the user wishes to perform the closeness check (BRW model)
#'@param B_check number of bootstrap samples for closeness check
#'@param DataPathTrade path of the Excel file containing the data (if any)
#'
#'
#'@examples
#'
#'ModelType <- "JLL original"
#'Eco <- c("China","Brazil","Mexico", "Uruguay")
#'DU <- "China"
#'Sig <- 1
#'NonOrtho <- 0
#'
#'
#'ListModelInputs(ModelType, Economies= Eco, DomUnit = DU, WishSigmas = Sig, SigmaNonOrtho= NonOrtho)
#'
#'
#'@export



ListModelInputs <-function(ModelType, Data= NULL, Economies, VARXtype= NULL, t_First_Wgvar= NULL, t_Last_Wgvar= NULL,
                           W_type = NULL, DomUnit= NULL, WishSigmas = NULL, SigmaNonOrtho= NULL, BiasCorrection = 0,
                           flag_mean = NULL, gamma= NULL, N_iter = NULL, N_burn= NULL, B= NULL, checkBRW = NULL,
                           B_check = NULL, DataPathTrade = NULL){

  if (any(ModelType == c("GVAR sepQ", "GVAR jointQ"))){
    GVARinputs <- list()
    GVARinputs$Economies <- Economies
    GVARinputs$GVARFactors <- Data$GVARFactors
    GVARinputs$VARXtype <- VARXtype

    GVARinputs$Wgvar <- Transition_Matrix(t_First_Wgvar, t_Last_Wgvar, Economies, W_type, DataPathTrade, Data$Wgvar)
    if (is.list(GVARinputs$Wgvar)){
      if (t_First_Wgvar !=t_Last_Wgvar){stop("For estimating GVAR models with time-varying interdependence,
                                             the start and ending dates of the transition matrix must coincide!")}

      GVARinputs$Wgvar <- GVARinputs$Wgvar[[t_First_Wgvar]]
    }
  } else { GVARinputs <- NULL }



  if (any(ModelType == c("JLL original", "JLL NoDomUnit","JLL jointSigma"))){
    JLLinputs <- list()
    JLLinputs$Economies <- Economies
    JLLinputs$DomUnit <- DomUnit
    JLLinputs$WishSigmas <- WishSigmas # Sigma matrix is estimated within the "InputsForMLEdensity" function
    JLLinputs$SigmaNonOrtho <- SigmaNonOrtho
    JLLinputs$JLLModelType <- ModelType
  }else{  JLLinputs <- NULL }


  if(BiasCorrection==1){
    BRWinputs <- list()
    BRWinputs$flag_mean <- flag_mean # 1 = compute the mean; 0 = compute the median
    BRWinputs$gamma <- gamma # Adjustment parameter
    BRWinputs$N_iter <- N_iter # Number of iteration to be conserved
    BRWinputs$N_burn <- N_burn  # Number of iteration to be discarded
    BRWinputs$B <- B # Number of bootstrap samples
    BRWinputs$checkBRW <- checkBRW
    BRWinputs$B_check <- B_check #
  } else{ BRWinputs <- NULL}


  outputs <- list(GVARinputs, JLLinputs, BRWinputs)
  names(outputs) <- c("GVARinputs","JLLinputs","BRWinputs")
  return(outputs)
}

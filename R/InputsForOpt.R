#' Generates inputs necessary to build the likelihood function for the ATSM model
#'
#' @param InitialSampleDate Start date of the sample period in the format "dd-mm-yyyy"
#' @param FinalSampleDate End date of the sample period in the format "dd-mm-yyyy"
#' @param ModelType A character vector indicating the model type to be estimated. Available options: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param Yields A numerical matrix with time series of yields (JxT or CJ x T)
#' @param GlobalMacro A numerical matrix with time series of the global risk factors (G x T)
#' @param DomMacro A numerical matrix with time series of the country-specific risk factors for all C countries (CM x T)
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param DataFrequency A character vector specifying the frequency of the data. Available options are: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", or "Annually".
#' @param GVARlist A list containing the necessary inputs for the estimation of GVAR-based models
#' @param JLLlist A list of necessary inputs for the estimation of JLL-based models. If the chosen model is "JLL original" or "JLL joint Sigma", then a dominant unit economy must be chosen. Otherwise, this list must be set as 'None'.
#' @param WishBRW Logical. Whether to estimate the physical parameter model with bias correction, based on the method by Bauer, Rudebusch and Wu (2012). Default is FALSE.
#' @param BRWlist List of necessary inputs for performing the bias-corrected estimation.
#' @param UnitYields A character string indicating the maturity unit of yields. Options are: "Month" for yields expressed in months, or "Year" for yields expressed in years. Default is "Month".
#' @param CheckInputs Logical. Whether to perform a prior check on the consistency of the provided input list. Default is TRUE.
#' @param BS_Adj Logical. Whether to adjust the global series for the sepQ models in the Bootstrap setting. Default is FALSE.
#'
#' @importFrom pracma null
#'
#' @return An object of class 'ATSMModelInputs' containing the necessary inputs for performing the model optimization.
#' @examples
#' \donttest{
#' # Example 1:
#' data(CM_GlobalMacroFactors)
#' data(CM_DomMacroFactors)
#' data(CM_Yields)
#'
#' ModelType <- "JPS original"
#' Economies <- "Mexico"
#' t0 <- "01-05-2007" # Initial Sample Date (Format: "dd-mm-yyyy")
#' tF <- "01-12-2018" # Final Sample Date (Format: "dd-mm-yyyy")
#' N <- 3
#' GlobalVar <- c("Gl_Eco_Act") # Global Variables
#' DomVar <- c("Eco_Act") # Domestic Variables
#' FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#'
#' DataFreq <- "Monthly"
#'
#' ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacroVar, DomesticMacroVar,
#'                              FactorLabels, Economies, DataFreq, CheckInputs = FALSE)
#'
#' # Example 2:
#' LoadData("CM_2024")
#'
#' ModelType <- "GVAR multi"
#'
#' Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#' t0 <- "01-05-2007" # InitialSampleDate (Format: "dd-mm-yyyy")
#' tF <- "01-12-2019" # FinalSampleDate (Format: "dd-mm-yyyy")
#' N <- 2
#' GlobalVar <- c("Gl_Eco_Act", "Gl_Inflation") # Global Variables
#' DomVar <- c("Inflation") # Domestic Variables
#' FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#'
#' DataFreq <- "Monthly"
#' GVARlist <- list(VARXtype = "unconstrained", W_type = "Sample Mean",
#'                  t_First_Wgvar = "2007", t_Last_Wgvar = "2019")
#'
#' ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacroVar, DomesticMacroVar,
#'                            FactorLabels, Economies, DataFreq, GVARlist, CheckInputs = FALSE)
#'
#' # Example 3:
#' if (requireNamespace('neldermead', quietly = TRUE)) {
#' LoadData("CM_2024")
#'
#' ModelType <- "JLL original"
#'
#' Economies <- c("China", "Brazil", "Uruguay")
#' t0 <- "01-05-2007" # InitialSampleDate (Format: "dd-mm-yyyy")
#' tF <- "01-12-2019" # FinalSampleDate (Format: "dd-mm-yyyy")
#' N <- 2
#' GlobalVar <- c("Gl_Eco_Act", "Gl_Inflation") # Global Variables
#' DomVar <- c("Eco_Act", "Inflation") # Domestic Variables
#' FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#'
#' JLLinputs <- list(DomUnit = "China")
#'
#' DataFrequency <- "Monthly"
#'
#' ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacroVar, DomesticMacroVar,
#'                            FactorLabels, Economies, DataFreq, JLLlist = JLLinputs,
#'                            CheckInputs = FALSE)
#' } else {
#'  message("skipping functionality due to missing Suggested dependency")
#' }
#' }
#' @export

InputsForOpt <- function(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
                         FactorLabels, Economies, DataFrequency, GVARlist = NULL, JLLlist = NULL,
                         WishBRW = FALSE, BRWlist = NULL, UnitYields = "Month", CheckInputs = TRUE, BS_Adj = FALSE) {

  # Print initial message with model type and sample period
  cat(paste("1) PREPARING INPUTS FOR THE ESTIMATION OF THE MODEL:", ModelType, ". SAMPLE PERIOD:", InitialSampleDate,
            "-", FinalSampleDate, "\n"))

  # Define labels for single and multi-country models
  Label_Single_Models <- c("JPS original", "JPS global", "GVAR single")
  Label_Multi_Models <- c("GVAR multi", "JPS multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma")

  # Check consistency of inputs across models if CheckInputs is TRUE
  if (isTRUE(CheckInputs)) {
  CheckInputsForMLE(InitialSampleDate, FinalSampleDate, Economies, DomMacro, GlobalMacro, UnitYields,
                    DataFrequency, Label_Single_Models, Label_Multi_Models, FactorLabels, GVARlist, ModelType) }


  # Construct the time-series of the risk factors
  cat("1.1) Constructing the time-series of the risk factors \n")
  RiskFactors <-  BuildATSM_RiskFactors(InitialSampleDate, FinalSampleDate, Yields, GlobalMacro, DomMacro,
                                       Economies, FactorLabels, ModelType, BS_Adj)

  Yields <- AdjustYieldsDates(Yields, RiskFactors, Economies)

  # Build the vector of common maturities across countries
  mat <- Maturities(Yields, Economies, UnitYields)

  # Collect general input lists
  ModelInputsGen <- GeneralMLEInputs(Yields, RiskFactors, FactorLabels, mat, DataFrequency, Label_Multi_Models,
                                     Economies, ModelType)

  # Collect model-specific input lists (GVARinputs, JLLinputs, and BRWlist)
  ModelInputsSpe <- SpecificMLEInputs(ModelType, Economies, RiskFactors, FactorLabels, GVARlist, JLLlist,
                                      WishBRW, BRWlist, DataPathTrade = NULL)

  # Estimate P-dynamic parameters
  cat("1.2) Estimating model P-dynamics parameters: \n")
  PdynPara <- GetPdynPara(RiskFactors, FactorLabels, Economies, ModelType, ModelInputsSpe$BRWinputs,
                          ModelInputsSpe$GVARinputs, ModelInputsSpe$JLLinputs, CheckInputs)

  # Gather outputs to export
  Outputs <- Outputs2exportMLE(Label_Multi_Models, Economies, RiskFactors, Yields, mat, ModelInputsGen,
                              ModelInputsSpe, PdynPara, ModelType)

  # Store metadata inside the class without explicitly exporting it
  attr(Outputs, "ModelInfo") <- list(
    ModelType = ModelType,
    Economies = Economies,
    InitialSampleDate = InitialSampleDate,
    FinalSampleDate = FinalSampleDate,
    RiskFactors = RiskFactors,
    Yields = Yields,
    DataFrequency = DataFrequency,
    Maturities = mat,
    G = length(FactorLabels$Global),
    NC = length(FactorLabels$Spanned)*length(Economies),
    MC = (length(FactorLabels$Domestic) - length(FactorLabels$Spanned))*length(Economies)
    )

  # Return the structured Outputs object
  return(structure(Outputs, class = "ATSMModelInputs"))
}
##########################################################################################################################
#' Print method for ATSMModelInputs objects
#' @param Object An object of class 'ATSMModelInputs'
#' @param ... Additional arguments (not used)
#'
#' @export

print.ATSMModelInputs <- function(Object, ...) {

  info <- attr(Object, "ModelInfo")

  cat("ATSM Model Inputs Object\n")
  cat("------------------------------------------\n")
  cat("Model Type:", info$ModelType, "\n")
  cat("Economic System:", paste(info$Economies, collapse = ", "), "\n")
  cat("Sample Period:", info$InitialSampleDate, "-", info$FinalSampleDate, "\n")
  cat("Data Frequency:", info$DataFrequency, "\n")
  cat("Common Maturities across Countries:", info$Maturities, "years \n")
  cat("------------------------------------------\n\n")

  cat("Key Structure of the Economic System: \n")
  cat("------------------------------------------\n")
  cat("Total amount of spanned factors in the system:", info$NC, "\n")
  cat("Total amount of global unspanned factors in the system:", info$G, "\n")
  cat("Total amount of country-specific unspanned factors in the system:", info$MC, "\n")
  cat("Total amount of risk factors in the system:", info$NC + info$G + info$MC, "\n")
  cat("------------------------------------------\n")

  invisible(Object)
}
##########################################################################################################################
#' Summary method for ATSMModelInputs objects
#' @param Object An object of class 'ATSMModelInputs'
#' @param ... Additional arguments (not used)
#'
#' @export

summary.ATSMModelInputs <- function(Object, ...) {


  info <- attr(Object, "ModelInfo")

  # Function to print the summary statistics
  summary_TS <- function(mat) {
    data.frame(
      Mean = round(rowMeans(mat, na.rm = TRUE), 3),
      Std_Dev = round(apply(mat, 1, sd, na.rm = TRUE), 3),
      Min = round(apply(mat, 1, min, na.rm = TRUE), 3),
      Max = round(apply(mat, 1, max, na.rm = TRUE), 3),
      row.names = rownames(mat)
    )
    }

  cat("Summary Statistics From the Time Series Components:\n")
  cat("------------------------------------------------------\n")

  cat("1) Risk Factors:\n")
  # Single-country estimation
  if(info$ModelType %in% c("JPS original", "JPS global", "GVAR single")){
    for(i in 1:length(info$Economies)){
    cat("\n Model", info$Economies[i] , "\n")
    summary_result <- summary_TS(info$RiskFactors[[i]])
    print(summary_result)
    }

    # Multicountry estimation
}else{
  summary_result <- summary_TS(info$RiskFactors)
  print(summary_result)
}


  cat("\n 2) Bond Yields:\n")
  # Single-country estimation
  if(info$ModelType %in% c("JPS original", "JPS global", "GVAR single")){

    J <- length(info$Maturities)

    for(i in 1:length(info$Economies)){
      cat("\n Model", info$Economies[i] , "\n")

      # Extract rows for the current country
      country_rows <- ((i - 1) * J + 1):(i * J)
      country_yields <- info$Yields[country_rows, ] * 100

      summary_result <- summary_TS(country_yields)
      print(summary_result)
    }

    # Multicountry estimation
}else{
  summary_result <- summary_TS((info$Yields)*100)
  print(summary_result)
}
  cat("------------------------------------------------------\n")

  invisible(Object)
}

##########################################################################################################################
#' Create a vector of numerical maturities in years
#'
#'@param DataYields matrix containing all yields of the system (JxT,if the model is single-country-based
#'                  or CJxT if the model is multy-country-based)
#'@param Economies vector containing names of all the economies of the system
#'@param UnitYields (i) "Month": if maturity of yields are expressed in months or
#'                  (ii) "Year": if maturity of yields are expressed in years
#'
#'@return Vector containing all observed maturities expressed in years
#'
#'@keywords internal


Maturities <- function(DataYields, Economies, UnitYields){


  if (UnitYields == "Month"){ fac <- 12}
  else if (UnitYields == "Year"){ fac <- 1}

  C <- length(Economies)
  s <- rownames(DataYields)
  AllMat <- readr::parse_number(s)
  J <- length(AllMat)/C
  mat <- AllMat[1:J]/fac

  return(mat)
}

#########################################################################################################################################
#' Concatenate the model-specific inputs in a list
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param Economies string-vector containing the names of the economies of the system
#'@param RiskFactors time series of risk factors (F x T)
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param GVARlist A list of required inputs to estimate the GVAR-based setups:
#'\enumerate{
#'      \item VARXtype  string-vector containing the VARX feature (see "GVAR" function) (GVAR-based models)
#'      \item t_First_Wgvar Sample starting date (year) (GVAR-based models)
#'      \item t_Last_Wgvar  Sample last date (year) (GVAR-based models)
#'      \item W_type  Criterion used in the computation of the star variables (see "Transition_Matrix" function)
#'               (GVAR-based models)
#'  }
#'@param JLLlist A list of required inputs to estimate the JLL-based setups:
#'\enumerate{
#'      \item DomUnit name of the economy which is assigned as the dominant unit (JLL-based models)
#'      \item WishSigmas equal to "1" if one wishes the variance-covariance matrices and the Cholesky factorizations (JLL-based models)
#'      \item SigmaNonOrtho NULL or some F x F matrix from the non-orthogonalized dynamics (JLL-based models)
#'}
#'@param WishBRW Whether the user wishes to estimate the physical parameter model with the Bias correction model from BRW (2012) (see "Bias_Correc_VAR" function).\cr
#'              Default is set to 0.
#'@param BRWlist A list of required inputs to estimate the bias corrected setups of the type of BRW:
#'\enumerate{
#'      \item BiasCorrection binary variable. it takes value equal to 1 if the user whishes the estimates to be bias-corrected
#'                      and 0, otherwise. (BRW model)
#'      \item flag_mean flag whether mean- (TRUE) or median- (FALSE) unbiased estimation is desired
#'      \item gamma adjustment parameter (BRW model)
#'      \item N_iter number of iterations (BRW model)
#'      \item N_burn number of burn-in iterations (BRW model)
#'      \item B  number of bootstrap samples (BRW model)
#'      \item checkBRW flag whether the user wishes to perform the closeness check (BRW model)
#'      \item B_check number of bootstrap samples for closeness check
#'}
#'@param DataPathTrade path of the Excel file containing the data (if any)
#'
#'@keywords internal


SpecificMLEInputs <-function(ModelType, Economies, RiskFactors, FactorLabels, GVARlist = NULL, JLLlist = NULL,
                             WishBRW=0, BRWlist= NULL, DataPathTrade = NULL){

  if (any(ModelType == c("GVAR single", "GVAR multi"))){
    # Compute the transition matrix and the
    # PRELIMINARY CHECK:  Makes sure that link matrices are correctly imputed in a model with time-varying interdependence
    if (GVARlist$W_type == "Time-varying" & GVARlist$t_First_Wgvar != GVARlist$t_Last_Wgvar){
      stop("For estimating GVAR models with time-varying interdependence, the start and ending dates
            of the transition matrix must coincide!")
    }

    t0 <- GVARlist$t_First_Wgvar
    tF <- GVARlist$t_Last_Wgvar
    LinkageMeasure <- GVARlist$DataConnectedness
    W_type <- GVARlist$W_type
    Wgvar <- Transition_Matrix(t0, tF, Economies, W_type, LinkageMeasure, DataPathTrade)
    if(ModelType == "GVAR single"){RiskFactors <- RiskFactors[[1]]}

    GVARFactors <- DataSet_BS(ModelType, RiskFactors, Wgvar, Economies, FactorLabels)

    # Build the list of the necessary inputs to estimate a GVAR model:
    GVARinputs <- list( Economies = Economies, GVARFactors = GVARFactors, VARXtype = GVARlist$VARXtype)

    if (W_type == "Time-varying"){
      GVARinputs$Wgvar <- Wgvar[[t0]]
    }else{  GVARinputs$Wgvar <- Wgvar} # constant interdependence

    } else { GVARinputs <- NULL }



  if (any(ModelType == c("JLL original", "JLL No DomUnit","JLL joint Sigma"))){
    JLLinputs <- list(
      Economies = Economies,
      DomUnit = JLLlist$DomUnit,
      WishSigmas = 1, # Sigma matrix is estimated within the "InputsForMLEdensity" function
      SigmaNonOrtho = NULL,
      JLLModelType = ModelType
    )
  }else{  JLLinputs <- NULL }


  if(WishBRW==1){BRWinputs <- BRWlist} else{ BRWinputs <- NULL}

  outputs <- list(GVARinputs = GVARinputs, JLLinputs = JLLinputs, BRWinputs = BRWinputs)
  return(outputs)
}

###########################################################################################################
#' Makes sure that the time series of yields and risk factors have coincident sample spans
#'
#'@param Yields time series of bond yields (CJ x T1 )
#'@param PdynamicsFactors time series of risk factors (F x T2)
#'@param Economies string-vector containing the names of the economies of the system.
#'
#'@keywords internal


AdjustYieldsDates <- function(Yields, PdynamicsFactors, Economies){

  # Restrict the sample span
  if(is.list(PdynamicsFactors)){
    t0 <- colnames(PdynamicsFactors[[1]])[1]
    tF <- colnames(PdynamicsFactors[[1]])[ncol(PdynamicsFactors[[1]])]
  }else{
    t0 <- colnames(PdynamicsFactors)[1]
    tF <- colnames(PdynamicsFactors)[ncol(PdynamicsFactors)]
  }

  Idx0 <- which(colnames(Yields)== t0)
  IdxF <- which(colnames(Yields)== tF)

  YieldsAdj <- Yields[, Idx0:IdxF]

  # Restrict the sample of countries
  IdxCountries <- sapply(Economies, function(economy) grep(economy, rownames(YieldsAdj)))
  YieldsAdj <- YieldsAdj[IdxCountries, ]

  return(YieldsAdj)
}



###############################################################################################################
#################################################################################################################
#' Loads data sets from several papers

#'@param DataPaper  Available options are \code{BR_2017} (Bauer and Rudebusch, 2017) , \code{CM_2023} (Candelon and Moura, 2023), \code{CM_2024} (Candelon and Moura, 2024)
#'
#'
#'@examples
#' #Example 1:
#' LoadData("BR_2017")
#'
#' #Example 2:
#' LoadData("CM_2023")
#'
#' #Example 3:
#' LoadData("CM_2024")
#'
#'
#'@returns
#' Complete set of data from several papers.
#'
#'@references
#'\enumerate{
#'\item Bauer and Rudebusch (2017). "Resolving the Spanning Puzzle in Macro-Finance Term Structure Models" (Review of Finance)
#'\item Candelon and Moura (2023). "Sovereign yield curves and the COVID-19 in emerging markets" (Economic Modelling)
#'\item Candelon and Moura (2024). "A Multicountry Model of the Term Structures of Interest Rates with a GVAR" (Journal of Financial Econometrics)
#'}
#'@export



LoadData <- function(DataPaper) {
  switch(DataPaper,
         "BR_2017" = utils::data("BR_jps_gro_R3"),
         "CM_2024" = {
           utils::data("CM_GlobalMacroFactors")
           utils::data('CM_DomMacroFactors')
           utils::data('CM_Trade')
           utils::data('CM_Yields')
         },
         "CM_2023" = {
           utils::data("CM_GlobalMacro_2023")
           utils::data('CM_DomMacro_2023')
           utils::data('CM_Trade_2023')
           utils::data('CM_Yields_2023')
         },
         stop("Database unavailable.")
  )
}


#################################################################################################################
#' Builds the time series of the risk factors that are used in the estimation of the ATSM
#'
#'@param InitialSampleDate Sample starting date
#'@param FinalSampleDate Sample last date
#'@param Yields matrix  (J x T), where  J - the number of maturities and  T - time series length
#'@param GlobalMacroFactors time series of the global macroeconomic risk factors (G x T)
#'@param DomMacroFactors time series of the country-specific macroeconomic risk factors for all C countries (CM x T)
#'@param Economies      string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels string-list based which contains the labels of all variables present in the model
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param BS_Adj adjustment of global series for sepQ model in the Bootstrap setting. Default is set to FALSE.
#'
#'@return
#'Generates the complete set of risk factors that are used in the estimation of the ATSM
#'
#'@keywords internal


BuildATSM_RiskFactors <- function(InitialSampleDate, FinalSampleDate, Yields, GlobalMacroFactors,
                                  DomMacroFactors, Economies, FactorLabels, ModelType, BS_Adj = FALSE){

  # Preliminary work
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)

  AllFactorLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)

  F <- C*(M+N) + G
  T <- ncol(DomMacroFactors)

  RiskFactors <- matrix(NA, nrow =  F, ncol = T)
  colnames(RiskFactors) <- colnames(DomMacroFactors)
  rownames(RiskFactors) <- AllFactorLabels

  # Input the global factors
  if (BS_Adj == FALSE || (any(ModelType == c("JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit",
                                             "JLL joint Sigma")))){
    IdxGlobal <- match(FactorLabels$Global, rownames(GlobalMacroFactors))
    IdxGlobalRF <- IdxGlobal[!is.na(IdxGlobal)]

    CommonLabelsGlob <- intersect(rownames(GlobalMacroFactors), AllFactorLabels)
    GlobalMacroFactorsClean <- GlobalMacroFactors[rownames(GlobalMacroFactors) %in% CommonLabelsGlob, ]

    RiskFactors[seq_len(G), ] <-  GlobalMacroFactorsClean
  }
  # Input the spanned factors
  PP <- Spanned_Factors(Yields, Economies, N)
  IdxSpa <- match(rownames(PP), AllFactorLabels)
  RiskFactors[IdxSpa, ] <- PP

  # Input the unspanned factors
  IdxUnspa <- match(rownames(DomMacroFactors), AllFactorLabels)
  IdxUnspaRF <- IdxUnspa[!is.na(IdxUnspa)]

  CommonLabelsUnsp <- intersect(rownames(DomMacroFactors), AllFactorLabels)
  DomMacroFactorsClean <- DomMacroFactors[rownames(DomMacroFactors) %in% CommonLabelsUnsp, ]

  RiskFactors[IdxUnspaRF, ] <- DomMacroFactorsClean

  # Adjust dates to the desired sample sample
  Idx0 <- which(colnames(RiskFactors)== InitialSampleDate)
  IdxF <- which(colnames(RiskFactors)== FinalSampleDate)

  RiskFactors <- RiskFactors[, Idx0:IdxF]

  # Adapt to the specifies of the single-country-based models
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    RiskFactorsList <- list()

    for (i in 1:C){
      if (BS_Adj){RiskFactors[seq_len(G), ] <- GlobalMacroFactors[((1+G*(i-1)):(G + G*(i-1))), ]   }

      if (ModelType == 'JPS original'){
        idxCountryFactors <- which(grepl(Economies[i], rownames(RiskFactors))) # index of the rows of the country-specific factors
        idxFactors <- c(seq_len(G), idxCountryFactors) # index global + country-specific factors
        RiskFactorsList[[Economies[i]]]  <-  RiskFactors[idxFactors, ]
      }  else {
        RiskFactorsList[[Economies[i]]] <- RiskFactors
      }
    }
    RiskFactors <- RiskFactorsList
  }

   return(RiskFactors)
}

###################################################################################################################
#' Get delta t
#'
#' @param DataFrequency single element character-based vector. Available options are: "Daily All Days", \cr
#'                      "Daily Business Days", "Weekly", "Monthly",  "Quarterly", "Annually"
#'
#'@keywords internal

Getdt <- function(DataFrequency){

  if (DataFrequency == "Daily All Days"){ dt <- 1/365}
  else if (DataFrequency == "Daily Business Days"){ dt <- 1/252}
  else if (DataFrequency == "Weekly"){ dt <- 1/52}
  else if (DataFrequency == "Monthly"){ dt <- 1/12}
  else if (DataFrequency == "Quarterly"){ dt <- 1/4}
  else if (DataFrequency == "Annually"){ dt <- 1}

  return(dt)
}

#############################################################################################################
#' Check consistence of inputs
#'
#'
#'@param t0 Sample starting date
#'@param tF Sample last date
#'@param Economies string-vector containing the names of the economies of the system.
#'@param DomesticMacroFac time series of the country-specific macroeconomic risk factors for all C countries (CM x T)
#'@param GlobalMacroFac time series of the global macroeconomic risk factors (G x T)
#'@param UnitYields (i) "Month": if maturity of yields are expressed in months or
#'                  (ii) "Year": if maturity of yields are expressed in years
#'@param DataFreq single element character-based vector. Available options are: "Daily All Days", \cr
#'                      "Daily Business Days", "Weekly", "Monthly",  "Quarterly", "Annually"
#'@param Label_Single_Models string-vector containing the names of the single country setups
#'@param Label_Multi_Models string-vector containing the names of the multicountry setups
#'@param FactorLabels string-list based which contains the labels of all variables present in the model
#'@param GVARlist list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'
#'@keywords internal

CheckInputsForMLE <- function(t0, tF, Economies, DomesticMacroFac, GlobalMacroFac, UnitYields, DataFreq,
                              Label_Single_Models, Label_Multi_Models, FactorLabels, GVARlist, ModelType){

  N <- length(FactorLabels$Spanned)
  DomVarLab <- utils::head(FactorLabels$Domestic, - N)
  GloVarLab <- FactorLabels$Global
  TimeSeriesLabels <- colnames(DomesticMacroFac)

  # CHECK 1: Check whether the model choice is compatible with the number of countries selected
  if (length(Economies) == 1 &  (any(ModelType == Label_Multi_Models))){
    stop("The models 'GVAR multi', 'JPS multi', 'JLL original', 'JLL No DomUnit', and 'JLL joint Sigma' are multicountry setups.
       Therefore, they require the estimation of several countries.")
  }

  # CHECK 2: Check whether the model type selected if available
  if (!(ModelType %in% c(Label_Single_Models, Label_Multi_Models))){stop("Model type is invalid.")}

  # CHECK 3: Label consistency in terms of economies, domestic and global sets
  Check_label_consistency(Economies, DomesticMacroFac, GlobalMacroFac, DomVarLab, GloVarLab)

  # CHECK 4: consistency of the chosen data sample
  if (!(all(c(t0, tF) %in% TimeSeriesLabels))){ stop("Initial and/or final sample dates are not part of the data sample span.") }

  # CHECK 5: data frequency chosen
  if(!(DataFreq %in% c("Daily All Days", "Daily Business Days", "Weekly", "Monthly",  "Quarterly", "Annually"))){
    stop(paste("Data frequency option", DataFreq, "is unavailable. Available options are: 'Daily All Days', 'Daily Business Days', 'Weekly', 'Monthly',  'Quarterly', 'Annually'."))}
  # CHECK 6: unit of bond yields chosen
  if(!(UnitYields %in% c("Month", "Year"))){
    stop(paste("Bond yield time-unit option", UnitYields,"is unavailable. Available options are 'Month' or 'Year'."))}

}
############################################################################################################
#'Gathers the general inputs for model estimation
#'
#'@param Yields matrix (CJ x T) or a list containing those matrices, where C is the number of countries,
#'              J - the number of maturities and T - time series length \cr
#'@param RiskFactors time series of risk factors (F x T). Could be stored in a list depending on the model
#'@param FactorLabels string-list based which contains the labels of all variables present in the model
#'@param mat vector of maturities (in years) used in the estimation
#'@param DataFrequency single element character-based vector. Available options are: "Daily All Days", \cr
#'                      "Daily Business Days", "Weekly", "Monthly",  "Quarterly", "Annually"
#'@param Label_Multi_Models string-vector containing the names of the multicountry setups
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'
#'@keywords internal


GeneralMLEInputs <- function(Yields, RiskFactors, FactorLabels, mat, DataFrequency, Label_Multi_Models,
                             Economies, ModelType){

  N <- length(FactorLabels$Spanned)
  dt <- Getdt(DataFrequency)
  J <- length(mat)

  # Compute the quantities of interest
  idxJ0 <- 0
  idxN0 <- 0

  ListInputs <- list()

  for (i in 1:length(Economies)){ # Country-specific inputs
    idxJ1 <- idxJ0 + J
    idxN1 <- idxN0 + N

    YCS <- Yields[(idxJ0+1):idxJ1,] # Yields (JxT)
    WpcaCS <- pca_weights_one_country(YCS, Economy= Economies[i])[1:N,] # matrix of weights for the portfolio without errors (N x J)
    WpcaCS <- 100*matrix(WpcaCS, nrow=N) # useful command for the case N = 1
    WeCS <- t(null(matrix(WpcaCS, nrow=N))) # matrix of weights for the yield portfolios priced with errors
    WpcaFullCS <- rbind(WpcaCS, WeCS)
    PPCS <- Spanned_Factors(YCS, Economies = Economies[i], N) # Spanned Factors (N x T)
    K1XQCS <-  Reg_K1Q(YCS, mat, PPCS, dt, type="Jordan")

    if ( any(ModelType == Label_Multi_Models)){
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
    } else{
      # Matrix of contemporaneous terms
      Gy.0 <- diag(nrow(RiskFactors[[1]]))
      # Store the outputs of interest for the single country model types
      ListInputs[[Economies[i]]] <- list( Wpca = WpcaCS, We = WeCS, WpcaFull = WpcaFullCS, Yields = YCS,
                                          SpaFact = PPCS, K1XQ = K1XQCS, Gy.0 = Gy.0)
    }

    idxJ0<- idxJ1
    idxN0 <- idxN1
  }

  # Matrix of contemporaneous terms
  if ( any(ModelType == Label_Multi_Models)){
    colnames(Y) <- colnames(RiskFactors)
    Gy.0 <- diag(nrow(RiskFactors))
  }

  # Store the outputs of interest for the multicountry model types
  if ( any(ModelType == Label_Multi_Models)){
    ListInputs <- list(Wpca = Wpca, We = We, WpcaFull = WpcaFull, Yields = Y, SpaFact = PP, K1XQ = K1XQ, Gy.0 = Gy.0)
  }

  return(ListInputs)
}


##############################################################################################################
#'Compute the parameters used in the P-dynamics of the model
#'
#'@param RiskFactors time series of risk factors (F x T). Could be stored in a list depending on the model
#'@param FactorLabels string-list based which contains the labels of all variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param BRWinputs list of necessary inputs for performing the bias-corrected estimation (see "Bias_Correc_VAR" function)
#'@param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#'@param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'
#'
#'@keywords internal


GetPdynPara <- function(RiskFactors, FactorLabels, Economies,  ModelType, BRWinputs, GVARinputs,
                        JLLinputs, CheckInputs = F){

  N <- length(FactorLabels$Spanned)

  # Parameters of the of the P-dynamics
  # Get bias-corrected parameters, if selected
  if (!is.null(BRWinputs)){
    cat("- With the bias-correction procedure from BRW (2012) \n")

    if(CheckInputs){
      if (any(ModelType == c("GVAR single", "GVAR multi"))){
        CheckInputsGVAR(GVARinputs, N)
      }else if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
        CheckJLLinputs(RiskFactors, JLLinputs)
      }
    }

    tryCatch({
    PdynPara <- GetPdynPara_BC(ModelType, BRWinputs, RiskFactors, Economies, FactorLabels, GVARinputs, JLLinputs)
    }, error = function(err) {
      stop("BRW procedure leads to a highly collinear system. Please choose another combination of BRW parameters.")})


  } else{

    cat("- Without the bias-correction procedure \n\n")
    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
      cat("JLL-based setup in progress: Estimating the variance-covariance matrix numerically.
          This may take some time. \n\n")
    }
    # Otherwise compute the non-biased correction model parameters
    PdynPara <- GetPdynPara_NoBC(ModelType, RiskFactors, Economies, N, GVARinputs, JLLinputs, CheckInputs)
  }


  return(PdynPara)

}

#####################################################################################################
#'Prepares inputs to export
#'
#'@param Label_Multi_Models string-vector containing the names of the multicountry setups
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param RiskFactors time series of risk factors (F x T). Could be stored in a list depending on the model
#'@param Yields matrix (CJ x T) or a list containing those matrices, where C is the number of countries,
#'              J - the number of maturities and T - time series length \cr
#'@param mat vector of maturities (in years) used in the estimation
#'@param ModelInputsGen List of generic inputs
#'@param ModelInputsSpec List of specific inputs
#'@param PdynPara Model parameters estimated in the P-dynamics the
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'@keywords internal


Outputs2exportMLE <- function(Label_Multi_Models, Economies, RiskFactors, Yields, mat, ModelInputsGen,
                              ModelInputsSpec, PdynPara, ModelType){



  if ( any(ModelType == Label_Multi_Models)){
    Output <- list(mat = mat,
                   Wpca = ModelInputsGen$Wpca,
                   We = ModelInputsGen$We,
                   WpcaFull = ModelInputsGen$WpcaFull,
                   Yields = ModelInputsGen$Y,
                   SpaFact = ModelInputsGen$SpaFact,
                   RiskFactors= RiskFactors,
                   Gy.0 = ModelInputsGen$Gy.0,
                   K1XQ = ModelInputsGen$K1XQ,
                   SSZ = PdynPara$SSZ,
                   K0Z = PdynPara$K0Z,
                   K1Z = PdynPara$K1Z,
                   JLLinputs = ModelInputsSpec$JLLinputs,
                   GVARinputs = ModelInputsSpec$GVARinputs)


  } else{

    C <- length(Economies)
    Output <- list()

    for (i in 1:C) {
      Output[[Economies[i]]] <- list(
        mat = mat,
        Wpca = ModelInputsGen[[Economies[i]]]$Wpca,
        We = ModelInputsGen[[Economies[i]]]$We,
        WpcaFull = ModelInputsGen[[Economies[i]]]$WpcaFull,
        Yields = ModelInputsGen[[Economies[i]]]$Y,
        SpaFact = ModelInputsGen[[Economies[i]]]$SpaFact,
        RiskFactors = RiskFactors[[Economies[i]]],
        Gy.0 = ModelInputsGen[[Economies[i]]]$Gy.0,
        K1XQ = ModelInputsGen[[Economies[i]]]$K1XQ,
        SSZ = PdynPara[[Economies[i]]]$SSZ,
        K0Z = PdynPara[[Economies[i]]]$K0Z,
        K1Z = PdynPara[[Economies[i]]]$K1Z,
        JLLinputs = ModelInputsSpec$JLLinputs,
        GVARinputs = ModelInputsSpec$GVARinputs
      )
    }

  }
  return(Output)
}
##############################################################################################################
#'Compute P-dynamics parameters using the bias correction method from BRW (2012)
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param BRWinputs list of necessary inputs for performing the bias-corrected estimation (see "Bias_Correc_VAR" function)
#'@param RiskFactors time series of risk factors (F x T). Could be stored in a list depending on the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels string-list based which contains the labels of all variables present in the model
#'@param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#'@param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'
#'
#'@keywords internal


GetPdynPara_BC <- function(ModelType, BRWinputs, RiskFactors, Economies, FactorLabels, GVARinputs, JLLinputs){

  N <- length(FactorLabels$Spanned)
  C <- length(Economies)

  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
    PdynPara <- list()
  for (i in 1:C){
    RiskFactors_CS <- RiskFactors[[Economies[i]]]

    cat(paste("Estimation for country:", Economies[i], "\n"))
    BiasCorrec <- Bias_Correc_VAR(ModelType, BRWinputs, t(RiskFactors_CS), N, Economies, FactorLabels,
                                  GVARinputs, JLLinputs)
    PdynPara[[Economies[i]]] <- list (K0Z = BiasCorrec$mu_tilde, K1Z = BiasCorrec$Phi_tilde,
                                      SSZ = BiasCorrec$V_tilde)
  }
} else{
  BiasCorrec <- Bias_Correc_VAR(ModelType, BRWinputs, t(RiskFactors), N, Economies, FactorLabels,
                                GVARinputs, JLLinputs)
  K0Z <- BiasCorrec$mu_tilde
  K1Z <- BiasCorrec$Phi_tilde
  SSZ <- BiasCorrec$V_tilde
}

# Export outputs
  if (any(ModelType == c('JPS original', 'JPS global', 'GVAR single'))){ Out <- PdynPara
}else{ Out <- list(K0Z = K0Z, K1Z= K1Z, SSZ=SSZ)}

  return(Out)
}


####################################################################################################################
#'Compute P-dynamics parameters without using the bias correction method from BRW (2012)
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param RiskFactors time series of risk factors (F x T). Could be stored in a list depending on the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param N number of country-specific spanned factors
#'@param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#'@param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'
#'@keywords internal

GetPdynPara_NoBC <- function(ModelType, RiskFactors, Economies, N, GVARinputs, JLLinputs, CheckInpts = F){

  # JPS-related models
if (any(ModelType == c('JPS original', 'JPS global'))){
  PdynPara <- list()
  for (i in 1:length(Economies)){
    RiskFactorsMat <- RiskFactors[[Economies[i]]]
    VARpara <- VAR(RiskFactorsMat, VARtype= 'unconstrained', Bcon = NULL)
    PdynPara[[Economies[i]]] <- VARpara
  }

} else if (ModelType== 'JPS multi'){
  VARpara <- VAR(RiskFactors, VARtype= 'unconstrained', Bcon = NULL)
  K0Z <- VARpara$K0Z
  K1Z <- VARpara$K1Z
  SSZ <- VARpara$SSZ
}
# GVAR-related models
else if (any(ModelType == c('GVAR single', 'GVAR multi'))){
  GVARpara <- GVAR(GVARinputs, N, CheckInpts)
  K0Z <- GVARpara$F0
  K1Z <- GVARpara$F1
  SSZ <- GVARpara$Sigma_y

  if (ModelType =='GVAR single'){
    PdynPara <- list()
    for (i in 1:length(Economies)){  PdynPara[[Economies[i]]] <- list(K0Z = K0Z, K1Z = K1Z, SSZ = SSZ)}

  }
}
# JLL-related models
else if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
  JLLinputs$WishSigmas <- 1
  JLLPara <- JLL(RiskFactors, N, JLLinputs, CheckInpts)
  K0Z <-JLLPara$k0
  K1Z <-JLLPara$k1
  JLLinputs$WishSigmas <- 0 # Ensures that the variance-covariance matrix will no longer be estimated within the JLL function

  if (any(ModelType == c("JLL original", "JLL No DomUnit"))){  SSZ <- JLLPara$Sigmas$VarCov_NonOrtho }
  else{ SSZ <- JLLPara$Sigmas$VarCov_Ortho
  # NOTE: the required vectorization that preceeds the numerical optimization of the variance-covariance matrix
  # is done on the orthogonalized variance-covariance matrix because we want to preserve the restrictions imposed in this matrix.
  # However, to compute the loadings A and B from Y= A + B*P we have to use the non-orthogonalized variance-covariance
  # matrix. We make this adjustment in the function 'A0N_MLEdensity_WOE...' by redefining SSZ before computing A and B.
  # (this avoids overcomplicating the overall structure of the code")
  }

}

  # Export outputs
  if (any(ModelType == c('JPS original', 'JPS global', 'GVAR single'))){ Out <- PdynPara
  }else{ Out <- list(K0Z = K0Z, K1Z= K1Z, SSZ=SSZ)}

  return(Out)
}
###################################################################################################################
#' Check consistency of labels (economies, domestic and global variables)
#'
#'
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param DomesticMacroFac time series of the country-specific domestic risk factors (C(M+N) x T)
#'@param GlobalMacroFac time series of the global risk factors (G X T)
#'@param DomVarLab string-vector containing the names of the desired domestic risk factors
#'@param GloVarLab string-vector containing the names of the desired global risk factors
#'
#'@keywords internal


Check_label_consistency <- function(Economies, DomesticMacroFac, GlobalMacroFac, DomVarLab, GloVarLab){

  DomVarLabels_ALL <- row.names(DomesticMacroFac)
  GlobalVarLabels_ALL <- row.names(GlobalMacroFac)

  # a) Check consistence in the economy set
  missing_economies <- Economies[!sapply(Economies, function(economy) any(grepl(economy, DomVarLabels_ALL)))]
  if (length(missing_economies) > 0) {
    print(paste("The following economy names are misspecified:", paste(missing_economies, collapse = ", "),
                ". Check whether (i) the label of the domestic variables is followed by the economy's name or (ii) the economy name is correctly spelled."))
  }

  # b) Check consistence in the domestic variable set
  missing_DomVariables <- DomVarLab[!sapply(DomVarLab, function(DomVarLab) any(grepl(DomVarLab, DomVarLabels_ALL)))]
  if (length(missing_DomVariables) > 0) {
    print(paste("Inconsitent domestic variable name(s):", paste(missing_DomVariables, collapse = ", ")))
  }

  # c) Check consistence in the global variable set
  missing_GlobalVariables <- GloVarLab[!sapply(GloVarLab, function(GloVarLab) any(grepl(GloVarLab, GlobalVarLabels_ALL)))]
  if (length(missing_GlobalVariables) > 0) {
    print(paste("Inconsitent global variable name(s):", paste(missing_GlobalVariables, collapse = ", ")))
  }

  if (any(lengths(list(missing_economies, missing_DomVariables, missing_GlobalVariables)) > 0)){
    stop("Inconsistent model input specifications. Check messages above.")
  }
}

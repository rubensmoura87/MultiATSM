#' Retrieve data from Excel and build the database used in the model estimation
#'
#' @param t0 Sample starting date (yyyy-mm-dd)
#' @param tF Sample last date (yyyy-mm-dd)
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param N   Number of country-specific spanned factor (scalar)
#' @param FactorLabels String-list based which contains the labels of all the variables present in the model
#' @param ModelType   String-vector containing the label of the model to be estimated
#' @param DataFrequency  Character-based-vector. Avaialable options are: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"
#' @param W_type Three possibilities:
#'\itemize{
#'      \item "Full Sample": if one wishes ALL weight matrices of each year from which data is available (it may extrapolate the sample period);
#'      \item "Sample Mean": if one wishes a SINGLE weight matrix containing the average of weights over of the entire sample period;
#'      \item Some year in particular (e.g. "1998", "2005" ...).
#'}
#' @param t_First_Wgvar Sample starting date (year)
#' @param t_Last_Wgvar Sample last date (year)
#' @param DataPathMacro  Path of the Excel file conating the macroeconomic data (if any). The default is linked to the excel file present in the package.
#' @param DataPathYields  Path of the Excel file conating the yields data (if any). The default is linked to the excel file present in the package.
#' @param DataPathTrade Path of the Excel file conating the trade data (if any). The default is linked to the excel file present in the package.
#'
#'
#'@return   A list containing the
#'\enumerate{
#'\item time series of the complete set of bond yields (matrix, JxT or CJxT);
#'\item time series of the complete set risk factors (matrix, KxT);
#'\item 'GVARFactors': list of all variables that are used in the estimation of the VARX \cr
#'                (see e.g. 'CM_Factors_GVAR' file). If the estimated model type is not GVAR-based, then returns NULL.
#'}
#'
#'@examples
#'DomVar <- c("Eco_Act", "Inflation")
#'GlobalVar <- c("GBC", "CPI_OECD")
#' t0 <- "2006-09-01"
#' tF <-  "2019-01-01"
#' Economies <- c("China", "Brazil", "Mexico", "Uruguay", "Russia")
#' N <- 2
#' ModelType <- "JPS"
#' FactorLabels <-  LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#' DataFrequency <- "Monthly"
#'
#'
#' DataModel <- DataForEstimation(t0, tF, Economies, N, FactorLabels, ModelType, DataFrequency)
#'
#' @export

DataForEstimation <- function(t0, tF, Economies, N, FactorLabels, ModelType, DataFrequency, W_type= NULL,
                              t_First_Wgvar= NULL, t_Last_Wgvar = NULL, DataPathMacro = NULL, DataPathYields = NULL,
                              DataPathTrade = NULL){


C <- length(Economies)

if (ModelType == 'GVAR sepQ' || ModelType == 'GVAR jointQ'){
  Wgvar <- Transition_Matrix(t_First_Wgvar, t_Last_Wgvar, Economies, W_type, DataPathTrade)
}
FactorSet <- DatabasePrep(t0, tF, Economies, N, FactorLabels, ModelType, Wgvar, DataPathMacro, DataPathYields)

# Gather all bond yields of all countries
Yields <- c()
for (i in 1:C){
  if (i== 1){
    Yields <- FactorSet[[Economies[i]]]$Yields
  }else{
    Yields <- rbind(Yields, FactorSet[[Economies[i]]]$Yields)
  }
}


# Gather the spanned and unspanned factors

RiskFactors <- RiskFactorsPrep(FactorSet, Economies, FactorLabels, t0, tF, DataFrequency)

# Gather the factors used in the GVAR estimation
if (ModelType == 'GVAR sepQ' || ModelType == 'GVAR jointQ'){
GVARFactors <-  DataSet_BS(ModelType, RiskFactors, Wgvar, Economies, FactorLabels)
} else { GVARFactors <- NULL}


Outputs <-list(Yields, RiskFactors, GVARFactors)
names(Outputs) <- c("Yields", "RiskFactors", "GVARFactors")


return(Outputs)
}


#' Collect the inputs that are used to construct the numerical and the graphical outputs
#'
#'@param ModelType String-vector containing the label of the model to be estimated
#'@param Horiz Single numerical vector conataining the desired horizon of analysis for the outputs
#'@param ListOutputWished List of desired graphical outputs. Available options are: "Fit","IRF", "FEVD", "GIRF", "GFEVD".
#'@param OutputLabel Name of the output label to be stored
#'@param WishStationarityQ User must set 1 is she whises to impose the largest eigenvalue under the Q to be strictly
#'                       smaller than 1, otherwise set 0.
#'@param WishGraphYields Binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise. Default is set as "0".
#'@param WishGraphRiskFactors Binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise. Default is set as "0".
#'@param WishOrthoJLLgraphs Binary variable: set 1, if the user wishes orthogonalized JLL-based graphs to be generated; or set 0, otherwise.
#'                         Default is set as "0"
#'@param WishBootstrap  Binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise. Default is set as "0".
#'@param ListBoot       List containing the four following elements:
#'\enumerate{
#'  \item "methodBS": Desired bootstrap method among (a) 'bs' (standard residual bootstrap), (b) 'wild' (wild bootstrap),
#'                   (c) 'block' (block  bootstrap);
#'  \item "BlockLength": if block bootstrap is chosen, then the user has to specify the lenght of the block (single numerical vector);
#'  \item "ndraws": number of draws;
#'  \item "pctg": level of confidence (single numerical vector expressed in basis points)
#'}
#'@param WishForecast Binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise. Default is set as "0".
#'@param ListForecast list containing the three following elements:
#'\enumerate{
#'  \item "ForHoriz": forecast horizon;
#'  \item "t0Sample": index of the first variable of the information set;
#'  \item "t0Forecast": index of the first forecast cut-off date.
#'}
#'
#'
#'
#'@examples
#'
#' ModelType <- "JPS"
#' Horiz <- 100
#' DesiredOutputGraphs <- c("Fit", "GIRF", "GFEVD")
#' OutputLabel <- "Test"
#' WishStationarityQ <- 1
#' WishGraphRiskFac <- 0
#' WishGraphYields <- 1
#'
#'
#'InputsList <- InputsForOutputs(ModelType, Horiz, DesiredOutputGraphs, OutputLabel,
#'                               WishStationarityQ, WishGraphYields, WishGraphRiskFac)
#'@returns
#' List of necessary inputs to generate the graphs of the outputs of the desired model
#'@export


InputsForOutputs <- function(ModelType, Horiz, ListOutputWished, OutputLabel, WishStationarityQ,
                             WishGraphYields = 0, WishGraphRiskFactors=0, WishOrthoJLLgraphs=0,
                             WishBootstrap =0, ListBoot = NULL, WishForecast = 0, ListForecast= NULL){


OutputTypeSet <- c("Fit","IRF", "FEVD", "GIRF", "GFEVD")
IdxWishOut <- match(ListOutputWished, OutputTypeSet)

if ('Fit' %in% ListOutputWished) {WishFitGraph <- 1}

# 1) Point estimate list
InputsForOutputs <- list()

InputsForOutputs[[1]] <- OutputLabel
InputsForOutputs[[2]] <- WishStationarityQ
names(InputsForOutputs) <- c("Label Outputs", "StationaryQ")


# Initialization
for (h in 1:length(OutputTypeSet)){

  if (h==1){
    InputsForOutputs[[ModelType]]$Fit$WishGraphs <- 0
  } else {

  InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$horiz <- Horiz

  InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactors <- 0
  InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$Yields <- 0
  InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactorsBootstrap <- 0
  InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$YieldsBootstrap <- 0

  if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma" ){
    InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactors <- 0
    InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$Yields <- 0
    InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactorsBootstrap <- 0
    InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$YieldsBootstrap <- 0
    }
}
}


# Filling up the desired outputs
for (h in IdxWishOut){

  if (h==1){
    InputsForOutputs[[ModelType]]$Fit$WishGraphs <- WishFitGraph
  } else {

InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$horiz <- Horiz

InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactors <- WishGraphRiskFactors
InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$Yields <- WishGraphYields
if (WishBootstrap ==1){
InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactorsBootstrap <- WishGraphRiskFactors
InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$YieldsBootstrap <- WishGraphYields
}

if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma" ){
if (WishOrthoJLLgraphs==1){
InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactors <- WishGraphRiskFactors
InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$Yields <- WishGraphYields
if (WishBootstrap ==1){
InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactorsBootstrap <- WishGraphRiskFactors
InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$YieldsBootstrap <- WishGraphYields
}
}
}
}
}


# 2) Bootstrap list
if (WishBootstrap==0){
  InputsForOutputs[[ModelType]]$Bootstrap$WishBoot <- 0
}else{
  InputsForOutputs[[ModelType]]$Bootstrap <- WishBootstrap
  names(InputsForOutputs[[ModelType]]$Bootstrap) <- "WishBoot"
  InputsForOutputs[[ModelType]]$Bootstrap <- append(InputsForOutputs[[ModelType]]$Bootstrap, ListBoot)
}


#3) Forecasting list
if (WishForecast==0){
  InputsForOutputs[[ModelType]]$Forecasting$WishForecast <- 0
}else{
  InputsForOutputs[[ModelType]]$Forecasting <- WishForecast
  names(InputsForOutputs[[ModelType]]$Forecasting) <- "WishForecast"
  InputsForOutputs[[ModelType]]$Forecasting <- append(InputsForOutputs[[ModelType]]$Forecasting, ListForecast)
}


return(InputsForOutputs)
}


###############################################################################################################



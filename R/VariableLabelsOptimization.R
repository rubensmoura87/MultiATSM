#' Create the variable labels used in the estimation
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param WishStationarityQ User must set "1" is she whises to impose the largest eigenvalue under the Q to be strictly
#'                       smaller than 1. Otherwise set "0"
#'
#'@export



ParaLabels <- function(ModelType, WishStationarityQ){

ParaLabelsList <- list()



ParaLabelsList[[ModelType]][["r0"]] <- "@r0: bounded"
ParaLabelsList[[ModelType]][["se"]] <- "@se: bounded"
ParaLabelsList[[ModelType]][["K0Z"]] <- "@K0Z: bounded"
ParaLabelsList[[ModelType]][["K1Z"]] <- "@K1Z: bounded"


# K1XQ
K1XQType <- K1XQStationary(WishStationarityQ)$SepQ

if (ModelType == "JPS" || ModelType == "JPS jointP" || ModelType == 'GVAR sepQ' ){
  ParaLabelsList[[ModelType]][["K1XQ"]] <- K1XQStationary(WishStationarityQ)$SepQ
} else{
  ParaLabelsList[[ModelType]][["K1XQ"]] <- K1XQStationary(WishStationarityQ)$JointQ
}


# SSZ
if (ModelType == "JPS" || ModelType == "JPS jointP" || ModelType == 'VAR jointQ' ){ ParaLabelsList[[ModelType]][["SSZ"]] <- "SSZ: psd" }
if (ModelType == "GVAR sepQ" || ModelType == 'GVAR jointQ'){  ParaLabelsList[[ModelType]][["SSZ"]] <- "SSZ: BlockDiag" }
if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit"){ ParaLabelsList[[ModelType]][["SSZ"]] <- "@SSZ: bounded" } # Variance-covariance matrix is not estimated under Q
if (ModelType == "JLL jointSigma"){ ParaLabelsList[[ModelType]][["SSZ"]] <- "SSZ: JLLstructure" }
# Ensures that the structure of the Variance-covariance matrix of the JLL is preserved

return(ParaLabelsList)
}

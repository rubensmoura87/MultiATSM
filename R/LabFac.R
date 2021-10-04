#' Generates the labels factors
#'
#'@param N number of spanned factors per country (scalar)
#'@param DomVar character-vector containing the names of the domestic variables
#'@param GlobalVar character-vector containing the names of the global variables
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'
#'@examples
#'N <- 2
#'DomVar <- c("inflation", "Economic growth")
#'GlobalVar <- "Commodity Prices"
#'Economies <- c("U.S.", "Canada", "Germany", "Japan")
#'ModelType <- "JPS"
#'
#'VarLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#'
#'@export

LabFac <- function(N, DomVar, GlobalVar, Economies, ModelType){

  M <- length(DomVar) # Number of domestic unspanned factors
  C <- length(Economies) # Number of countries of the economic system

  FactorLabels <- list()
  FactorLabels$Spanned <- LabelsSpanned(N)
  FactorLabels$Domestic <- c(DomVar, FactorLabels$Spanned)
  FactorLabels$Star <-  LabelsStar(FactorLabels$Domestic)
  FactorLabels$Global <- GlobalVar


  labelsDomVar <- c()
  labelsDomVarJLL <- c()
  j <-0
  for (i in 1:C){
    count0 <- (i-1)*j
    for (j in 1:(N+M)){
      count1 <- count0 + j
      labelsDomVar[count1] <- paste(FactorLabels$Domestic[j], Economies[i])
      if ("JLL original" %in% ModelType || "JLL NoDomUnit" %in% ModelType || "JLL jointSigma" %in% ModelType){
        labelsDomVarJLL[count1] <- paste(FactorLabels$Domestic[j], Economies[i], "JLL")
      }
    }
    FactorLabels$Tables[[Economies[i]]] <- labelsDomVar[(count0+1):count1]
    FactorLabels$TablesJLL[[Economies[i]]] <- labelsDomVarJLL[(count0+1):count1]
  }


  idx0 <- 0
  for (a in 1:C){
    idx1 <- idx0 + N + M
    FactorLabels$Tables$AllCountries[(idx0+1):idx1] <- do.call(cbind,lapply(FactorLabels$Tables[[Economies[a]]],matrix,nrow=1))
    if ("JLL original" %in% ModelType || "JLL NoDomUnit" %in% ModelType || "JLL jointSigma" %in% ModelType){
      FactorLabels$Tables$AllCountriesJLL[(idx0+1):idx1] <- do.call(cbind,lapply(FactorLabels$TablesJLL[[Economies[a]]],matrix,nrow=1))
    }
    idx0 <- idx1
  }



  FactorLabels$TablesJLL <- NULL # Remove the country-specific JLL factors list


  return(FactorLabels)
}


########################################################################################################
#' Generate the labels of the spanned factors
#'
#' @param N number of spanned factors
#'

LabelsSpanned <-function(N){

  # Goal: assign the labels of the spanned factors depending on the number of spanned factors

  LabelsSpannedALL <- c("Level", "Slope", "Curvature", "Fourth PC", "Fifth PC", "Sixth PC", "Seventh PC", "Eighth PC")

  LabelsSpan <- LabelsSpannedALL[1:N]


  return(LabelsSpan)
}

###################################################################################################
#' Generate the labels of the star variables
#'
#'@param FactorLabels Factor labels


LabelsStar <-function(FactorLabels){

  L <- length(FactorLabels)

  Label <- c()
  for (j in 1:L){
    Label[j] <- paste(FactorLabels[j], ".Star", sep="")
  }

  return(Label)
}



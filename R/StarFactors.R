#' Generates the star variables necessary for the GVAR estimation
#'
#'@param RiskFactors time series of the risk factors (F x T)
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'@param W GVAR transition matrix (C x C)

#'
#'@return List of star factors
#'@examples
#'data(CM_Factors)

#' Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#' Wgvar <- matrix( c(0, 0.83, 0.86, 0.38, 0.65, 0, 0.13, 0.55,
#'          0.32, 0.12, 0, 0.07, 0.03, 0.05, 0.01, 0), nrow = 4, ncol = 4)
#' rownames(Wgvar) <- Economies
#' colnames(Wgvar) <- Economies

#' StarFactors(RiskFactors, Economies, Wgvar)
#'
#'@export



StarFactors <- function(RiskFactors, Economies, W){

C <- length(Economies)
T <- ncol(RiskFactors)

StarLabel <- c()
ListFactors <- list()


# Re-arrange country-specific factors per country
for (i in 1:C){
IdxRF <- grepl(Economies[i], rownames(RiskFactors))
ListFactors[[Economies[i]]]$Factors <- RiskFactors[IdxRF,] # Country-specific factors

VarLabel <- rownames(RiskFactors)[IdxRF] # Label of the country-specific factors
if (i==1){
  StarLabels <- paste(VarLabel, ".Star", sep="")
}else{
  StarPrep<- paste(VarLabel, ".Star", sep="")
  StarLabels <- append(StarLabels, StarPrep) # Label of the star-variables
}
}


NumDomFac <- length(IdxRF[IdxRF== TRUE])
Z <- list()

for (j in 1:NumDomFac){
  X <- matrix(NA, nrow= C, ncol=T)
  for (i in 1:C){
    X[i,] <- ListFactors[[Economies[i]]]$Factors[j,]
    Z[[j]] <- X # Each element of the list contains the same country-specific variable of all countries
  }
}


# Compute the star variables
StarVariables <- matrix(NA, nrow = C*NumDomFac, ncol = T)
  for (i in 1:C){
    idxx <- (i-1)*NumDomFac
    for (j in 1:NumDomFac){
      StarVariables[idxx + j,]   <- t(W[i,]%*%Z[[j]])
    }
  }

rownames(StarVariables) <- StarLabels
colnames(StarVariables) <- colnames(RiskFactors)


return(StarVariables)
}




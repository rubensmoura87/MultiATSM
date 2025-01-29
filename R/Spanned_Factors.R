#' Computes the country-specific spanned factors
#
#'@param Yields     matrix  (J x T), where J is the number of maturities and T is the length of the time series.
#'@param Economies  A character vector containing the names of the economies included in the system.
#'@param N     Scalar representing the desired number of country-specific spanned factors (maximum allowed is N = J).
#'
#'@return Matrix containing the N spanned for all the countries of the system  (CJ xT)
#'@examples
#' data(CM_Yields)
#'Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#'N <- 3
#'SpaFact_TS <- Spanned_Factors(Yields, Economies, N)
#'
#'@export


Spanned_Factors <- function(Yields, Economies, N){

  C <- length(Economies)   # number of the economies of the system

  #Label of the spanned factors
  LabelsSpannedALL <- c("Level", "Slope", "Curvature", "Fourth PC", "Fifth PC", "Sixth PC", "Seventh PC", "Eighth PC")
  LabelsSpanned <- LabelsSpannedALL[1:N]

  # Spanned factors:
  P <- c()

  for (i in 1:C){
    Idx <- grepl(Economies[i], rownames(Yields)) # Extract the yields of a single country
    Y_CS <- Yields[Idx,] # Country-specific yields
    W <- pca_weights_one_country(Y_CS, Economies[i]) # Weight matrix
    if(N > nrow(Y_CS)){stop("The Number of country-specific spanned factors cannot exceed the number of country-specific bond yields.")}
    W <- W[1:N,]*100 # Select only the set of weights that will be used to compute the first N PCs
    if (i == 1){
      P <- W%*%Y_CS # Country-specific spanned factor
      row.names(P) <- paste(LabelsSpanned, Economies[i])
      }else{
    Ptemp <- W%*%Y_CS # Country-specific spanned factor
    row.names(Ptemp) <- paste(LabelsSpanned, Economies[i])
    P <- rbind(P, Ptemp)
  }
  }


  return(P)
}

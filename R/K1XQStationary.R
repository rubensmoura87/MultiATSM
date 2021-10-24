#' Impose stationarity under the Q-measure
#'
#'@param StationaryEigenvalues Binary variable: set "1" if the user whises the largest eigenvalue
#'                            to be strictly smaller than 1. Set "0", otherwise
#'
#'@examples
#' stat <- 1 # Takes values 1 and 0
#'K1XQStationary(stat)
#'
#'
#'@returns
#'list
#'
#'@export

K1XQStationary<- function(StationaryEigenvalues){

  K1Type <- list()

  if (StationaryEigenvalues == 1){
    K1Type$SepQ <- paste("K1XQ: ", "Jordan", "; stationary", sep="")
    K1Type$JointQ <- paste("K1XQ: ", "Jordan", " MultiCountry", "; stationary", sep="")
  }else{
    K1Type$SepQ <- paste("K1XQ: ", "Jordan", sep="")
    K1Type$JointQ <- paste("K1XQ: ", "Jordan", " MultiCountry", sep="")
  }

  return(K1Type)
}

#' Create a vector of numerical maturities in years
#'
#'@param DataYields matrix containing all yields of the system (JxT,if the model is single-country-based
#'                  or CJxT if the model is multy-country-based)
#'@param Economies vector containing names of all the economies of the system
#'@param UnitYields (i) "Month": if maturity of yields are expressed in months or
#'                  (ii) "Year": if maturity of yields are expressed in years
#'@importFrom readr parse_number
#'
#'@return Vector containing all observed maturities expressed in years
#'
#'@examples
#'data('CM_Yields')
#'Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#'Maturities(Yields, Economies, "Month")
#'@export



Maturities <- function(DataYields, Economies, UnitYields){


  if (UnitYields == "Month"){ fac <- 12}
  if (UnitYields == "Year"){ fac <- 1}

  C <- length(Economies)
  s <- rownames(DataYields)
  AllMat <- parse_number(s)
  J <- length(AllMat)/C
  mat <- AllMat[1:J]/fac


  return(mat)
}

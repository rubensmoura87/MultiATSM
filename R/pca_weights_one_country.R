#' Weight matrix from principal components
#
#'@param Y        A numeric matrix (T x J) representing the time series of bond yields of a specific country
#'@param Economy  String-vector containing the name of a single economy.
#'
#'
#'@return A matrix (J x J) that corresponds to the eigenvectors of the variance-covariance matrix of yields
#'@examples
#'data("CM_Yields")
#'pca_weights_one_country(Yields, Economy= "Brazil")
#'@export


pca_weights_one_country <- function(Y, Economy) {


  # Extract the yields of a single country
  Idx <- grepl(Economy, rownames(Y))
  Y_CS <- Y[Idx,] # Country-specific yields

    # Store the eigenvectors of the matrix of yields
  H <- eigen(stats::cov(t(Y_CS)))
  W <- t(H$vectors)

  J <- nrow(Y_CS) # Total number of maturities
  M <- round(stats::median(1:J)) # Median of the maturity range


  # Adjust the weights to ease the interpretability of the spanned factors
  # a) LEVEL: if the weights of the level are negative, then invert the sign
  # (ideally, high values of level should be interpreted as a high values for all yields)
  if (all(W[1,] < 0)){   W[1,] <- W[1,]*(-1)   }
  # b) SLOPE: if the slope is downward sloping, then invert the sign
  # (ideally, high values of slope should be interpreted as a steep high curve)
  if (W[2,1] > W[2,J]){   W[2,] <- W[2,]*(-1)   }
  # c) Curvature: if the curvature is U-shapped, then invert the sign
  # (ideally, the curvature factor should capture an effect that is linked to mid-term maturities)
  if (W[3,1] > W[3,M] & W[3,J] > W[3,M]){ W[3,] <- W[3,]*(-1) }


  return(W)
}

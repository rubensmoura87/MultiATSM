#' Computes the transition matrix required in the estimation of the GVAR model
#'
#' @param t_First Sample starting date (in the format: yyyy).
#' @param t_Last  Sample ending date (in the format: yyyy).
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param type A character string indicating the method for computing interdependence. Possible options include:
#' \itemize{
#'      \item \code{Time-varying}: Computes time-varying interdependence and returns the weight matrices for each year based on available data (may extrapolate the sample period).
#'      \item \code{Sample Mean}: Returns a single weight matrix containing the average weights over the entire sample period, suitable for time-invariant interdependence.
#'      \item A specific year (e.g., "1998", "2005"): Used to compute time-invariant interdependence for the specified year.
#' }
#' @param DataConnectedness Data used to compute the transition matrix. Default is set to NULL.
#' @param DataPath Path to the Excel file containing the data (if applicable). The default is linked to the Excel file available in the package.
#'
#' @return matrix or list of matrices
#'
#' @details
#' If there is missing data for any country of the system for that particularly year,
#' then the transition matrix will include only NAs.
#'
#' @examples
#' data(CM_Trade)
#'
#' t_First <- "2006"
#' t_Last <-  "2019"
#' Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#' type <- "Sample Mean"
#' W_mat <- Transition_Matrix(t_First, t_Last, Economies, type, DataConnectedness = TradeFlows)
#'
#' @export

Transition_Matrix <- function(t_First, t_Last, Economies, type, DataConnectedness = NULL, DataPath = NULL) {

  # 1) Load data if Data Connectedness is not provided
  if (is.null(DataConnectedness) || length(DataConnectedness) == 0) {
    if (is.null(DataPath)) {
      DataPath <- system.file("extdata", "TradeData.xlsx", package = "MultiATSM")
    }

    tab_names_Trade <- readxl::excel_sheets(DataPath)
    list_all_Trade <- suppressMessages(lapply(tab_names_Trade, function(x) readxl::read_excel(path = DataPath, sheet = x)))
    names(list_all_Trade) <- tab_names_Trade

    L <- length(list_all_Trade)

    for (i in 1:L) {
      Countries <- list_all_Trade[[i]][[1]]
      list_all_Trade[[i]] <- as.data.frame(list_all_Trade[[i]][,-1])
      rownames(list_all_Trade[[i]]) <- Countries
    }

    DataConnectedness <- list_all_Trade
  }

  # 2) Pre-allocation of variables
  DataAdj <- lapply(DataConnectedness, function(x) as.data.frame(t(x)))
  C <- length(Economies)
  T <- ncol(DataConnectedness[[1]])
  WgvarAllYears <- vector(mode='list', length = T)
  names(WgvarAllYears) <- colnames(DataConnectedness[[1]])
  num <- matrix(NA, nrow = C, ncol= C)
  dem <- c()
  Wyear <- matrix(NA, nrow = C, ncol= C)

  # 3) Generate the matrix of weights year-by-year
  for (k in 1:T) {
    for (h in 1:C) {
      for (j in 1:C) {
        num[h,j] <- DataAdj[[Economies[h]]][[Economies[j]]][k]
      }
      dem[h] <- sum(num[h,])
      Wyear[h,] <- num[h,] / dem[h]
    }
    if (any(is.na(Wyear))) { # If there is missing data for any country of the system for that particularly year, then return a matrix of NAs.
      Wyear <- matrix(NA, nrow = C, ncol= C)
    }
    WgvarAllYears[[k]] <- Wyear
  }

  # Desired output
  # a) Time-varying: Gathers the interdependence matrices from the full Sample
  if (type == "Time-varying") {
    Wgvar <- Filter(function(x) all(stats::complete.cases(x)), WgvarAllYears)
  }
  # b) Sample Mean
  else if (type == "Sample Mean") {
    Y_First <- substring(t_First, 1,4)
    Y_Last <- substring(t_Last, 1,4)
    TimeLable <- colnames(DataConnectedness[[1]])
    idx0 <- which(TimeLable == Y_First)
    idx1 <- which(TimeLable == Y_Last)

    WgvarSubSample <- WgvarAllYears[idx0:idx1]
    Wgvar <- apply(simplify2array(WgvarSubSample), 1:2, mean)

    colnames(Wgvar) <- Economies
    rownames(Wgvar) <- Economies
  }
  # c) Specific Year
  else {
    Wgvar <- WgvarAllYears[[type]]
    colnames(Wgvar) <- Economies
    rownames(Wgvar) <- Economies
  }

  return(Wgvar)
}

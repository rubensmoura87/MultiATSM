#' Compute the transition matrix required in the estimation of the GVAR model
#'
#'@param t_First      Sample starting date (year)
#'@param t_Last       Sample last date (year)
#'@param Economies    Vector containing the names of all the economies of the system.
#'@param type         Three possibilities:
#'\itemize{
#'      \item "Full Sample": if one wishes ALL weight matrices of each year from which data is available (it may extrapolate the sample period);
#'      \item "Sample Mean": if one wishes a SINGLE weight matrix containing the average of weights over of the entire sample period;
#'      \item Some year in particular (e.g. "1998", "2005" ...).
#'}
#'
#'@param DataPath     path of the Excel file containing the data (if any). The default is linked to the Excel file available in the package.
#'@param Data         Data for computing the transition matrix. Default is set to NULL.
#'
#
#'@return matrix or list of matrices
#'
#'@details
#' NOTE: if there is missing data for any country of the system for that particularly year,
#' then the transition matrix will include only NAs.
#'
#'@examples
#'data(CM_Trade)
#'
#' t_First <- "2006"
#' t_Last <-  "2019"
#' Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#' type <- "Sample Mean"
#' Transition_Matrix(t_First, t_Last, Economies, type, DataPath = NULL, Data = TradeFlows)
#'
#'
#'@export


Transition_Matrix <- function(t_First, t_Last, Economies, type, DataPath = NULL, Data = NULL){


  if (sjmisc::is_empty(Data)){

  if (sjmisc::is_empty(DataPath)){ DataPath <- system.file("extdata", "TradeData.xlsx", package = "MultiATSM") }

  tab_names_Trade <- readxl::excel_sheets(DataPath)
  list_all_Trade <- suppressMessages(lapply(tab_names_Trade, function(x) readxl::read_excel(path = DataPath, sheet = x)))
  names(list_all_Trade) <- tab_names_Trade

  L <- length(list_all_Trade)


  for (i in 1:L){
    Countries <- list_all_Trade[[i]][[1]]
    list_all_Trade[[i]] <- list_all_Trade[[i]][,-1]
    row.names(list_all_Trade[[i]]) <- suppressMessages(Countries)
  }


  Data <- list_all_Trade
  }

# a) Pre-allocation of variables
DataAdj <- lapply(Data, function(x) as.data.frame(t(x)))
C <- length(Economies)
T <- ncol(Data[[1]])
WgvarAllYears <- vector(mode='list', length = T)
names(WgvarAllYears) <- colnames(Data[[1]])
num <- matrix(NA, nrow = C, ncol= C)
dem <- c()
Wyear <- matrix(NA, nrow = C, ncol= C)

# b) Generate the matrix of weights year-by-year
for(k in 1:T){
  for (h in 1:C){
    for (j in 1:C){
      num[h,j] <- DataAdj[[Economies[h]]][[Economies[j]]][k]
    }
    dem[h] <- sum(num[h,])
    Wyear[h,] <- num[h,]/dem[h]
      }
  if(any(is.na(Wyear))){ # If there is missing data for any country of the system for that particularly year, then return a matrix of NAs.
    Wyear <- matrix(NA, nrow = C, ncol= C)
  }
  WgvarAllYears[[k]] <- Wyear
}


# c) desired output
# c.1) "Full Sample"
if (type == "Full Sample"){
  Wgvar <- WgvarAllYears
}

# c.2) "Sample Mean"
if (type == "Sample Mean"){

  Y_First <- substring(t_First, 1,4)
  Y_Last <- substring(t_Last, 1,4)
  TimeLable <- colnames(Data[[1]])
  idx0 <- which(TimeLable == Y_First)
  idx1 <- which(TimeLable == Y_Last)

  WgvarSubSample <- WgvarAllYears[idx0:idx1]
  Wgvar <- apply(simplify2array(WgvarSubSample), 1:2, mean)

  colnames(Wgvar) <- Economies
  rownames(Wgvar) <- Economies
  }

# c.3) Specific Year
if (type !="Sample Mean" & type !="Full Sample" ){
  Wgvar <- WgvarAllYears[[type]]
  colnames(Wgvar) <- Economies
  rownames(Wgvar) <- Economies

  }


return(Wgvar)
}

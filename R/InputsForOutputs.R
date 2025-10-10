#' Collects the inputs that are used to construct the numerical and the graphical outputs
#'
#' @param ModelType A character vector indicating the model type to be estimated.
#' @param Horiz A numeric scalar specifying the desired analysis horizon for the outputs.
#' @param ListOutputWished A list of desired graphical outputs. Available options are: "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremia".
#' @param OutputLabel A string for the name of the output label to be stored.
#' @param WishStationarityQ A binary variable (1 or 0) indicating whether to impose that the largest eigenvalue under Q is strictly
#'                          smaller than 1. Set to 1 to impose the restriction, or 0 otherwise.
#' @param DataFrequency A character vector specifying the data frequency. Available options: "Daily All Days", "Daily Business Days",
#'                      "Weekly", "Monthly", "Quarterly", "Annually".
#' @param WishGraphYields A binary variable (1 or 0) indicating whether the user wishes to generate graphs for yields. Default is 0.
#' @param WishGraphRiskFactors A binary variable (1 or 0) indicating whether the user wishes to generate graphs for risk factors. Default is 0.
#' @param WishOrthoJLLgraphs A binary variable (1 or 0) indicating whether the user wishes to generate orthogonalized JLL-based graphs.
#'                           Default is 0.
#' @param WishForwardPremia A binary variable (1 or 0) indicating whether the user wishes to generate forward premia graphs. Default is 0.
#' @param LimFP A numeric vector containing the maturities associated with the start and end dates of the loan.
#' @param WishBootstrap A binary variable (1 or 0) indicating whether the user wishes to perform bootstrap-based estimation. Default is 0.
#' @param ListBoot A List containing the following four elements:
#' \enumerate{
#'  \item \code{methodBS}: Desired bootstrap method: (a) 'bs' for standard residual bootstrap, (b) 'wild' for wild bootstrap,
#'                    or (c) 'block' for block bootstrap.
#'  \item \code{BlockLength}: If block bootstrap is chosen, specify the block length (numeric scalar).
#'  \item \code{ndraws}: Number of bootstrap draws.
#'  \item \code{pctg}: Confidence level expressed in basis points (numeric vector).
#' }
#' @param WishForecast A binary variable (1 or 0) indicating whether the user wishes to generate forecasts. Default is 0.
#' @param ListForecast A list containing the following three elements:
#' \enumerate{
#'  \item \code{ForHoriz}: forecast horizon;
#'  \item \code{t0Sample}: Index of the first variable in the information set.
#'  \item \code{t0Forecast}: Index of the first forecast cut-off date.
#'  \item \code{ForType}: A string specifying the desired forecast type. Available options are: "Rolling" or "Expanding".
#' }
#' @param UnitYields A character string indicating the maturity unit of yields. Options are: (i) "Month" for yields expressed in months, or (ii) "Year" for yields expressed in years. Default is "Month".
#'
#' @examples
#'
#' ModelType <- "JPS original"
#' Horiz <- 100
#' DesiredOutputGraphs <- c("Fit", "GIRF", "GFEVD")
#' OutputLabel <- "Test"
#' WishStationarityQ <- 1
#' WishGraphRiskFac <- 0
#' WishGraphYields <- 1
#'
#' InputsList <- InputsForOutputs(
#'   ModelType, Horiz, DesiredOutputGraphs, OutputLabel,
#'   WishStationarityQ, WishGraphYields, WishGraphRiskFac
#' )
#' @returns
#' List of necessary inputs to generate the graphs of the outputs of the desired model
#' @export

InputsForOutputs <- function(ModelType, Horiz, ListOutputWished, OutputLabel, WishStationarityQ, DataFrequency,
                             WishGraphYields = 0, WishGraphRiskFactors = 0, WishOrthoJLLgraphs = 0,
                             WishForwardPremia = 0, LimFP = NULL, WishBootstrap = 0, ListBoot = NULL,
                             WishForecast = 0, ListForecast = NULL, UnitYields = "Month") {
  OutputTypeSet <- c("RiskFactors", "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremia", "ForwardPremia")
  IdxWishOut <- match(ListOutputWished, OutputTypeSet)

  if ("RiskFactors" %in% ListOutputWished) {
    WishRFGraph <- 1
  }
  if ("Fit" %in% ListOutputWished) {
    WishFitGraph <- 1
  }
  if ("TermPremia" %in% ListOutputWished) {
    WishTPGraph <- 1
  }

  InputsForOutputs <- list(
    "Label Outputs" = OutputLabel,
    StationaryQ = WishStationarityQ,
    UnitMatYields = UnitYields,
    DataFreq = DataFrequency,
    ForwardPremia = WishForwardPremia
  )

  # Initialization
  for (h in 1:length(OutputTypeSet)) {
    if (h == 1) {
      InputsForOutputs[[ModelType]]$RiskFactors$WishGraphs <- 0
    }
    if (h == 2) {
      InputsForOutputs[[ModelType]]$Fit$WishGraphs <- 0
    }
    if (h >= 3 & h <= 6) {
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$horiz <- Horiz

      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactors <- 0
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$Yields <- 0
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactorsBootstrap <- 0
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$YieldsBootstrap <- 0

      if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactors <- 0
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$Yields <- 0
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactorsBootstrap <- 0
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$YieldsBootstrap <- 0
      }
    }
    if (h > 6) {
      InputsForOutputs[[ModelType]]$RiskPremia$WishGraphs <- 0
      InputsForOutputs[[ModelType]]$ForwardPremia$Limits <- LimFP
    }
  }


  # Filling up the desired outputs
  for (h in IdxWishOut) {
    if (h == 1) {
      InputsForOutputs[[ModelType]]$RiskFactors$WishGraphs <- WishRFGraph
    }
    if (h == 2) {
      InputsForOutputs[[ModelType]]$Fit$WishGraphs <- WishFitGraph
    }

    if (h >= 3 & h <= 6) {
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$horiz <- Horiz

      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactors <- WishGraphRiskFactors
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$Yields <- WishGraphYields
      if (WishBootstrap == 1) {
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactorsBootstrap <- WishGraphRiskFactors
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$YieldsBootstrap <- WishGraphYields
      }

      if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
        if (WishOrthoJLLgraphs == 1) {
          InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactors <- WishGraphRiskFactors
          InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$Yields <- WishGraphYields
          if (WishBootstrap == 1) {
            InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactorsBootstrap <- WishGraphRiskFactors
            InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$YieldsBootstrap <- WishGraphYields
          }
        }
      }
    }

    if (h == 7) {
      InputsForOutputs[[ModelType]]$RiskPremia$WishGraphs <- WishTPGraph
    }
  }

  # 2) Bootstrap list
  if (WishBootstrap == 0) {
    InputsForOutputs[[ModelType]]$Bootstrap$WishBoot <- 0
  } else {
    InputsForOutputs[[ModelType]]$Bootstrap <- WishBootstrap
    names(InputsForOutputs[[ModelType]]$Bootstrap) <- "WishBoot"
    InputsForOutputs[[ModelType]]$Bootstrap <- append(InputsForOutputs[[ModelType]]$Bootstrap, ListBoot)
  }


  # 3) Forecasting list
  if (WishForecast == 0) {
    InputsForOutputs[[ModelType]]$Forecasting$WishForecast <- 0
  } else {
    InputsForOutputs[[ModelType]]$Forecasting <- WishForecast
    names(InputsForOutputs[[ModelType]]$Forecasting) <- "WishForecast"
    InputsForOutputs[[ModelType]]$Forecasting <- append(InputsForOutputs[[ModelType]]$Forecasting, ListForecast)
  }


  return(InputsForOutputs)
}

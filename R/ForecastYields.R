#' Gather bond yields forecasts for all the model types
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param ModelPara  List of model parameter estimates (See the "Optimization" function)
#'@param InputsForOutputs list conataining the desired horizon of analysis for the IRFs, GIRFs, FEVDs, and GFEVDs
#'@param FactorLabels  a string-list based which contains all the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param DataFrequency text: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"
#'@param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'@param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#'

#'@examples
#' # See examples in the vignette file of this package (Section 4).
#'
#'
#'@returns
#' List containg the following elements
#' \enumerate{
#' \item Out-of-sample forecasts of bond yields per forecast horizon
#' \item Out-of-sample forecast errors of bond yields per forecast horizon
#' \item Root mean square errors per forecast horizon
#'}
#'
#'@export

ForecastYields <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies, DataFrequency,
                               JLLinputs, GVARinputs){


  WishForecast <- InputsForOutputs[[ModelType]]$Forecasting$WishForecast

  if (WishForecast ==0){ print( "No bond yields forecasts were generated")
    }else{
      Jmisc::tic()
  # If one chooseS models in which the estimation is done country-by-country
  if ( "JPS" %in% ModelType || "JPS jointP" %in% ModelType ||  "GVAR sepQ" %in% ModelType){
    ForeYie <- ForecastYieldsSepQ(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies, DataFrequency,
                                  JLLinputs, GVARinputs)

  }
  # If one chooseS models in which the estimation is done jointly for all countries
  if ( "GVAR jointQ" %in% ModelType || "VAR jointQ" %in% ModelType || "JLL original" %in% ModelType
       || "JLL NoDomUnit" %in% ModelType || "JLL jointSigma" %in% ModelType){

    ForeYie <- ForecastYieldsJointQ(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies, DataFrequency,
                                  JLLinputs, GVARinputs)

  }
      Jmisc::toc()
      return(ForeYie)
    }

}

################################################################################################################
#' Bond yields forecasts  ("sep Q" models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param ModelPara  List of model parameter estimates (See the "Optimization" function)
#'@param InputsForOutputs list conataining the desired horizon of analysis for the IRFs, GIRFs, FEVDs, and GFEVDs
#'@param FactorLabels  a string-list based which contains all the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param DataFrequency character-based vector: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"
#'@param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'@param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#'
#'


ForecastYieldsSepQ <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies, DataFrequency,
                               JLLinputs, GVARinputs){



  print('#########################################################################################################')
  print( paste('#################################', 'Forecasting', ModelType, '#################################' ))
  print('#########################################################################################################')

  # 1) Redefine some general model outputs
  StationarityUnderQ <- InputsForOutputs$StationaryQ
  t0Sample <- InputsForOutputs[[ModelType]]$Forecasting$t0Sample
  t0Forecast <- InputsForOutputs[[ModelType]]$Forecasting$t0Forecast
  H <- InputsForOutputs[[ModelType]]$Forecasting$ForHoriz


  if (t0Forecast < t0Sample ){stop("The first forecast cut-off date is earlier than the start of the sample.")}

  FullModelParaList <- list()
  OutofSampleForecast <- list()


  T <- ncol(ModelPara[[ModelType]][[Economies[1]]]$inputs$Y)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  G <- length(FactorLabels$Global)
  C <- length(Economies)

  dt <- ModelPara[[ModelType]][[Economies[1]]]$inputs$dt # peridiocity of the data in years (i.e. if monthly dt =1/12)
  mat <- ModelPara[[ModelType]][[Economies[1]]]$inputs$mat
  J <- length(mat)

  nloops <- T- t0Forecast- H + 1 # Number of times that the model will be estimated


  ###################################################################################################
  if (t0Forecast > T ){stop("The first forecast cut-off date is longer than the sample length.")}

  for (tt in 1:nloops){

    if (nloops <= 0){stop("Impossible to generate forecast errors: sample period is extrapolated!")}

    if (tt ==1){
      tlastObserved <- t0Forecast
    }else{
      tlastObserved <- tlastObserved + 1
    }

    Ttemp <- tlastObserved  # Time dimension of the model (updated after each iteration)


    # 2) Model Estimation

    for (i in 1:C){

      # Redefine the dataset used in the estimation
      ZZfull <- ModelPara[[ModelType]][[Economies[i]]]$inputs$AllFactors
      YieldsFull <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y

      i <<- i # Re-define i in the global envirnonment
    # Yields
    YYtemp <- YieldsFull[, t0Sample:tlastObserved]
    # Spanned factors
    PPall <- SpannedFactorsSepQ(ModelType, ModelPara, Economies, t0Sample, tlastObserved)
    # All risk factors
    IdxSpa <- IdxAllSpanned(ModelType, FactorLabels, Economies)
    ZZtemp <- ZZfull[ , t0Sample:tlastObserved]
    ZZtemp[IdxSpa, ] <- PPall

    # For the GVAR-based models
    if (ModelType == 'GVAR sepQ'){
      GVARinputs$GVARFactors <- DataSet_BS(ModelType, ZZtemp, GVARinputs$Wgvar, Economies, FactorLabels)
      #) NOTE:  To avoid over complicating the code, we keep the transition matrix as the one
      #  estimated conditionally on the full information set
    }
    # Compute the inputs that go directly into the log-likelihood function
    ATSMInputs <- InputsForMLEdensity(ModelType, YYtemp, ZZtemp, FactorLabels, mat, Economies, DataFrequency,
                                        JLLinputs, GVARinputs)

  # Initial guesses for Variables that will be concentrared out of from the log-likelihood function
      K1XQ <- ATSMInputs$K1XQ
      SSZ <- ATSMInputs$SSZ

  # Build the objective function
  f <- Functionf(ATSMInputs, Economies, mat, DataFrequency, FactorLabels, ModelType)

      # Choose the optimization settings
      VarLab <- ParaLabels(ModelType, StationarityUnderQ)

      varargin <- list()
      varargin$K1XQ <-list(K1XQ, VarLab[[ModelType]][["K1XQ"]] , NULL , NULL)
      varargin$SSZ <- list(SSZ, VarLab[[ModelType]][["SSZ"]], NULL, NULL)
      varargin$r0 <- list(NULL, VarLab[[ModelType]][["r0"]], NULL, NULL)
      varargin$se <- list(NULL, VarLab[[ModelType]][["se"]], 1e-6, NULL)
      varargin$K0Z <- list(NULL, VarLab[[ModelType]][["K0Z"]], NULL, NULL)
      varargin$K1Z <- list(NULL, VarLab[[ModelType]][["K1Z"]], NULL, NULL)
      varargin$OptRun <-  c("iter off")

      LabelVar<- c('Value', 'Label', 'LB', 'UB') # Elements of each parameter
      for (d in 1:(length(varargin)-1)){ names(varargin[[d]]) <-  LabelVar}

      tol <- 1e-4

      FullModelParaList[[ModelType]][[Economies[i]]] <- Optimization(f, tol, varargin, FactorLabels, Economies, ModelType)$Summary


      # 3) Forecasting

      # Define general inputs
      A <- FullModelParaList[[ModelType]][[Economies[i]]]$rot$P$A
      Bfull <- BUnspannedAdapSep(G,M, FullModelParaList, Economies, Economy = Economies[i], ModelType)
      K0Z <- FullModelParaList[[ModelType]][[Economies[i]]]$ests$K0Z
      K1Z <- FullModelParaList[[ModelType]][[Economies[i]]]$ests$K1Z

      # Forecast of yields
      ForecastYields <- matrix(NA, nrow = J, ncol = H)
      LabelForecastPeriod <- colnames(YieldsFull[,(Ttemp+1):(Ttemp+H)])
      rownames(ForecastYields) <- rownames(YieldsFull)
      colnames(ForecastYields) <- LabelForecastPeriod

      ZZtt <- ZZtemp[,Ttemp]
      K1ZsumOld <- 0

      for (hh in 1:H){
        K1Znew <- powerplus::Matpow(K1Z, numer = hh -1 )
        VARforecast <- (K1ZsumOld + K1Znew)%*%K0Z + powerplus::Matpow(K1Z, numer = hh)%*%ZZtt
        ForecastYields[,hh] <- A + Bfull%*%(VARforecast)

        K1ZsumOld <- K1ZsumOld + K1Znew
      }

      # Actual yields for the period of the forecasting
      YieldsObsForPer <- YieldsFull[,(Ttemp+1):(Ttemp+H)]

      # 4) Forecast error
      ForecastError <- YieldsObsForPer - ForecastYields

      ForecastDate <- colnames(ZZfull)[Ttemp]
      OutofSampleForecast[[ModelType]][[Economies[i]]][[ForecastDate]]$Forcast <- ForecastYields
      OutofSampleForecast[[ModelType]][[Economies[i]]][[ForecastDate]]$Error <- ForecastError



      print(paste(ModelType, Economies[i], ": Out-of-sample forecast for information set available until", ForecastDate))

      saveRDS(OutofSampleForecast, paste(tempdir(),"/Forecast_", InputsForOutputs$'Label Outputs','.rds',sep=""))
    }
  }



  # 5) RMSE
  RMSE <- list(RMSEsep(OutofSampleForecast))
  names(RMSE) <- "RMSE"
  OutofSampleForecast <- append(OutofSampleForecast[[ModelType]], RMSE)


  saveRDS(OutofSampleForecast, paste(tempdir(),"/Forecast_", InputsForOutputs$'Label Outputs','.rds',sep=""))

  return(OutofSampleForecast)

  }

###############################################################################################################
#' Bond yields forecasts ("joint Q" models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param ModelPara  List of model parameter estimates (See the "Optimization" function)
#'@param InputsForOutputs list conataining the desired horizon of analysis for the IRFs, GIRFs, FEVDs, and GFEVDs
#'@param FactorLabels  a string-list based which contains all the labels of all the variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param DataFrequency character-based vector: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"
#'@param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'@param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#'
#'


ForecastYieldsJointQ <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies, DataFrequency,
                               JLLinputs, GVARinputs){

  print('#########################################################################################################')
  print( paste('#################################', 'Forecasting', ModelType, '#################################' ))
  print('#########################################################################################################')

  # 1) Redefine some general model outputs
  StationarityUnderQ <- InputsForOutputs$StationaryQ
  t0Sample <- InputsForOutputs[[ModelType]]$Forecasting$t0Sample
  t0Forecast <- InputsForOutputs[[ModelType]]$Forecasting$t0Forecast
  H <- InputsForOutputs[[ModelType]]$Forecasting$ForHoriz

  if (t0Forecast < t0Sample ){stop("The first forecast cut-off date is earlier than the start of the sample.")}

  FullModelParaList <- list()
  OutofSampleForecast <- list()

  T <- ncol(ModelPara[[ModelType]]$inputs$Y)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  G <- length(FactorLabels$Global)
  C <- length(Economies)

  dt <- ModelPara[[ModelType]]$inputs$dt # peridiocity of the data in years (i.e. if monthly dt =1/12)
  mat <- ModelPara[[ModelType]]$inputs$mat
  J <- length(mat)

  nloops <- T- t0Forecast- H + 1 # Number of times that the model will be re-estimated


  ###################################################################################################
  if (t0Forecast > T ){stop("The first forecast cut-off date is longer than the sample length.")}

  for (tt in 1:nloops){

    if (nloops <= 0){stop("Impossible to generate forecast errors: sample period is extrapolated!")}

    if (tt ==1){
      tlastObserved <- t0Forecast
    }else{
      tlastObserved <- tlastObserved + 1
    }

    Ttemp <- tlastObserved  # Time dimension of the model (updated after each iteration)


    # 2) Model Estimation

      # Redefine the dataset used in the estimation
      ZZfull <- ModelPara[[ModelType]]$inputs$AllFactors
      YieldsFull <- ModelPara[[ModelType]]$inputs$Y

      # Yields
      YYtemp <- YieldsFull[, t0Sample:tlastObserved]
      # Spanned factors
      PPall <- SpannedFactorsjointQ(ModelType, ModelPara, Economies, t0Sample, tlastObserved)
      # All risk factors
      IdxSpa <- IdxAllSpanned(ModelType, FactorLabels, Economies)
      ZZtemp <- ZZfull[ , t0Sample:tlastObserved]
      ZZtemp[IdxSpa, ] <- PPall

      # For the GVAR-based models
      if (ModelType == 'GVAR jointQ'){
        GVARinputs$GVARFactors <- DataSet_BS(ModelType, ZZtemp, GVARinputs$Wgvar, Economies, FactorLabels)
        #) NOTE:  To avoid over complicating the code, we keep the transition matrix as the one
        #  estimated conditionally on the full information set
      }
      # Compute the inputs that go directly into the log-likelihood function
      ATSMInputs <- InputsForMLEdensity(ModelType, YYtemp, ZZtemp, FactorLabels, mat, Economies, DataFrequency,
                                        JLLinputs, GVARinputs)

      # Initial guesses for Variables that will be concentrared out of from the log-likelihood function
      K1XQ <- ATSMInputs$K1XQ
      if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit" ){ SSZ <- NULL} else{SSZ <- ATSMInputs$SSZ}

      # Build the objective function
      f <- Functionf(ATSMInputs, Economies, mat, DataFrequency, FactorLabels, ModelType)

      # Choose the optimization settings
      VarLab <- ParaLabels(ModelType, StationarityUnderQ)

      varargin <- list()
      varargin$K1XQ <-list(K1XQ, VarLab[[ModelType]][["K1XQ"]] , NULL , NULL)
      varargin$SSZ <- list(SSZ, VarLab[[ModelType]][["SSZ"]], NULL, NULL)
      varargin$r0 <- list(NULL, VarLab[[ModelType]][["r0"]], NULL, NULL)
      varargin$se <- list(NULL, VarLab[[ModelType]][["se"]], 1e-6, NULL)
      varargin$K0Z <- list(NULL, VarLab[[ModelType]][["K0Z"]], NULL, NULL)
      varargin$K1Z <- list(NULL, VarLab[[ModelType]][["K1Z"]], NULL, NULL)
      varargin$OptRun <-  c("iter off")

      LabelVar<- c('Value', 'Label', 'LB', 'UB') # Elements of each parameter
      for (d in 1:(length(varargin)-1)){ names(varargin[[d]]) <-  LabelVar}

      tol <- 1e-4


      invisible(utils::capture.output(FullModelParaList[[ModelType]] <- Optimization(f, tol, varargin, FactorLabels,
                                                     Economies, ModelType, JLLinputs)$Summary))


      # 3) Forecasting

      # Define general inputs
      A <- FullModelParaList[[ModelType]]$rot$P$A
      K0Z <- FullModelParaList[[ModelType]]$ests$K0Z
      K1Z <- FullModelParaList[[ModelType]]$ests$K1Z

      Bspanned <- FullModelParaList[[ModelType]]$rot$P$B
      Bfull <- BUnspannedAdapJoint(G,M, N, C, J, Bspanned)

      # Forecast of yields
      ForecastYields <- matrix(NA, nrow = C*J, ncol = H)
      LabelForecastPeriod <- colnames(YieldsFull[,(Ttemp+1):(Ttemp+H)])
      rownames(ForecastYields) <- rownames(YieldsFull)
      colnames(ForecastYields) <- LabelForecastPeriod

      ZZtt <- ZZtemp[,Ttemp]
      K1ZsumOld <- 0

      for (hh in 1:H){
        K1Znew <- powerplus::Matpow(K1Z, numer = hh -1 )
        VARforecast <- (K1ZsumOld + K1Znew)%*%K0Z + powerplus::Matpow(K1Z, numer = hh)%*%ZZtt
        ForecastYields[,hh] <- A + Bfull%*%(VARforecast)

        K1ZsumOld <- K1ZsumOld + K1Znew
      }

      # Actual yields for the period of the forecasting
      YieldsObsForPer <- YieldsFull[,(Ttemp+1):(Ttemp+H)]

      # 4) Forecast error
      ForecastError <- YieldsObsForPer - ForecastYields

      ForecastDate <- colnames(ZZfull)[Ttemp]
      OutofSampleForecast[[ModelType]][[ForecastDate]]$Forcast <- ForecastYields
      OutofSampleForecast[[ModelType]][[ForecastDate]]$Error <- ForecastError



      print(paste(ModelType, ": Out-of-sample forecast for information set available until", ForecastDate))

      saveRDS(OutofSampleForecast, paste(tempdir(),"/Forecast_", InputsForOutputs$'Label Outputs','.rds',sep=""))
  }



  # 5) RMSE
  RMSE <- list(RMSEjoint(OutofSampleForecast))
  names(RMSE) <- "RMSE"
  OutofSampleForecast <- append(OutofSampleForecast[[ModelType]], RMSE)


  saveRDS(OutofSampleForecast, paste(tempdir(),"/Forecast_", InputsForOutputs$'Label Outputs','.rds',sep=""))

  return(OutofSampleForecast)

}


################################################################################################################
#' Compute the root mean square error ("sep Q" models)
#'
#'@param ForecastOutputs  List of country-specific forecasts (see "ForecastYieldsSepQ" function)
#'


RMSEsep <- function(ForecastOutputs){

  nfor <- length(ForecastOutputs[[1]][[1]])
  Economies <- names(ForecastOutputs[[1]])
  ModelType <- names(ForecastOutputs)

  H <- ncol(ForecastOutputs[[1]][[1]][[1]][[1]])
  C <- length(Economies)


  FElist <- list()
    for (i in 1:C){
      for (h in 1:nfor){
        FElist[[Economies[i]]][[h]] <- ForecastOutputs[[ModelType]][[Economies[i]]][[h]]$Error
      }
    }



  rmse <- list()
    for (i in 1:C){
      rmse[[Economies[i]]] <- sqrt(Reduce("+", lapply(FElist[[Economies[i]]],
                                                                     function(x, N= nfor) x^2/N)))
      colnames(rmse[[Economies[i]]]) <- 1:H
    }



  return(rmse)
}


##############################################################################################################
#' Compute the root mean square error ("joint Q" models)
#'
#'@param ForecastOutputs  List of country-specific forecasts (see "ForecastYieldsjointQ" function)
#'


RMSEjoint <- function(ForecastOutputs){

  nfor <- length(ForecastOutputs[[1]])
  Modeljoint <- names(ForecastOutputs)

  H <- ncol(ForecastOutputs[[1]][[1]][[1]])

  FElist <- list()
    for (h in 1:nfor){
      FElist[[Modeljoint]][[h]] <- ForecastOutputs[[Modeljoint]][[h]]$Error
    }

    rmse <- list()
    rmse[[Modeljoint]] <- sqrt(Reduce("+", lapply(FElist[[Modeljoint]], function(x, N= nfor) x^2/N)))
    colnames(rmse[[Modeljoint]]) <- 1:H


  return(rmse)
}


################################################################################################################
#' Find the indexes of the spanned factors
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'

IdxAllSpanned <- function(ModelType, FactorLabels, Economies){

  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)

  IdxSpanned <- c()

  if (ModelType== "JPS"){ IdxSpanned <- (G+M+1):(G+M+N) } else{
    idxSpa0 <- G + M
    for (j in 1:C){
      idxSpa1 <- idxSpa0 + N

      if (j ==1){ IdxSpanned <- (idxSpa0+1):idxSpa1
      }  else{
        IdxSpanned <- c(IdxSpanned, (idxSpa0+1):idxSpa1)
      }
      idxSpa0 <- idxSpa1 + M
    }
  }

  return(IdxSpanned)
}

#################################################################################################################
#' Gather all spanned factors ("sep Q" models)
#'
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelPara  set of model parameters
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'@param t0Sample   index for the initial sample date
#'@param tlastObserved index for the last observation of the information set
#'

SpannedFactorsSepQ <- function(ModelType, ModelPara, Economies, t0Sample, tlastObserved){

  C <- length(Economies)
  N <- ModelPara[[ModelType]][[Economies[1]]]$inputs$N

  PPALL <- matrix()

  if (ModelType== "JPS"){
    i <- get("i", globalenv())
    YCS <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y[, t0Sample:tlastObserved]
    PPALL <- Spanned_Factors(YCS, Economies = Economies[i], N)
  } else {
    for (i in 1:C){
      YCS <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y[,t0Sample:tlastObserved]
      if (i==1 ){PPALL <- Spanned_Factors(YCS, Economies = Economies[i], N)
      } else{
        PPtemp <- Spanned_Factors(YCS, Economies = Economies[i], N)
        PPALL <- rbind(PPALL, PPtemp)
      }
    }
  }

  return(PPALL)
}

#######################################################################################################
#' Gather all spanned factors ("joint Q" models)
#'
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelPara  set of model parameters
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param t0Sample   index for the initial sample date
#'@param tlastObserved index for the last observation of the information set
#'

SpannedFactorsjointQ <- function(ModelType, ModelPara, Economies, t0Sample, tlastObserved){

  C <- length(Economies)
  N <- ModelPara[[ModelType]]$inputs$N
  J <- length(ModelPara[[ModelType]]$inputs$mat)

  PPALL <- matrix()

  idxJ0 <- 0

  for (i in 1:C){ # Country-specific inputs
    idxJ1 <- idxJ0 + J

    YCS <- ModelPara[[ModelType]]$inputs$Y[(idxJ0+1):idxJ1, t0Sample:tlastObserved][,]

  if (i==1 ){PPALL <- Spanned_Factors(YCS, Economies = Economies[i], N)
    }else{
        PPtemp <- Spanned_Factors(YCS, Economies = Economies[i], N)
        PPALL <- rbind(PPALL, PPtemp)
      }

    idxJ0<- idxJ1

    }


  return(PPALL)
}


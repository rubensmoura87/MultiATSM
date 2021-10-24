#' Set up the vector-valued objective function (Point estimate)

#'@param MLEinputs Set of inputs that are necessary to the log-likelihood function
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param mat vector of maturities (in years) of yields used in estimation (J x 1)
#'@param DataFrequency  character-based vector: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'
#'
#'@examples
#'\dontrun{
#' # See examples in the vignette file of this package (Section 4).
#'}
#'
#'@returns
#'objective function
#'
#'@export



Functionf<- function(MLEinputs, Economies, mat, DataFrequency, FactorLabels, ModelType){


  if (DataFrequency == "Daily All Days"){ dt <- 1/365}
  if (DataFrequency == "Daily Business Days"){ dt <- 1/252}
  if (DataFrequency == "Weekly"){ dt <- 1/52}
  if (DataFrequency == "Monthly"){ dt <- 1/12}
  if (DataFrequency == "Quarterly"){ dt <- 1/4}
  if (DataFrequency == "Annually"){ dt <- 1}


# 1) If one choose models in which the estimation is done country-by-country
if (ModelType == 'JPS' || ModelType == 'JPS jointP' || ModelType == "GVAR sepQ"){

  i <- get("i", globalenv())

  f <-functional::Curry(MLEdensity_sepQ, r0 = NULL, MLEinputs$K0Z, MLEinputs$K1Z, se= NULL, MLEinputs$Gy.0, mat, MLEinputs$Y,
            MLEinputs$ZZ, MLEinputs$PP, MLEinputs$Wpca, MLEinputs$We, MLEinputs$WpcaFull, dt, Economies[i], FactorLabels,
            ModelType, MLEinputs$GVARinputs)
}



# 2) If one choose models in which the estimation is done jointly for all countries of the system
if (ModelType == 'GVAR jointQ' || ModelType == 'VAR jointQ' ||  ModelType == "JLL jointSigma") {

  if(length(Economies) == 1){stop('The chosen model type has to be estimated for multiple countries at the time.
                                The input "Economies" should contain several members.')}

    f <-functional::Curry(MLEdensity_jointQ, r0 = NULL, MLEinputs$K0Z, MLEinputs$K1Z, se= NULL, MLEinputs$Gy.0, mat, MLEinputs$Y,
              MLEinputs$ZZ, MLEinputs$PP, MLEinputs$Wpca, MLEinputs$We, MLEinputs$WpcaFull, dt, Economies,
              FactorLabels, ModelType, MLEinputs$GVARinputs, MLEinputs$JLLinputs)
}



# 3) If one choose models in which the estimation is done jointly for all countries of the system with
  # the Sigma matrix estimated exclusively under the P-dynamics
if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit") {
  if(length(Economies) == 1){stop('The chosen model type has to be estimated for multiple countries at the time.
                                The input "Economies" should contain several members.')}

  f <- functional::Curry(MLEdensity_jointQ_sepSigma, r0 = NULL, MLEinputs$SSZ, MLEinputs$K0Z, MLEinputs$K1Z, se= NULL, MLEinputs$Gy.0,
            mat, MLEinputs$Y, MLEinputs$ZZ, MLEinputs$PP, MLEinputs$Wpca, MLEinputs$We, MLEinputs$WpcaFull, dt,
            Economies, FactorLabels, ModelType, MLEinputs$JLLinputs)
}




return(f)
}


############################################################################################################
#' Set up the vector-valued objective function (Bootstrap)

#'@param ModelType string-vector containing the label of the model to be estimated
#'@param MLEinputsBS Set of inputs that are necessary to the log-likelihood function
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param mat vector of maturities (in years) of yields used in estimation (J x 1)
#'@param dt   adjusted yearly frequency of the data
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param residBS   indexes of the re-ordered bootstrap residuals
#'@param MaxEigen   largest eigenvalue under the P-dynamics
#'@param JLLinputs   necessary inputs for the estimation of JLL-based models
#'@param GVARinputs  necessary inputs for the estimation of GVAR-based models
#'


Functionf_Boot <- function(ModelType, MLEinputsBS, Economies, mat, dt, FactorLabels, residBS, MaxEigen,
                            JLLinputs, GVARinputs){



  Y_artificial <- MLEinputsBS$Y
  ZZ_artificial <- MLEinputsBS$ZZ
  PP_artificial <- MLEinputsBS$PP
  Wpca_artificial <- MLEinputsBS$Wpca
  We_artificial <- MLEinputsBS$We
  WpcaFull_artificial <- MLEinputsBS$WpcaFull
  K1XQ <- MLEinputsBS$K1XQ
  SSZ <- MLEinputsBS$SSZ
  K0Z <- MLEinputsBS$K0Z
  K1Z <- MLEinputsBS$K1Z
  Gy.0 <- MLEinputsBS$Gy.0

# 1) If one choose models in which the estimation is done country-by-country
if (ModelType == 'JPS' || ModelType == 'JPS jointP' || ModelType == "GVAR sepQ"){
  i <- get("i", globalenv())
   f <-functional::Curry(A0N_MLEdensity_WOE__sepQ_Bootstrap, r0 = NULL, K0Z, K1Z, se= NULL, Gy.0, mat, Y_artificial,
            ZZ_artificial, PP_artificial, Wpca_artificial, We_artificial, WpcaFull_artificial, dt, Economies[i],
            FactorLabels, ModelType, residBS,  MaxEigen, GVARinputs)
}

  # 2) If one choose models in which the estimation is done jointly for all countries of the system
if (ModelType == 'GVAR jointQ' || ModelType == 'VAR jointQ' ||  ModelType == "JLL jointSigma") {
  f <-functional::Curry(A0N_MLEdensity_WOE__jointQ_Bootstrap, r0= NULL, K0Z, K1Z, se= NULL, Gy.0,  mat, Y_artificial,
            ZZ_artificial, PP_artificial, Wpca_artificial, We_artificial, WpcaFull_artificial, dt, Economies,
            FactorLabels, ModelType, residBS, MaxEigen, GVARinputs, JLLinputs)
}


# 3) If one choose models in which the estimation is done jointly for all countries of the system
if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit") {
  f <-functional::Curry(A0N_MLEdensity_WOE__jointQ_sepSigma_Bootstrap, r0 = NULL, SSZ, K0Z, K1Z, se= NULL, Gy.0, mat,
            Y_artificial, ZZ_artificial, PP_artificial, Wpca_artificial, We_artificial, WpcaFull_artificial,
            dt, Economies, FactorLabels, ModelType, residBS, MaxEigen, GVARinputs, JLLinputs)
}

  return(f)
}

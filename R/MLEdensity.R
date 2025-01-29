#' Compute the maximum likelihood function of all models
#'
#'@param K1XQ risk-neutral feedback matrix (N x N or CN x CN)
#'@param r0   long-run interest rate  (scalar or vector with length C)
#'@param SSZ  variance-covariance matrix (F x F)
#'@param K0Z  intercept from the P-dynamics (F x 1)
#'@param K1Z  feedback matrix from the P-dynamics (F x F)
#'@param se   Variance of the portfolio of yields observed with error (scalar). Default is set to NULL.
#'@param Gy.0 matrix of contemporaneous terms from the P-dynamics (F x F)
#'@param mat  vector of maturities (in years) of yields used in estimation (J x 1)
#'@param Y    matrix of yields used in estimation  (J x T or CJ x T)
#'@param Z    complete set of spanned and unspanned factors (F x T)
#'@param P    complete set of spanned factors (N x T or CN x T)
#'@param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#'@param We   matrix of weights of the portfolios observed with errors ((J-N) x J or C(J-N) x CJ)
#'@param WpcaFull composite matrix of weights the portfolios observed with and without errors
#'@param dt  time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param GVARinputs if the model chosen is the "GVAR single" or "GVAR multi", the "GVARinputs" should be specified (see "GVAR" function)
#'@param JLLinputs if the model chosen is JLL-based. "JLLinputs" should contain (i) DomUnit, (ii) WishSigmas, (iii) SigmaNonOrtho, (iv) JLLModelType (See JLL function)
#'@param BS_outputs Generates simplified output list in the bootstrap setting. Default is set to FALSE.
#'@param nargout if nargout== 1: provides only the values of the likelihood; if nargout== 2: complete ATSM outputs
#'
#'@importFrom pracma mrdivide
#'
#'@references
#' This function is modified version of the "A0N_MLEdensity_WOE" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal


MLEdensity <- function(K1XQ, r0, SSZ, K0Z, K1Z, se, Gy.0, mat, Y, Z, P, Wpca, We, WpcaFull, dt, Economies,
                       FactorLabels, ModelType, GVARinputs = NULL, JLLinputs = NULL, BS_outputs = F, nargout){

  # 0) Initialize some variables, if necessary
  if (!exists("r0")){ r0 <- as.numeric()}
  if (!exists("se")){ se <- as.numeric()}

  N <- length(FactorLabels$Spanned) # number of country-specific spanned factors
  if (ModelType == 'JLL joint Sigma'){ SSZ <- Update_SSZ_JLL(SSZ, Z, N, JLLinputs)}
  # x is a column vector that contains all the elements of each parameter of the model.
  x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, t(t(as.vector(K0Z))), t(t(as.vector(K1Z))) )

  if ( any(is.nan(x)) ||any(is.infinite(x)) || any(Im(x)!= 0) ){
    browser()
    y <- 1e6*rep(1, times = T-1)
    out <- c()
  }else{

    # 1) Compute loadings:
    # Slope coefficients (Bn's)
    LoadBs <- Get_Bs(mat, dt, K1XQ, SSZ, Wpca, FactorLabels, Economies, ModelType)
    # Intercepts (An's):
    if (!exists("r0")|| sjmisc::is_empty(r0) ){ r0 <- Get_r0(Y, P, N, mat, dt, LoadBs, Wpca, We, Economies, ModelType)  }
    LoadAs <- Get_As(LoadBs, Wpca, r0, dt, Economies, ModelType)

    # 2) Build Log-likelihood density function
    llkOut <- Get_llk(P, Y, Z, N, mat, We, Wpca, K0Z, K1Z, SSZ, LoadBs, LoadAs, ModelType)
    y <- llkOut$y
    se <- llkOut$se

    # 3) Variance-yields
    VarYields <- Get_SigmaYields(Y, N, mat, WpcaFull, se, ModelType)

    # 4) Output to export
    if (nargout>1){
        Out <- OptOutputs(Y, Z , mat, N, dt, Wpca, K1XQ, SSZ, LoadAs, LoadBs, r0, se, K0Z, K1Z, Gy.0, VarYields,
                        y, GVARinputs, JLLinputs, Economies, ModelType, BS_outputs)
    }
  }

    if (nargout ==1){ return(y) } else{ return(Out)}
}

##########################################################################################################
##########################################################################################################
#' Extract the indexes related to the spanned factors in the variance-covariance matrix
#'
#'@param G number of global unspanned factors (scalar)
#'@param M number of domestic unspanned factors per country (scalar)
#'@param N number of domestic spanned factors per country (scalar)
#'@param C number of countries of the economic system (scalar)
#'
#'@keywords internal

IdxSpanned <- function(G,M,N,C){

  K <- C*(N+M) + G

  vector <- 1:K

  idxA <- 0
  idxB <- G + M
  IDXX <- c()
  for (i in 1:C){
    idxAA <- idxA + N
    idxBB <- idxB + N
    IDXX[(idxA+1):idxAA] <- vector[(idxB+1):idxBB]
    idxA <- idxAA
    idxB <- idxBB +M
  }

  return(IDXX)
}
##################################################################################################################
##################################################################################################################
#' Compute r0 for the various models
#'
#'@param Y matrix of yields used in estimation  (J x T or CJ x T)
#'@param P complete set of spanned factors (N x T or CN x T)
#'@param N number of country-specific spanned factors
#'@param mat vector of maturities (in years) of yields used in estimation (J x 1)
#'@param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#'@param B_list list containing the B loadings
#'@param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#'@param We matrix of weights of the portfolios observed with errors ((J-N) x J or C(J-N) x CJ)
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'@keywords internal


Get_r0 <- function(Y, P, N, mat, dt, B_list, Wpca, We, Economies, ModelType){

  # General procedure:
  # (i) A0 = (I - Bx(W*Bx)^(-1)*W)*Ax0;
  # (ii) A1 = (I - Bx(W*Bx)^(-1)*W)*Axr;
  # (iii) APer= We * A1
  # (iv) vt = We*(Yt - Bx(W*Bx)^(-1)*Pt - A0)
  # Recall that BnP= Bx(W*Bx)^(-1)

  J <- numel(mat) # number of country-specific yields used in estimation;
  T <- dim(Y)[2]
  t <- 2:T

  BnP <- B_list$BnP
  betan <- B_list$betan

  # 1) Get r0 for the models estimated on a country-by-country basis
  if (any(ModelType == c("JPS original","JPS global", 'GVAR single'))){

    A0 <- (diag(J)-BnP%*%Wpca)%*%betan
    A1 <- (diag(J)-BnP%*%Wpca)%*%matrix(1,J,1)/dt

    #  r0 = APer'*vt/(APer'*APer)
    r0 <- solve((t(A1)%*%(t(We)%*%We)%*%A1), # Numerator from equation r0 equation
                (t(A1)%*%(t(We)%*%We))%*%(rowMeans(Y[,t] - BnP%*%P[,t]) - A0)) # Denominator from r0 equation

  } else{
    # 2) Get r0 for the models estimated on a jointly basis
    C <- length(Economies)

    r0 <- c()

    idxA <- 0
    idxB <- 0
    idxC <- 0

    for (i in 1:C){
      idxAA <- idxA + N
      idxBB <- idxB + J
      idxCC <- idxC + J-N

      BnPCS <- matrix(BnP[(idxB+1):idxBB , (idxA+1):idxAA], nrow = J)
      WpcaCS <- Wpca[(idxA+1):idxAA,(idxB+1):idxBB]
      betanCS <- betan[(idxB+1):idxBB]
      WeCS <- We[(idxC+1):(idxCC),(idxB+1):idxBB]
      YCS <- Y[(idxB+1):idxBB,]
      PCS <- P[(idxA+1):idxAA,]

      A0 <- (diag(J)-BnPCS%*%WpcaCS)%*%betanCS
      A1 <- (diag(J)-BnPCS%*%WpcaCS)%*%matrix(1,J,1)/dt

      # r0 = APer'*vt/(APer'*APer)
      if (any(ModelType == c("JPS multi", "GVAR multi"))){
        PCS_adj <- if (N==1) PCS[t] else PCS[,t]
        r0[i] <- solve( (t(A1)%*%(t(WeCS)%*%WeCS)%*%A1), # Numerator from r0 equation
                        (t(A1)%*%(t(WeCS)%*%WeCS))%*%(rowMeans(YCS[,t] - BnPCS%*%PCS_adj) - A0) ) # Denominator from r0 equation
      } else {
        r0[i] <- solve( (t(A1)%*%(t(WeCS)%*%WeCS)%*%A1), # Numerator from r0 equation
                        (t(A1)%*%(t(WeCS)%*%WeCS))%*%(rowMeans(YCS[,t] - BnPCS%*%PCS[,t]) - A0) ) # Denominator from  r0 equation
      }

      idxA <- idxAA
      idxB <- idxBB
      idxC <- idxCC
    }
  }

  return(r0)
}
##############################################################################################################################
#' BUild the B loadings
#'
#'@param mat vector of maturities (in years) of yields used in estimation (J x 1)
#'@param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1
#'@param K1XQ risk-neutral feedback matrix (N x N or CN x CN)
#'@param SSZ variance-covariance matrix (F x F)
#'@param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param Economy string-vector containing the names of the economies which are part of the economic system
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'@keywords internal

Get_Bs <- function(mat, dt, K1XQ, SSZ, Wpca, FactorLabels, Economy, ModelType){

  Lab_SingleQ <- c("JPS original","JPS global", "GVAR single")

  N <- length(FactorLabels$Spanned) # number of country-specific spanned factors
  M <- length(FactorLabels$Domestic) - N  # number of country-specific unspanned factors

  # 1) Get BnX and BnP
    # Yields are affine function of the vector Pt,i.e., Y(t) = AnP +BnP*P(t).
    # Further, we define Z(t) as an affine function of P(t) such that: Z(t) = phi0+ phi1*P(t).
    # As such, we can write P(t) = phi1^(-1)*(Z(t) - phi0).
    BnX <- A0N__BnAn(round(mat/dt),K1XQ, ModelType, Economies = Economy)$BnX/dt
    # NOTE: the function "A0N__computeBnAn" generates outputs for interest rates per unit of time interval.
    # Hence, multiplying by "1/dt" produces annualized results.
    if (any(ModelType == Lab_SingleQ)){
    if (N==1){ BnP <- mrdivide(matrix(BnX),Wpca%*%BnX)}
    else{BnP <- mrdivide(BnX,Wpca%*%BnX) } #  BnP = BnX*(W*BnX)^(-1)

    } else{
    BnP <- mrdivide(BnX,Wpca%*%BnX)
      }

    # 2) Covariance matrix of latent states X:
  if (any(ModelType == Lab_SingleQ)){
  AllLabels <- GetLabels_sepQ(Economy, ModelType, FactorLabels)
    dimnames(SSZ) <- list(AllLabels, AllLabels)
    LabelSpannedCS <- c(FactorLabels$Tables[[Economy]][-(1:M)])
    idxSpa <- match(LabelSpannedCS, AllLabels)
  } else{
    G <- length(FactorLabels$Global)
    C <- length(Economy)
    idxSpa <- IdxSpanned(G,M,N,C)
  }
    WBX <- Wpca%*%BnX # WBX = W * BnX
    SSP <- SSZ[idxSpa,idxSpa]

    SSX <- matrix(Inf, nrow = nrow(SSP), ncol = ncol(SSP)) # Initialize with Inf

    SSX <- mrdivide(solve(WBX, SSP, tol = 1e-50), t(WBX)) # Try to assign the correct value


    # 3) Optimal estimate of r0. NOTE: We set r0=0, because when r0=0, An = Betan.
    if (any(ModelType == Lab_SingleQ)){r0 <- 0}else{r0 <- rep(0, times=C)}

    betan <- A0N__BnAn(round(mat/dt), K1XQ, ModelType, r0 = r0, SSX=SSX, Economies =Economy)$betan
    betan <- t(t(betan/dt))


  out <- list(BnP = BnP, BnX = BnX, betan = betan, SSP = SSP, SSX = SSX, idxSpa = idxSpa)

  return(out)
}
############################################################################################################
#' Compute the A loadings
#'
#'@param LoadBs list containing the B loadings
#'@param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#'@param r0 long-run interest rate  (scalar or vector with length C)
#'@param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'@keywords internal

Get_As <-  function(LoadBs, Wpca, r0, dt, Economies, ModelType){

  betan <- LoadBs$betan
  BnP <- LoadBs$BnP

  if (any(ModelType == c("JPS original","JPS global", "GVAR single"))){
  AnX <- betan + as.numeric(r0/dt)
  } else{

    C <- length(Economies)
    J <- nrow(betan)/C

    for (i in 1:C){
      Qtemp <- t(t(rep(r0[i], times= J)/dt))
      if (i==1){
        Q <- Qtemp
      }else{
        Q <- rbind(Q, Qtemp)
      }
    }
    AnX <- betan + Q
  }

  AnP <- AnX - LoadBs$BnP%*%(Wpca%*%AnX) # AnP = AnX - BnX*(W*BnX)^(-1)*W*AnX

  return(list(AnX = AnX, AnP = AnP))
}
############################################################################################################
#'Generate the factor labels for models estimated on a country-by-country bases
#'
#'@param Economy  string containing the names of the economy to be estimated
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param FactorLabels list containing the factor labels
#'
#'@keywords internal

GetLabels_sepQ <- function(Economy, ModelType, FactorLabels){
  if (ModelType == "JPS original"){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]]) }
  else if (any(ModelType == c("JPS global", 'GVAR single'))){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

  return(AllLabels)
}


#########################################################################################################
#' Compute the log-likelihood function
#'
#'@param P time-series of spanned factors (N x T or CN x T)
#'@param Y time-series of yields (J x T or CJ x T)
#'@param Z time-series of risk factors (F x T)
#'@param N number of country-specific spanned factors
#'@param mat vector of maturities (in years) of yields used in estimation (J x 1)
#'@param We matrix of weights of the portfolios observed with errors ((J-N) x J or C(J-N) x CJ)
#'@param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x CJ)
#'@param K0Z matrix of intercepts (P-dynamics)
#'@param K1Z feedback matrix (P-dynamics)
#'@param SSZ variance-covariance matrix (P-dynamics)
#'@param LoadBs list containing the B loadings
#'@param LoadAs list containing the A loadings
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'@keywords internal

Get_llk <- function(P, Y, Z, N, mat, We, Wpca, K0Z, K1Z, SSZ, LoadBs, LoadAs, ModelType){

  T <- dim(P)[2]
  t <- 2:T
  J <- length(mat)

  Peo <- We%*%Y # portfolio observed WITH errors
  BnP <- LoadBs$BnP
  AnP <- LoadAs$AnP


  # 1) density of yields pricing errors:
  Pe <- matrix(0, nrow= nrow(We) , ncol= T)
  for( f in 1:T){
    Pe[,f] <- We%*% (BnP%*%P[,f] + AnP) # Pe_t = We*(AP+BP*P_t): "model-implied" yields for the portfolio of yields observed with errors
  }
  eQ <- Peo[,t] - Pe[,t]

  # 2) Standard deviation of the measurement error.
  if (!exists("se")|| sjmisc::is_empty(se) ){
    if (any(ModelType == c("JPS original","JPS global", "GVAR single"))){
      se <- sqrt(mean(as.vector(eQ^2))) # Scalar
    } else{
      C <- nrow(P)/N
      se <- rep(NA, times= C)
      idx0 <- 0
      for (i in 1:C){
        idx1 <- idx0+ J-N
        se[i] <- sqrt(mean(as.vector(eQ[(idx0+1):idx1,]^2))) # Scalar - standard deviation of the measurement error.
        idx0 <- idx1
      }
}
}

  # 3) The log-likelihood function:
  if (any(is.nan(se) ==1)){
    y <- 1e6*matrix(1,T-1,1)
  }else{

    # Cross-sectional density (i.e. density for the portfolios observed with measurement error)
    if (any(ModelType == c("JPS original","JPS global", "GVAR single"))){
      y <- GaussianDensity(eQ, se^2*diag(J-N))
    } else{
            aa <- se
      idx0 <- 0
      for (h in 1:C){
        idx1 <- idx0+ J-N
        se[(idx0+1):idx1] <- rep(aa[h], times=J-N)
        idx0 <- idx1
      }

      se <- se[seq(1, C*J-C*N, by=J-N)] # Recast se
      y <- GaussianDensity(eQ, se^2*diag(C*J-C*N))

    }

    # time series density:
    eP <- matrix(0, nrow(K0Z), T-1)
    for (f in 1:length(t)){ eP[,f] <-  Z[,f+1] - K1Z%*%Z[,f] - K0Z }

    y <- y + GaussianDensity(eP, SSZ) # Cross-sectional density + time-series density (final likelihood function, except for the Jacobian term)

    # Jacobian: (so that the density is for the observed yields (and the non-yields variables Z)
    #           and not the yields portfolios which depend on W).

    y <-  y + 0.5*log(abs(det(Wpca%*%t(Wpca))))
  }

  if ( any(is.nan(y)) ||any(is.infinite(y)) || any(Im(y)!= 0) ){
    y <- 1e6*rep(1, times = T-1)
  }else{
    y <- -t(as.vector(y)) # necessary for minimization
  }

  return(list(y = y, se = se))
}
########################################################################################################
#' computes the density function of a gaussian process

#'@param res    matrix of residuals (N x T)
#'@param SS     covariance matrice or array of covariance matrices \cr
#'               (If dim(SS) > 3, then the model has a stochastic volatility) (N x N) or (N x N x T)
#'
#'@param invSS  Inverse of SS (N x N) or (N x N x T) - optional input
#'@param logabsdetSS   log(abs(|SS|)) (1 x T) - optional input
#'
#'
#'@keywords internal
#'
#'@return y   - vector of density (1 x T)
#'
#'@references
#' This function is based on the "Gaussian" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'



GaussianDensity <- function(res,SS, invSS, logabsdetSS){


  nargin <- nargs()


  N <- dim(res)[1]
  T <- dim(res)[2]


  if (is.na(dim(SS)[3])){ # If SS is a matrix and NOT a 3-dimension array.
    if (!missing("invSS")){
      SSres <-  invSS%*%res
    }else{
      SSres <- solve(SS,res, tol = 1e-50)} # SSres <- (SS)^(-1)*res


    if (missing("logabsdetSS")){
      logabsdetSS <- 0.5*logdet(SS%*%t(SS),opt="chol")} # why do we multiply by 0.5? Because in "logdet", we do
    # 2 *sum(log(abs)).

  }else{
    if(nargin<3){
      invSS <- mult__inv(SS, whichoutput= NULL , 2)$inva
      logabsdetSS <- mult__inv(SS, whichoutput= NULL , 2)$logabsdet
    }else if (nargin==3){
      logabsdetSS <- mult__inv(SS, whichoutput= NULL,'logabsdet')$logabsdet
    }

    a <- invSS[[1]]
    b <- array(res, c(N, 1, T))
    SSres <- array(mult__prod(a, b, c=NULL), c(N, T))
    logabsdetSS <- logabsdetSS[[1]]
  }
  y <- -0.5*N*log(2*pi)-0.5*logabsdetSS-0.5*abs(colSums(res*SSres)) # Last term is equal to "t(res)*(SS)^(-1)*res"

  return(y)

}

#############################################################################################################
#' Compute the variance-covariance matrix of the bond yields
#'
#'@param YieldsTS matrix of yields used in estimation  (J x T or CJ x T)
#'@param N number of country-specific spanned factors
#'@param mat vector of maturities (in years) of yields used in estimation (J x 1)
#'@param WpcaFull composite matrix of weights the portfolios observed with and without errors
#'@param se Variance of the portfolio of yields observed with error (scalar). Default is set to NULL
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'@keywords internal

Get_SigmaYields <- function(YieldsTS, N, mat, WpcaFull, se, ModelType){
  J <- length(mat)

  # Single-country models
  if (any(ModelType == c("JPS original","JPS global", "GVAR single"))){
  Sigma_e_Row  <- c(rep(0, times = N), rep(se,times= J-N)^2) # Variances of portfolio WITHOUT and WITH errors
  Sigma_e <- Sigma_e_Row*diag(J)
  SIGMA_yields <- solve(WpcaFull)%*%Sigma_e%*%t(solve(WpcaFull))
  VarYields <- diag(SIGMA_yields)
  } else{

  # Multi-country models
    Sigma_e_CS <- c()
    Sigma_e_AllRows <- c()

    C <- length(se)

    for (i in 1:C){
      Sigma_e_CS  <- c(rep(0, times = N), rep(se[i],times= J-N)^2) # Variances of portfolio WITHOUT and WITH errors
      Sigma_e_AllRows <- append(Sigma_e_AllRows, Sigma_e_CS)
    }

    Sigma_e_AllCountries <- Sigma_e_AllRows*diag(C*J)
    SIGMA_yields <- solve(WpcaFull)%*%Sigma_e_AllCountries%*%t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)

    names(VarYields) <- rownames(YieldsTS)
    VarYields <- t(t(VarYields))

}
  return(VarYields)
}

################################################################################################################
#' Update the variance-covariance matrix from the "JLL joint Sigma" model. Necessary for optimization
#'
#'@param SSZ Variance-covariance matrix from JLL model
#'@param Z complete set of spanned and unspanned factors (F x T)
#'@param N number of country-specific spanned factors
#'@param JLLinputs List of inputs from JLL models
#'
#'@keywords internal

Update_SSZ_JLL <- function(SSZ, Z, N, JLLinputs){
  Sigma_Ye <- t(chol(SSZ)) # in this case, SSZ is the Variance-covariance matrix of the ORTHOGONALIZED dynamics
  JLLinputs$WishSigmas <- 0 # Avoid recomputing the variance-covariance matrices
  JLLPara <- JLL(Z, N, JLLinputs)
  Sigma_Y <- JLLPara$PI%*%Sigma_Ye
  SSZ <- Sigma_Y%*%t(Sigma_Y) # Variance-covariance matrix non-orthogonalized

  return(SSZ)
}

#####################################################################################################################
#' Prepare outputs to export after the model optimization
#'
#'@param Y matrix of yields used in estimation  (J x T or CJ x T)
#'@param Z complete set of spanned and unspanned factors (F x T)
#'@param mat vector of maturities (in years) of yields used in estimation (J x 1)
#'@param N number of country-specific spanned factors
#'@param dt time interval unit of the model (scalar)
#'@param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#'@param K1XQ risk-neutral feedback matrix (N x N or CN x CN)
#'@param SSZ variance-covariance matrix (F x F)
#'@param LoadAs list containing the A loadings
#'@param LoadBs list containing the B loadings
#'@param r0 long-run interest rate  (scalar or vector with length C)
#'@param se Variance of the portfolio of yields observed with error (scalar).
#'@param K0Z intercept from the P-dynamics (F x 1)
#'@param K1Z feedback matrix from the P-dynamics (F x F)
#'@param Gy.0 matrix of contemporaneous terms from the P-dynamics (F x F)
#'@param VarYields variance-covariance matrix of the bond yields
#'@param y likelihood of each time series (Tx1)
#'@param GVARinputs List of inputs from GVAR models
#'@param JLLinputs List of inputs from JLL models
#'@param Economies string containing the names of the economy to be estimated
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param BS_out
#'
#'@keywords internal

OptOutputs <- function(Y, Z , mat, N, dt, Wpca, K1XQ, SSZ, LoadAs, LoadBs, r0, se, K0Z, K1Z, Gy.0, VarYields,
                       y, GVARinputs, JLLinputs, Economies, ModelType, BS_out = FALSE){

    # a) Build LIST 1
    # a.1) List "input"
    inputs <- list(Y = Y, AllFactors = Z , mat = mat, N= N, dt= dt, Wpca = Wpca)
    if (any(ModelType == c("GVAR single", "GVAR multi"))) {  inputs$Wgvar <-  GVARinputs$Wgvar}
    else if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
      inputs$DomUnit <- JLLinputs$DomUnit }

    # a.2) List "ests"
    ests <- list(K1XQ = K1XQ, SSZ = SSZ, SSP = LoadBs$SSP, r0 = r0, se = se, K0Z = K0Z, K1Z = K1Z,
                 Gy.0 = Gy.0,  VarYields = VarYields)

    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
      JLLinputs$WishSigmas <- 1
      JLLinputs$SigmaNonOrtho <- SSZ
      JLLoutcomesOrtho <- JLL(Z, N, JLLinputs)
      # Remove the non-orthogonalized outputs
      JLLoutcomesOrtho$k0 <- NULL
      JLLoutcomesOrtho$k1 <- NULL
      if (BS_out){
        JLLoutcomesOrtho$a_W <- NULL
        JLLoutcomesOrtho$b <- NULL
        JLLoutcomesOrtho$c <- NULL
        JLLoutcomesOrtho$a_DU_CS <- NULL
        JLLoutcomesOrtho$PIb <- NULL
        JLLoutcomesOrtho$PIac <- NULL
        JLLoutcomesOrtho$Ye <- NULL
      }

      #Add the zeros to the variance-covariance matrix (to avoid rounding problems)
      ZeroIdxVarCov <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$VarCovOrtho
      ZeroIdxSigma_Ye <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$Sigma_Ye
      JLLoutcomesOrtho$Sigmas$VarCov_Ortho[ZeroIdxVarCov] <- 0
      JLLoutcomesOrtho$Sigmas$Sigma_Ye[ZeroIdxSigma_Ye] <- 0

      PI <- JLLoutcomesOrtho$PI
      JLLoutcomesOrtho$Sigmas$Sigma_Y <- PI%*%JLLoutcomesOrtho$Sigmas$Sigma_Ye

      ests$JLLoutcomes <- JLLoutcomesOrtho
    }

    # a.3) List "llk"
    llk <- list(-t(y))

    # a.4) List Q dynamics: X as risk factors:
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      Q <- list(K0 = matrix(0,N,1), K1 = K1XQ, SS = LoadBs$SSX)
    }else { Q <- list(K0 = matrix(0, N*length(Economies),1), K1 = K1XQ, SS = LoadBs$SSX) }
    X <- list(B = LoadBs$BnX, A = LoadAs$AnX, Q = Q)
    rot <- list(X = X)

    # Summary: List 1
    Out <- list(inputs = inputs, ests = ests, llk= llk, rot= rot)

    # b) Build LIST 2
    # Q dynamics: PCN as risk factors: % PCN = U0 + U1*X, where U0 = W*AnX and U1 = W*BnX
    U1 <- Wpca%*%LoadBs$BnX
    U0 <- Wpca%*%LoadAs$AnX
    out.rot.P <-  FMN__Rotate(Out$rot$X, U1, U0)

    # updating the generated outputs
    Out$rot$P <- out.rot.P
    if (BS_out){
    Out$rot$X <- NULL # Delete the parameters related to the latent factors
    #(Not needded for the Bootstrap, save some time and pc memory)
    } else{
    # P-dynamics: Z as risk factors:
    out.rot.P.P <- list( K0= K0Z, K1 = K1Z, SS = SSZ)
    #names(out.rot.P.P) <- c("K0", "K1", "SS" )
    Out$rot$P$P <- out.rot.P.P

    # P-dynamics: PCN as risk factors:
    idxSpanned <- LoadBs$idxSpa
    K0P <- K0Z[idxSpanned]
    K1P <- K1Z[idxSpanned, idxSpanned]
    SSP <- SSZ[idxSpanned, idxSpanned]
    out.rot.P.P2 <- list(K0 = K0P, K1 = K1P, SS = SSP)

    out.rot.X.P <- FMN__Rotate(out.rot.P.P2, solve(U1, tol = 1e-50), -solve(U1,U0, tol = 1e-50))

    Out$rot$X$P <- out.rot.X.P

    # PC as risk factors:
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      PC1NW <- pca_weights_one_country(Y, Economies)
      PC1NW <- PC1NW[1:N,]*100
      U1 <- PC1NW%*%LoadBs$BnX
      U0 <- PC1NW%*%LoadAs$AnX
    }else{
      U1 <- Wpca%*%LoadBs$BnX
      U0 <- Wpca%*%LoadAs$AnX
    }

    out.rot.PC <- FMN__Rotate(Out$rot$X, U1, U0)
    Out$rot$PC <- out.rot.PC

    # P-dynamics: PC as risk factors:
    out.rot.PC.P <- list( K0 = K0Z, K1 = K1Z, SS = SSZ)
    Out$rot$PC$P <- out.rot.PC.P
    }
return(Out)
}

#################################################################################################################
#' computes the logarithm of determinant of a matrix A

#'@param A squared matrix
#'@param opt  "chol" or "NULL" (text)
#'
#'
#'@importFrom pracma lu


#'@keywords internal
#'@details
#' Theoretically, this function should be functionally equivalent to log(det(A)). However, it avoids the
#' overflow/underflow problems that are likely to happen when applying det to large matrices.

#' The key idea is based on the mathematical fact that the determinant of a triangular matrix equals the
#' product of its diagonal elements. Hence, the matrix's log-determinant is equal to the sum of their logarithm
#' values. By keeping all computations in log-scale, the problem of underflow/overflow caused by product of
#' many numbers can be effectively circumvented. The implementation is based on LU factorization.

#'@references
#' This function is based on the "logdet" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'



logdet <-function(A, opt){


  if (opt == "chol"){
    v <- 2 * sum(log(abs(diag(chol(A)))))
  }else{
    # lu function: Triangular (LU) Decomposition Of A Matrix
    L <- lu(A)$L # generates a lower triangular matrix
    U <- lu(A)$U # generates an upper triangular matrix
    P <- A%*% solve(L%*%U) # As A = P*L*U, P= A(LU)^(-1).
    du <- t(t(diag(U)))
    c <- det(P) * prod(sign(du))
    v <- log(c) + sum(log(abs(du)))
  }

  return(v)
}

######################################################################################################
###############################################################################################################
#' Inverts an array of matrices so that:  inva[,,i] = inv(a[,,i])

#'@param a        matrix array (N x N x T)
#'@param whichoutput      if = 'lobabsdet' computes the log(abs(det(a))) only (text).
#'@param nargout    "nargout == 1" or "nargout == 2"(scalar)


#'@importFrom pracma numel size
#'
#'@keywords internal
#'@return "nargout == 1" returns inva -> matrix array: a^(-1) (N x N x T) \cr
#'        "nargout == 2" returns inva -> matrix array: a^(-1) (N x N x T) and  logabsdet ->  vector of log(abs(det(a)))  (1 x T)

#'@references
#' This function is modified version of the "mult__inv" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."\cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029



mult__inv <- function(a, whichoutput, nargout){


  if (!exists("whichoutput")){
    whichoutput = ""
  }

  ss <- size(a)

  if (numel(ss)<=2){
    if (nargout>1){
      inva <- solve(a)
      logabsdet <- log(abs(det(a)))
    } else {
      if (whichoutput =="logabsdet"){
        inva <- log(abs(det(a)))
      } else{
        inva <- solve(a)
      }
    }
  }

  if (numel(ss)>2){

    if (dim(a)[1]<5){
      inva <- mult_inv_small(a)$y
      dd <-  mult_inv_small(a)$dd

      if (nargout>1) {
        logabsdet <- log(abs(dd))
      } else{
        if (whichoutput == 'logabsdet'){
          inva <- log(abs(dd))
        }
      }
    } else{
      if (nargout>1){
        inva <- mult_inv_large(a)
        logabsdet <- mult_logabsdet(a)
      } else{
        if (whichoutput == 'logabsdet') {
          inva <- mult_logabsdet(a)
        } else{
          inva <- mult_inv_large(a)
        }
      }
    }
  }



  # Printing the outputs:

  if (nargout==1){
    return(inva)
  }else{
    output <- list(inva,logabsdet)
    names(output) <- c("inva", "logabsdet")
    return(output)
  }


}

#############################################################################################################
#############################################################################################################
#' Efficient computation of matrix product for arrays

#'@param a array
#'@param b array
#'@param c array

#'@keywords internal
#'
#'@details
#'#' Efficiently computes matrix product for arrays a and b and c:
#'  d[,,i] = a[,,i] b[,,i] c[,,i].
#' (efficiently = without using an inefficient loop)
#'@references
#' This function is modified version of the "mult__prod" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."\cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029



mult__prod <- function(a,b,c){


  if (!exists("c", inherits = FALSE) ||is.array(c) || is.null(c) ) { # "inherits = FALSE" ensures that we serach for variables defined with that name only.
    d <-  multiprod_2terms(a,b)
  } else {
    sa <- dim(a)
    sb <- dim(b)
    sc <- dim(c)
    if (numel(a)==1){
      d <- a*multiprod_2terms(b,c)
    } else if (numel(b)==1) {
      d <- b*multiprod_2terms(a,c)
    } else if (numel(c)==1) {
      d <- c*multiprod_2terms(a,b)
    } else if (length(sa)==2 && length(sb)==2 && length(sc)==2){
      d <- a*b*c
    } else {

      if (sa[2]!= sb[1] || sb[2]!= sc[1]) {
        stop("arrays do not commute")
      } else {
        k <- array(a, c(sa[1:2], 1, 1, sa[seq_len(max(0, length(sa) - 2)) +2] )) #
        l <- array(b, c(1, sb[1:2], 1, sb[seq_len(max(0, length(sb) - 2)) +2]) )
        d <- array(0, c(dim(k),T))
        for (j in 1:T){ d[ , , , ,j] <- k*l[ , , , , j] }
        d <- colSums(d)
        m <- array(c, c(1, 1, sc))
        o <- array(0, c(dim(d)[1], 1, 1, dim(m)[length(dim(m))-1], T) )
        for (j in 1:T){
          for (p in 1:dim(d)[1]){
            for (q in 1:dim(m)[length(dim(m))-1]){
              o[ p, , , q,j] <- d[p,,,j]*m[,,,q,j] }}}
        d <- o
        sd <- dim(d)
        d <- array(d, c(sd[1], sd[seq_len(max(0, length(sd) - 3)) +3 ] )) #
      }
    }
  }

  return(d)
}

#############################################################################################################
#' Inverse each 2D slice of an array (M) with arbitrary dimensions support

#'@param M multi-dimension array
#'
#'
#'@keywords internal


#'@return n_D array (m x m x [p x q x  ...]), with same size as M
#'@details
#'  Inverse every 2D slice (the first two dimensions of M) for multi-dimension array M.
#'   M[,,p,q,...] * X[,,p,q,...] <- repmat(diag(m),[1,1,p,q,...])

#'@references
#' This function is based on the "mult_inv_large" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#' Available at: https://cepr.org/40029



mult_inv_large <- function(M){


  sn <- dim(M)
  m <- sn[1]
  n <- sn[2]

  if (m!=n){
    stop("mult_inv_large: The first two dimensions of M must be m x m slices.")
  }

  p <- prod(seq_len(max(0, length(sn) - 2)) +2)


  RHS <- kronecker(array(1, c(p,1)), diag(m))
  X <- solve(M, RHS)


  X <- array(X, c(n, p, m))
  X <- aperm(X, c(1,3,2))
  X <- array(X, c(n,m, seq_len(max(0, length(sn) - 2)) +2 ))


  return(X)
}

#############################################################################################################
#' Inverse the (m,m,T) array of matrices for m<=4

#'@param A multi-dimensional array
#'
#'

#'@keywords internal


#'@return dd  vector of determinants (Tx1)
#'
#'@details
#' equivalent to the following for loop:
#'
#' \code{ for (i in 1:T){ y[ , ,i] = inv(A[,,i]) } }
#'
#'@references
#' This function is based on the "mult_inv_small" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029



mult_inv_small <- function(A){


  y <- as.numeric()
  dd <- as.numeric()
  m <- dim(A)[1]
  mm <- dim(A)


  if (m==1){
    dd <- t(t(as.vector(A)))
    y <-  1/ A
  }


  if (m==2){
    dd <- A[1,1,]*A[2,2,] - A[1,2,]*A[2,1,]
    y <- -A

    y <-  y/dd

  }

  if (m==3){
    a <- A[1,1,]
    b <- A[1,2,]
    c <- A[1,3,]
    d <- A[2,1,]
    e <- A[2,2,]
    f <- A[2,3,]
    g <- A[3,1,]
    h <- A[3,2,]
    k <- A[3,3,]
    ekfh <- e*k-f*h
    fgdk <- f*g-d*k
    dheg <- d*h-e*g
    dd <- b*fgdk + a*ekfh + c*dheg

    y <- A
    y[1,1,] <- ekfh
    y[1,2,] <- c*h-b*k
    y[1,3,] <- b*f-c*e
    y[2,1,] <- fgdk
    y[2,2,] <- a*k-c*g
    y[2,3,] <- c*d-a*f
    y[3,1,] <- dheg
    y[3,2,] <- b*g-a*h
    y[3,3,] <- a*e-b*d


    y <-  y/dd

  }

  outputs <- list(y = y, dd = dd)
  return(outputs)
}

#############################################################################################################
#' Inverse each 2D slice of an array (M) with arbitrary dimensions support

#'@param M multi-dimensional array
#'
#'@importFrom pracma lu

#'@keywords internal


#'@return X  : (n-2)_D array ([p x q x  ...])
#'
#'@details
#'  Inverse every 2D slice (the first two dimensions of M) for multi-dimension array M. In Matlab language:
#'   M[ , , p,q,...] * X[ , ,p,q,...] = repmat(diag(m),[1,1,p,q,...])
#'
#'
#'@references
#' This function is based on the "mult_logabsdet" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029



mult_logabsdet <- function(M){

  sn <- dim(M)
  m <- sn[1]
  n <- sn[2]

  if (m!=n){
    stop("mult_logabsdet: The first two dimensions of M must be m x m slices.")
  }

  p <- prod(seq_len(max(0, length(sn) - 2)) +2) #


  u <-lu(M)$U
  X <- colSums(array(log(abs(diag(u))), c(m, p)))

  return(X)
}


#############################################################################################################
#' computes matrix product for arrays a and b:   c[,,i] = a[,,i] b[,,i]

#'@param a array (M x N x T)
#'@param b array (N x K x T)
#'@return c  array (M x K x T)
#'
#'@keywords internal
#'@references
#' This function is based on the "multiprod_2terms" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029



multiprod_2terms <- function(a,b){



  if (numel(a)==1||numel(b)==1||(length(dim(a))==2 && length(dim(b))==2)){ # either a or b is scalar or matrices
    c <- a%*%b
  }else{
    sa <- dim(a)
    sb <- dim(b)

    if ((sa[1]==1 && sa[2]==1)||(sb[1]==1 && sb[2]==1)) { # either a or b is scalar
      c <- array(0, c(nrow(a), ncol(a) ,T))

      for (f in 1:T) {  c[,,f] <- a*b[,,f]   }
    } else{


      if (sa[2]!=sb[1]){
        stop('arrays do not commute')
      }else{
        d <- array(a, c(sa[1:2], 1, seq_len(max(0, length(sa) - 2)) +2))
        g <- array(b, c(1, sb))
        h <- array(0, c(dim(d),T))
        for (f in 1:T){    h[,,,f] <- d*g[,,,f]     } # Note: this loop does not accommodate arrays of other dimensions
        c <- colSums(h)
        sc <- dim(c)
        c <- array(c, c(sc[1:length(sc)]))
      }

    }
  }

  return(c)

}

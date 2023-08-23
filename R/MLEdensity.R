#' Compute the maximum likelihood function ("sep Q" models)
#'
#'@param K1XQ risk-neutral feedback matrix (NxN)
#'@param r0   long-run interest rate  (scalar)
#'@param SSZ  variance-covariance matrix (KxK)
#'@param K0Z  intercept from the P-dynamics (Kx1)
#'@param K1Z  feedback matrix from the P-dynamics (KxK)
#'@param se   Variance of the portfolio of yields observed with error (scalar)
#'@param Gy.0 matrix of contemporaneous terms from the P-dynamics (KxK)
#'@param mat  vector of maturities (in years) of yields used in estimation (J x 1)
#'@param Y    matrix of yields used in estimation  (J x T)
#'@param Z    complete set of spanned and unspanned factors (KxT)
#'@param P    complete set of spanned factors (NxT)
#'@param Wpca matrix of weights of the portfolios observed without errors (NxJ)
#'@param We   matrix of weights of the portfolios observed with errors ((J-N)xJ)
#'@param WpcaFull composite matrix of weights the portfolios observed with and without errors
#'@param dt  time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#'@param Economy name of the economy under study
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param ModelType Feasible options are: (i) "JPS", (ii) "JPS jointP" or (iii) "GVAR sepQ"
#'@param GVARinputs if the model chosen is the "GVAR sepQ", the "GVARinputs" should be specified (see "GVAR" function)
#'@param nargout if nargout== 1: provides only the values of the likelihood; if nargout== 2: complete ATSM outputs
#'
#'@importFrom pracma mldivide mrdivide
#'
#'@references
#' This function is modified version of the "A0N_MLEdensity_WOE" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal




MLEdensity_sepQ <- function(K1XQ, r0, SSZ, K0Z, K1Z, se, Gy.0, mat, Y, Z, P, Wpca, We, WpcaFull,
                            dt, Economy, FactorLabels, ModelType, GVARinputs = NULL, nargout){

  N <- length(FactorLabels$Spanned) # number of country-specific spanned factors
  M <- length(FactorLabels$Domestic) - N  # number of country-specific unspanned factors
  J <- numel(mat) # number of country-specific yields used in estimation;
  T <- dim(Y)[2] # time-dimension of the model


  if (ModelType == "JPS"){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]]) }
  if (ModelType == "JPS jointP" || ModelType == 'GVAR sepQ'){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

  rownames(SSZ) <- AllLabels
  colnames(SSZ) <- AllLabels
  LabelSpannedCS <- c(FactorLabels$Tables[[Economy]][-(1:M)])
  idxSpanned <- match(LabelSpannedCS, AllLabels)

  Peo <- We%*%Y # portfolio observed WITH errors

  t <- 2:T

  # Intialize some variables, if necessary
  if (!exists("r0")){ r0 <- as.numeric()}
  if (!exists("se")){ se <- as.numeric()}
  if (!exists("K0Z")){ K0Z <- matrix(, nrow=N, ncol=0)}
  if (!exists("K1Z")){ K1Z <- matrix(, nrow=N*N, ncol=0)}

  if (sjmisc::is_empty(K0Z) || sjmisc::is_empty(K1Z) ){
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, K0Z, K1Z )
  } else{
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, t(t(as.vector(K0Z))), t(t(as.vector(K1Z))) )
  }
  # x is a column vector that contains all the elements of each parameter of the model.

  if ( any(is.nan(x)) ||any(is.infinite(x)) || any(Im(x)!= 0) ){
    y <- 1e6*rep(1, times = T-1)
    out <- as.numeric()
  }else{

    # Loadings:
    # Yields can be affine function of Pt's,i.e., Y(t) = AnP +BnP*P(t).
    # Further, we define Z(t) as an affine function of P(t) such that: Z(t) = phi0+ phi1*P(t).
    # As such, we can write P(t) = phi1^(-1)*(Z(t) - phi0).
    BnX <- A0N__computeBnAn_sepQ(round(mat/dt),K1XQ, dX= NULL, r0= NULL, SSX= NULL)[[1]]/dt
    # NOTE: the function "A0N__computeBnAn" generates outputs for interest rates per unit of time interval.
    # Hence, by multiplying by "1/dt" we obtain annualized results.
    BnP <- mrdivide(BnX,Wpca%*%BnX) #  BnP = BnX*(W*BnX)^(-1)


    # Covariance matrix of latent states X:
    WBX <- Wpca%*%BnX # WBX = W * BnX
    SSP <- SSZ[idxSpanned,idxSpanned]
    SSX <- mrdivide(solve(WBX,SSP, tol = 1e-50),t(WBX)) # SSX = (W*BnX)^(-1)*SSP*t((*W*BnX)^(-1))
    # Optimal estimate of r0:
    betan <- A0N__computeBnAn_sepQ(round(mat/dt), K1XQ, dX=NULL , r0=0, SSX)[[3]]  # NOTE: We set r0=0, because when r0=0, An = Betan.
    betan <- t(t(betan/dt))


    if (!exists("r0")|| sjmisc::is_empty(r0) ){
      # (i) A0 = (I - Bx(W*Bx)^(-1)*W)*Ax0;
      # (ii) A1 = (I - Bx(W*Bx)^(-1)*W)*Axr;
      # (iii) APer= We * A1
      # (iv) vt = We*(Yt - Bx(W*Bx)^(-1)*Pt - A0)
      # Recall that BnP= Bx(W*Bx)^(-1)
      A0 <- (diag(J)-BnP%*%Wpca)%*%betan
      A1 <- (diag(J)-BnP%*%Wpca)%*%matrix(1,J,1)/dt

      #  r0 = APer'*vt/(APer'*APer)
      r0 <- solve((t(A1)%*%(t(We)%*%We)%*%A1), # Numerator from equation r0 equation
                  (t(A1)%*%(t(We)%*%We))%*%(rowMeans(Y[,t] - BnP%*%P[,t]) - A0)) # Denominator from r0 equation
    }


    # compute the intercepts:
    AnX <- betan + as.numeric(r0/dt)
    AnP <- AnX -BnP%*%(Wpca%*%AnX) # AnP = AnX - BnX*(W*BnX)^(-1)*W*AnX

    # density of yields pricing errors:
    Pe <- matrix(0, nrow= nrow(We) , ncol= T)
    for( f in 1:T){
      Pe[,f] <- We%*% (BnP%*%P[,f] + AnP) # Pe_t = We*(AP+BP*P_t): "model-implied" yields for the portfolio of yields observed with errors
    }
    eQ <- Peo[,t] - Pe[,t]

    if (!exists("se")|| sjmisc::is_empty(se) ){
      se <- sqrt(mean(as.vector(eQ^2))) # Scalar - standard deviation of the measurament error.
    }

    ## The log-likelihood function:
    if (is.nan(se)){
      y <- 1e6*matrix(1,T-1,1)
    }else{
      # Cross-sectional density (i.e. density for the portfolios observed with measurament error)
      y <- GaussianDensity(eQ, se^2*diag(J-N))

      # time series density:
      if ((!exists("K0Z")||sjmisc::is_empty(K0Z)) & (ModelType == 'JPS' || ModelType == 'JPS jointP' )){
        VARpara <- VAR(Z, VARtype= "unconstrained", Bcon = NULL)
        K0Z <- VARpara$K0Z  # Column vector (Kx1)
        K1Z <-  VARpara$K1Z # matrix (KxK)
        Gy.0 <- diag(length(K0Z))
      }
      if ( (!exists("K0Z")||sjmisc::is_empty(K0Z)) & (ModelType == 'GVAR sepQ') )  {
       # Estimate GVAR(1)
        GVARPara <- GVAR(GVARinputs, N)
        K0Z <- GVARPara$F0 # Column vector (Kx1)
        K1Z <- GVARPara$F1 # matrix (KxK)
        Gy.0 <- GVARPara$Gy.0
      }

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

    # Variance-yields
    Sigma_e_Row  <- c(rep(0, times = N), rep(se,times= J-N)^2) # Variances of portfolio WITHOUT and WITH errors

    Sigma_e <- Sigma_e_Row*diag(J)
    SIGMA_yields <- solve(WpcaFull)%*%Sigma_e%*%t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)



    if (nargout>1){

      if (ModelType == 'JPS' || ModelType == 'JPS jointP' ){ inputs <- list(Y, Z, mat, N, dt, Wpca)}
      if (ModelType == 'GVAR sepQ' ){ inputs <- list(Y, Z, mat, N, dt, Wpca, GVARinputs$Wgvar)}

      ests <- list(K1XQ, SSZ, SSP, r0, se, K0Z, K1Z, Gy.0, VarYields)
      llk <- list(-t(y)) # [1 x T-p] actual densities (not negative of densities)


      # Q dynamics: X as risk factors:
      Q <- list(zeros(N,1), K1XQ, SSX)
      X <- list(BnX, AnX, Q)
      rot <- list(X)

      Out <- list(inputs, ests, llk, rot)


      names(Out) <- c("inputs", "ests", "llk", "rot")

      if (ModelType == 'JPS' || ModelType == 'JPS jointP' ){ names(Out$inputs)  <-c('Y','AllFactors' , 'mat', 'N', 'dt', 'Wpca')}
      if (ModelType == 'GVAR sepQ' ){ names(Out$inputs)  <-c('Y','AllFactors' , 'mat', 'N', 'dt', 'Wpca', 'Wgvar')}

      names(Out$ests) <- c('K1XQ', 'SSZ', 'SSP', 'r0', 'se', 'K0Z', 'K1Z', 'Gy.0', 'VarYields')
      names(Out$rot) <- c("X")

      names(Out$rot$X) <-  c("B","A","Q")
      names(Out$rot$X$Q) <- c("K0", "K1", "SS" )


      # Q dynamics: PCN as risk factors: % PCN = U0 + U1*X, where U0 = W*AnX and U1 = W*BnX
      U1 <- Wpca%*%BnX
      U0 <- Wpca%*%AnX
      out.rot.P <-  FMN__Rotate(Out$rot$X, U1, U0)

      # updating the generated outputs
      Out$rot[["P"]] <- out.rot.P

      # P-dynamics: Z as risk factors:
      out.rot.P.P <- list( K0Z, K1Z, SSZ)
      names(out.rot.P.P) <- c("K0", "K1", "SS" )
      Out$rot$P[["P"]] <- out.rot.P.P


      # P-dynamics: PCN as risk factors:
      K0P <- K0Z[idxSpanned]
      K1P <- K1Z[idxSpanned, idxSpanned]
      SSP <- SSZ[idxSpanned, idxSpanned]
      out.rot.P.P2 <- list(K0P, K1P, SSP)
      names(out.rot.P.P2) <- c("K0", "K1", "SS" )
      out.rot.X.P <- FMN__Rotate(out.rot.P.P2, solve(U1), -solve(U1,U0))
      # updating the generated outputs, once again...
      Out$rot$X[["P"]] <- out.rot.X.P


      # PC as risk factors:
      PC1NW <- pca_weights_one_country(Y, Economy)
      PC1NW <- PC1NW[1:N,]*100

      U1 <- PC1NW%*%BnX
      U0 <- PC1NW%*%AnX
      out.rot.PC <- FMN__Rotate(Out$rot$X, U1, U0)

      Out$rot[["PC"]] <- out.rot.PC

      # P-dynamics: PC as risk factors:
      out.rot.PC.P <- list( K0Z, K1Z, SSZ)
      names(out.rot.PC.P) <- c("K0", "K1", "SS" )
      Out$rot$PC[["P"]] <- out.rot.PC.P
    }

  }

  if (nargout ==1){
    return(y)
  } else{
    return(Out)
  }


}

##############################################################################################################
#################################################################################################################
#' Compute the maximum likelihood function ("joint Q" models)
#'
#'@param K1XQ risk-neutral feedback matrix (NxN)
#'@param r0   long-run interest rate  (scalar)
#'@param SSZ  variance-covariance matrix (KxK)
#'@param K0Z  intercept from the P-dynamics (Kx1)
#'@param K1Z  feedback matrix from the P-dynamics (KxK)
#'@param se   Variance of the portfolio of yields observed with error (scalar)
#'@param Gy.0 matrix of contemporaneous terms from the P-dynamics (KxK)
#'@param mat  vector of maturities (in years) of yields used in estimation (J x 1)
#'@param Y    matrix of yields used in estimation  (J x T)
#'@param Z    complete set of spanned and unspanned factors (KxT)
#'@param P    complete set of spanned factors (NxT)
#'@param Wpca matrix of weights of the portfolios observed without errors (NxJ)
#'@param We   matrix of weights of the portfolios observed with errors ((J-N)xJ)
#'@param WpcaFull composite matrix of weights the portfolios observed with and without errors
#'@param dt  time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#'@param Economies set of economies that are part of the economic system (vector of text)
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param ModelType feasible options are (i) "VAR jointQ", (ii) "GVAR jointQ" or (iii) "JLL jointSigma"
#'@param GVARinputs if the model chosen is the "GVAR sepQ", the "GVARinputs" should be specified (see "GVAR" function)
#'@param JLLinputs if the model chosen is the "JLL jointSigma". "JLLinputs" should contain (i) DomUnit, (ii) WishSigmas, (iii) SigmaNonOrtho, (iv) JLLModelType (See JLL function)
#'@param nargout if nargout== 1: provides only the values of the likelihood; if nargout== 2: complete ATSM outputs
#'
#'@importFrom pracma mldivide mrdivide
#'
#'@references
#' This function is an extended version of the "A0N_MLEdensity_WOE" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal



MLEdensity_jointQ <- function(K1XQ, r0, SSZ, K0Z, K1Z, se, Gy.0, mat, Y, Z, P, Wpca, We, WpcaFull,
                              dt, Economies, FactorLabels, ModelType, GVARinputs, JLLinputs, nargout){


  # Country-specific inputs
  C <- length(Economies) # Number of economies
  G <- length(FactorLabels$Global)
  N <- nrow(K1XQ)/C
  J <- nrow(Y)/C
  T <- dim(Y)[2]
  M <- (nrow(Z)-G)/C-N # number of country-specific unspanned factors

  # Country-specific inputs
  CJ <- C*J # total number of all yields used in estimation;
  CN <- C*N # total number of country-specific spanned factors of the entire system

  Peo <- We%*%Y # portfolio observed WITH errors
  t <- 2:T

  if (ModelType == 'JLL jointSigma'){
    Sigma_Ye <- t(chol(SSZ)) # in this case, SSZ is the Variance-covariance matrix of the ORTHOGONALIZED dynamics
    JLLinputs$WishSigmas == 0 # Avoid recomputing the variance-covariance matrices
    JLLPara <- JLL(Z, N, JLLinputs)
    PI <-  JLLPara$PI
    Sigma_Y <- PI%*%Sigma_Ye
    SSZ <- Sigma_Y%*%t(Sigma_Y) # Variance-covariance matrix non-orthogonalized
  }



  # Intialize some variables, if necessary
  if (!exists("r0")){ r0 <- as.numeric()}
  if (!exists("se")){ se <- as.numeric()}
  if (!exists("K0Z")){ K0Z <- matrix(, nrow=N, ncol=0)}
  if (!exists("K1Z")){ K1Z <- matrix(, nrow=N*N, ncol=0)}

  if (sjmisc::is_empty(K0Z) || sjmisc::is_empty(K1Z) ){
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, K0Z, K1Z )
  } else{
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, t(t(as.vector(K0Z))), t(t(as.vector(K1Z))) )
  }
  # x is a column vector that contains all the elements of each parameter of the model.
  if ( any(is.nan(x)) ||any(is.infinite(x)) || any(Im(x)!= 0) ){
    y <- 1e6*rep(1, times= T-1)
    out <- as.numeric()
  }else{

    # Loadings:
    # Yields can be affine function of Pt's,i.e., Y(t) = AnP +BnP*P(t).
    # Further, we define Z(t) as an affine function of P(t) such that: Z(t) = phi0+ phi1*P(t).
    # As such, we can write P(t) = phi1^(-1)*(Z(t) - phi0).
    BnX <- A0N__computeBnAn_jointQ(round(mat/dt),K1XQ, dX= NULL, r0= NULL, SSX= NULL, Economies)[[1]]/dt
    # NOTE: the function "A0N__computeBnAn" generates outputs for interest rates per unit of time interval.
    # Hence by multiplying by "1/dt" we obtain annualized results.

    BnP <- mrdivide(BnX,Wpca%*%BnX) #  BnP = BnX*(W*BnX)^(-1)


    # Covariance matrix of latent states X:
    WBX <- Wpca%*%BnX # WBX = W * BnX
    b <- IdxSpanned(G,M,N,C)
    SSP <- SSZ[b, b]
    SSX <- mrdivide(solve(WBX,SSP, tol = 1e-50),t(WBX)) # SSX = (W*BnX)^(-1)*SSP*t((*W*BnX)^(-1))

    # Optimal estimate of r0:
    betan <- A0N__computeBnAn_jointQ(round(mat/dt), K1XQ, dX=NULL , r0=rep(0, times=C), SSX, Economies)[[3]]  # Note: r0=0, because when r0=0, An = Betan
    betan <- t(t(betan/dt))


    if (!exists("r0")|| sjmisc::is_empty(r0) ){
      # (i) A0 = (I - Bx(W*Bx)^(-1)*W)*Ax0;
      # (ii) A1 = (I - Bx(W*Bx)^(-1)*W)*Axr;
      # (iii) APer= We * A1
      # (iv) vt = We*(Yt - Bx(W*Bx)^(-1)*Pt - A0)
      # Recall that BnP= Bx(W*Bx)^(-1)

      r0 <- c()

      idxA <- 0
      idxB <- 0
      idxC <- 0

      for (i in 1:C){
        idxAA <- idxA + N
        idxBB <- idxB + J
        idxCC <- idxC + J-N

        BnPCS <- BnP[(idxB+1):idxBB , (idxA+1):idxAA]
        WpcaCS <- Wpca[(idxA+1):idxAA,(idxB+1):idxBB]
        betanCS <- betan[(idxB+1):idxBB]
        WeCS <- We[(idxC+1):(idxCC),(idxB+1):idxBB]
        YCS <- Y[(idxB+1):idxBB,]
        PCS <- P[(idxA+1):idxAA,]

        A0 <- (diag(J)-BnPCS%*%WpcaCS)%*%betanCS
        A1 <- (diag(J)-BnPCS%*%WpcaCS)%*%matrix(1,J,1)/dt

        # r0 = APer'*vt/(APer'*APer)
        r0[i] <- solve( (t(A1)%*%(t(WeCS)%*%WeCS)%*%A1), # Numerator from r0 equation
                        (t(A1)%*%(t(WeCS)%*%WeCS))%*%(rowMeans(YCS[,t] - BnPCS%*%PCS[,t]) - A0) ) # Denominator from r0 equation

        idxA <- idxAA
        idxB <- idxBB
        idxC <- idxCC
      }
    }


    # compute the intercepts:
    for (i in 1:C){
      Qtemp <- t(t(rep(r0[i], times= J)/dt))
      if (i==1){
        Q <- Qtemp
      }else{
        Q <- rbind(Q, Qtemp)
      }
    }

    AnX <- betan + Q
    AnP <- AnX -BnP%*%(Wpca%*%AnX) # AnP = AnX - BnX*(W*BnX)^(-1)*W*AnX

    # density of yields pricing errors:
    Pe <- matrix(0, nrow= nrow(We) , ncol= T)
    for( f in 1:T){
      Pe[,f] <- We%*% (BnP%*%P[,f] + AnP) # Pe_t = We*(AP+BP*P_t): "model-implied" yields for the portfolio of yields observed with errors
    }
    eQ <- Peo[,t] - Pe[,t]

    if (!exists("se")|| sjmisc::is_empty(se) ){
      se<- rep(NA, times= C)
      idx0 <-0
      for (i in 1:C){
        idx1 <- idx0+ J-N
        se[i] <- sqrt(mean(as.vector(eQ[(idx0+1):idx1,]^2))) # Scalar - standard deviation of the measurament error.
        idx0 <- idx1
      }
    }


    ## The log-likelihood function:
    if (any(is.nan(se) ==1)){
      y <- 1e6*matrix(1,T-1,1);
    }else{
      # Cross-sectional density (i.e. density for the portfolios observed with measurament error)
      aa <- se
      idx0 <- 0
      for (h in 1:C){
        idx1 <- idx0+ J-N
        se[(idx0+1):idx1] <- rep(aa[h], times=J-N)
        idx0 <- idx1
      }

      y <- GaussianDensity(eQ, se^2*diag(CJ-CN))

      se <- se[seq(1, CJ-CN, by=J-N)] # Recast se

      # time series density:
      if ((!exists("K0Z")|| sjmisc::is_empty(K0Z)) & (ModelType == 'VAR jointQ')){
        VARpara <- VAR(Z, VARtype= "unconstrained", Bcon = NULL)
        K0Z <- VARpara$K0Z  # Column vector (Kx1)
        K1Z <-  VARpara$K1Z # matrix (KxK)
      }
      if ( (!exists("K0Z")|| sjmisc::is_empty(K0Z)) & (ModelType == 'GVAR sepQ' || ModelType == 'GVAR jointQ') ){
        # Estimate GVAR(1)
          GVARPara <- GVAR(GVARinputs, N)
          K0Z <- GVARPara$F0 # Column vector (Kx1)
          K1Z <- GVARPara$F1 # matrix (KxK)
          }
      if ((!exists("K0Z")|| sjmisc::is_empty(K0Z)) & (ModelType == 'JLL jointSigma')){
        K0Z <-  JLLPara$k0 # Column vector (Kx1)
        K1Z <-  JLLPara$k1 # matrix (KxK)
        }


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
      y = -t(as.vector(y)) # for minimization;
    }


    # Variance of yields
    Sigma_e_CS <- c()
    Sigma_e_AllRows <- c()

    for (i in 1:C){
      Sigma_e_CS  <- c(rep(0, times = N), rep(se[i],times= J-N)^2) # Variances of portfolio WITHOUT and WITH errors
      Sigma_e_AllRows <- append(Sigma_e_AllRows, Sigma_e_CS)
    }

    Sigma_e_AllCountries <- Sigma_e_AllRows*diag(CJ)
    SIGMA_yields <- solve(WpcaFull)%*%Sigma_e_AllCountries%*%t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)

    names(VarYields) <- rownames(Y)
    VarYields <- t(t(VarYields)) # make a column vector

    if (nargout>1){


      if (ModelType == 'GVAR jointQ') { inputs <- list(Y, Z, mat, N, dt, Wpca, GVARinputs$Wgvar) }
      if (ModelType == 'VAR jointQ') { inputs <- list(Y, Z, mat, N, dt, Wpca) }
      if (ModelType == "JLL jointSigma") { inputs <- list(Y, Z, mat, N, dt, Wpca, JLLinputs$DomUnit) }

      if (ModelType == "JLL jointSigma" ) {
        JLLinputs$WishSigmas <- 1
        JLLinputs$SigmaNonOrtho <- SSZ
        JLLoutcomesOrtho <- JLL(Z, N, JLLinputs)
        # Remove the non-orthogonalized outputs
        JLLoutcomesOrtho$k0 <- NULL
        JLLoutcomesOrtho$k1 <- NULL
        #Add the zeros to the variance-covariance matrix (to avoid rounding problems)
        ZeroIdxVarCov <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$VarCovOrtho
        ZeroIdxSigma_Ye <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$Sigma_Ye

        JLLoutcomesOrtho$Sigmas$VarCov_Ortho[ZeroIdxVarCov] <- 0
        JLLoutcomesOrtho$Sigmas$Sigma_Ye[ZeroIdxSigma_Ye] <- 0
        PI<- JLLoutcomesOrtho$PI
        JLLoutcomesOrtho$Sigmas$Sigma_Y <- PI%*%JLLoutcomesOrtho$Sigmas$Sigma_Ye

        ests <- list(K1XQ, SSZ, SSP, r0, se, K0Z, K1Z, Gy.0, VarYields, JLLoutcomesOrtho)
      }else{
        ests <- list(K1XQ, SSZ, SSP, r0, se, K0Z, K1Z, Gy.0, VarYields)
      }

      llk <- list(-t(y)) # [1 x T-p] actual densities (not negative of densities)


      # Q dynamics: X as risk factors:
      Q <- list(zeros(CN,1), K1XQ, SSX)
      X <- list(BnX, AnX, Q)
      rot <- list(X)


      Out <- list(inputs, ests, llk, rot)
      names(Out) <- c("inputs", "ests", "llk", "rot")

      if (ModelType == 'GVAR jointQ') { names(Out$inputs)  <-c('Y','AllFactors' , 'mat', 'N', 'dt', 'Wpca', 'Wgvar')}
      if (ModelType == 'VAR jointQ') { names(Out$inputs)  <-c('Y','AllFactors' , 'mat', 'N', 'dt', 'Wpca')}
      if (ModelType == "JLL jointSigma"){ names(Out$inputs) <- c('Y','AllFactors' , 'mat', 'N', 'dt', 'Wpca', 'DomUnit') }


      if (ModelType == "JLL jointSigma" ) {
        names(Out$ests) <- c('K1XQ', 'SSZ', 'SSP', 'r0', 'se', 'K0Z', 'K1Z', 'Gy.0', 'VarYields', 'JLLoutcomes')
      }else{
        names(Out$ests) <- c('K1XQ', 'SSZ', 'SSP', 'r0', 'se', 'K0Z', 'K1Z', 'Gy.0', 'VarYields')
      }


      names(Out$rot) <- c("X")
      names(Out$rot$X) <-  c("B","A","Q")
      names(Out$rot$X$Q) <- c("K0", "K1", "SS" )


      # Q dynamics: PCN as risk factors: % PCN = U0 + U1*X, where U0 = W*AnX and U1 = W*BnX
      U1 <- Wpca%*%BnX
      U0 <- Wpca%*%AnX
      out.rot.P <-  FMN__Rotate(Out$rot$X, U1, U0)

      # updating the generated outputs...
      Out$rot[["P"]] <- out.rot.P


      # P dynamics: Z as risk factors:
      out.rot.P.P <- list( K0Z, K1Z, SSZ)
      names(out.rot.P.P) <- c("K0", "K1", "SS" )
      Out$rot$P[["P"]] <- out.rot.P.P


      # P dynamics: PCN as risk factors:
      K0P <- K0Z[b]
      K1P <- K1Z[b, b]
      SSP <- SSZ[b, b]
      out.rot.P.P2 <- list(K0P, K1P, SSP)
      names(out.rot.P.P2) <- c("K0", "K1", "SS" )
      out.rot.X.P <- FMN__Rotate(out.rot.P.P2, solve(U1,tol = 1e-50), -solve(U1,U0, tol = 1e-50))
      # updating the generated outputs... Once again...
      Out$rot$X[["P"]] <- out.rot.X.P


      U1 <- Wpca%*%BnX
      U0 <- Wpca%*%AnX
      out.rot.PC <- FMN__Rotate(Out$rot$X, U1, U0)

      Out$rot[["PC"]] <- out.rot.PC

      # P dynamics: PC as risk factors:
      out.rot.PC.P <- list( K0Z, K1Z, SSZ)
      names(out.rot.PC.P) <- c("K0", "K1", "SS" )
      Out$rot$PC[["P"]] <- out.rot.PC.P
    }

  }

  if (nargout ==1){
    return(y)
  } else{
    return(Out)
  }


}

#################################################################################################################
#################################################################################################################
#' Compute the maximum likelihood function ("joint Q" models for separate Sigma estimation)
#'
#'@param K1XQ risk-neutral feedback matrix (NxN)
#'@param r0   long-run interest rate  (scalar)
#'@param SSZ  variance-covariance matrix (KxK)
#'@param K0Z  intercept from the P-dynamics (Kx1)
#'@param K1Z  feedback matrix from the P-dynamics (KxK)
#'@param se   Variance of the portfolio of yields observed with error (scalar)
#'@param Gy.0 matrix of contemporaneous terms from the P-dynamics (KxK)
#'@param mat  vector of maturities (in years) of yields used in estimation (J x 1)
#'@param Y    matrix of yields used in estimation  (J x T)
#'@param Z    complete set of spanned and unspanned factors (KxT)
#'@param P    complete set of spanned factors (NxT)
#'@param Wpca matrix of weights of the portfolios observed without errors (NxJ)
#'@param We   matrix of weights of the portfolios observed with errors ((J-N)xJ)
#'@param WpcaFull composite matrix of weights the portfolios observed with and without errors
#'@param dt  time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#'@param Economies Set of economies that are part of the economic system (vector of text)
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param ModelType feasible options are (i) "JLL original" or (ii) "JLL NoDomUnit"
#'@param JLLinputs if the model chosen is the "JLL jointSigma", "JLLinputs" should be specified (see "JLL" function)
#'@param nargout if nargout== 1: provides only the values of the likelihood; if nargout== 2: complete ATSM outputs
#'
#'@importFrom pracma mldivide mrdivide
#'
#'@references
#' This function is an extended version of the "A0N_MLEdensity_WOE" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal



MLEdensity_jointQ_sepSigma <- function(K1XQ, r0, SSZ, K0Z, K1Z, se, Gy.0, mat, Y, Z, P, Wpca, We, WpcaFull,
                                       dt, Economies, FactorLabels, ModelType, JLLinputs, nargout){


  # Country-specific inputs
  C <- length(Economies) # Number of economies
  G <- length(FactorLabels$Global)
  N <- nrow(K1XQ)/C
  J <- nrow(Y)/C
  T <- dim(Y)[2]
  M <- (nrow(Z)-G)/C-N # number of country-specific unspanned factors


  # Country-specific inputs
  CJ <- C*J # number of all yields used in estimation;
  CN <- C*N # number of country-specific spanned factors


  Peo <- We%*%Y # portfolio observed WITH errors
  t <- 2:T


  # Intialize variables, if necessary
  if (!exists("SSZ")){SSZ <- matrix(, nrow=N*N, ncol=0)}
  if (!exists("r0")){ r0 <- as.numeric()}
  if (!exists("se")){ se <- as.numeric()}
  if (!exists("K0Z")){ K0Z <- matrix(, nrow=N, ncol=0)}
  if (!exists("K1Z")){ K1Z <- matrix(, nrow=N*N, ncol=0)}


  if (sjmisc::is_empty(K0Z) || sjmisc::is_empty(K1Z) ){
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, K0Z, K1Z )
  } else{
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, t(t(as.vector(K0Z))), t(t(as.vector(K1Z))) )
  }
  # x is a column vector that contains all the elements of each parameter of the model.

  if ( any(is.nan(x)) ||any(is.infinite(x)) || any(Im(x)!= 0) ){
    y <- 1e6*rep(1, times= T-1)
    out <- as.numeric()
  }else{

    # Loadings:
    # Yields can be affine function of Pt's,i.e., Y(t) = AnP +BnP*P(t).
    # Further, we define Z(t) as an affine function of P(t) such that: Z(t) = phi0+ phi1*P(t).
    # As such, we can write P(t) = phi1^(-1)*(Z(t) - phi0).
    BnX <- A0N__computeBnAn_jointQ(round(mat/dt),K1XQ, dX= NULL, r0= NULL, SSX= NULL, Economies)[[1]]/dt
    # NOTE: the function "A0N__computeBnAn" generates outputs for interest rates per unit of time interval.
    # Hence by multiplying by "1/dt" we obtain annualized results.
    BnP <- mrdivide(BnX,Wpca%*%BnX) #  BnP = BnX*(W*BnX)^(-1)


    WBX <- Wpca%*%BnX # WBX = W * BnX
    b <- IdxSpanned(G,M,N,C)
    SSP <- SSZ[b, b]
    SSX <- mrdivide(solve(WBX,SSP, tol = 1e-50),t(WBX)) # SSX = (W*BnX)^(-1)*SSP*t((*W*BnX)^(-1))

    # Optimal estimate of r0:
    betan <- A0N__computeBnAn_jointQ(round(mat/dt), K1XQ, dX=NULL , r0=rep(0, times=C), SSX, Economies)[[3]]  # Note: r0=0, because when r0=0, An = Betan
    betan <- t(t(betan/dt))


    if (!exists("r0")|| sjmisc::is_empty(r0) ){
      # (i) A0 = (I - Bx(W*Bx)^(-1)*W)*Ax0;
      # (ii) A1 = (I - Bx(W*Bx)^(-1)*W)*Axr;
      # (iii) APer= We * A1
      # (iv) vt = We*(Yt - Bx(W*Bx)^(-1)*Pt - A0)
      # Recall that BnP= Bx(W*Bx)^(-1)

      r0 <- c()

      idxA <- 0
      idxB <- 0
      idxC <- 0
      for (i in 1:C){
        idxAA <- idxA + N
        idxBB <- idxB + J
        idxCC <- idxC + J-N

        BnPCS <- BnP[(idxB+1):idxBB , (idxA+1):idxAA]
        WpcaCS <- Wpca[(idxA+1):idxAA,(idxB+1):idxBB]
        betanCS <- betan[(idxB+1):idxBB]
        WeCS <- We[(idxC+1):idxCC,(idxB+1):idxBB]
        YCS <- Y[(idxB+1):idxBB,]
        PCS <- P[(idxA+1):idxAA,]

        A0 <- (diag(J)-BnPCS%*%WpcaCS)%*%betanCS
        A1 <- (diag(J)-BnPCS%*%WpcaCS)%*%matrix(1,J,1)/dt

        #r0 = APer'*vt/(APer'*APer)
        r0[i] <- solve( (t(A1)%*%(t(WeCS)%*%WeCS)%*%A1), # Numerator from r0 equation
                        (t(A1)%*%(t(WeCS)%*%WeCS))%*%(rowMeans(YCS[,t] - BnPCS%*%PCS[,t]) - A0) ) # Denominator from  r0 equation


        idxA <- idxAA
        idxB <- idxBB
        idxC <- idxCC
      }
    }


    # compute the intercepts:
    for (i in 1:C){
      Qtemp <- t(t(rep(r0[i], times= J)/dt))
      if (i==1){
        Q <- Qtemp
      }else{
        Q <- rbind(Q, Qtemp)
      }
    }

    AnX <- betan + Q
    AnP <- AnX -BnP%*%(Wpca%*%AnX) # AnP = AnX - BnX*(W*BnX)^(-1)*W*AnX


    # density of yields pricing errors:
    Pe <- matrix(0, nrow= nrow(We) , ncol= T)
    for( f in 1:T){
      Pe[,f] <- We%*% (BnP%*%P[,f] + AnP) # Pe_t = We*(AP+BP*P_t): "model-implied" yields for the portfolio of yields observed with errors
    }
    eQ <- Peo[,t] - Pe[,t]

    if (!exists("se")|| sjmisc::is_empty(se) ){
      se<- rep(NA, times= C)
      idx0 <-0
      for (i in 1:C){
        idx1 <- idx0+ J-N
        se[i] <- sqrt(mean(as.vector(eQ[(idx0+1):idx1,]^2))) # Scalar - standard deviation of the measurament error.
        idx0 <- idx1
      }
    }


    ## The log-likelihood function:
    if (any(is.nan(se) ==1)){
      y <- 1e6*matrix(1,T-1,1);
    }else{
      # Cross-sectional density (i.e. density for the portfolios observed with measurament error)
      aa <- se
      idx0 <- 0
      for (h in 1:C){
        idx1 <- idx0+ J-N
        se[(idx0+1):idx1] <- rep(aa[h], times=J-N)
        idx0 <- idx1
      }

      y <- GaussianDensity(eQ, se^2*diag(CJ-CN))

      se <- se[seq(1, CJ-CN, by=J-N)] # Recast se

      # time series density:
      if (!exists("K0Z")|| sjmisc::is_empty(K0Z)){
      JLLinputs$WishSigmas <- 1
      JLLinputs$SigmaNonOrtho <- SSZ

      JLLfac <- JLL(Z, N, JLLinputs)
      K0Z <- JLLfac$k0
      K1Z <- JLLfac$k1
      Gy.0 <- diag(C*(M+N) + G)
}

      eP <- matrix(0, nrow(K0Z), T-1)
      for (f in 1:length(t)){ eP[,f] <-  Z[,f+1] - K1Z%*%Z[,f] - K0Z }

      y <- y + GaussianDensity(eP, SSZ) # Cross-sectional density + time-series density (final likelihood function, except for the Jacobian term)

      # Jacobian: (so that the density is for the observed yields (and the non-yields variables Z)
      #           and not the yields portfolios which depend on W).

      y <-  y + 0.5*log(abs(det(Wpca%*%t(Wpca))))
    }

    if ( any(is.nan(y)) ||any(is.infinite(y)) || any(Im(y)!= 0) ){
      y <- 1e6*rep(1, times=T-1)
    }else{
      y = -t(as.vector(y)) # for minimization;
    }


    # Variance of yields
    Sigma_e_CS <- c()
    Sigma_e_AllRows <- c()

    for (i in 1:C){
      Sigma_e_CS  <- c(rep(0, times = N), rep(se[i],times= J-N)^2) # Variances of portfolio WITHOUT and WITH errors
      Sigma_e_AllRows <- append(Sigma_e_AllRows, Sigma_e_CS)
    }

    Sigma_e_AllCountries <- Sigma_e_AllRows*diag(CJ)
    SIGMA_yields <- solve(WpcaFull)%*%Sigma_e_AllCountries%*%t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)
    names(VarYields) <- rownames(Y)
    VarYields <- t(t(VarYields))

    if (nargout>1){

      #Inputs:
      inputs <- list(Y, Z, mat, N, dt, Wpca, JLLinputs$DomUnit)

      #Estimates:
      JLLinputs$WishSigmas <- 1
      JLLinputs$SigmaNonOrtho <- SSZ
      JLLoutcomesOrtho <- JLL(Z, N, JLLinputs)
      ZeroIdxVarCov <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$VarCovOrtho
      ZeroIdxSigma_Ye <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$Sigma_Ye
      PI<- JLLoutcomesOrtho$PI

      # Remove the non-orthogonalized outputs
      JLLoutcomesOrtho$k0 <- NULL
      JLLoutcomesOrtho$k1 <- NULL
      #Add the orthogonalized Variance-covariance matrix (set zeros to avoid rounding prolems)
      JLLoutcomesOrtho$Sigmas$VarCov_Ortho[ZeroIdxVarCov] <- 0
      JLLoutcomesOrtho$Sigmas$Sigma_Ye[ZeroIdxSigma_Ye] <- 0
      JLLoutcomesOrtho$Sigmas$Sigma_Y <- PI%*%JLLoutcomesOrtho$Sigmas$Sigma_Ye

      ests <- list(K1XQ, SSZ, SSP, r0, se, K0Z, K1Z, Gy.0, VarYields, JLLoutcomesOrtho)

      llk <- list(-t(y)) # [1 x T-p] actual densities (not negative of densities)

      # Q dynamics: X as risk factors:
      Q <- list(zeros(CN,1), K1XQ, SSX)
      X <- list(BnX, AnX, Q)
      rot <- list(X)

      Out <- list(inputs, ests, llk, rot)
      names(Out) <- c("inputs", "ests", "llk", "rot")
      names(Out$ests) <- c('K1XQ', 'SSZ', 'SSP', 'r0', 'se', 'K0Z', 'K1Z', 'Gy.0', 'VarYields', 'JLLoutcomes')
      names(Out$inputs)  <-c('Y','AllFactors' , 'mat', 'N', 'dt', 'Wpca', 'DomUnit')


      names(Out$rot) <- c("X")
      names(Out$rot$X) <-  c("B","A","Q")
      names(Out$rot$X$Q) <- c("K0", "K1", "SS" )


      # Q dynamics: PCN as risk factors: % PCN = U0 + U1*X, where U0 = W*AnX and U1 = W*BnX
      U1 <- Wpca%*%BnX
      U0 <- Wpca%*%AnX
      out.rot.P <-  FMN__Rotate(Out$rot$X, U1, U0)

      # updating the generated outputs.
      Out$rot[["P"]] <- out.rot.P


      # P dynamics: Z as risk factors:
      out.rot.P.P <- list( K0Z, K1Z, SSZ)
      names(out.rot.P.P) <- c("K0", "K1", "SS" )
      Out$rot$P[["P"]] <- out.rot.P.P


      # P dynamics: PCN as risk factors:
      K0P <- K0Z[b]
      K1P <- K1Z[b, b]
      SSP <- SSZ[b, b]
      out.rot.P.P2 <- list(K0P, K1P, SSP)
      names(out.rot.P.P2) <- c("K0", "K1", "SS" )
      out.rot.X.P <- FMN__Rotate(out.rot.P.P2, solve(U1), -solve(U1,U0))
      # updating the generated outputs, once again...
      Out$rot$X[["P"]] <- out.rot.X.P


      U1 <- Wpca%*%BnX
      U0 <- Wpca%*%AnX
      out.rot.PC <- FMN__Rotate(Out$rot$X, U1, U0)

      Out$rot[["PC"]] <- out.rot.PC

      # P dynamics: PC as risk factors:
      out.rot.PC.P <- list( K0Z, K1Z, SSZ)
      names(out.rot.PC.P) <- c("K0", "K1", "SS" )
      Out$rot$PC[["P"]] <- out.rot.PC.P
    }

  }

  if (nargout ==1){
    return(y)
  } else{
    return(Out)
  }


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


  K <- C*(N+M) +G

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


#################################################################################################################
#################################################################################################################
################################ BOOTSTRAP FUNCTIONS ############################################################
#################################################################################################################
#################################################################################################################
#' Compute the maximum likelihood function ("sep Q" models) - Bootstrap version
#'
#'@param K1XQ risk-neutral feedback matrix (NxN)
#'@param r0   long-run interest rate  (scalar)
#'@param SSZ  variance-covariance matrix (KxK)
#'@param K0Z  intercept from the P-dynamics (Kx1)
#'@param K1Z  feedback matrix from the P-dynamics (KxK)
#'@param se   Variance of the portfolio of yields observed with error (scalar)
#'@param Gy.0 matrix of contemporaneous terms from the P-dynamics (KxK)
#'@param mat  vector of maturities (in years) of yields used in estimation (J x 1)
#'@param Y    matrix of yields used in estimation  (J x T)
#'@param Z    complete set of spanned and unspanned factors (KxT)
#'@param P    complete set of spanned factors (NxT)
#'@param Wpca matrix of weights of the portfolios observed without errors (NxJ)
#'@param We   matrix of weights of the portfolios observed with errors ((J-N)xJ)
#'@param WpcaFull composite matrix of weights the portfolios observed with and without errors
#'@param dt  time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#'@param Economy Name of the economies under study
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param ModelType Feasible options are: (i) "JPS", (ii) "JPS jointP" or (iii) "GVAR sepQ"
#'@param residBS  index of the re-ordered bootstrap residuals
#'@param MaxEigen largest eigenvalue under the P-dynamics
#'@param GVARinputs if the model chosen is the "GVAR sepQ", "GVARinputs" should be specified (see "GVAR" function)
#'@param nargout if nargout== 1: provides only the values of the likelihood; if nargout== 2: complete ATSM outputs
#'
#'
#'@importFrom pracma mldivide mrdivide
#'
#'@references
#' This function is modified version of the "A0N_MLEdensity_WOE" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal



A0N_MLEdensity_WOE__sepQ_Bootstrap <- function(K1XQ, r0, SSZ, K0Z, K1Z, se, Gy.0, mat, Y, Z, P, Wpca, We, WpcaFull,
                                               dt, Economy, FactorLabels, ModelType, residBS, MaxEigen, GVARinputs,
                                               nargout){



  N <- length(FactorLabels$Spanned) # number of country-specific spanned factors
  M <- length(FactorLabels$Domestic) - N # number of country-specific unspanned factors
  J <- numel(mat) # number of country-specific yields used in estimation;
  T <- dim(Y)[2]

  Peo <- We%*%Y # portfolio observed WITH errors
  t <- 2:T


  if (ModelType == "JPS"){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]]) }
  if (ModelType == "JPS jointP" || ModelType == 'GVAR sepQ'){ AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)}

  rownames(SSZ) <- AllLabels
  colnames(SSZ) <- AllLabels
  LabelSpannedCS <- c(FactorLabels$Tables[[Economy]][-(1:M)])
  idxSpanned <- match(LabelSpannedCS, AllLabels)

  ##################################################################################
  if (!exists("r0")){ r0 <- as.numeric()}
  if (!exists("se")){ se <- as.numeric()}
  if (!exists("K0Z")){ K0Z <- matrix(, nrow=N, ncol=0)}
  if (!exists("K1Z")){ K1Z <- matrix(, nrow=N*N, ncol=0)}
  #############################################################################
  if (sjmisc::is_empty(K0Z) || sjmisc::is_empty(K1Z) ){
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, K0Z, K1Z )
  } else{
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, t(t(as.vector(K0Z))), t(t(as.vector(K1Z))) )
  }
  # x is a column vector that contains all the elements of each parameter of the model.

  if ( any(is.nan(x)) ||any(is.infinite(x)) || any(Im(x)!= 0) ){
    y <- 1e6*rep(1, times = T-1)
    out <- as.numeric()
  }else{

    # Loadings:
    # Yields can be affine function of Pt's,i.e., Y(t) = AnP +BnP*P(t).
    # Further, we define Z(t) as an affine function of P(t) such that: Z(t) = phi0+ phi1*P(t).
    # As such, we can write P(t) = phi1^(-1)*(Z(t) - phi0).7
    BnX <- A0N__computeBnAn_sepQ(round(mat/dt),K1XQ, dX= NULL, r0= NULL, SSX= NULL)[[1]]/dt
    # NOTE: the function "A0N__computeBnAn" generates outputs for interest rates per unit of time interval.
    # Hence by multiplying by "1/dt" we obtain annualized results.
    BnP <- mrdivide(BnX,Wpca%*%BnX) #  BnP = BnX*(W*BnX)^(-1)


    rownames(BnX) <- rownames(Y)

    # Covariance matrix of latent states X:
    WBX <- Wpca%*%BnX # WBX = W * BnX
    SSP <- SSZ[idxSpanned , idxSpanned]
    SSX <- mrdivide(solve(WBX,SSP,tol = 1e-50),t(WBX)) # SSX = (W*BnX)^(-1)*SSP*t((*W*BnX)^(-1))
    # Optimal estimate of r0:
    betan <- A0N__computeBnAn_sepQ(round(mat/dt), K1XQ, dX=NULL , r0=0, SSX)[[3]]
    betan <- t(t(betan/dt))


    if (!exists("r0")|| sjmisc::is_empty(r0) ){
      # (i) A0 = (I - Bx(W*Bx)^(-1)*W)*Ax0;
      # (ii) A1 = (I - Bx(W*Bx)^(-1)*W)*Axr;
      # (iii) APer= We * A1
      # (iv) vt = We*(Yt - Bx(W*Bx)^(-1)*Pt - A0)
      # Recall that BnP= Bx(W*Bx)^(-1)
      A0 <- (diag(J)-BnP%*%Wpca)%*%betan
      A1 <- (diag(J)-BnP%*%Wpca)%*%matrix(1,J,1)/dt

      r0 <- solve((t(A1)%*%(t(We)%*%We)%*%A1), # Numerator of r0
                  (t(A1)%*%(t(We)%*%We))%*%(rowMeans(Y[,t] - BnP%*%P[,t]) - A0)) # Denominator of r0
      # i.e: r0 = APer'*vt/(APer'*APer)

    }


    # compute the intercepts:
    AnX <- betan + as.numeric(r0/dt)
    AnP <- AnX -BnP%*%(Wpca%*%AnX) # AnP = AnX - BnX*(W*BnX)^(-1)*W*AnX

    rownames(AnX) <- rownames(Y)

    # density of yields pricing errors:
    Pe <- matrix(0, nrow= nrow(We) , ncol= T)
    for( f in 1:T){
      Pe[,f] <- We%*% (BnP%*%P[,f] + AnP) # Pe_t = We*(AP+BP*P_t): "model-implied" yields for the portfolio of yields observed with errors
    }
    eQ <- Peo[,t] - Pe[,t]


    if (!exists("se")|| sjmisc::is_empty(se) ){
      se <- sqrt(mean(as.vector(eQ^2))) # Scalar - standard deviation of the measurament error.
    }


    ## The log-likelihood function:
    if (is.nan(se)){
      y <- 1e6*matrix(1,T-1,1)
    }else{
      # Cross-sectional density (i.e. density for the portfolios observed with measurament error)
      y <- GaussianDensity(eQ, se^2*diag(J-N))

      # time series density:
      if ((!exists("K0Z")|| sjmisc::is_empty(K0Z)) & (ModelType == 'JPS' || ModelType == 'JPS jointP' )){
        VARpara <- VAR(Z, VARtype= "unconstrained", Bcon = NULL)
        K0Z <- VARpara$K0Z  # Column vector (Kx1)
        K1Z <-  VARpara$K1Z # matrix (KxK)
        Gy.0 <- diag(length(K0Z))
      }
      if ( (!exists("K0Z")|| sjmisc::is_empty(K0Z)) || (ModelType == 'GVAR sepQ') )  {
        # Estimate GVAR(1)
        GVARPara <- GVAR(GVARinputs, N)
        K0Z <- GVARPara$F0 # Column vector (Kx1)
        K1Z <- GVARPara$F1 # matrix (KxK)
        Gy.0 <- GVARPara$Gy.0
      }

      eP <- matrix(0, nrow(K0Z), T-1)
      for (f in 1:length(t)){ eP[,f] <-  Z[,f+1] - K1Z%*%Z[,f] - K0Z }

      y <- y + GaussianDensity(eP, SSZ) # Cross-sectional density + time-series density (final likelihood function, except for the Jacobian Term)

      # Jacobian: (so that the density is for the observed yields (and the non-yields variables Z)
      #           and not the yields portfolios which depend on W).

      y <-  y + 0.5*log(abs(det(Wpca%*%t(Wpca))))
    }

    if ( any(is.nan(y)) ||any(is.infinite(y)) || any(Im(y)!= 0) ){
      y <- 1e6*rep(1, times = T-1)
    }else{
      y = -t(as.vector(y)) # for minimization;
    }


    # Variance-yields
    Sigma_e_Row  <- c(rep(0, times = N), rep(se,times= J-N)^2) # Variances of portfolio WITHOUT and WITH errors

    Sigma_e <- Sigma_e_Row*diag(J)
    SIGMA_yields <- solve(WpcaFull)%*%Sigma_e%*%t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)


    if (nargout>1){

      ests <- list(K1XQ, SSZ, r0, se, K0Z, K1Z, Gy.0, MaxEigen, residBS, VarYields)
      llk <- mean(y) # Value of objective function at the maximum


      # Q dynamics: X as risk factors:
      Q <- list(zeros(N,1), K1XQ, SSX)
      X <- list(BnX, AnX, Q)
      rot <- list(X)

      Out <- list(ests, llk, rot)

      # Label list
      names(Out) <- c("ests", "llkMax", "rot")

      names(Out$ests) <- c('K1XQ', 'SSZ', 'r0', 'se', 'K0Z', 'K1Z', 'Gy.0' , 'MaxEigen', 'residBS', 'VarYields')

      names(Out$rot) <- c("X")

      names(Out$rot$X) <-  c("B","A","Q")
      names(Out$rot$X$Q) <- c("K0", "K1", "SS" )


      # Q dynamics: PCN as risk factors: % PCN = U0 + U1*X, where U0 = W*AnX and U1 = W*BnX
      U1 <- Wpca%*%BnX
      U0 <- Wpca%*%AnX
      out.rot.P <-  FMN__Rotate(Out$rot$X, U1, U0)

      # updating the generated outputs...
      Out$rot[["P"]] <- out.rot.P



      Out$rot$X <- NULL # Delete the parameters related to the latent factors
      #(Not needded for the Bootstrap, save some time and pc memory)

    }

  }

  if (nargout ==1){
    return(y)
  } else{
    return(Out)
  }


}

##################################################################################################################
#################################################################################################################
#' Compute the maximum likelihood function (joint Q models) - Bootstrap version
#'
#'@param K1XQ risk-neutral feedback matrix (NxN)
#'@param r0   long-run interest rate  (scalar)
#'@param SSZ  variance-covariance matrix (KxK)
#'@param K0Z  intercept from the P-dynamics (Kx1)
#'@param K1Z  feedback matrix from the P-dynamics (KxK)
#'@param se   Variance of the portfolio of yields observed with error (scalar)
#'@param Gy.0 matrix of contemporaneous terms from the P-dynamics (KxK)
#'@param mat  vector of maturities (in years) of yields used in estimation (J x 1)
#'@param Y    matrix of yields used in estimation  (J x T)
#'@param Z    complete set of spanned and unspanned factors (KxT)
#'@param P    complete set of spanned factors (NxT)
#'@param Wpca matrix of weights of the portfolios observed without errors (NxJ)
#'@param We   matrix of weights of the portfolios observed with errors ((J-N)xJ)
#'@param WpcaFull composite matrix of weights the portfolios observed with and without errors
#'@param dt  time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param ModelType feasible options are (i) "VAR jointQ", (ii) "GVAR jointQ" or (iii) "JLL jointSigma"
#'@param residBS  index of the re-ordered bootstrap residuals
#'@param MaxEigen largest eigenvalue under the P-dynamics
#'@param GVARinputs if the model chosen is the "GVAR sepQ", "GVARinputs" should be specified (see "GVAR" function )
#'@param JLLinputs if the model chosen is the "JLL jointSigma". "JLLinputs" should be specified (see "JLL" function)
#'@param nargout if nargout== 1: provides only the values of the likelihood; if nargout== 2: complete ATSM outputs
#'
#'@importFrom pracma mldivide mrdivide
#'
#'@references
#' This function is modified version of the "A0N_MLEdensity_WOE" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal




A0N_MLEdensity_WOE__jointQ_Bootstrap <- function(K1XQ, r0, SSZ, K0Z, K1Z, se, Gy.0, mat, Y, Z, P, Wpca, We, WpcaFull,
                                                dt, Economies, FactorLabels, ModelType, residBS, MaxEigen,
                                                GVARinputs, JLLinputs, nargout){


  # Country-specific inputs
  C <- length(Economies) # Number of economies
  G <- length(FactorLabels$Global)
  N <- nrow(K1XQ)/C
  J <- nrow(Y)/C
  T <- dim(Y)[2]
  M <- (nrow(Z)-G)/C-N # number of country-specific unspanned factors

  # Country-specific inputs
  CJ <- C*J # total number of all yields used in estimation;
  CN <- C*N # total number of country-specific spanned factors of the entire system

  Peo <- We%*%Y # portfolio observed WITH errors
  t <- 2:T
  ###########################################################################################
  if (!exists("r0")){ r0 <- as.numeric()}
  if (!exists("se")){ se <- as.numeric()}
  if (!exists("K0Z")){ K0Z <- matrix(, nrow=N, ncol=0)}
  if (!exists("K1Z")){ K1Z <- matrix(, nrow=N*N, ncol=0)}
  #############################################################################################
  if (sjmisc::is_empty(K0Z) || sjmisc::is_empty(K1Z) ){
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, K0Z, K1Z )
  } else{
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, t(t(as.vector(K0Z))), t(t(as.vector(K1Z))) )
  }
  # x is a column vector that contains all the elements of each parameter of the model.
  if ( any(is.nan(x)) ||any(is.infinite(x)) || any(Im(x)!= 0) ){
    y <- 1e6*rep(1, times = T-1)
    out <- as.numeric()
  }else{

    # Loadings:
    # Yields can be affine function of Pt's,i.e., Y(t) = AnP +BnP*P(t).
    # Further, we define Z(t) as an affine function of P(t) such that: Z(t) = phi0+ phi1*P(t).
    # As such, we can write P(t) = phi1^(-1)*(Z(t) - phi0).
    BnX <- A0N__computeBnAn_jointQ(round(mat/dt),K1XQ, dX= NULL, r0= NULL, SSX= NULL, Economies)[[1]]/dt # NOTE: the function "A0N__computeBnAn" generates outputs for interest rates per unit of time interval. Hence by multiplying by "1/dt" we obtain annualized results.

    BnP <- mrdivide(BnX,Wpca%*%BnX) #  BnP = BnX*(W*BnX)^(-1)
    rownames(BnX) <- rownames(Y)

    # Covariance matrix of latent states X:
    WBX <- Wpca%*%BnX # WBX = W * BnX
    b<- IdxSpanned(G,M,N,C)
    SSP <- SSZ[b, b]
    SSX <- mrdivide(solve(WBX,SSP,tol = 1e-50),t(WBX)) # SSX = (W*BnX)^(-1)*SSP*t((*W*BnX)^(-1))


    # Optimal estimate of r0:
    betan <- A0N__computeBnAn_jointQ(round(mat/dt), K1XQ, dX=NULL , r0=rep(0, times=C), SSX, Economies)[[3]]
    betan <- t(t(betan/dt))


    if (!exists("r0")|| sjmisc::is_empty(r0) ){
      # (i) A0 = (I - Bx(W*Bx)^(-1)*W)*Ax0;
      # (ii) A1 = (I - Bx(W*Bx)^(-1)*W)*Axr;
      # (iii) APer= We * A1
      # (iv) vt = We*(Yt - Bx(W*Bx)^(-1)*Pt - A0)
      # Recall that BnP= Bx(W*Bx)^(-1)

      r0 <- c()

      idxA <- 0
      idxB <- 0
      idxC <- 0

      for (i in 1:C){
        idxAA <- idxA + N
        idxBB <- idxB + J
        idxCC <- idxC + J-N

        BnPCS <- BnP[(idxB+1):idxBB , (idxA+1):idxAA]
        WpcaCS <- Wpca[(idxA+1):idxAA,(idxB+1):idxBB]
        betanCS <- betan[(idxB+1):idxBB]
        WeCS <- We[(idxC+1):(idxCC),(idxB+1):idxBB]
        YCS <- Y[(idxB+1):idxBB,]
        PCS <- P[(idxA+1):idxAA,]

        A0 <- (diag(J)-BnPCS%*%WpcaCS)%*%betanCS
        A1 <- (diag(J)-BnPCS%*%WpcaCS)%*%matrix(1,J,1)/dt

        r0[i] <- solve( (t(A1)%*%(t(WeCS)%*%WeCS)%*%A1), # Numerator from r0
                        (t(A1)%*%(t(WeCS)%*%WeCS))%*%(rowMeans(YCS[,t] - BnPCS%*%PCS[,t]) - A0) ) # Denominator from r0
        # i.e: r0 = APer'*vt/(APer'*APer)

        idxA <- idxAA
        idxB <- idxBB
        idxC <- idxCC
      }
    }


    # compute the intercepts:
    for (i in 1:C){
      Qtemp <- t(t(rep(r0[i], times= J)/dt))
      if (i==1){
        Q <- Qtemp
      }else{
        Q <- rbind(Q, Qtemp)
      }
    }

    AnX <- betan + Q
    AnP <- AnX -BnP%*%(Wpca%*%AnX) # AnP = AnX - BnX*(W*BnX)^(-1)*W*AnX

    rownames(AnX) <- rownames(Y)

    # density of yields pricing errors:
    Pe <- matrix(0, nrow= nrow(We) , ncol= T)
    for( f in 1:T){
      Pe[,f] <- We%*% (BnP%*%P[,f] + AnP) # Pe_t = We*(AP+BP*P_t): "model-implied" yields for the portfolio of yields observed with errors
    }
    eQ <- Peo[,t] - Pe[,t]

    if (!exists("se")|| sjmisc::is_empty(se) ){
      se<- rep(NA, times= C)
      idx0 <-0
      for (i in 1:C){
        idx1 <- idx0+ J-N
        se[i] <- sqrt(mean(as.vector(eQ[(idx0+1):idx1,]^2))) # Scalar - standard deviation of the measurament error.
        idx0 <- idx1
      }
    }


    ## The log-likelihood function:
    if (any(is.nan(se) ==1)){
      y <- 1e6*matrix(1,T-1,1)
    }else{
      # Cross-sectional density (i.e. density for the portfolios observed with measurament error)
      aa <- se
      idx0 <- 0
      for (h in 1:C){
        idx1 <- idx0+ J-N
        se[(idx0+1):idx1] <- rep(aa[h], times=J-N)
        idx0 <- idx1
      }

      y <- GaussianDensity(eQ, se^2*diag(CJ-CN))

      se <- se[seq(1, CJ-CN, by=J-N)] # Recast se

      # time series density:
      if ((!exists("K0Z")|| sjmisc::is_empty(K0Z)) & (ModelType == 'VAR jointQ')){
        VARpara <- VAR(Z, VARtype= "unconstrained", Bcon = NULL)
        K0Z <- VARpara$K0Z  # Column vector (Kx1)
        K1Z <-  VARpara$K1Z # matrix (KxK)
      }
      if ( (!exists("K0Z")|| sjmisc::is_empty(K0Z)) & (ModelType == 'GVAR sepQ' || ModelType == 'GVAR jointQ') ){
        # Estimate GVAR(1)
        GVARPara <- GVAR(GVARinputs, N)
        K0Z <- GVARPara$F0 # Column vector (Kx1)
        K1Z <- GVARPara$F1 # matrix (KxK)
      }
      if ((!exists("K0Z")|| sjmisc::is_empty(K0Z)) & (ModelType == 'JLL jointSigma')){
        JLLinputs$WishSigmas == 0 # Avoid recomputing the variance-covariance matrices
        JLLPara <- JLL(Z, N, JLLinputs)
        K0Z <-  JLLPara$k0 # Column vector (Kx1)
        K1Z <-  JLLPara$k1 # matrix (KxK)
      }



      eP <- matrix(0, nrow(K0Z), T-1)
      for (f in 1:length(t)){ eP[,f] <-  Z[,f+1] - K1Z%*%Z[,f] - K0Z }

      y <- y + GaussianDensity(eP, SSZ) # Cross-sectional density + time-series density (final likelihood function, except for the Jacobian Term

      # Jacobian: (so that the density is for the observed yields (and the non-yields variables Z)
      #           and not the yields portfolios which depend on W).

      y <-  y + 0.5*log(abs(det(Wpca%*%t(Wpca))))
    }

    if ( any(is.nan(y)) ||any(is.infinite(y)) || any(Im(y)!= 0) ){
      y <- 1e6*rep(1, times = T-1)
    }else{
      y = -t(as.vector(y)) # for minimization;
    }

    # Variance of yields
    Sigma_e_CS <- c()
    Sigma_e_AllRows <- c()

    for (i in 1:C){
      Sigma_e_CS  <- c(rep(0, times = N), rep(se[i],times= J-N)^2) # Variances of portfolio WITHOUT and WITH errors
      Sigma_e_AllRows <- append(Sigma_e_AllRows, Sigma_e_CS)
    }

    Sigma_e_AllCountries <- Sigma_e_AllRows*diag(CJ)
    SIGMA_yields <- solve(WpcaFull)%*%Sigma_e_AllCountries%*%t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)

    names(VarYields) <- rownames(Y)
    VarYields <- t(t(VarYields)) # make a column vector


    if (nargout>1){


      if (ModelType == "JLL jointSigma" ) {
        JLLinputs$WishSigmas <- 1
        JLLinputs$SigmaNonOrtho <- SSZ
        JLLoutcomesOrtho <- JLL(Z, N, JLLinputs)
        # Remove the non-orthogonalized outputs
        JLLoutcomesOrtho$k0 <- NULL
        JLLoutcomesOrtho$k1 <- NULL
        JLLoutcomesOrtho$a_W <- NULL
        JLLoutcomesOrtho$b <- NULL
        JLLoutcomesOrtho$c <- NULL
        JLLoutcomesOrtho$a_DU_CS <- NULL
        JLLoutcomesOrtho$PIb <- NULL
        JLLoutcomesOrtho$PIac <- NULL
        JLLoutcomesOrtho$Ye <- NULL
        #Add the orthogonalized Variance-covariance matrix (to avoid rounding problems)
        ZeroIdxVarCov <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$VarCovOrtho
        ZeroIdxSigma_Ye <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$Sigma_Ye
        JLLoutcomesOrtho$Sigmas$VarCov_Ortho[ZeroIdxVarCov] <- 0
        JLLoutcomesOrtho$Sigmas$Sigma_Ye[ZeroIdxSigma_Ye] <- 0

        PI<- JLLoutcomesOrtho$PI
        JLLoutcomesOrtho$Sigmas$Sigma_Y <- PI%*%JLLoutcomesOrtho$Sigmas$Sigma_Ye

        ests <- list(K1XQ, SSZ, r0, se, K0Z, K1Z, Gy.0, MaxEigen, VarYields, residBS, JLLoutcomesOrtho)
      }else{
        ests <- list(K1XQ, SSZ, r0, se, K0Z, K1Z, Gy.0, MaxEigen, residBS, VarYields)
      }

      llk <- mean(y) # [1 x T-p] actual densities (not negative of densities)


      # Q dynamics: X as risk factors:
      Q <- list(zeros(CN,1), K1XQ, SSX)
      X <- list(BnX, AnX, Q)
      rot <- list(X)

      Out <- list(ests, llk, rot)

      # Label list
      names(Out) <- c("ests", "llkMax", "rot")

      if (ModelType == "JLL jointSigma" ) {
        names(Out$ests) <- c('K1XQ', 'SSZ', 'r0', 'se', 'K0Z', 'K1Z', 'Gy.0', 'MaxEigen', 'residBS',
                             'VarYields', 'JLLoutcomes')
      }else{
        names(Out$ests) <- c('K1XQ', 'SSZ', 'r0', 'se', 'K0Z', 'K1Z', 'Gy.0', 'MaxEigen', 'residBS', 'VarYields')
      }

      names(Out$rot) <- c("X")
      names(Out$rot$X) <-  c("B","A","Q")
      names(Out$rot$X$Q) <- c("K0", "K1", "SS" )


      # Q dynamics: PCN as risk factors: % PCN = U0 + U1*X, where U0 = W*AnX and U1 = W*BnX
      U1 <- Wpca%*%BnX
      U0 <- Wpca%*%AnX
      out.rot.P <-  FMN__Rotate(Out$rot$X, U1, U0)

      # updating the generated outputs.
      Out$rot[["P"]] <- out.rot.P


      Out$rot$X <- NULL # Delete the parameters related to the latent factors
      #(Not needded for the Bootstrap, save some time and pc memory)

    }

  }

  if (nargout ==1){
    return(y)
  } else{
    return(Out)
  }


}

#################################################################################################################
#' Compute the maximum likelihood function ("joint Q" models for separate Sigma estimation) - Bootstrap version
#'
#'@param K1XQ risk-neutral feedback matrix (NxN)
#'@param r0   long-run interest rate  (scalar)
#'@param SSZ  variance-covariance matrix (KxK)
#'@param K0Z  intercept from the P-dynamics (Kx1)
#'@param K1Z  feedback matrix from the P-dynamics (KxK)
#'@param se   Variance of the portfolio of yields observed with error (scalar)
#'@param Gy.0 matrix of contemporaneous terms from the P-dynamics (KxK)
#'@param mat  vector of maturities (in years) of yields used in estimation (J x 1)
#'@param Y    matrix of yields used in estimation  (J x T)
#'@param Z    complete set of spanned and unspanned factors (KxT)
#'@param P    complete set of spanned factors (NxT)
#'@param Wpca matrix of weights of the portfolios observed without errors (NxJ)
#'@param We   matrix of weights of the portfolios observed with errors ((J-N)xJ)
#'@param WpcaFull composite matrix of weights the portfolios observed with and without errors
#'@param dt  time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'@param ModelType feasile options are (i) "JLL original" or (ii) "JLL NoDomUnit"
#'@param residBS  indexes of the re-ordered bootstrap residuals
#'@param MaxEigen largest eigenvalue under the P-dynamics
#'@param GVARinputs if the model chosen is the "GVAR sepQ", "GVARinputs" must be specified (see "GVAR" function )
#'@param JLLinputs if the model chosen is the "JLL jointSigma", "JLLinputs" must be specified (see "JLL" function)
#'@param nargout if nargout== 1: provides only the values of the likelihood; if nargout== 2: complete ATSM outputs
#'
#'@importFrom pracma mldivide mrdivide
#'
#'@references
#' This function is modified version of the "A0N_MLEdensity_WOE" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@keywords internal


A0N_MLEdensity_WOE__jointQ_sepSigma_Bootstrap <- function(K1XQ, r0, SSZ, K0Z, K1Z, se, Gy.0, mat, Y, Z, P, Wpca, We,
                                                          WpcaFull, dt, Economies, FactorLabels, ModelType, residBS,
                                                          MaxEigen, GVARinputs, JLLinputs, nargout){


  # Country-specific inputs
  C <- length(Economies) # Number of economies
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  J <- nrow(Y)/C
  T <- dim(Y)[2]
  M <- (nrow(Z)-G)/C-N # number of country-specific unspanned factors

  # Country-specific inputs
  CJ <- C*J # number of all yields used in estimation;
  CN <- C*N # number of country-specific spanned factors

  Peo <- We%*%Y # portfolio observed WITH errors
  t <- 2:T


  ###########################################################################
  if (!exists("SSZ")){SSZ <- matrix(, nrow=N*N, ncol=0)}
  if (!exists("r0")){ r0 <- as.numeric()}
  if (!exists("se")){ se <- as.numeric()}
  if (!exists("K0Z")){ K0Z <- matrix(, nrow=N, ncol=0)}
  if (!exists("K1Z")){ K1Z <- matrix(, nrow=N*N, ncol=0)}
  #############################################################################
  if (sjmisc::is_empty(K0Z) || sjmisc::is_empty(K1Z) ){
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, K0Z, K1Z )
  } else{
    x <- rbind(t(t(as.vector(K1XQ))), t(t(as.vector(SSZ))), r0, se, t(t(as.vector(K0Z))), t(t(as.vector(K1Z))) )
  }
  # x is a column vector that contains all the elements of each parameter of the model.

  if ( any(is.nan(x)) ||any(is.infinite(x)) || any(Im(x)!= 0) ){
    browser()
    y <- 1e6*rep(1, times = T-1)
    out <- as.numeric()
  }else{

    # Loadings:
    # Yields can be affine function of Pt's,i.e., Y(t) = AnP +BnP*P(t).
    # Further, we define Z(t) as an affine function of P(t) such that: Z(t) = phi0+ phi1*P(t).
    # As such, we can write P(t) = phi1^(-1)*(Z(t) - phi0).
    BnX <- A0N__computeBnAn_jointQ(round(mat/dt),K1XQ, dX= NULL, r0= NULL, SSX= NULL, Economies)[[1]]/dt
    # NOTE: the function "A0N__computeBnAn" generates outputs for interest rates per unit of time interval.
    # Hence by multiplying by "1/dt" we obtain annualized results.
    BnP <- mrdivide(BnX,Wpca%*%BnX) #  BnP = BnX*(W*BnX)^(-1)


    rownames(BnX) <- rownames(Y)


    WBX <- Wpca%*%BnX # WBX = W * BnX
    b<- IdxSpanned(G,M,N,C)
    SSP <- SSZ[b, b]
    SSX <- mrdivide(solve(WBX,SSP, tol = 1e-50),t(WBX)) # SSX = (W*BnX)^(-1)*SSP*t((*W*BnX)^(-1))

    # Optimal estimate of r0:
    betan <- A0N__computeBnAn_jointQ(round(mat/dt), K1XQ, dX=NULL , r0=rep(0, times=C), SSX, Economies)[[3]]
    betan <- t(t(betan/dt))


    if (!exists("r0")|| sjmisc::is_empty(r0) ){
      # (i) A0 = (I - Bx(W*Bx)^(-1)*W)*Ax0;
      # (ii) A1 = (I - Bx(W*Bx)^(-1)*W)*Axr;
      # (iii) APer= We * A1
      # (iv) vt = We*(Yt - Bx(W*Bx)^(-1)*Pt - A0)
      # Recall that BnP= Bx(W*Bx)^(-1)

      r0 <- c()

      idxA <- 0
      idxB <- 0
      idxC <- 0
      for (i in 1:C){
        idxAA <- idxA + N
        idxBB <- idxB + J
        idxCC <- idxC + J-N

        BnPCS <- BnP[(idxB+1):idxBB , (idxA+1):idxAA]
        WpcaCS <- Wpca[(idxA+1):idxAA,(idxB+1):idxBB]
        betanCS <- betan[(idxB+1):idxBB]
        WeCS <- We[(idxC+1):idxCC,(idxB+1):idxBB]
        YCS <- Y[(idxB+1):idxBB,]
        PCS <- P[(idxA+1):idxAA,]

        A0 <- (diag(J)-BnPCS%*%WpcaCS)%*%betanCS
        A1 <- (diag(J)-BnPCS%*%WpcaCS)%*%matrix(1,J,1)/dt

        r0[i] <- solve( (t(A1)%*%(t(WeCS)%*%WeCS)%*%A1), # Numerator of r0
                        (t(A1)%*%(t(WeCS)%*%WeCS))%*%(rowMeans(YCS[,t] - BnPCS%*%PCS[,t]) - A0) ) # Denominator of r0
        # i.e: r0 = APer'*vt/(APer'*APer)

        idxA <- idxAA
        idxB <- idxBB
        idxC <- idxCC
      }
    }


    # compute the intercepts:
    for (i in 1:C){
      Qtemp <- t(t(rep(r0[i], times= J)/dt))
      if (i==1){
        Q <- Qtemp
      }else{
        Q <- rbind(Q, Qtemp)
      }
    }

    AnX <- betan + Q
    AnP <- AnX -BnP%*%(Wpca%*%AnX) # AnP = AnX - BnX*(W*BnX)^(-1)*W*AnX

    rownames(AnX) <- rownames(Y)

    # density of yields pricing errors:
    Pe <- matrix(0, nrow= nrow(We) , ncol= T)
    for( f in 1:T){
      Pe[,f] <- We%*% (BnP%*%P[,f] + AnP) # Pe_t = We*(AP+BP*P_t): "model-implied" yields for the portfolio of yields observed with errors
    }
    eQ <- Peo[,t] - Pe[,t]

    if (!exists("se")|| sjmisc::is_empty(se) ){
      se<- rep(NA, times= C)
      idx0 <-0
      for (i in 1:C){
        idx1 <- idx0+ J-N
        se[i] <- sqrt(mean(as.vector(eQ[(idx0+1):idx1,]^2))) # Scalar - standard deviation of the measurament error.
        idx0 <- idx1
      }
    }


    ## The log-likelihood function:
    if (any(is.nan(se) ==1)){
      y <- 1e6*matrix(1,T-1,1)
    }else{
      # Cross-sectional density (i.e. density for the portfolios observed with measurament error)
      aa <- se
      idx0 <- 0
      for (h in 1:C){
        idx1 <- idx0+ J-N
        se[(idx0+1):idx1] <- rep(aa[h], times=J-N)
        idx0 <- idx1
      }

      y <- GaussianDensity(eQ, se^2*diag(CJ-CN))

      se <- se[seq(1, CJ-CN, by=J-N)] # Recast se


      # time series density:

      if (!exists("K0Z")|| sjmisc::is_empty(K0Z)){
      JLLinputs$WishSigmas <- 1
      JLLinputs$SigmaNonOrtho <- SSZ

      JLLfac <- JLL(Z, N, JLLinputs)
      K0Z <- JLLfac$k0
      K1Z <- JLLfac$k1
      Gy.0 <- diag(C*(M+N) + G)
      }


      eP <- matrix(0, nrow(K0Z), T-1)
      for (f in 1:length(t)){ eP[,f] <-  Z[,f+1] - K1Z%*%Z[,f] - K0Z }

      y <- y + GaussianDensity(eP, SSZ) # Cross-sectional density + time-series density (final likelihood function, except for the Jacobian Term

      # Jacobian: (so that the density is for the observed yields (and the non-yields variables Z)
      #           and not the yields portfolios which depend on W).

      y <-  y + 0.5*log(abs(det(Wpca%*%t(Wpca))))
    }

    if ( any(is.nan(y)) ||any(is.infinite(y)) || any(Im(y)!= 0) ){
      y <- 1e6*rep(1, times = T-1)
    }else{
      y = -t(as.vector(y)) # for minimization;
    }


    # Variance of yields
    Sigma_e_CS <- c()
    Sigma_e_AllRows <- c()

    for (i in 1:C){
      Sigma_e_CS  <- c(rep(0, times = N), rep(se[i],times= J-N)^2) # Variances of portfolio WITHOUT and WITH errors
      Sigma_e_AllRows <- append(Sigma_e_AllRows, Sigma_e_CS)
    }

    Sigma_e_AllCountries <- Sigma_e_AllRows*diag(CJ)
    SIGMA_yields <- solve(WpcaFull)%*%Sigma_e_AllCountries%*%t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)

    names(VarYields) <- rownames(Y)
    VarYields <- t(t(VarYields)) # make a column vector


    if (nargout>1){

          #Estimates:


        JLLinputs$WishSigmas <- 1
        JLLinputs$SigmaNonOrtho <- SSZ
        JLLoutcomesOrtho <- JLL(Z, N, JLLinputs)

        # Remove the non-orthogonalized outputs
        JLLoutcomesOrtho$k0 <- NULL
        JLLoutcomesOrtho$k1 <- NULL
        JLLoutcomesOrtho$a_W <- NULL
        JLLoutcomesOrtho$b <- NULL
        JLLoutcomesOrtho$c <- NULL
        JLLoutcomesOrtho$a_DU_CS <- NULL
        JLLoutcomesOrtho$PIb <- NULL
        JLLoutcomesOrtho$PIac <- NULL
        JLLoutcomesOrtho$Ye <- NULL
        #Add the orthogonalized Variance-covariance matrix
        ZeroIdxVarCov <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$VarCovOrtho
        ZeroIdxSigma_Ye <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$Sigma_Ye
        JLLoutcomesOrtho$Sigmas$VarCov_Ortho[ZeroIdxVarCov] <- 0
        JLLoutcomesOrtho$Sigmas$Sigma_Ye[ZeroIdxSigma_Ye] <- 0

        PI <- JLLoutcomesOrtho$PI
        JLLoutcomesOrtho$Sigmas$Sigma_Y <- PI%*%JLLoutcomesOrtho$Sigmas$Sigma_Ye

        ests <- list(K1XQ, SSZ, SSP, r0, se, K0Z, K1Z, Gy.0, VarYields, residBS, JLLoutcomesOrtho)

      llk <- mean(y) # [1 x T-p] actual densities (not negative of densities)


      # Q dynamics: X as risk factors:
      Q <- list(zeros(CN,1), K1XQ, SSX)
      X <- list(BnX, AnX, Q)
      rot <- list(X)

      Out <- list(ests, llk, rot)

      names(Out) <- c("ests", "llkMax", "rot")

      names(Out$ests) <- c('K1XQ', 'SSZ', 'SSP', 'r0', 'se', 'K0Z', 'K1Z', 'Gy.0', 'VarYields', 'residBS', 'JLLoutcomes')
      names(Out$rot) <- c("X")
      names(Out$rot$X) <-  c("B","A","Q")
      names(Out$rot$X$Q) <- c("K0", "K1", "SS" )


      # Q dynamics: PCN as risk factors: % PCN = U0 + U1*X, where U0 = W*AnX and U1 = W*BnX
      U1 <- Wpca%*%BnX
      U0 <- Wpca%*%AnX
      out.rot.P <-  FMN__Rotate(Out$rot$X, U1, U0)

      # updating the generated outputs...
      Out$rot[["P"]] <- out.rot.P


      Out$rot$X <- NULL # Delete the parameters related to the latent factors
      #(Not needded for the Bootstrap, save some time and pc memory)

    }

  }

  if (nargout ==1){
    return(y)
  } else{
    return(Out)
  }


}

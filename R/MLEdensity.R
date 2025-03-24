#' Compute the maximum likelihood function of all models
#'
#' @param K1XQ risk-neutral feedback matrix (N x N or CN x CN)
#' @param r0 long-run interest rate (scalar or vector with length C)
#' @param SSZ variance-covariance matrix (F x F)
#' @param K0Z intercept from the P-dynamics (F x 1)
#' @param K1Z feedback matrix from the P-dynamics (F x F)
#' @param se Variance of the portfolio of yields observed with error (scalar). Default is set to NULL.
#' @param Gy.0 matrix of contemporaneous terms from the P-dynamics (F x F)
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param Y matrix of yields used in estimation (J x T or CJ x T)
#' @param Z complete set of spanned and unspanned factors (F x T)
#' @param P complete set of spanned factors (N x T or CN x T)
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#' @param We matrix of weights of the portfolios observed with errors ((J-N) x J or C(J-N) x CJ)
#' @param WpcaFull composite matrix of weights the portfolios observed with and without errors
#' @param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param GVARinputs if the model chosen is the "GVAR single" or "GVAR multi", the "GVARinputs" should be specified (see "GVAR" function)
#' @param JLLinputs if the model chosen is JLL-based. "JLLinputs" should contain (i) DomUnit, (ii) WishSigmas, (iii) SigmaNonOrtho, (iv) JLLModelType (See JLL function)
#' @param BS_outputs Generates simplified output list in the bootstrap setting. Default is set to FALSE.
#' @param nargout if nargout == 1: provides only the values of the likelihood; if nargout == 2: complete ATSM outputs
#'
#' @references
#' This function is modified version of the "A0N_MLEdensity_WOE" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#' @keywords internal

MLEdensity <- function(K1XQ, r0, SSZ, K0Z, K1Z, se, Gy.0, mat, Y, Z, P, Wpca, We, WpcaFull, dt, Economies,
                       FactorLabels, ModelType, GVARinputs = NULL, JLLinputs = NULL, BS_outputs = FALSE, nargout) {

  # 0) Initialize some variables, if necessary
  if (!exists("r0")) r0 <- numeric()
  if (!exists("se")) se <- numeric()

  N <- length(FactorLabels$Spanned)
  if (ModelType == 'JLL joint Sigma') SSZ <- Update_SSZ_JLL(SSZ, Z, N, JLLinputs)

  # x is a column vector that contains all model parameters.
  x <- rbind(matrix(K1XQ), matrix(SSZ), r0, se, matrix(K0Z), matrix(K1Z))

  # Test for the optimization feasibility
  if (any(is.nan(x)) || any(is.infinite(x)) || any(Im(x) != 0)) {
    T <- ncol(Z)
    y <- 1e6 * rep(1, times = T - 1)
    out <- list()
  } else {

  # 1) Compute loadings:
  # Slope coefficients (Bn's)
    LoadBs <- Get_Bs(mat, dt, K1XQ, SSZ, Wpca, FactorLabels, Economies, ModelType)
  # Intercepts (An's):
    if (is.null(r0) || length(r0) == 0) r0 <- Get_r0(Y, P, N, mat, dt, LoadBs, Wpca, We, Economies, ModelType)
    LoadAs <- Get_As(LoadBs, Wpca, r0, dt, Economies, ModelType)

  # 2) Build Log-likelihood density function
    llkOut <- Get_llk(P, Y, Z, N, mat, We, Wpca, K0Z, K1Z, SSZ, LoadBs, LoadAs, ModelType)
    y <- llkOut$y
    se <- llkOut$se

  # 3) Variance-yields
    VarYields <- Get_SigmaYields(Y, N, mat, WpcaFull, se, ModelType)

  # 4) Output to export
    if (nargout > 1) {
      out <- OptOutputs(Y, Z, mat, N, dt, Wpca, K1XQ, SSZ, LoadAs, LoadBs, r0, se, K0Z, K1Z, Gy.0, VarYields,
                        y, GVARinputs, JLLinputs, Economies, ModelType, BS_outputs)
    }
  }

  if (nargout == 1) {
    return(y)
  } else {
    return(out)
  }
}
##############################################################################################################
#' Extract the indexes related to the spanned factors in the variance-covariance matrix
#'
#' @param G number of global unspanned factors (scalar)
#' @param M number of domestic unspanned factors per country (scalar)
#' @param N number of domestic spanned factors per country (scalar)
#' @param C number of countries of the economic system (scalar)
#'
#' @keywords internal

IdxSpanned <- function(G, M, N, C) {

  K <- C * (N + M) + G
  vector <- 1:K

  idxA <- seq(from = G + M, by = N + M, length.out = C)
  IDX_F <- unlist(lapply(1:C, function(i) vector[(idxA[i] + 1):(idxA[i] + N)]))

  return(IDX_F)
}
##############################################################################################################
#' Compute r0 for the various models
#'
#' @param Y matrix of yields used in estimation (J x T or CJ x T)
#' @param P complete set of spanned factors (N x T or CN x T)
#' @param N number of country-specific spanned factors
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#' @param B_list list containing the B loadings
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#' @param We matrix of weights of the portfolios observed with errors ((J-N) x J or C(J-N) x CJ)
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

Get_r0 <- function(Y, P, N, mat, dt, B_list, Wpca, We, Economies, ModelType){

  # General procedure:
  # (i) A0 = (I - Bx(W*Bx)^(-1)*W)*Ax0;
  # (ii) A1 = (I - Bx(W*Bx)^(-1)*W)*Axr;
  # (iii) APer= We * A1
  # (iv) vt = We*(Yt - Bx(W*Bx)^(-1)*Pt - A0)
  # Recall that BnP= Bx(W*Bx)^(-1)

  J <- length(mat) # number of country-specific yields used in estimation;
  T <- dim(Y)[2]
  tt <- 2:T

  BnP <- B_list$BnP
  betan <- B_list$betan

  # 1) Get r0 for the models estimated on a country-by-country basis
  if (any(ModelType == c("JPS original", "JPS global", 'GVAR single'))) {
    A0 <- (diag(J) - BnP %*% Wpca) %*% betan
    A1 <- (diag(J) - BnP %*% Wpca) %*% matrix(1, J, 1) / dt

    #  r0 = APer'*vt/(APer'*APer)
    A_per <- t(A1) %*% (t(We) %*% We)
    r0 <- solve(A_per %*% A1, # Numerator from equation r0 equation
                A_per %*% (rowMeans(Y[, tt] - BnP %*% P[, tt]) - A0)) # Denominator from r0 equation
  } else {

  # 2) Get r0 for the models estimated on a jointly basis
    C <- length(Economies)
    r0 <- numeric(C)

    idxA <- 0
    idxB <- 0
    idxC <- 0

    for (i in 1:C) {
      idxAA <- idxA + N
      idxBB <- idxB + J
      idxCC <- idxC + J - N

      BnPCS <- matrix(BnP[(idxB + 1):idxBB, (idxA + 1):idxAA], nrow = J)
      WpcaCS <- Wpca[(idxA + 1):idxAA, (idxB + 1):idxBB]
      betanCS <- betan[(idxB + 1):idxBB]
      WeCS <- We[(idxC + 1):idxCC, (idxB + 1):idxBB]
      YCS <- Y[(idxB + 1):idxBB, ]
      PCS <- P[(idxA + 1):idxAA, ]

      A0 <- (diag(J) - BnPCS %*% WpcaCS) %*% betanCS
      A1 <- (diag(J) - BnPCS %*% WpcaCS) %*% matrix(1, J, 1) / dt

      # r0 = APer'*vt/(APer'*APer)
      A_per <- t(A1) %*% (t(WeCS) %*% WeCS)
      PCS_adj <- if (N == 1) PCS[tt] else PCS[, tt]

      r0[i] <- solve(A_per %*% A1, # Numerator from r0 equation
                     A_per %*% (rowMeans(YCS[, tt] - BnPCS %*% PCS_adj) - A0))  # Denominator from r0 equation

      idxA <- idxAA
      idxB <- idxBB
      idxC <- idxCC
    }
  }

  return(r0)
}
##############################################################################################################
#' Build the B loadings
#'
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1
#' @param K1XQ risk-neutral feedback matrix (N x N or CN x CN)
#' @param SSZ variance-covariance matrix (F x F)
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economy string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

Get_Bs <- function(mat, dt, K1XQ, SSZ, Wpca, FactorLabels, Economy, ModelType) {

  Lab_SingleQ <- c("JPS original", "JPS global", "GVAR single")

  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N

  # 1) Get BnX and BnP
  # Yields are affine function of the vector Pt,i.e., Y(t) = AnP +BnP*P(t).
  # Further, we define Z(t) as an affine function of P(t) such that: Z(t) = phi0+ phi1*P(t).
  # As such, we can write P(t) = phi1^(-1)*(Z(t) - phi0).
  BnX <- A0N__BnAn(round(mat / dt), K1XQ, ModelType, Economies = Economy)$BnX / dt
  # We multiply by "1/dt" to produce annualized results.

  if (any(ModelType == Lab_SingleQ)) {
    if (N==1){ BnP <- matrix(BnX) %*% solve(Wpca %*% BnX)}
    else{ BnP <- BnX %*% solve(Wpca %*% BnX) }
  } else {
    BnP <- BnX %*% solve(Wpca %*% BnX)
  }

  # 2) Covariance matrix of latent states X:
  if (any(ModelType == Lab_SingleQ)) {
    AllLabels <- GetLabels_sepQ(Economy, ModelType, FactorLabels)
    dimnames(SSZ) <- list(AllLabels, AllLabels)
    LabelSpannedCS <- c(FactorLabels$Tables[[Economy]][-(1:M)])
    idxSpa <- match(LabelSpannedCS, AllLabels)
  } else {
    G <- length(FactorLabels$Global)
    C <- length(Economy)
    idxSpa <- IdxSpanned(G, M, N, C)
  }

  WBX <- Wpca %*% BnX # WBX = W * BnX
  SSP <- SSZ[idxSpa, idxSpa]

  if (is.null(dim(SSP))) {
    SSX <- matrix(Inf, nrow = 1, ncol = 1) # Initialize with Inf
  } else {
    SSX <- matrix(Inf, nrow = nrow(SSP), ncol = ncol(SSP)) # Try to assign the correct value
  }

  SSX <- solve(WBX, SSP, tol = 1e-50) %*% solve(t(WBX))

  # 3) Optimal estimate of r0. NOTE: We set r0=0, because when r0=0, An = Betan.
  if (any(ModelType == Lab_SingleQ)) {
    r0 <- 0
  } else {
    r0 <- rep(0, times = C)
  }

  betan <- A0N__BnAn(round(mat / dt), K1XQ, ModelType, r0 = r0, SSX = SSX, Economies = Economy)$betan
  betan <- t(t(betan / dt))

  out <- list(BnP = BnP, BnX = BnX, betan = betan, SSP = SSP, SSX = SSX, idxSpa = idxSpa)

  return(out)
}
##############################################################################################################
#' Compute the A loadings
#'
#' @param LoadBs list containing the B loadings
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#' @param r0 long-run interest rate (scalar or vector with length C)
#' @param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

Get_As <- function(LoadBs, Wpca, r0, dt, Economies, ModelType) {

  betan <- LoadBs$betan
  BnP <- LoadBs$BnP

  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    AnX <- betan + as.numeric(r0 / dt)
  } else {
    C <- length(Economies)
    J <- nrow(betan) / C

    Q <- do.call(rbind, lapply(1:C, function(i) matrix(rep(r0[i], times = J) / dt, nrow = J, ncol = 1)))
    AnX <- betan + Q
  }

  AnP <- AnX - BnP %*% (Wpca %*% AnX) # AnP = AnX - BnX*(W*BnX)^(-1)*W*AnX

  return(list(AnX = AnX, AnP = AnP))
}
##############################################################################################################
#' Generate the factor labels for models estimated on a country-by-country basis
#'
#' @param Economy string containing the names of the economy to be estimated
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param FactorLabels list containing the factor labels
#'
#' @keywords internal

GetLabels_sepQ <- function(Economy, ModelType, FactorLabels) {

  if (ModelType == "JPS original") {
    AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]])
  } else if (any(ModelType == c("JPS global", 'GVAR single'))) {
    AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  }

  return(AllLabels)
}
##############################################################################################################
#' Compute the log-likelihood function
#'
#' @param P time-series of spanned factors (N x T or CN x T)
#' @param Y time-series of yields (J x T or CJ x T)
#' @param Z time-series of risk factors (F x T)
#' @param N number of country-specific spanned factors
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param We matrix of weights of the portfolios observed with errors ((J-N) x J or C(J-N) x CJ)
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x CJ)
#' @param K0Z matrix of intercepts (P-dynamics)
#' @param K1Z feedback matrix (P-dynamics)
#' @param SSZ variance-covariance matrix (P-dynamics)
#' @param LoadBs list containing the B loadings
#' @param LoadAs list containing the A loadings
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

Get_llk <- function(P, Y, Z, N, mat, We, Wpca, K0Z, K1Z, SSZ, LoadBs, LoadAs, ModelType) {

  T <- dim(P)[2]
  t <- 2:T
  J <- length(mat)

  Peo <- We %*% Y # portfolio observed WITH errors
  BnP <- LoadBs$BnP
  AnP <- LoadAs$AnP

  # 1) density of yields pricing errors:
  MatOnes <- matrix(1, ncol = ncol(P), nrow = 1)
  Pe <- We %*% (BnP %*% P + AnP%*%MatOnes)

  eQ <- Peo[, t] - Pe[, t]

  # 2) Standard deviation of the measurement error.
  if (!exists("se") || is.null(se) || length(se) == 0) {
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      se <- sqrt(mean(as.vector(eQ^2)))
    } else {
      C <- nrow(P) / N
      se <- numeric(C)
      idx0 <- 0
      for (i in 1:C) {
        idx1 <- idx0 + J - N
        se[i] <- sqrt(mean(as.vector(eQ[(idx0 + 1):idx1, ]^2))) # Standard deviation of the measurement error (scalar).
        idx0 <- idx1
      }
    }
  }

  # 3) The log-likelihood function:
  if (any(is.nan(se))) {
    y <- 1e6 * matrix(1, T - 1, 1)
  } else {

    # Cross-sectional density (i.e. density for the portfolios observed with measurement error)
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      y <- GaussianDensity(eQ, se^2 * diag(J - N))
    } else {
      aa <- se
      idx0 <- 0
      for (h in 1:C) {
        idx1 <- idx0 + J - N
        se[(idx0 + 1):idx1] <- rep(aa[h], times = J - N)
        idx0 <- idx1
      }

      se <- se[seq(1, C * J - C * N, by = J - N)] # Recast se
      y <- GaussianDensity(eQ, se^2 * diag(C * J - C * N))
    }

    # Time series density:
    MatOne_Z <- matrix(1, nrow = 1, ncol = ncol(Z) - 1)
    eP <- Z[, -1] - K1Z %*% Z[ , 1:(T - 1)] - K0Z%*%MatOne_Z

    y <- y + GaussianDensity(eP, SSZ) # Cross-sectional density + time-series density (final likelihood function, except for the Jacobian term)

    # Jacobian: (so that the density is for the observed yields (and the non-yields variables Z)
    #           and not the yields portfolios which depend on W).
    y <- y + 0.5 * log(abs(det(Wpca %*% t(Wpca))))
  }

  if (any(is.nan(y)) || any(is.infinite(y)) || any(Im(y) != 0)) {
    y <- 1e6 * rep(1, times = T - 1)
  } else {
    y <- -t(as.vector(y)) # necessary for minimization
  }

  return(list(y = y, se = se))
}
##############################################################################################################
#' computes the density function of a gaussian process
#'
#' @param res matrix of residuals (N x T)
#' @param SS covariance matrice or array of covariance matrices (N x N) or (N x N x T)
#' @param invSS Inverse of SS (N x N) or (N x N x T) - optional input
#' @param logabsdetSS log(abs(|SS|)) (1 x T) - optional input
#'
#' @keywords internal
#'
#' @return y vector of density (1 x T)
#'
#' @references
#' This function is based on the "Gaussian" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029

GaussianDensity <- function(res, SS, invSS, logabsdetSS) {

  N <- dim(res)[1]
  T <- dim(res)[2]

  if (!missing("invSS")) {
    SSres <- invSS %*% res
  } else {
    SSres <- solve(SS, res, tol = 1e-50)
  }

  if (missing("logabsdetSS")) {
    LodDet <- 2 * sum(log(abs(diag(chol(SS %*% t(SS))))))
    logabsdetSS <- 0.5 * LodDet
  }

  y <- -0.5 * N * log(2 * pi) - 0.5 * logabsdetSS - 0.5 * abs(colSums(res * SSres))

  return(y)
}
##############################################################################################################
#' Compute the variance-covariance matrix of the bond yields
#'
#' @param YieldsTS matrix of yields used in estimation (J x T or CJ x T)
#' @param N number of country-specific spanned factors
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param WpcaFull composite matrix of weights the portfolios observed with and without errors
#' @param se Variance of the portfolio of yields observed with error (scalar). Default is set to NULL
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

Get_SigmaYields <- function(YieldsTS, N, mat, WpcaFull, se, ModelType) {

  J <- length(mat)

  # Single-country models
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    Sigma_e_Row <- c(rep(0, times = N), rep(se, times = J - N)^2) # Variances of portfolio WITHOUT and WITH errors
    Sigma_e <- Sigma_e_Row * diag(J)
    SIGMA_yields <- solve(WpcaFull) %*% Sigma_e %*% t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)
  } else {

    # Multicountry models
    C <- length(se)

    Sigma_e_AllRows <- rep(0, times = N * C)  # Pre-allocate the vector with N*C zeros
    for (i in 1:C) {
      Sigma_e_CS <- c(rep(0, N), rep(se[i]^2, J - N))  # Variances of portfolio WITHOUT and WITH errors
      Sigma_e_AllRows[((i - 1) * J + 1):(i * J)] <- Sigma_e_CS  # Update the correct segment
    }

    Sigma_e_AllCountries <- diag(C * J) * Sigma_e_AllRows
    SIGMA_yields <- solve(WpcaFull) %*% Sigma_e_AllCountries %*% t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)

    names(VarYields) <- rownames(YieldsTS)
    VarYields <- t(t(VarYields))
  }

  return(VarYields)
}
##############################################################################################################
#' Update the variance-covariance matrix from the "JLL joint Sigma" model. Necessary for optimization
#'
#' @param SSZ Variance-covariance matrix from JLL model
#' @param Z complete set of spanned and unspanned factors (F x T)
#' @param N number of country-specific spanned factors
#' @param JLLinputs List of inputs from JLL models
#'
#' @keywords internal

Update_SSZ_JLL <- function(SSZ, Z, N, JLLinputs) {

  Sigma_Ye <- t(chol(SSZ)) # in this case, SSZ is the Variance-covariance matrix of the ORTHOGONALIZED dynamics
  JLLinputs$WishSigmas <- 0 # Avoid recomputing the variance-covariance matrices
  JLLPara <- JLL(Z, N, JLLinputs)
  Sigma_Y <- JLLPara$PI %*% Sigma_Ye
  SSZ <- Sigma_Y %*% t(Sigma_Y) # Variance-covariance matrix non-orthogonalized

  return(SSZ)
}
##############################################################################################################
#' Prepare outputs to export after the model optimization
#'
#' @param Y matrix of yields used in estimation (J x T or CJ x T)
#' @param Z complete set of spanned and unspanned factors (F x T)
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param N number of country-specific spanned factors
#' @param dt time interval unit of the model (scalar)
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#' @param K1XQ risk-neutral feedback matrix (N x N or CN x CN)
#' @param SSZ variance-covariance matrix (F x F)
#' @param LoadAs list containing the A loadings
#' @param LoadBs list containing the B loadings
#' @param r0 long-run interest rate (scalar or vector with length C)
#' @param se Variance of the portfolio of yields observed with error (scalar).
#' @param K0Z intercept from the P-dynamics (F x 1)
#' @param K1Z feedback matrix from the P-dynamics (F x F)
#' @param Gy.0 matrix of contemporaneous terms from the P-dynamics (F x F)
#' @param VarYields variance-covariance matrix of the bond yields
#' @param y likelihood of each time series (Tx1)
#' @param GVARinputs List of inputs from GVAR models
#' @param JLLinputs List of inputs from JLL models
#' @param Economies string containing the names of the economy to be estimated
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param BS_out Bootstrap output. Default is FALSE.
#'
#' @keywords internal

OptOutputs <- function(Y, Z, mat, N, dt, Wpca, K1XQ, SSZ, LoadAs, LoadBs, r0, se, K0Z, K1Z, Gy.0, VarYields,
                       y, GVARinputs, JLLinputs, Economies, ModelType, BS_out = FALSE) {

  # a) Build LIST 1
  # a.1) List "input"
  inputs <- list(Y = Y, AllFactors = Z, mat = mat, N = N, dt = dt, Wpca = Wpca)
  if (ModelType  %in% c("GVAR single", "GVAR multi")) {
    inputs$Wgvar <- GVARinputs$Wgvar
  } else if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    inputs$DomUnit <- JLLinputs$DomUnit
  }

  # a.2) List "ests"
  ests <- list(K1XQ = K1XQ, SSZ = SSZ, SSP = LoadBs$SSP, r0 = r0, se = se, K0Z = K0Z, K1Z = K1Z,
               Gy.0 = Gy.0, VarYields = VarYields)

  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    JLLinputs$WishSigmas <- 1
    JLLinputs$SigmaNonOrtho <- SSZ
    JLLoutcomesOrtho <- JLL(Z, N, JLLinputs)
    # Remove the non-orthogonalized outputs
    JLLoutcomesOrtho$k0 <- NULL
    JLLoutcomesOrtho$k1 <- NULL
    if (BS_out) {
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
    JLLoutcomesOrtho$Sigmas$Sigma_Y <- PI %*% JLLoutcomesOrtho$Sigmas$Sigma_Ye

    ests$JLLoutcomes <- JLLoutcomesOrtho
  }

  # a.3) List "llk"
  llk <- list(-t(y))

  # a.4) List Q dynamics: X as risk factors:
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    Q <- list(K0 = matrix(0, N, 1), K1 = K1XQ, SS = LoadBs$SSX)
  } else {
    Q <- list(K0 = matrix(0, N * length(Economies), 1), K1 = K1XQ, SS = LoadBs$SSX)
  }
  X <- list(B = LoadBs$BnX, A = LoadAs$AnX, Q = Q)
  rot <- list(X = X)

  # Summary: List 1
  Out <- list(inputs = inputs, ests = ests, llk = llk, rot = rot)

  # b) Build LIST 2
  # Q dynamics: PCN as risk factors: % PCN = U0 + U1*X, where U0 = W*AnX and U1 = W*BnX
  U1 <- Wpca %*% LoadBs$BnX
  U0 <- Wpca %*% LoadAs$AnX
  Out$rot$P <- FMN__Rotate(Out$rot$X, U1, U0)

  if (BS_out) {
    Out$rot$X <- NULL # Delete the parameters related to the latent factors
    #(Not needed for the Bootstrap, save some time and pc memory)
  } else {
    # P-dynamics: Z as risk factors:
    Out$rot$P$P <- list(K0 = K0Z, K1 = K1Z, SS = SSZ)

    # P-dynamics: PCN as risk factors:
    idxSpanned <- LoadBs$idxSpa
    K0P <- K0Z[idxSpanned]
    K1P <- K1Z[idxSpanned, idxSpanned]
    SSP <- SSZ[idxSpanned, idxSpanned]
    out_rot_P_P2 <- list(K0 = K0P, K1 = K1P, SS = SSP)

    Out$rot$X$P <- FMN__Rotate(out_rot_P_P2, solve(U1, tol = 1e-50), -solve(U1, U0, tol = 1e-50))

    row.names(Out$rot$P$A) <- row.names(Y)
    row.names(Out$rot$P$B) <- row.names(Y)
    colnames(Out$rot$P$B) <- row.names(Z)[idxSpanned]

    # PC as risk factors:
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      PC1NW <- pca_weights_one_country(Y, Economies)
      PC1NW <- PC1NW[1:N, ] * 100
      U1 <- PC1NW %*% LoadBs$BnX
      U0 <- PC1NW %*% LoadAs$AnX
    } else {
      U1 <- Wpca %*% LoadBs$BnX
      U0 <- Wpca %*% LoadAs$AnX
    }

    Out$rot$PC <- FMN__Rotate(Out$rot$X, U1, U0)

    # P-dynamics: PC as risk factors:
    Out$rot$PC$P <- list(K0 = K0Z, K1 = K1Z, SS = SSZ)
  }

  return(Out)
}

##############################################################################################################
#' Efficient computation of matrix product for arrays
#'
#' @param a array
#' @param b array
#'
#' @keywords internal
#'
#' @references
#' This function is a simplified version of the "mult__prod" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."\cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029

mult__prod <- function(a, b) {
  d <- a %*% b
  return(d)
}


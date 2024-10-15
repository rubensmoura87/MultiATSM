#' Estimates an unbiased VAR(1) using stochastic approximation (Bauer, Rudebusch and Wu, 2012)
#'
#'
#'@param ModelType A character vector indicating the model type to be estimated.
#'@param BRWinputs A list containing the necessary inputs for the BRW model estimation:
#'\enumerate{
#'\item \code{flag_mean}: Logical. Determines whether mean- (TRUE) or median- (FALSE) unbiased estimation is desired. Default is TRUE.
#'\item \code{gamma}: Numeric. Adjustment parameter between 0 and 1. Default is 0.5.
#'\item \code{N_iter}: Integer. Number of iterations for the stochastic approximation algorithm after burn-in. Default is 5,000.
#'\item \code{N_burn}: Integer. Number of burn-in iterations. Default is 15% of \code{N_iter}.
#'\item \code{B}: Integer. Number of bootstrap samples per iteration for calculating the noisy measure of the OLS estimator's mean or median. Default is 50.
#'\item \code{check}: Logical. Indicates whether to perform a closeness check. Default is TRUE.
#'\item \code{B_check}: Integer. Number of bootstrap samples for the closeness check. Default is 100,000.
#'}
#'@param RiskFactors A numeric matrix (T x F) representing the time series of risk factors.
#'@param N Integer. Number of country-specific spanned factors.
#'@param Economies A character vector containing the names of the economies included in the system.
#'@param FactorLabels A list of character vectors with labels for all variables in the model.
#'@param GVARinputs List. Inputs for GVAR model estimation (see \code{\link{GVAR}} function). Default is NULL.
#'@param JLLinputs List. Inputs for JLL model estimation (see \code{\link{JLL}} function). Default is NULL.
#'@param ev_restr Numeric. Restriction on the largest eigenvalue under the P-measure. Default is 1.
#'@param nargout Integer. Number of elements in the output list. Default is 4.
#'
#'
#'
#'@examples
#'\donttest{
#' data(CM_Factors)
#' Factors <- t(RiskFactors[1:7,])
#'
#' BRWinputs <- list(flag_mean = TRUE, gamma = 0.4, N_iter = 1000, N_burn = 100,
#'                   B = 10, check = 1, B_check = 5000)
#'
#' Economies <- "China"
#' N <- 3
#' ModelType <- "JPS original"
#' FactorLabels <- NULL
#'
#'
#' BRWpara <- Bias_Correc_VAR(ModelType, BRWinputs, Factors, N, Economies, FactorLabels)
#'}
#'
#'@return Bias-corrected VAR paramaters based on the framework of Bauer, Rudebusch and Wu (2012). The list contains:
#'\enumerate{
#'\item \code{Phi_tilde}: estimated coefficient matrix (F x F);
#'\item \code{mu_tilde}: estimated intercept (F x 1);
#'\item \code{V_tilde}: estimated variance-covariance matrix (F x F);
#'\item \code{dist}: root mean square distance (scalar);
#'\item \code{Phi_sample}: sample estimated variance-covariance matrix used in the checks (F x F x B_check) - this output is
#'      reported if nargout is 5.
#'}
#'
#'
#'@references
#' Bauer, Rudebusch and, Wu (2012). "Correcting Estimation Bias in Dynamic Term Structure Models" \cr
#' This function is based on the \code{est_unb_var} Matlab function available at Cynthia Wu's website
#' (https://sites.google.com/view/jingcynthiawu/).
#'@export



Bias_Correc_VAR  <- function(ModelType, BRWinputs, RiskFactors, N, Economies, FactorLabels, GVARinputs = NULL,
                             JLLinputs = NULL, ev_restr = 1, nargout= 4){

  flag_mean <- ifelse(is.null(BRWinputs$flag_mean), TRUE, BRWinputs$flag_mean)
  gamma <- ifelse(is.null(BRWinputs$gamma), 0.5, BRWinputs$gamma)
  N_iter <- ifelse(is.null(BRWinputs$N_iter), 5000, BRWinputs$N_iter)
  N_burn <- ifelse(is.null(BRWinputs$N_burn), N_iter * 0.15, BRWinputs$N_burn)
  B <- ifelse(is.null(BRWinputs$B), 50, BRWinputs$B)
  check <- ifelse(is.null(BRWinputs$check), FALSE, BRWinputs$check)
  B_check <- ifelse(is.null(BRWinputs$B_check), 100000, BRWinputs$B_check)
  checkSigma <- ifelse(is.null(BRWinputs$checkSigma), TRUE, BRWinputs$checkSigma)

  T <- dim(RiskFactors)[1]
  k <- dim(RiskFactors)[2]

  # 1) OLS
  Phi_hat <- estVARbrw(RiskFactors, ModelType, N, GVARinputs, JLLinputs, FactorLabels, Economies,
                       demean = TRUE, intercept = FALSE)$Gamma_hat

  # 2) initialization for SA algorithm
  theta <- matrix(0, nrow = k^2, ncol=N_burn+N_iter)
  theta_hat <- Phi_hat[1:(ncol(Phi_hat)*nrow(Phi_hat))]
  theta[ ,1] <- theta_hat # starting value

  # 3) SA algorithm
  for (j in 1:(N_burn+N_iter-1)){
    Phi_new <- m_var(theta[ ,j], B, RiskFactors, N, GVARinputs, JLLinputs, FactorLabels, Economies,
                     ModelType, flag_mean)$Phi_new
    theta_new <- c(Phi_new)
    d <- theta_hat - theta_new
    theta[ ,j+1] <- theta[ ,j] + gamma*d

    if ( (j>N_burn) & ( ( (j-N_burn) %% 100) == 0) ){
      # print some diagnostics
      theta_tilde <- rowMeans(theta[ ,(N_burn+1):j])
      Phi_tilde <- matrix(theta_tilde,k,k)
      cat(paste('  *** Iteration:', j-N_burn, '. Largest absolute eigenvalue: ',
                round(max(abs(eigen(Phi_tilde)$value)), 6),"\n"))
    }
  }

  theta_tilde <- rowMeans(theta[ ,(N_burn+1):(N_burn+N_iter)])

  # 4) check whether mean/median of OLS is close to actual OLS estimates
  if (check){
    dist <- Check_comparison__OLS(theta_tilde, theta_hat, B_check, RiskFactors, N, GVARinputs, JLLinputs,
                                  FactorLabels, Economies, ModelType)
  }else{
    Phi_sample <- NaN
  }

  # 5) bias_adjusted estimates
  Phi_tilde <- matrix(theta_tilde,k,k)

  # 6) impose restriction on eigenvalue
  Phi_tilde <- shrink_Phi(Phi_tilde, Phi_hat, ev_restr)

  # 7) choose intercept to match sample mean
  mu_tilde <- (diag(k) - Phi_tilde)%*% t(t(colMeans(RiskFactors)))

  # 8) residuals and their variance-covariance matrix
  if(checkSigma){
  V_tilde <- Get_V_tilde_BC(Phi_tilde, N, RiskFactors, GVARinputs, JLLinputs, FactorLabels, ModelType)
  }else{ V_tilde <- NaN}

  # 9) Compilled desired outputs
  if (nargout==5){
    output <- list( Phi_tilde = Phi_tilde, mu_tilde = mu_tilde, V_tilde = V_tilde, dist = dist, Phi_sample= Phi_sample)
    }else{output <- list(Phi_tilde = Phi_tilde, mu_tilde = mu_tilde, V_tilde = V_tilde, dist = dist) }

  return(output)
}
###############################################################################################################
###############################################################################################################
#' Estimate a VAR(1) - suited to Bauer, Rudebusch and Wu (2012) methodology
#'
#'@param RiskFactors time series of the risk factors (T x F)
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param N number of country-specific spanned factors (scalar)
#'@param GVARinputs inputs used in the estimation of the GVAR-based models (see "GVAR" function)
#'@param JLLinputs inputs used in the estimation of the JLL-based models (see "JLL" function)
#'@param FactorLabels string-list based which contains the labels of all variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param demean demean the data before estimation. Default is set to FALSE
#'@param intercept Include intercept in the VAR model. Default is set to TRUE
#'
#'@keywords internal
#'
#'@return list containing VAR(1) parameters
#'#'\enumerate{
#'\item Gamma_hat: feedback matrix (F X F)
#'\item alpha_hat: intercept (F x 1)
#'}
#'
#'#'@references
#' Bauer, Rudebusch and, Wu (2012). "Correcting Estimation Bias in Dynamic Term Structure Models". \cr
#' This function is similar to the "estVAR" Matlab function available at Cynthia Wu's website
#' (https://sites.google.com/view/jingcynthiawu/).



estVARbrw <- function(RiskFactors, ModelType, N, GVARinputs, JLLinputs, FactorLabels, Economies,
                      demean = FALSE, intercept = TRUE){


  T <- dim(RiskFactors)[1]; k <- dim(RiskFactors)[2]

  # 1) Data prep
  if (demean == TRUE){  RiskFactors <- RiskFactors - matrix(1, T, 1)%*%colMeans(RiskFactors) }

  Yt <- RiskFactors[1:(T-1), ]  # (T-1)*k
  Ytp1 <- RiskFactors[2:T, ]  # (T-1)*k
  Y <- t(Ytp1)  # k*(T-1)

    if (intercept== TRUE){Z <- t(cbind(matrix(1, T-1, 1), Yt)) # (k+1)*(T-1)
    }else{  Z <- t(Yt)} # k*(T-1)

  # 2) Generate feedback matrix of a VAR(1) with demeaned data
    # a) JPS-related models
    if(any(ModelType == c("JPS original", "JPS global", "JPS multi"))){
    A <- Y%*%t(Z)%*%solve(Z%*%t(Z), tol = 1e-50) # k*(k+1) / k*k
  if (intercept== TRUE){
    alpha_hat = A[ ,1]
    Gamma_hat = A[,2:ncol(A)]
  }else{
    alpha_hat <- NaN
    Gamma_hat <- A
  }


  # b) GVAR-related models
  }else if(any(ModelType == c("GVAR single", "GVAR multi")))
  {
  Economies <- GVARinputs$Economies # necessary to ensure that the code works for the "GVAR sepQ" model
  GVARinputs$GVARFactors <- DataSet_BS(ModelType, t(RiskFactors), GVARinputs$Wgvar, Economies, FactorLabels)
  GVARpara <- GVAR(GVARinputs, N)
  if (intercept== TRUE){
    alpha_hat <- GVARpara$F0
    Gamma_hat <-GVARpara$F1
  }else{
    alpha_hat <- NaN
    Gamma_hat <- GVARpara$F1
  }
  # c) JLL-related models
  } else if(any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma")))
  {JLLinputs$WishSigmas <- 0
  JLLPara <- JLL(t(RiskFactors), N, JLLinputs)
    if (intercept== TRUE){
    alpha_hat <- JLLPara$k0
    Gamma_hat <-JLLPara$k1
  }else{
    alpha_hat <- NaN
    Gamma_hat <- JLLPara$k1
  }
  }


  Output <- list(Gamma_hat = Gamma_hat, alpha_hat = alpha_hat)

  return(Output)
}

###################################################################################################################
###################################################################################################################
#' Generate M data sets from VAR(1) model
#'
#'@param Phi feedback matrix (F x F)
#'@param M number of Monte Carlo replications
#'@param RiskFactors time series of the risk factors (T x F)
#'
#'@keywords internal
#'
#'@references
#' Bauer, Rudebusch and, Wu (2012). "Correcting Estimation Bias in Dynamic Term Structure Models". \cr
#' This function is based on to the"genVAR" Matlab function available at Cynthia Wu's website
#' (https://sites.google.com/view/jingcynthiawu/).


genVARbrw <- function(Phi, M, RiskFactors){


  T <- dim(RiskFactors)[1];  k <- dim(RiskFactors)[2]

  Y_mean <- colMeans(RiskFactors)
  Y_mean_rep <- t(t(Y_mean))%*%matrix(1,1,M)
  X_sim <- array(0, c(k, T, M))

    # VAR(1)
    # 1. obtain residuals
    # use bootstrapped residuals
    resid <- matrix(0,k, T-1)
    for (t in 2:T){ resid[,(t-1)] <- t(t(RiskFactors[t, ])) -t(t(Y_mean)) - (Phi %*% t(t(RiskFactors[t-1, ] - Y_mean)))}

    # 2. generate series
    # randomly select initial values from the data Y
    ind_start <- sample(T, M, replace = TRUE)
    X_sim[ ,1, ] <- t(RiskFactors[ind_start, ])
    for (t in 2:T){
      u_sim <- resid[ , sample(T-1, M, replace = TRUE)]
      X_sim[ ,t, ] <- Y_mean_rep + Phi%*% (drop(X_sim[ ,t-1,]) - Y_mean_rep) + u_sim
    }


    rownames(X_sim) <- colnames(RiskFactors)

  return(X_sim)
}

###################################################################################################################
###################################################################################################################
#' Find mean or median of OLS when DGP is VAR(1)
#'
#'@param theta parameters from the feedback matrix in vector form
#'@param M number of Monte Carlo replications
#'@param RiskFactors time series of the risk factors (T x F)
#'@param N number of country-specific spanned factors (scalar)
#'@param GVARinputs inputs used in the estimation of the GVAR-based models (see "GVAR" function). Default is set to NULL
#'@param JLLinputs inputs used in the estimation of the JLL-based models (see "JLL" function). Default is set to NULL
#'@param FactorLabels string-list based which contains the labels of all variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param flag_mean flag whether mean- (TRUE) or median- (FALSE) unbiased estimation is desired. Default is set to TRUE
#'
#'@keywords internal
#'
#'@references
#' Bauer, Rudebusch and, Wu (2012). "Correcting Estimation Bias in Dynamic Term Structure Models". \cr
#' This function is similar to the "m_var" Matlab function available at Cynthia Wu's website
#' (https://sites.google.com/view/jingcynthiawu/).
#'

m_var <- function(theta, M, RiskFactors, N, GVARinputs, JLLinputs, FactorLabels, Economies, ModelType, flag_mean = TRUE){

  T <- dim(RiskFactors)[1]; k <- dim(RiskFactors)[2]

    if (flag_mean== TRUE){  Phi_sample <- array(0,c(k,k,M))}else{ Phi_sample <- NaN }

    # 1) get coefficient matrix
    Phi_tilde <-  matrix(theta, k, k)

    # 2) simulate M datasets
    X_sim <- genVARbrw(Phi_tilde, M, RiskFactors)

    # 3) estimate VAR(1) on each series
    theta_new_i <- matrix(0, M, k^2)
    for (m in 1:M){
      Phi_new <- estVARbrw(t(X_sim[ , ,m]), ModelType, N, GVARinputs, JLLinputs, FactorLabels, Economies,
                           demean = TRUE, intercept = FALSE)$Gamma_hat

      if (flag_mean == TRUE){ Phi_sample[ , , m] <- Phi_new }
      theta_new_i[m, ] <- c(Phi_new) # All elements as a vector
      }

    # 4) Find mean/median of OLS estimates
    if (flag_mean== TRUE){
      # Mean
      Phi_new <- matrix(colMeans(theta_new_i), k, k )
    } else{
      # Median
      Phi_new <- matrix(apply(theta_new_i, 2 , stats::median), k, k )
    }


  Outputs <- list(Phi_new = Phi_new, Phi_sample = Phi_sample)

  return(Outputs)
}

###################################################################################################################
###################################################################################################################
#' Killan's VAR stationarity adjustment
#'
#'@param Phi_tilde VAR (1) bias-corrected feedback matrix from Bauer, Rudebusch and, Wu (2012)
#'@param Phi_hat  unrestricted VAR(1) feedback matrix
#'@param ev_restr maximum eigenvalue desired in the feedback matrix after the adjustement
#'
#'@keywords internal
#'
#'@return
#'stationary VAR(1)
#'
#'@references
#' Bauer, Rudebusch and, Wu (2012). "Correcting Estimation Bias in Dynamic Term Structure Models". \cr
#' This function is an adapted version of the"shrink_Phi" Matlab function available at Cynthia Wu's website
#' (https://sites.google.com/view/jingcynthiawu/).
#'

shrink_Phi <- function(Phi_tilde, Phi_hat, ev_restr){


  ev_bc <- max(abs(eigen(Phi_tilde)$value))   # get scalar largest absolute eigenvalues
  ev_ols <- max(abs(eigen(Phi_hat)$value))

  if (ev_bc>ev_restr){
    if (ev_ols<ev_restr){ # if Phi_hat actually satisfies eigenvalue restriction, then then shrink Phi_tilde to Phi_hat
  Phi_diff0 <- Phi_hat - Phi_tilde
  delta <- 1
  while ((ev_bc > ev_restr) & delta > 0){ # delta > 0 is necessary to ensure that Phi_tilde won't get explosive
    delta <- delta - 0.01
    Phi_diff <- delta*Phi_diff0
    Phi_tilde <- Phi_hat - Phi_diff
    ev_bc <- max(abs(eigen(Phi_tilde)$value))
    }
    } else{
    #  # OLS estimates do not satisfy restriction
    Phi_tilde <- Phi_hat
    }
}


  return(Phi_tilde)
}

###################################################################################################################
#' Compute the variance-covariance matrix after the bias correction procedure
#'
#'@param Phi_tilde Feedback matrix resulting from the bias-correction procedure
#'@param N number of country-specific spanned factors (scalar)
#'@param RiskFactors time series of the risk factors (T x F)
#'@param GVARinputs inputs used in the estimation of the GVAR-based models (see "GVAR" function). Default is set to NULL
#'@param JLLinputs inputs used in the estimation of the JLL-based models (see "JLL" function). Default is set to NULL
#'@param FactorLabels string-list based which contains the labels of all variables present in the model
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'
#'@keywords internal


Get_V_tilde_BC <- function(Phi_tilde, N, RiskFactors, GVARinputs, JLLinputs, FactorLabels, ModelType){

  if(any(ModelType == c("JPS original", "JPS global", "JPS multi"))){
    T <- nrow(RiskFactors)
    Xdem <- RiskFactors - matrix(1,T,1)%*%colMeans(RiskFactors)
    resid_tilde <- t(Xdem[2:T, ]) - Phi_tilde %*% t(Xdem[1:(T-1), ])
    V_tilde <- resid_tilde%*%t(resid_tilde)/(T-1)
  }
  else if(any(ModelType == c("GVAR single", "GVAR multi"))){
    Economies <- GVARinputs$Economies # necessary to ensure that the code works for the "GVAR single" model
    GVARinputs$GVARFactors <- DataSet_BS(ModelType, t(RiskFactors), GVARinputs$Wgvar, Economies, FactorLabels)
    V_tilde <- GVAR(GVARinputs, N)$Sigma_y
  }
  else if(any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma")))
  {JLLinputs$WishSigmas <- 1
  JLLPara <- JLL(t(RiskFactors), N, JLLinputs)
  if (any(ModelType == c("JLL original", "JLL No DomUnit"))){  V_tilde <- JLLPara$Sigmas$VarCov_NonOrtho }
  if (ModelType == "JLL joint Sigma"){ V_tilde <- JLLPara$Sigmas$VarCov_Ortho}
  }

  return(V_tilde)
}



###############################################################################################################
#' check whether mean/median of OLS is close to actual OLS estimates
#'
#'@param theta_tilde feedback matrix after the bias correction procedure
#'@param theta_OLS feedback matrix obtained by OLS estimation (without bias correction)
#'@param B_check number of bootstrap samples used in the closeness check
#'@param RiskFactors time series of the risk factors (F x T)
#'@param N number of country-specific spanned factors (scalar)
#'@param GVARinputs inputs used in the estimation of the GVAR-based models (see "GVAR" function). Default is set to NULL
#'@param JLLinputs inputs used in the estimation of the JLL-based models (see "JLL" function). Default is set to NULL
#'@param FactorLabels string-list based which contains the labels of all variables present in the model
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param ModelType string-vector containing the label of the model to be estimated
#'
#'
#'@keywords internal

Check_comparison__OLS <- function(theta_tilde, theta_OLS, B_check, RiskFactors, N, GVARinputs, JLLinputs,
                                  FactorLabels, Economies, ModelType){


  cat('... checking closeness of mean/median to actual estimates ... \n')
  Phi_set <- m_var(theta_tilde, B_check, RiskFactors, N, GVARinputs, JLLinputs, FactorLabels, Economies,
                   ModelType, flag_mean = TRUE)
  Phi_new <- Phi_set$Phi_new
  Phi_sample <- Phi_set$Phi_sample  # return sample
  theta_new <- c(Phi_new)

  dist <- sqrt(sum((theta_new - theta_OLS)^2)/length(theta_new))
  cat(paste('Root mean square distance: ', round(dist,6), "\n"))
  cat('... Done! \n\n')

  return(dist)
}

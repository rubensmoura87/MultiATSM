#' Estimate an unbiased VAR(1) using stochastic approximation (Bauer, Rudebusch and Wu, 2012)
#'
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param BRWinputs List containing the following necessary inputs for the estimation of the BRW model:
#'\enumerate{
#'\item flag_mean: flag whether mean- (TRUE) or median- (FALSE) unbiased estimation is desired. Default is set to TRUE;
#'\item gamma: adjustment parameter. Value parameters should vary between 0 and 1. Default is set to 0.5;
#'\item N_iter: number of iterations used in the stochatic approximation algorithm after burn-in. Default is set to 5,000;
#'\item N_burn: number of burn-in iterations used in the stochatic approximation algorithm. Default is set to 0.15*N_iter;
#'\item B: number of bootstrap samples per iteration to calculate noisy measure of mean/median of the OLS estimator.
#'        Default is set to 50;
#'\item check: flag whether the user wishes to perform the closeness check. Default is set to TRUE;
#'\item B_check: number of bootstrap samples used in the closeness check. Default is set to 100,000.
#'}
#'@param RiskFactors time series of the risk factors (T x F)
#'@param N number of country-specific spanned factors (scalar)
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels string-list based which contains the labels of all variables present in the model
#'@param GVARinputs inputs used in the estimation of the GVAR-based models (see "GVAR" function). Default is set to NULL
#'@param JLLinputs inputs used in the estimation of the JLL-based models (see "JLL" function). Default is set to NULL
#'@param ev_restr largest eigenvalue restriction under the P-measure. Default is set to 1
#'@param nargout number of elements present in the list of outputs. Default is set to 4
#'
#'
#'
#'@examples
#'\donttest{
#' data(CM_Factors)
#' Factors <- t(RiskFactors[1:7,])
#'
#' BRWinputs <- list()
#' BRWinputs$flag_mean <- TRUE
#' BRWinputs$gamma <- 0.4
#' BRWinputs$N_iter <- 1000
#' BRWinputs$N_burn <- 100
#' BRWinputs$B <- 10
#' BRWinputs$check <- 1
#' BRWinputs$B_check <- 5000
#'
#' Economies <- "China"
#' N <- 3
#' ModelType <- "JPS"
#' FactorLabels <- NULL
#'
#'
#' BRWpara <- Bias_Correc_VAR(ModelType, BRWinputs, Factors, N, Economies, FactorLabels)
#'}
#'
#'@return Bias-corrected VAR paramaters based on the framework of Bauer, Rudebusch and Wu (2012). The list contains:
#'\enumerate{
#'\item Phi_tilde estimated coefficient matrix (F x F);
#'\item mu_tilde: estimated intercept (F x 1);
#'\item V_tilde: estimated variance-covariance matrix (F x F);
#'\item dist: root mean square distance (scalar);
#'\item Phi_sample: sample estimated variance-covariance matrix used in the checks (F x F x B_check) - this output is
#'      reported if nargout is set to 5.
#'}
#'
#'
#'@references
#' Bauer, Rudebusch and, Wu (2012). "Correcting Estimation Bias in Dynamic Term Structure Models" \cr
#' This function is based on the "est_unb_var" Matlab function available at Cynthia Wu's website
#' (https://sites.google.com/view/jingcynthiawu/).
#'@export




Bias_Correc_VAR  <- function(ModelType, BRWinputs, RiskFactors, N, Economies, FactorLabels, GVARinputs = NULL,
                             JLLinputs = NULL, ev_restr = 1, nargout= 4){


  print('#########################################################################################################')
  print('############################ Bias correction estimation (BRW, 2012) - ###################################')
  print('#########################################################################################################')


  flag_mean <- BRWinputs$flag_mean
  gamma <- BRWinputs$gamma
  N_iter <- BRWinputs$N_iter
  N_burn <- BRWinputs$N_burn
  B <- BRWinputs$B
  check <- BRWinputs$check
  B_check <- BRWinputs$B_check


  if(is.null(flag_mean)){ flag_mean <- TRUE}
  if(is.null(gamma)){ gamma <- 0.5}
  if(is.null(N_iter)){ N_iter <- 5000}
  if(is.null(N_burn)){N_burn <- N_iter*0.15 }
  if(is.null(B)){B <- 50}

  if(is.null(check)){ check <- FALSE }
  if(is.null(B_check)){B_check <- 100000}



  T <- dim(RiskFactors)[1]
  k <- dim(RiskFactors)[2]

  # 1) OLS
  Phi_hat <- estVARbrw(RiskFactors, ModelType, N, GVARinputs, JLLinputs, FactorLabels, Economies,
                       demean = TRUE, intercept = FALSE)$Gamma_hat
  #print(paste('largest absolute eigenvalue OLS:', round(max(abs(eigen(Phi_hat)$value)), 6)))

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
      print(paste('****** iteration ******', j-N_burn))
      theta_tilde <- rowMeans(theta[ ,(N_burn+1):j])
      Phi_tilde <- matrix(theta_tilde,k,k)
      print(paste('largest absolute eigenvalue:  ', round(max(abs(eigen(Phi_tilde)$value)), 6)))
    }
  }

  theta_tilde <- rowMeans(theta[ ,(N_burn+1):(N_burn+N)])

  # 4) check whether mean/median of OLS is close to actual OLS estimates
  if (check ==TRUE){
    print('... checking closeness of mean/median to actual estimates ...')
    Phi_set <- m_var(theta_tilde, B_check, RiskFactors, N, GVARinputs, JLLinputs, FactorLabels, Economies,
                     ModelType, flag_mean = TRUE)
    Phi_new <- Phi_set$Phi_new
    Phi_sample <- Phi_set$Phi_sample  # return sample
    theta_new <- c(Phi_new)

    dist <- sqrt(sum((theta_new - theta_hat)^2)/length(theta_new))
  print(paste('root mean square distance: ', round(dist,6)))
  print('... Done!')
  }else{
    Phi_sample <- NaN
  }


  # 5) bias_adjusted estimates
  Phi_tilde <- matrix(theta_tilde,k,k)
  #ev_bc <- max(abs(eigen(Phi_tilde)$value))   # get scalar largest absolute eigenvalues
#  print(paste('largest eigenvalue after BC: ', round(ev_bc,6)))

  # impose restriction on eigenvalue
 Phi_tilde <- shrink_Phi(Phi_tilde, Phi_hat, ev_restr)


  # 7) choose intercept to match sample mean
  mu_tilde <- (diag(k) - Phi_tilde)%*% t(t(colMeans(RiskFactors)))

  # 8) residuals and their variance-covariance matrix
  Xdem <- RiskFactors - matrix(1,T,1)%*%colMeans(RiskFactors)

  if(any(ModelType == c("JPS", "JPS jointP", "VAR jointQ"))){
    resid_tilde <- t(Xdem[2:T,]) - Phi_tilde %*% t(Xdem[1:(T-1), ])
    V_tilde <- resid_tilde%*%t(resid_tilde)/(T-1)
  }

  if(any(ModelType == c("GVAR sepQ", "GVAR jointQ"))){
    Economies <- GVARinputs$Economies # necessary to ensure that the code works for the "GVAR sepQ" model
    GVARinputs$GVARFactors <- DataSet_BS(ModelType, t(RiskFactors), GVARinputs$Wgvar, Economies, FactorLabels)
    V_tilde <- GVAR(GVARinputs, N)$Sigma_y
  }

  if(any(ModelType == c("JLL original", "JLL NoDomUnit", "JLL jointSigma")))
  {JLLinputs$WishSigmas <- 1
  JLLPara <- JLL(t(RiskFactors), N, JLLinputs)
  if (any(ModelType == c("JLL original", "JLL NoDomUnit"))){  V_tilde <- JLLPara$Sigmas$VarCov_NonOrtho }
  if (ModelType == "JLL jointSigma"){ V_tilde <- JLLPara$Sigmas$VarCov_Ortho}
  }

  # 9) Compilled desired outputs
  if (nargout==5){
    output <- list(Phi_tilde, mu_tilde, V_tilde, dist, Phi_sample)
    names(output) <- c("Phi_tilde", "mu_tilde", "V_tilde", "dist", "Phi_sample")
  }else{
    output <- list(Phi_tilde, mu_tilde, V_tilde, dist)
    names(output) <- c("Phi_tilde", "mu_tilde", "V_tilde", "dist")
  }

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
    if(any(ModelType == c("JPS", "JPS jointP", "VAR jointQ"))){
    A <- Y%*%t(Z)%*%solve(Z%*%t(Z), tol = 1e-50) # k*(k+1) / k*k
  if (intercept== TRUE){
    alpha_hat = A[ ,1]
    Gamma_hat = A[,2:ncol(A)]
  }else{
    alpha_hat <- NaN
    Gamma_hat <- A
  }
    }

  # b) GVAR-related models
  if(any(ModelType == c("GVAR sepQ", "GVAR jointQ")))
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
  }


  # c) JLL-related models
  if(any(ModelType == c("JLL original", "JLL NoDomUnit", "JLL jointSigma")))
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




  Output <- list(Gamma_hat, alpha_hat)
  names(Output) <- c("Gamma_hat", "alpha_hat")

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
#'@references
#' Bauer, Rudebusch and, Wu (2012). "Correcting Estimation Bias in Dynamic Term Structure Models". \cr
#' This function is similar to the"genVAR" Matlab function available at Cynthia Wu's website
#' (https://sites.google.com/view/jingcynthiawu/).
#'

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



  Outputs <- list(Phi_new, Phi_sample)
  names(Outputs) <- c("Phi_new", "Phi_sample")

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

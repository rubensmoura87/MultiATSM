#' Estimate a GVAR(1) and a VARX(1,1,1)
#'
#
#'@param GVARinputs List containing the following necessary inputs for the estimation of the GVAR:
#' \enumerate{
#'        \item  Economies:  string-vector containing the names of the economies which are part of the economic system
#'        \item  'GVARFactors': list of all variables that are used in the estimation of the VARX \cr
#'                (see e.g. 'CM_Factors_GVAR' file);
#'        \item 'VARXtype': character-vector containing two possibilities:
#'  \itemize{
#'        \item 'unconstrained': model is estimated without any constrained (each equation is estimated individually by OLS);
#'        \item "constrained": model is estimated taking into account the fact that foreign-pricing-factors
#'                              do NOT impirge on
#'                               (i) domestic economic variables and (ii) domestic pricing factors.
#'                               (equations are estimated by restricted least squares)
#'          }
#'          \item 'Wgvar':  GVAR transition matrix (C x C) - see the output from 'Transition_Matrix' function
#' }
#'@param N    number of country-specific spanned factors (scalar)
#'
#'@importFrom pracma numel
#'
#'@return   A list containing
#'\enumerate{
#'\item parameters of the country-specific VARX(1,1,1)
#'\itemize{
#'\item intercept (M+Nx1);
#'\item phi_1   (M+N x M+N);
#'\item phi_1^star (M+N x M+N);
#'\item phi_g (M+N x M+N);
#'\item Sigma (M+N x G)
#'}
#'\item parameters of the GVAR.
#'\itemize{
#' \item F0 (F X 1);
#' \item F1 (F x F);
#' \item Sigma_y (F x F)
#'}
#'}
#'
#'
#'@examples
#' data(CM_Factors_GVAR)
#'
#' N <- 3
#'
#' GVARinputs <- list()
#' GVARinputs$Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#' GVARinputs$GVARFactors <- FactorsGVAR
#' GVARinputs$VARXtype <- "unconstrained"
#' GVARinputs$Wgvar <- matrix( c(0, 0.83, 0.86, 0.38,
#'                               0.65, 0, 0.13, 0.55,
#'                               0.32, 0.12, 0, 0.07,
#'                               0.03, 0.05, 0.01, 0), nrow = 4, ncol = 4)
#'
#' GVAR(GVARinputs, N)

#'@references
#' Chudik and Pesaran, (2016). "Theory and Practice of GVAR modelling" (Journal of Economic Surveys)
#'@export



GVAR <- function(GVARinputs, N){



  # 0. Preliminary work
# Time series dimension
T <- numel(GVARinputs$GVARFactors[[GVARinputs$Economies[1]]]$Factors[[1]])

# Labels of group of variables
DomAndStarLabels <- names(GVARinputs$GVARFactors[[GVARinputs$Economies[1]]]$Factors)
L <- length(DomAndStarLabels)

DomLabels <- DomAndStarLabels[1:(L/2)]
StarLabels <- DomAndStarLabels[(L/2+1):L]
GlobalLabels <- names(GVARinputs$GVARFactors$Global)

# Dimension of group of variables
M <- length(DomLabels) - N # Number of country-specific macro variables
G <- length(GVARinputs$GVARFactors$Global) # Number of economies in the system
C <- length(GVARinputs$Economies) # Number of economies in the system



#######################################################################################################
##################################### VARX(1,1,1) #####################################################
#######################################################################################################


# 1) Prepare variables to be used in the estimation
  # a) Z.t:
  Z.t <- list()
  for (i in 1:C){
    X <- matrix(NA, nrow= length(DomLabels), ncol=(T-1))
    for (j in 1:(M+N)){
      X[j,] <- GVARinputs$GVARFactors[[GVARinputs$Economies[i]]]$Factors[[j]][2:T]
    }
    Z.t[[i]] <- X
  }
  names(Z.t) <- GVARinputs$Economies

  # b) Z lagged (Z.Lt)
  Z.Lt <- list()
  for (i in 1:C){
    X <- matrix(NA, nrow= length(DomLabels), ncol=(T-1))
    for (j in 1:(M+N)){
      X[j,] <- GVARinputs$GVARFactors[[GVARinputs$Economies[i]]]$Factors[[j]][1:(T-1)]
    }
    Z.Lt[[i]] <- X
  }
  names(Z.Lt) <- GVARinputs$Economies

  # c) Z star lagged (Zstar.Lt)
  idx1 <- M + N
  Zstar.Lt <- list()
  for (i in 1:C){
    X <- matrix(NA, nrow= length(StarLabels), ncol=(T-1))
    for (j in 1:(M+N)){
      X[j,] <- GVARinputs$GVARFactors[[GVARinputs$Economies[i]]]$Factors[[idx1+j]][1:(T-1)]
    }
    Zstar.Lt[[i]] <- X
  }
  names(Zstar.Lt) <- GVARinputs$Economies

  # d) Global lagged (Global.Lt)
  Global.Lt <- list()
  X <- matrix(NA, nrow= G, ncol=(T-1))
  for (j in seqi(1,G)){
    X[j,] <- GVARinputs$GVARFactors$Global[[j]][1:(T-1)]
  }
  Global.Lt <- X
  rownames(Global.Lt) <- GlobalLabels



  # 2) Estimate the VARX(1,1)
  # Prepare regressor set in LHS and RHS
  LHS <- vector(mode='list', length = C)
  RHS <- vector(mode='list', length = C)
  names(LHS) <- GVARinputs$Economies
  names(RHS) <- GVARinputs$Economies
  for (i in 1:C){
    LHS[[i]] <- Z.t[[i]]
    RHS[[i]] <- rbind(rep(1, times=T-1), Z.Lt[[i]], Zstar.Lt[[i]], Global.Lt)
    rownames(RHS[[i]]) <- c("Intercept", DomLabels, StarLabels, GlobalLabels)
  }

  # Set the foreign contemporaneous matrices to be equal to zero
  phi0_star <- list()
  X <- matrix(0, nrow= M + N, ncol= M + N)
  for (i in 1:C){
    phi0_star[[i]]<- X
  }
  names(phi0_star) <- GVARinputs$Economies

  # a) unconstrained system:
  Coeff <- list()

  if (GVARinputs$VARXtype == 'unconstrained'){
    ModelEstimates <- list()
    et <- list()
    Sigma <- list()
    for (i in 1:C){
      ModelEstimates[[i]] <- stats::lm( t(LHS[[i]]) ~ t(RHS[[i]])-1) # RHS and LHS are Tx(M+N)
      Coeff[[i]] <- t(ModelEstimates[[i]]$coefficients)
      et[[i]] <-t(ModelEstimates[[i]]$residual)
      Sigma[[i]] <- (et[[i]]%*%t(et[[i]]))/T
      rownames(Coeff[[i]]) <- DomLabels
    }
  }

  # b) constrained system:
  if (GVARinputs$VARXtype == 'constrained'){
    idxIntercept<- 1
    idxM <- idxIntercept + M
    idxP <- idxM + N
    idxM.star <- idxP + M
    idxP.star <- idxM.star + N
    idx.global <- idxP.star + G

    Bcon <- list()
    Coeff <- list()
    eT <- list()
    Sigma <- list()

    for (i in 1:C){
      Bcon[[i]] <- matrix(NaN, nrow= nrow(LHS[[i]]), ncol=nrow(RHS[[i]])) # (M+N)x (2*(M+N)+G+1)
      # b.1) if constrains is avoid the effect of P* on P and M:
      Bcon[[i]][,(idxM.star+1):idxP.star] <- 0

      Coeff[[i]] <- Reg__OLSconstrained(Y= LHS[[i]], X= RHS[[i]],Bcon[[i]], G=NULL)
      colnames(Coeff[[i]]) <- c("Intercept", DomLabels, StarLabels, GlobalLabels)
      rownames(Coeff[[i]]) <- DomLabels

      eT[[i]] <- LHS[[i]] - Coeff[[i]]%*%RHS[[i]]
      Sigma[[i]] <- (eT[[i]]%*%t(eT[[i]]))/T
    }
    names(Coeff) <- GVARinputs$Economies
  }

  # 3) Prepare outputs:
  ParaVARX <- list()
  for (i in 1:C){
    ParaVARX[[i]] <- vector(mode = 'list', length = 6)
    names(ParaVARX[[i]]) <- c('Phi0', 'Phi1', 'Phi1.star', 'Phi.global', 'Sigma', 'Phi0.star')

    idxPhi0 <- 1
    idxPhi1 <- idxPhi0 + (N+M)
    idxPhi1.star <- idxPhi1 + (N+M)
    idxPhi.global <- idxPhi1.star+G

    ParaVARX[[i]][[1]] <- Coeff[[i]][,idxPhi0]
    ParaVARX[[i]][[2]] <- Coeff[[i]][,(idxPhi0+1):idxPhi1]
    ParaVARX[[i]][[3]] <- Coeff[[i]][,(idxPhi1+1):idxPhi1.star]
    ParaVARX[[i]][[4]] <- Coeff[[i]][, seqi(idxPhi1.star+1,idxPhi.global)]
    ParaVARX[[i]][[5]] <- Sigma[[i]]
    ParaVARX[[i]][[6]] <- phi0_star[[i]]
  }
  names(ParaVARX) <- GVARinputs$Economies


#######################################################################################################
##################################### GVAR(1) #########################################################
#######################################################################################################


  #  1) Build a0, Ai0, Ai1 and Wi
  # a) a0
  a0 <- matrix(NA, nrow=C*(N+M), ncol=1)

  count0 <- 0
  for (j in 1:C){
    count1 <- count0+(N+M)
    a0[(count0+1):count1] <- do.call(rbind,lapply(ParaVARX[[GVARinputs$Economies[j]]]$Phi0,matrix,ncol=1))
    count0  <- count1
  }

  # b) Ai0:
  Ai0 <- list()
  for (i in 1:C){
    Ai0[[i]] <- cbind(diag(N+M), ParaVARX[[GVARinputs$Economies[i]]]$Phi0.star)
  }
  names(Ai0) <- GVARinputs$Economies

  # c) Ai1
  Ai1 <- list()
  for (i in 1:C){
    Ai1[[i]] <- cbind(ParaVARX[[GVARinputs$Economies[i]]]$Phi1, ParaVARX[[GVARinputs$Economies[i]]]$Phi1.star)
  }
  names(Ai1) <- GVARinputs$Economies

  # d) Wi
  # Bottom part of Wi
  bottomWi <- list()
  for(i in 1:C){
    b <- matrix(NA, nrow=N+M, ncol= C*(N+M))
    bottomWi[[i]] <- b
  }

  a <- matrix(NA, nrow=N+M, ncol= N+M)
  for (j in 1:C){
    count0 <- 0
    for (i in 1:C){
      count1 <- count0 +  (N+M)
      a <- GVARinputs$Wgvar[j,i]*diag(N+M)
      bottomWi[[j]][,(count0+1):count1] <- a
      count0 <- count1
    }
  }
  names(bottomWi) <- GVARinputs$Economies

  # Top part of Wi
  topWi <- list()
  for(i in 1:C){
    c <- matrix(NA, nrow= N+M, ncol= C*(N+M))
    topWi[[i]] <- c
  }

  d <- matrix(NA, nrow=N+M, ncol= N+M)
  IndexPos <- diag(C)
  for (j in 1:C){
    count0 <- 0
    for (i in 1:C){
      count1 <- count0 +  (N+M)
      d <- IndexPos[j,i]*diag(N+M)
      topWi[[j]][,(count0+1):count1] <- d
      count0 <- count1
    }
  }
  names(topWi) <- GVARinputs$Economies

  # Concatenate TopWi and bottomWi in Wi
  Wi <- list()
  for (i in 1:C){
    Wi[[i]] <- rbind(topWi[[i]], bottomWi[[i]])
  }
  names(Wi) <- GVARinputs$Economies


  # 2) Compute G0, G1 and Sigma
  # a) G0 and G1
  G0prep <- list()
  G1prep <- list()
  for (i in 1:C){
    G0prep[[i]] <- Ai0[[i]]%*%Wi[[i]]
    G1prep[[i]] <- Ai1[[i]]%*%Wi[[i]]
  }
  G0 <- do.call(rbind,lapply(G0prep,matrix,ncol=C*(N+M)))
  G1 <- do.call(rbind,lapply(G1prep,matrix,ncol=C*(N+M)))

  # b) Sigma
  Sigma <- matrix(0, ncol=C*(N+M), nrow=C*(N+M) )
  count0 <- 0
  for (i in 1:C){
    count1 <- count0+ (N+M)
    Sigma[(count0+1):count1,(count0+1):count1] <- ParaVARX[[GVARinputs$Economies[i]]]$Sigma
    count0 <- count1
  }


  # 3) Compute Gy.0,  Gy.1
  # a) Estimate the dynamics of the global variables: VAR(1)
  X <- do.call(rbind,lapply(GVARinputs$GVARFactors$Global,matrix,ncol=T))
  if (length(X) != 0 ){
  X <- t(X) # useful command for the case in which G <- 1
  RHS <- as.matrix(X[2:T,])
  LHS <- as.matrix(X[1:(T-1),])

  Phi.w1 <- matrix(NA, nrow=G, ncol= G)
  Phi.w0 <- matrix(NA, nrow=G, ncol= 1)
  for (i in seqi(1,G)){
    Phi.w0[i,] <- stats::lm( LHS[,i] ~ RHS)$coefficients[1]
    Phi.w1[i,] <- stats::lm( LHS[,i] ~ RHS)$coefficients[2:(G+1)]
  }

  eta.t <- matrix(NA, nrow= G, ncol=T-1)
  for (j in 1:(T-1)){
    eta.t[,j] <- LHS[j,] - Phi.w0 - Phi.w1%*%RHS[j,]
  }

  Sigma_w <- (eta.t%*%t(eta.t))/T
  } else{
    Phi.w0 <- c()
    Phi.w1 <- matrix(,ncol = 0, nrow=0)
    Sigma_w <- matrix(,ncol = 0, nrow=0)
}


  # b) Gy.0:
  Gy.0 <- magic::adiag(diag(G), G0)
  # c) Gy.1:
  if (G != 0 ){
  D1 <- list()
  for (i in 1:C){
    D1[[i]] <- ParaVARX[[GVARinputs$Economies[[i]]]]$Phi.global
  }
  D1 <- do.call(rbind,lapply(D1,matrix,ncol=G))
} else { D1 <- c()}

  topGy.1 <- matrix(0, nrow =G, ncol= C*(N+M) +G)
  topGy.1[seqi(1,G),seqi(1,G)] <- Phi.w1

  bottomGy.1 <- cbind(D1,G1)
  Gy.1 <- rbind(topGy.1, bottomGy.1)

  # 4) Build the GVAR(1): y_t = F0 + F1* y_{t-1} + (Gy.0)^(-1)ey_t ( equation 19)
  # a) F0 and F1:
  Phi0VARX <- list()
  for (i in 1:C){
    Phi0VARX[[i]] <-  ParaVARX[[GVARinputs$Economies[[i]]]]$Phi0
  }
  Phi0VARX <- do.call(rbind,lapply(Phi0VARX,matrix,ncol=1))

  F0 <- solve(Gy.0)%*%rbind(Phi.w0, Phi0VARX)
  F1 <- solve(Gy.0)%*%Gy.1

  # b) Sigma_y:
  Sigma_y <- magic::adiag(Sigma_w, Sigma)


  # 5) Prepare labels of the tables
  labelsDomVar <- c()
  for (i in 1:C){
    labelsDomVarCS <- paste(DomLabels, GVARinputs$Economies[i])
    labelsDomVar <- append(labelsDomVar,labelsDomVarCS)
  }

  labelsTables <- c(GlobalLabels,labelsDomVar)

  colnames(F1) <- labelsTables
  rownames(F1) <- labelsTables

  colnames(Sigma_y) <- labelsTables
  rownames(Sigma_y) <- labelsTables

  rownames(F0) <- labelsTables


  GVARoutputs<-list(ParaVARX, Gy.0,F0,F1, Sigma_y)
  names(GVARoutputs) <- c("VARX","Gy.0","F0", "F1", "Sigma_y")

  return(GVARoutputs)

}



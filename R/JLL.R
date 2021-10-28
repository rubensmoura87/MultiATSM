#' Set of inputs present at JLL's P-dynamics
#'
#'@param NonOrthoFactors Risk factors before the orthogonalization (FxT)
#'@param N Number of country-specific spanned factors
#'@param JLLinputs List of necessary inputs to estimate JLL outputs:
#'  \enumerate{
#'      \item Economies:  set of economies that are part of the economic system (string-vector)
#'      \item "DomUnit": name of the economy which is assigned as the dominant unit. \cr
#'                  If no dominant unit is assigned, then this variable is defined as "None"
#'      \item WishSigmas: equal to "1" if one wishes the variance-covariance matrices and the Cholesky factorizations
#'                  (can take long if they need to be estimated). Set "0", otherwise.
#'      \item SigmaNonOrtho: NULL or some K x K matrix from the non-orthogonalized dynamics
#'      \item JLLModelType: available options are "JLL original", "JLL jointSigma"  or "JLL NoDomUnit"
#'  }

#'@examples
#'\donttest{
#'if (requireNamespace('neldermead', quietly = TRUE)) {
#'
#'data(CM_Factors)
#' ZZ <- RiskFactors
#' N <- 3
#'
#' JLLinputs <- list()
#' JLLinputs$Economies <- c( "China", "Brazil", "Mexico", "Uruguay")
#' JLLinputs$DomUnit <- "China"
#' JLLinputs$WishSigmas <- 1
#' JLLinputs$SigmaNonOrtho <- NULL
#' JLLinputs$JLLModelType <- "JLL original"

#' JLL(ZZ, N, JLLinputs)
#'
#' #'} else {
#'  message("skipping functionality due to missing Suggested dependency")
#'}
#'}

#'@details
#' For the models 'JLL original' or "JLL jointSigma" the name of one dominant economy must assigned.\cr
#' For the model 'JLL  NoDomUnit', the name of one dominant economy must be set as "None".
#'@references
#' Jotiskhatira, Le and Lundblad (2015). "Why do interest rates in different currencies co-move?" (Journal of Financial Economics)
#'@return      List of model parameters from both the orthogonalized and non-orthogonalized versions of the JLL's based models
#
#'@export


JLL <- function(NonOrthoFactors, N, JLLinputs){


# 0. Preliminary works/checks

 if ((grepl("JLL original", JLLinputs$JLLModelType) ||
      grepl("JLL jointSigma", JLLinputs$JLLModelType)) & JLLinputs$DomUnit == "None"){
    stop("In 'JLL original' DomUnit cannot be 'None'. One dominant country is required! " )
 }

  if ( grepl("JLL NoDomUnit", JLLinputs$JLLModelType) & JLLinputs$DomUnit != "None"){
    stop("In 'JLL NoDomUnit' DomUnit cannot cannot contain a name of a country. DomUnit must be set as 'None'! " )
  }


  if (JLLinputs$DomUnit != "None"){
    IdxDomUnit <- which(JLLinputs$DomUnit== JLLinputs$Economies) # Index of the dominant country
  }

# System dimension
T <- ncol(NonOrthoFactors)
K <- nrow(NonOrthoFactors)
C <- length(JLLinputs$Economies)

# Extract the number of global factors
G <- c()
for (h in 1:K){
  G[h] <- all(sapply(JLLinputs$Economies, grepl, rownames(NonOrthoFactors))[h,] == 0)
}

G <- length(G[G==TRUE]) # Number of global unspnned factors
M <- (K-G)/C - N # Number of domestic unspnned factors

# Factor labels
FactorsJLL <- c()
FactorLabels <- list()

for (i in 1:C){
  Idx <- which(grepl(JLLinputs$Economies[i], rownames(NonOrthoFactors)))
  FactorLabels[[JLLinputs$Economies[i]]] <-  rownames(NonOrthoFactors)[Idx]

  JLLlab <- paste(FactorLabels[[JLLinputs$Economies[i]]], "JLL")
  FactorsJLL <- append(FactorsJLL, JLLlab)
}

FactorLabels$Global <- rownames(NonOrthoFactors)[seqi(1,G)]
LabelsJLL <- c(FactorLabels$Global, FactorsJLL)


# 1) Pre-allocation of the factors set
MacroGlobal <- NonOrthoFactors[c(FactorLabels$Global),]
FullFactorsSet <- list()

for (i in 1:C){
  FullFactorsSet[[JLLinputs$Economies[i]]]$Macro <- NonOrthoFactors[FactorLabels[[JLLinputs$Economies[i]]], ][1:M,]
  FullFactorsSet[[JLLinputs$Economies[i]]]$Pricing <- NonOrthoFactors[FactorLabels[[JLLinputs$Economies[i]]][(M+1):(N+M)],]
}


# 2) Orthogonalization of pricing factors
PricingRegressEQ6 <- list()
b <- list()
P_e <- list()
c <- list()

# Equation 6
for (i in 1:C){
  PricingRegressEQ6[[JLLinputs$Economies[i]]] <- stats::lm( t(FullFactorsSet[[JLLinputs$Economies[i]]]$Pricing)~ t(FullFactorsSet[[JLLinputs$Economies[i]]]$Macro) -1)
  b[[JLLinputs$Economies[i]]] <- t(PricingRegressEQ6[[JLLinputs$Economies[i]]]$coefficients)
  P_e[[JLLinputs$Economies[i]]] <- t(PricingRegressEQ6[[JLLinputs$Economies[i]]]$residuals)
}

# Equation 10
if (JLLinputs$DomUnit == "None"){
  for (i in 1:C){
    c[[JLLinputs$Economies[i]]] <- matrix(0, N, N)
  }
}else{
  PricingRegressEQ10 <- list()
  P_e_star <- list()
  c <- list()

  for (j in 1:(C-1)){
    PricingRegressEQ10[[JLLinputs$Economies[-IdxDomUnit][j]]] <- stats::lm( t(P_e[[JLLinputs$Economies[-IdxDomUnit][j]]]) ~ t(P_e[[JLLinputs$DomUnit]]) -1)
    P_e_star[[JLLinputs$Economies[-IdxDomUnit][j]]] <- t(PricingRegressEQ10[[JLLinputs$Economies[-IdxDomUnit][j]]]$residuals)
    c[[JLLinputs$Economies[-IdxDomUnit][j]]] <- t(PricingRegressEQ10[[JLLinputs$Economies[-IdxDomUnit][j]]]$coefficients)
  }
}


# 3) Orthogonalization of macro factors
MacroRegressEQ8 <- list()
a_W <- list()
M_e <- list()
a_DU_CS <- list() # Dominant unit country
if (JLLinputs$DomUnit == "None"){
  # Equation 8
  for (i in 1:C){
    MacroRegressEQ8[[JLLinputs$Economies[i]]] <- stats::lm( t(FullFactorsSet[[JLLinputs$Economies[i]]]$Macro)~ t(MacroGlobal) -1)
    a_W[[JLLinputs$Economies[i]]] <- t(MacroRegressEQ8[[JLLinputs$Economies[i]]]$coefficients)
    M_e[[JLLinputs$Economies[i]]] <- t(MacroRegressEQ8[[JLLinputs$Economies[i]]]$residuals)
    a_DU_CS[[JLLinputs$Economies[i]]] <- matrix(0, M, G)
  }

}else{
  MacroRegressEQ9 <- list()
  M_e_CS <- list()
  # Equation 8
  MacroRegressEQ8[[JLLinputs$DomUnit]] <- stats::lm( t(FullFactorsSet[[JLLinputs$DomUnit]]$Macro)~ t(MacroGlobal) -1)
  a_W[[JLLinputs$DomUnit]] <- t(MacroRegressEQ8[[JLLinputs$DomUnit]]$coefficients)
  M_e[[JLLinputs$DomUnit]] <- t(MacroRegressEQ8[[JLLinputs$DomUnit]]$residuals)

  # Equation 9
  for (j in 1:(C-1)){
    MacroRegressEQ9[[JLLinputs$Economies[-IdxDomUnit][j]]] <- stats::lm( t(FullFactorsSet[[JLLinputs$Economies[-IdxDomUnit][j]]]$Macro)~ t(MacroGlobal) + t(M_e[[JLLinputs$DomUnit]]) -1)
    a_W[[JLLinputs$Economies[-IdxDomUnit][j]]] <- t(MacroRegressEQ9[[JLLinputs$Economies[-IdxDomUnit][j]]]$coefficients)[,seqi(1,G)]
    a_DU_CS[[JLLinputs$Economies[-IdxDomUnit][j]]] <-t(MacroRegressEQ9[[JLLinputs$Economies[-IdxDomUnit][j]]]$coefficients)[,(G+1):(G+M)]
    M_e_CS[[JLLinputs$Economies[-IdxDomUnit][j]]] <- t(MacroRegressEQ9[[JLLinputs$Economies[-IdxDomUnit][j]]]$residuals)
  }
}


# 4) Build the Pi matrices:
# PIb
PIb <- diag(K)
idxRow0 <- G+M
idxCol0 <- G
for (i in 1:C){
  idxRow1 <- idxRow0+N
  idxCol1 <- idxCol0 + M
  PIb[(idxRow0+1):idxRow1, (idxCol0+1):idxCol1] <- b[[JLLinputs$Economies[i]]]
  idxRow0 <- idxRow1+M
  idxCol0 <- idxCol1 +N
}

rownames(PIb) <- rownames(NonOrthoFactors)
colnames(PIb) <- LabelsJLL

# PIac
PIac <- diag(K)
idxRow00 <- G
for (i in 1:C){
  idxRow11 <- idxRow00+ M
  PIac[(idxRow00+1):idxRow11, 1:M] <- a_W[[JLLinputs$Economies[i]]]
  idxRow00 <- idxRow11+N
}

if (JLLinputs$DomUnit !='None' ){
  # Place the orthogonalization of the pricing factors with respect to the dominant unit
  idxRowaDU0 <- G + M +N
  idxColaDU0 <- G
  idxColaDU1 <- idxColaDU0 + M
  for (j in 1:(C-1)){
    idxRowaDU1 <- idxRowaDU0 + M
    PIac[(idxRowaDU0+1):idxRowaDU1, (idxColaDU0+1):idxColaDU1] <- a_DU_CS[[JLLinputs$Economies[-IdxDomUnit][j]]]
    idxRowaDU0 <- idxRowaDU1 + N
  }


  # Place the orthogonalization of the pricing factors with respect to the dominant unit
  # (c coefficients from equation 10)
  idxRowc0 <- G + M + N + M
  idxColc0 <- G + M
  idxColc1 <- idxColc0 + N
  for (j in 1:(C-1)){
    idxRowc1 <- idxRowc0 + N
    PIac[(idxRowc0+1):idxRowc1, (idxColc0+1):idxColc1] <- c[[JLLinputs$Economies[-IdxDomUnit][j]]]
    idxRowc0 <- idxRowc1 +M
  }

}


rownames(PIac) <- rownames(NonOrthoFactors)
colnames(PIac) <- LabelsJLL


# PI
PI <- PIb%*%PIac

rownames(PI) <- rownames(NonOrthoFactors)
colnames(PI) <- LabelsJLL

# 5) VAR(1) with orthogonalized factors
AllDomFactorsOrtho <- matrix(NA, nrow= C*(N+M), ncol = T)
count0 <- 0
if (JLLinputs$DomUnit== 'None'){
  for(i in 1:C){
    count1 <- N+M + count0
    AllDomFactorsOrtho[(count0+1):count1,] <- rbind(M_e[[JLLinputs$Economies[i]]],P_e[[JLLinputs$Economies[i]]])
    count0 <- count1
  }
}else{
  DomUnitOrtho <- rbind(M_e[[JLLinputs$DomUnit]], P_e[[JLLinputs$DomUnit]])
  NoDomUnitOrtho <- matrix(NA, nrow= (C-1)*(N+M), ncol = T)
  for(i in 1:(C-1)){
    count1 <- N+M + count0
    NoDomUnitOrtho[(count0+1):count1,] <- rbind(M_e_CS[[JLLinputs$Economies[-IdxDomUnit][i]]], P_e_star[[JLLinputs$Economies[-IdxDomUnit][i]]])
    count0 <- count1
  }
  AllDomFactorsOrtho <- rbind(DomUnitOrtho, NoDomUnitOrtho)
}

Ye <- rbind(MacroGlobal, AllDomFactorsOrtho)
rownames(Ye) <- LabelsJLL

# Build the vector of non-orthogonalized factors
# (to make sure that the dominant unit country will be placed right after the global factors)
if (JLLinputs$DomUnit == "None"){
  Y <- NonOrthoFactors
}else{
  Y <- matrix(NA, nrow= K, ncol = T)
  rownames(Y) <- rownames(NonOrthoFactors)
  Y[seqi(1,G), ] <- NonOrthoFactors[seqi(1,G),] # Global factors
  Y[(G+1):(G+M+N),] <- NonOrthoFactors[FactorLabels[[IdxDomUnit]],]  # Dominant country

  COUNTER0 <- G+M+N
  for(i in 1:(C-1)){ # Non-dominant countries
    COUNTER1 <- N+M + COUNTER0
    Y[(COUNTER0+1):COUNTER1,] <- NonOrthoFactors[FactorLabels[[JLLinputs$Economies[-IdxDomUnit][i]]],]
    COUNTER0 <- COUNTER1
  }
}



# 5.a) Set the constrains on the feedback matrix
Bcon<-  FeedbackMatrixRestrictionsJLL(JLLinputs$DomUnit, K,G,M,N)

# 5.b) Estimate the VAR(1) with the orthogonalized variables
intercept <- rep(1, times= T-1)
RHS <- rbind(intercept, Ye[,1:(T-1)])
LHS <- Ye[,2:T]

Coeff <- Reg__OLSconstrained(Y= LHS, X= RHS, Bcon, G=NULL)
k0_e <- Coeff[,1]
k1_e <- Coeff[,2:(K+1)]

# Add Labels to k1_e
rownames(k1_e) <- LabelsJLL
colnames(k1_e) <- LabelsJLL

# 6) Obtain the non-orthogonalized factors:
k0 <- PI%*%k0_e
k1 <- PI%*%k1_e%*%solve(PI)

# Ensures that the almost zero elements of k1 are actually zeros (this procedure avoids weired IRFs)
idxZEROS <- which(Bcon[,2:(K+1)] == 0)
k1[idxZEROS] <- 0



# 6) Sigmas/Cholesky factorizations

if (JLLinputs$WishSigmas == 1){
# If the Variance-covariance matrix of the orthogonalized factors are NOT provided
if (is.null(JLLinputs$SigmaNonOrtho)){

  T <- ncol(Ye)

  LHS <- Ye[,2:T]
  RHS <- Ye[,1:(T-1)]


  et <- LHS - k0_e - k1_e%*%RHS
  SIGMA_Unres <- et%*%t(et)/dim(et)[2]

  #Labels
  rownames(SIGMA_Unres) <- LabelsJLL
  colnames(SIGMA_Unres) <- LabelsJLL

  # If the estimation of SIGMA_Ye is necessary
  Sigma_Ye <- EstimationSigma_Ye(SIGMA_Unres, et, M, G, JLLinputs$Economies, JLLinputs$DomUnit)

  # Cholesky term (non-orthogonalized factors)
  Sigma_Y <- PI%*%Sigma_Ye

  # Variance-covariance matrices
  Sigma_Res_Ortho <- Sigma_Ye%*%t(Sigma_Ye) #  Orthogonalized dynamics
  Sigma_Res_NonOrtho <- Sigma_Y%*%t(Sigma_Y) # Non-orthogonalized dynamics


}else{
  Sigma_Res_NonOrtho <- JLLinputs$SigmaNonOrtho
  Sigma_Y <- t(chol(JLLinputs$SigmaNonOrtho))
  Sigma_Ye <- solve(PI)%*%Sigma_Y
  Sigma_Res_Ortho <- Sigma_Ye%*%t(Sigma_Ye)
}

  ZeroIdxSigmaJLL <- IDXZeroRestrictionsJLLVarCovOrtho(M, N, G, JLLinputs$Economies, JLLinputs$DomUnit) # Identify the zero elements of the orthogonalized variance-covariance matrix
  # (usuful for distinguishing real zeros from nearly zero elements, later on in the code)
  Sigmas <- list(Sigma_Res_Ortho, Sigma_Res_NonOrtho, Sigma_Y, Sigma_Ye, ZeroIdxSigmaJLL)

  names(Sigmas) <- c("VarCov_Ortho", "VarCov_NonOrtho", "Sigma_Y", "Sigma_Ye", "ZeroIdxSigmaJLLOrtho")
} else{
  Sigmas <- NULL
}


# 7) Prepare the outputs

outputs <- list(a_W, a_DU_CS, b, c, PIb, PIac, PI, Ye, k0_e, k1_e, k0, k1, Sigmas)

names(outputs) <- c('a_W','a_DU_CS', 'b', 'c', 'PIb', 'PIac', 'PI', 'Ye', 'k0_e','k1_e', 'k0', 'k1', 'Sigmas')


return(outputs)
}


###############################################################################################################
#' Find the indexes of zero-restrictions from the orthogonalized variance-covariance matrix from the JLL-based models
#

#'@param M number of country-specific unspanned factors (scalar)
#'@param N number of country-specific spanned factors (scalar)
#'@param G number of global unspanned factors (scalar)
#'@param Economies Set of economies that are part of the economic system (string-vector)
#'@param DomUnit Name of the economy which is assigned as the dominant unit. \cr
#'               If no dominant unit is assigned, then this variable is defined as "None"

#
#'@keywords internal
#'@return      restricted version of the JLL of the Cholesky factorization (F x F)
#



IDXZeroRestrictionsJLLVarCovOrtho <- function(M, N, G, Economies, DomUnit){

  C <- length(Economies)
  K <- (M+N)*C + G

  MatOnes <- matrix(1, nrow = K, ncol = K)
  # Transform the matrix to be the Cholesky form
  MatOnes[upper.tri(MatOnes)] <- 0
  CholOrtho <- MatOnes

  if (DomUnit != "None"){
    IdxDomUnit <- which(DomUnit== Economies) # Index of the dominant country
  }


  # Zero restrictions of global variables on spanned factors
  idx0Global <- G + M
  for (h in 1:C){
    idx1Global <- idx0Global+ N
    CholOrtho[(idx0Global+1):idx1Global, seqi(1,G)] <- 0
    idx0Global <- idx1Global + M
  }

  # Zero restrictions of macro domestic variables on spanned factors
  for (i in 1:C){
    idx0RowMacroSpanned <- G + M
    idx0ColMacroSpanned <- G +(i-1)*(M+N)
    idx1ColMacroSpanned <- idx0ColMacroSpanned + M
    for (h in 1:C){ # Fix the columns and loop through the rows
      idx1RowMacroSpanned <- idx0RowMacroSpanned + N
      CholOrtho[(idx0RowMacroSpanned+1):idx1RowMacroSpanned, (idx0ColMacroSpanned+1):idx1ColMacroSpanned] <- 0
      idx0RowMacroSpanned <- idx1RowMacroSpanned + M
    }
  }


  # Zero restrictions of spanned factors on macro domestic variables
  for (i in 1:C){
    idx0RowSpannedMacro <- G
    idx0ColSpannedMacro <- G + M +(i-1)*(M+N)
    idx1ColSpannedMacro <- idx0ColSpannedMacro + N
    for (h in 1:C){ # Fix the columns and loop through the rows
      idx1RowSpannedMacro <- idx0RowSpannedMacro + M
      CholOrtho[(idx0RowSpannedMacro+1):idx1RowSpannedMacro, (idx0ColSpannedMacro+1):idx1ColSpannedMacro] <- 0
      idx0RowSpannedMacro <- idx1RowSpannedMacro + N
    }
  }


  # Zero restrictions of Macro country i on Macro country j
  if ( DomUnit != "None"){
    for (i in 1:C){
      if (i!=IdxDomUnit){
        idx0RowMacroMacro <- G
        idx0ColMacroMacro <- G +(i-1)*(M+N)
        idx1ColMacroMacro <- idx0ColMacroMacro + M
        for (h in 1:C){ # Fix the columns and loop through the rows
          idx1RowMacroMacro <- idx0RowMacroMacro + M
          if (i != h){
            CholOrtho[(idx0RowMacroMacro+1):idx1RowMacroMacro, (idx0ColMacroMacro+1):idx1ColMacroMacro] <- 0
          }
          idx0RowMacroMacro <- idx1RowMacroMacro + N
        }
      }
    }
  }
  # Zero restrictions of Spanned factors of country i on Spanned factors country j
  if (DomUnit != "None"){
    for (i in 1:C){
      if ( i!=IdxDomUnit){
        idx0RowSpannedSpanned <- G + M
        idx0ColSpannedSpanned <- G + M +(i-1)*(M+N)
        idx1ColSpannedSpanned <- idx0ColSpannedSpanned + N
        for (h in 1:C){ # Fix the columns and loop through the rows
          idx1RowSpannedSpanned <- idx0RowSpannedSpanned + N
          if (i != h){
            CholOrtho[(idx0RowSpannedSpanned+1):idx1RowSpannedSpanned, (idx0ColSpannedSpanned+1):idx1ColSpannedSpanned] <- 0
          }
          idx0RowSpannedSpanned <- idx1RowSpannedSpanned + M
        }
      }
    }
  }


  VarCovOrtho <- CholOrtho%*%t(CholOrtho)
  IdxZerosVarCovOrtho  <- which(VarCovOrtho ==0)
  IdxZerosSigma_Ye<- which(CholOrtho == 0)

  IDXzerosJLL <- list(IdxZerosSigma_Ye, IdxZerosVarCovOrtho)
  names(IDXzerosJLL) <- c("Sigma_Ye","VarCovOrtho")

  return(IDXzerosJLL)

}


#' Generates the bootstrap-related outputs
#'
#'@param ModelType   string-vector containing the label of the model to be estimated
#'@param ModelParaPE point estimate from the model parameters (see the outputs of the "Optimization" function)
#'@param NumOutPE    point estimate from the numerical outputs (see the outputs of the "NumOutputs" function)
#'@param mat         vector of maturities (in years) used in the estimation
#'@param Economies   string-vector containing the names of the economies which are part of the economic system
#'@param InputsForOutputs list containing the desired inputs for the construction of IRFs, GIRFs, FEVDs, and GFEVDs.
#'@param FactorLabels   string-list based which contains the labels of all the variables present in the model
#'@param DataFrequency    character-based vector: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"
#'@param vararginPE    list containg starting values and constraints (see arguments of the "Optimization" function)
#'@param JLLinputs    list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'@param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#'
#'
#'@importFrom pracma ceil rand strcmp num2str tril
#'
#'
#'
#'@examples
#' # See examples in the vignette file of this package (Section 4).
#'
#'@returns
#'list containing the following elements:
#'\itemize{
#' \item list of model parameters for one each one the draws;
#' \item list of numerical outputs (IRFs, GIRFs, FEVDs, GFEVDs) for each one of the draws;
#' \item Confidence bands for the chosen level of significance.
#' }

#'@references
#' This function is a modified and extended version of the "VARirbound" function from "A toolbox for VAR analysis"
#' by Ambrogio Cesa-Bianchi (https://github.com/ambropo/VAR-Toolbox)
#' @export


Bootstrap <- function(ModelType, ModelParaPE, NumOutPE, mat, Economies, InputsForOutputs,
                      FactorLabels, DataFrequency, vararginPE, JLLinputs = NULL, GVARinputs = NULL){


  WishBoot<- InputsForOutputs[[ModelType]]$Bootstrap$WishBoot

  if (WishBoot ==0){ print("No Bootstrap analysis was generated")
    }
  else{


    methodBS <- InputsForOutputs[[ModelType]]$Bootstrap$methodBS
    BlockLength <- InputsForOutputs[[ModelType]]$Bootstrap$BlockLength
    ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
    pctg <- InputsForOutputs[[ModelType]]$Bootstrap$pctg

    N <- length(FactorLabels$Spanned)
    M <- length(FactorLabels$Domestic) - N
    G <- length(FactorLabels$Global)

    if (DataFrequency == "Daily All Days"){ dt <- 1/365}
    if (DataFrequency == "Daily Business Days"){ dt <- 1/252}
    if (DataFrequency == "Weekly"){ dt <- 1/52}
    if (DataFrequency == "Monthly"){ dt <- 1/12}
    if (DataFrequency == "Quarterly"){ dt <- 1/4}
    if (DataFrequency == "Annually"){ dt <- 1}

    ###################################### PRE-ALLOCATION ###########################################################
    #################################################################################################################
    # Pre-allocation of list of outputs
    ModelBootstrap <- list()
    ModelParaPE_artificial <- list()
    NumericalOutputs_artificial <- list()
    Bounds <- list()

    GeneralInputs <- list(mat, dt, N)
    names(GeneralInputs) <- c("mat", "dt", "N")

    ModelBootstrap <- list(GeneralInputs, ModelParaPE_artificial, NumericalOutputs_artificial, Bounds)
    names(ModelBootstrap) <- c("GeneralInputs","ParaDraws", "NumOutDraws", "ConfBounds")

    C <- length(Economies)

    ## Loop over the number of draws
    tt <- 1 # numbers of accepted draws
    ww <- 1 # index for printing on screen

    # P-dynamics settings
    nlag  <- 1
    const <- 1


    for (i in 1:C){

      if (( ModelType == "GVAR jointQ" || ModelType == "VAR jointQ" || ModelType == "JLL original"
            || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma" )   & i >1 ){break} # the break avoids the models with joint estimation under the Q to be estimated C times


      # Select the original data and labels (preliminary work)

      if (ModelType == 'JPS' || ModelType == 'JPS jointP' || ModelType == "GVAR sepQ"){  # country-by-country estimation

        print('#########################################################################################################')
        print( paste('#################################', 'Bootstrap for model', ModelType, '-', Economies[i], '#################################' ))
        print('#########################################################################################################')

        RiskFactorLabels <- rownames(ModelParaPE[[ModelType]][[Economies[i]]]$inputs$AllFactors)
        YieldsLabels  <- rownames(ModelParaPE[[ModelType]][[Economies[i]]]$inputs$Y)

        YY <- ModelParaPE[[ModelType]][[Economies[i]]]$inputs$Y
        ZZ <- ModelParaPE[[ModelType]][[Economies[i]]]$inputs$AllFactors # K x T

        A <- ModelParaPE[[ModelType]][[Economies[i]]]$rot$P$A
        B <- ModelParaPE[[ModelType]][[Economies[i]]]$rot$P$B
        D0Z <-  ModelParaPE[[ModelType]][[Economies[i]]]$ests$K0Z
        D1Z <-  ModelParaPE[[ModelType]][[Economies[i]]]$ests$K1Z

      }else{ # Jointly estimation of countries

        print('#########################################################################################################')
        print( paste('#################################', 'Bootstrap for model', ModelType,  '#################################' ))
        print('#########################################################################################################')

        YieldsLabels  <- rownames(ModelParaPE[[ModelType]]$inputs$Y)

        YY <- ModelParaPE[[ModelType]]$inputs$Y
        ZZ <- ModelParaPE[[ModelType]]$inputs$AllFactors # K x T

        A <- ModelParaPE[[ModelType]]$rot$P$A
        B <- ModelParaPE[[ModelType]]$rot$P$B
        D0Z <-  ModelParaPE[[ModelType]]$ests$K0Z
        D1Z <-  ModelParaPE[[ModelType]]$ests$K1Z
      }

      T <- ncol(ZZ)
      K  <- nrow(ZZ) # NUmber of risk factors
      ZZ_row  <- t(ZZ) # T x K
      Y_row <- t(YY) # TxJ or TxCJ
      J <- length(mat)

      # Residuals from the original model
      # a) P-dynamics
      residPdyn <- PdynamicsSet_BS(ModelType, ZZ, FactorLabels, Economies, JLLinputs, GVARinputs)$eZ
      # b) observational equation
      # For models estimated on a country-by-country basis
      if (ModelType == 'JPS' || ModelType == 'JPS jointP' || ModelType == "GVAR sepQ"){
        BFull <- matrix(0, nrow = J, ncol= K)
        LabelSpannedCS <- c(FactorLabels$Tables[[Economies[i]]][-(1:M)])
        idxSpanned <- match(LabelSpannedCS, RiskFactorLabels)
        BFull[,idxSpanned] <- B
      } else{
        # For models estimated jointly
        BFull <- BUnspannedAdapJoint(G,M,N,C, J, B)
      }

      YYhat <- matrix(A, nrow= nrow(YY), ncol=ncol(YY)) + BFull%*%ZZ # Model-implied yields
      residYie <- t(YY - YYhat)


      ###################################### OPTIMIZATION #############################################################
      #################################################################################################################
      Jmisc::tic()
      while (tt<=ndraws){

        ## Create the matrices for the loop
        ZZ_artificial <- matrix(0, nrow= T -1 +nlag, ncol = K)
        rownames(ZZ_artificial) <- rownames(ZZ_row)
        colnames(ZZ_artificial) <- colnames(ZZ_row)

        Y_artificial <- matrix(0, nrow= T  , ncol = nrow(YY))
        rownames(Y_artificial) <- rownames(Y_row)
        colnames(Y_artificial) <- colnames(Y_row)

        # Display number of loops

        if (tt==10*ww){
          print(paste('Loop ', num2str(tt, fmt=0) , ' / ', num2str(ndraws, fmt=0), ' draws'))
          ww <- ww +1

        }


        ## STEP 1: Choose the method and generate the residuals
        if (strcmp(methodBS,'bs')){
          # Use the residuals to bootstrap: generate a random number bounded
          # between 0 and the number of residuals, then use the ceil function to select
          # that row of the residuals (this is equivalent to sampling with replacement)
          rr <- ceil((T-nlag)*rand(T,1)) # T x 1
          uPdyn <- residPdyn[rr[1:(T-nlag)], ] # (T-1) x K
          uYiel <- residYie[rr, ] # TxJ  or T x CJ
        } else if (strcmp(methodBS,'wild')){

          # Wild bootstrap based on simple distribution (~Rademacher)
          rr <- 1-2*(rand(T,1)>0.5)
          uPdyn<- residPdyn*(rr[1:(T-nlag)]%*%matrix(1,nrow=1, ncol=K))
          uYiel <- residYie*(rr%*%matrix(1,nrow=1, ncol= ncol(residYie)) )

        } else if (strcmp(methodBS,'block')){

          # Blocks overlap and are drawn with replacement
          FullBlocksSet <- dim(residPdyn)[1] - BlockLength +1 # all possible blocks that can be drawn
          SampleBlock <- ceil((T-nlag)/BlockLength) #

          bb <- ceil(SampleBlock*rand(SampleBlock,1))
          IdxBlocks <- matrix(NA, nrow= BlockLength, ncol= FullBlocksSet)
          for (mm in 1:FullBlocksSet){
            IdxBlocks[,mm] <- mm:(mm + BlockLength-1)
          }
          rr <- as.vector(IdxBlocks[,bb])[1:T]
          uPdyn <- residPdyn[rr[1:(T-nlag)],]
          uYiel <- residYie[rr, ]
        }else{
          stop(paste('The method ', methodBS, ' is not available'))
        }

        ## STEP 2: generate the artificial data
        # 2.1) initial values for the artificial data
        # Intialize the first nlag observations with real data
        LAG <- c()
        for (jj in 1:nlag){
          ZZ_artificial[jj,] <- ZZ_row[jj,]
          LAG <- rbind(ZZ_artificial[jj,], LAG)
        }
        # Initialize the artificial series and the LAGplus vector
        if (const==0){
          LAGplus <- LAG
        }else if (const==1){
          LAGplus <- cbind(1, LAG)
        }


        # 2.2) generate artificial series
        Ft <- rbind(t(D0Z), t(D1Z))
        # From observation nlag+1 to nobs, compute the artificial data
        for (jj in (nlag+1):(T-1+nlag)){
          for (mm in 1:K){
            # Compute the value for time=jj
            ZZ_artificial[jj,mm] = LAGplus %*% as.matrix(Ft[,mm]) + uPdyn[jj-nlag,mm]
          }
          # now update the LAG matrix
          if (jj<T-1+nlag){
            LAG <- rbind(ZZ_artificial[jj, ], LAG[1,seqi(1,((nlag-1)*K)) ] )
            if (const==0){LAGplus <- LAG} else if (const==1){
              LAGplus <- cbind(1, LAG)
            }
          }
        }


        # Yields
        Y_artificial <-  uYiel + matrix(A, nrow= T, ncol=nrow(YY), byrow=T) + ZZ_artificial%*%t(BFull)

        ## STEP 3: generate the artificial parameters of the ATSM

        i <<- i # Re-define i in the global envirnonment

        # 3.1) BUild the "Factor Set" list, which is necessary for the estimation of the GVAR
        FactorSet_BS <- DataSet_BS(ModelType, t(ZZ_artificial), GVARinputs$Wgvar, Economies, FactorLabels)
        GVARinputs$GVARFactors <- FactorSet_BS

        # 3.2) Test whether the VAR is stationary (if not, drop the draw)
        K1Z_artificial <- PdynamicsSet_BS(ModelType, t(ZZ_artificial), FactorLabels, Economies,
                                          JLLinputs, GVARinputs)$K1Z
        MaxEigen <-max(abs(eigen(K1Z_artificial)$value))


      if (MaxEigen < 0.99999 ) {
      # 3.3) Prepare the inputs used the in the llk
      InputsMLE_artificial <- InputsForMLEdensity_BS(ModelType, Y_artificial, ZZ_artificial, FactorLabels, mat,
                                                        Economies, DataFrequency, JLLinputs, GVARinputs)

      # 3.4) Variables that will be concentrared out of from the log-likelihood function
      K1XQ <- InputsMLE_artificial$K1XQ

    if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit" ){ SSZ <- NULL} else{SSZ <- InputsMLE_artificial$SSZ}


    ## STEP 4: collect the model inpits to build the llk
      f <- Functionf_Boot(ModelType, InputsMLE_artificial, Economies, mat, dt, FactorLabels, rr,
                               MaxEigen, JLLinputs, GVARinputs)

    ## STEP 5: prepare the optimization settings
    varargin <- list()

    varargin$K1XQ <- list(K1XQ, vararginPE$K1XQ$Label, vararginPE$K1XQ$LB, vararginPE$K1XQ$UB)
    varargin$SSZ <- list(SSZ, vararginPE$SSZ$Label, vararginPE$SSZ$LB, vararginPE$SSZ$UB)
    varargin$r0 <- list(NULL, vararginPE$r0$Label, vararginPE$r0$LB, vararginPE$r0$UB)
    varargin$se <- list(NULL, vararginPE$se$Label, vararginPE$se$LB, vararginPE$se$UB)
    varargin$K0Z <- list(NULL, vararginPE$K0Z$Label, vararginPE$K0Z$LB, vararginPE$K0Z$UB)
    varargin$K1Z <- list(NULL, vararginPE$K1Z$Label, vararginPE$K1Z$LB, vararginPE$K1Z$UB)
    varargin$OptRun <- c("iter off", "fminsearch only") # makes the estimation faster

    LabelVar <- c('Value', 'Label', 'LB', 'UB') # Elements of each parameter

    for (d in 1:(length(varargin)-1)){ names(varargin[[d]]) <-  LabelVar}

    tol <- 1e-1 # To avoid spending unecessary amount of time on flat regions of the llk functions
    ## STEP 6: Run the optimization
  if (ModelType == 'JPS' || ModelType == 'JPS jointP' || ModelType == "GVAR sepQ"){
  ModelBootstrap$ParaDraws[[ModelType]][[Economies[i]]][[tt]] <- Optimization_Boot(f, tol, varargin, FactorLabels,
                                                                                   Economies, ModelType, JLLinputs,
                                                                                   GVARinputs)
  }else{
  ModelBootstrap$ParaDraws[[ModelType]][[tt]] <- Optimization_Boot(f, tol, varargin, FactorLabels,
                                                                  Economies, ModelType, JLLinputs, GVARinputs)
  }

saveRDS(ModelBootstrap, paste(tempdir(),"/Bootstrap_", InputsForOutputs$'Label Outputs','.rds',sep="")) # if the optmization crashs after some draws,
          # we can keep the ouputs of draws before
      tt<-tt+1

        }
      }

      tt <- 1
      ww <- 1

    }

    print('-- Done!')

    Jmisc::toc()


    ###################################### NUMERICAL OUTPUTS ########################################################
    ModelBootstrap$NumOutDraws <- NumOutputs_Bootstrap(ModelType,ModelBootstrap, InputsForOutputs, FactorLabels,
                                                       Economies)


    ###################################### CONFIDENCE BOUNDS ########################################################
    ModelBootstrap$ConfBounds <- BootstrapBoundsSet(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies)


    ## To save space, clean the repeated outputs from the JLL outputs
    if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit"  || ModelType == "JLL jointSigma"){
        for (tt in 1:ndraws){
          ModelBootstrap$NumOutDraws$IRF[[ModelType]][[tt]]$Yields$Ortho <- NULL
          ModelBootstrap$NumOutDraws$FEVD[[ModelType]][[tt]]$Yields$Ortho <- NULL
        }
        ModelBootstrap$ConfBounds$IRF[[ModelType]]$Yields$Ortho <- NULL
        ModelBootstrap$ConfBounds$FEVD[[ModelType]]$Yields$Ortho <- NULL

    }

    saveRDS(ModelBootstrap, paste(tempdir(),"/Bootstrap_", InputsForOutputs$'Label Outputs','.rds',sep=""))

    return(ModelBootstrap)

  }
}


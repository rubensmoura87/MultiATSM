#' Builds the confidence bounds and graphs (Bootstrap set)
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelBootstrap list containing the complete set of model parameters after the bootstrap estimation procedure
#'@param NumOutPE      point estimate from the numerical outputs (see the outputs of the "NumOutputs" function)
#'@param InputsForOutputs list conataining the desired inputs for the construction of IRFs, GIRFs, FEVDs, and GFEVDs
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@keywords internal


BootstrapBoundsSet <- function(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies){

  # Generate the graph paths and the graph folders
  dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap", sep=""))
  PathsGraphs <- FolderCreationBoot(ModelType, Economies)


  # If one chooseS models in which the estimation is done country-by-country
  if ( any(ModelType ==c("JPS original", "JPS global", "GVAR single"))){

    for (i in 1:length(Economies)){dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/Model ", Economies[i], sep=""))}


    IRFandGIRFsep <- IRFandGIRFbs_sepQ(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs)
    FEVDandGFEVDsep <- FEVDandGFEVDbs_sepQ(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs)


    NumOutSep <- append(IRFandGIRFsep, FEVDandGFEVDsep)
  }

  # If one chooseS models in which the estimation is done jointly for all countries
  else {


    IRFandGIRFjoint <- IRFandGIRFbs_jointQ(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs,
                                           Economies, PathsGraphs)
    FEVDandGFEVDjoint <- FEVDandGFEVDbs_jointQ(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs,
                                               Economies, PathsGraphs)

    # Orthogonalized Outputs for JLL
    if ( any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){

      IRFandGIRFjoint_Ortho <- IRFandGIRFbs_jointQ_Ortho(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs,
                                                         Economies, PathsGraphs)

      FEVDandGFEVDjoint_Ortho <- FEVDandGFEVDbs_jointQ_Ortho(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs,
                                                             Economies, PathsGraphs)


        # Merge the lists of orthogonalized and non-orthogonalized factors
        for ( d in 1:length(IRFandGIRFjoint_Ortho)){ # IRFs abd GIRFs
          IRFandGIRFjoint[[d]][[ModelType]]$Factors <- list(IRFandGIRFjoint[[d]][[ModelType]]$Factors, IRFandGIRFjoint_Ortho[[d]][[ModelType]]$Factors)
         IRFandGIRFjoint[[d]][[ModelType]]$Yields <- list(IRFandGIRFjoint[[d]][[ModelType]]$Yields, IRFandGIRFjoint_Ortho[[d]][[ModelType]]$Yields)
          names(IRFandGIRFjoint[[d]][[ModelType]]$Factors) <- c("NonOrtho", "Ortho")
          names(IRFandGIRFjoint[[d]][[ModelType]]$Yields) <- c("NonOrtho", "Ortho")
        }

        for ( d in 1:length(FEVDandGFEVDjoint_Ortho)){ # FEVDs abd GFEVDs
          FEVDandGFEVDjoint[[d]][[ModelType]]$Factors <- list(FEVDandGFEVDjoint[[d]][[ModelType]]$Factors, FEVDandGFEVDjoint_Ortho[[d]][[ModelType]]$Factors)
          FEVDandGFEVDjoint[[d]][[ModelType]]$Yields <- list(FEVDandGFEVDjoint[[d]][[ModelType]]$Yields, FEVDandGFEVDjoint_Ortho[[d]][[ModelType]]$Yields)
          names(FEVDandGFEVDjoint[[d]][[ModelType]]$Factors) <- c("NonOrtho", "Ortho")
          names(FEVDandGFEVDjoint[[d]][[ModelType]]$Yields) <- c("NonOrtho", "Ortho")
        }

    }

    NumOutJoint <- append(IRFandGIRFjoint, FEVDandGFEVDjoint)

  }

  AllNumOutputs <- list()

  # Prepare final list of outputs
  if (!exists("NumOutJoint")){ AllNumOutputs <- NumOutSep}
  if (!exists("NumOutSep")){ AllNumOutputs <- NumOutJoint}
  if ( exists("NumOutSep") & exists("NumOutJoint")){

    for (i in 1:length(NumOutSep)){AllNumOutputs[[i]] <- append(NumOutSep[[i]], NumOutJoint[[i]]) }
    names(AllNumOutputs) <- names(NumOutSep)
 }

  cat(paste("Desired graphs are saved in your temporary directory. Please, check:",tempdir(), "\n\n"))
  return(AllNumOutputs)

}

######################################################################################################
######################################################################################################
####################### OUTPUTS FOR MODELS IN WHICH THE ESTIMATION ###################################
########################       IS DONE COUNTRY-BY-COUNTRY      #######################################
######################################################################################################
######################################################################################################
#' Creates the confidence bounds and the graphs of IRFs and GIRFs after bootstrap ("sep Q" models)
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelBootstrap list containing the complete set of model parameters after bootstrap estimation procedure
#'@param NumOutPE  list of model parameter point estimates
#'@param InputsForOutputs list conataining the desired inputs for the construction of the outputs of interest
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param PathsGraphs path of the folder in which the graphs will be saved
#'
#'@keywords internal


IRFandGIRFbs_sepQ <- function(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs){

  ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
  pctg <-   InputsForOutputs[[ModelType]]$Bootstrap$pctg

  #Define the percentiles
  pctg_inf <- (100-pctg)/2
  pctg_sup <- 100 - (100-pctg)/2
  quants <- c(pctg_inf, 50, pctg_sup)/100 # Desired quantiles

  # initializarion
  LabIRF <- c("IRF","GIRF")
  OutNames <- names(ModelBootstrap$NumOutDraws)
  C <- length(Economies)
  J <- length(ModelBootstrap$GeneralInputs$mat)

  HorizNumOut <- c(InputsForOutputs[[ModelType]]$IRF$horiz, InputsForOutputs[[ModelType]]$FEVD$horiz,
                   InputsForOutputs[[ModelType]]$GIRF$horiz, InputsForOutputs[[ModelType]]$GFEVD$horiz)

  NumOutBounds <- list()


    for (nn in match(LabIRF, OutNames) ){
      for (i in 1:C){

        K <- nrow(ModelBootstrap$ParaDraws[[ModelType]][[Economies[i]]][[1]]$ests$K1Z)

        ############################################# Factors ######################################################################

        INFfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # lower bound
        MEDfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # Median
        SUPfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # upper bound

        DimLabelsFac <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[Economies[i]]][[1]]$Factors)
        dimnames(INFfacs) <- DimLabelsFac
        dimnames(MEDfacs) <- DimLabelsFac
        dimnames(SUPfacs) <- DimLabelsFac


        # Allocation
        AllShocksOnePeriodFacs <- array(NA, c(ndraws, K, K))
        Facs <- matrix(NA, nrow = ndraws, ncol = K)

        for (thor in 1:HorizNumOut[[nn]]){
          for (h in 1:K){ # loop through the shocks
            for (g in 1:ndraws){ # # loop through the draws
              Facs[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[Economies[i]]][[g]]$Factors[thor,,h] # All responses to one shock in one horizon
            }

            AllShocksOnePeriodFacs[,,h]  <- apply(Facs,2, sort) # Ensures that each column is in ascending order

            INFfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[ , ,h], 2, stats::quantile, probs = quants[1])
            MEDfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[ , ,h], 2, stats::quantile, probs = quants[2])
            SUPfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[ , ,h], 2, stats::quantile, probs = quants[3])
          }
        }

        NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Factors$INF <- INFfacs
        NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Factors$MED <- MEDfacs
        NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Factors$SUP <- SUPfacs




        ############################################# Yields ######################################################################
        INFyies <- array(NA, c(HorizNumOut[[nn]], J, K)) # lower bound
        MEDyies <- array(NA, c(HorizNumOut[[nn]], J, K)) # Median
        SUPyies <- array(NA, c(HorizNumOut[[nn]], J, K)) # upper bound


        DimLabelsYies <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[Economies[i]]][[1]]$Yields)
        dimnames(INFyies) <- DimLabelsYies
        dimnames(MEDyies) <- DimLabelsYies
        dimnames(SUPyies) <- DimLabelsYies


        #Allocation
        AllShocksOnePeriodyies <- array(NA, c(ndraws, J, K))
        yies <- matrix(NA, nrow = ndraws, ncol = J)

        for (thor in 1:HorizNumOut[[nn]] ){
          for (h in 1:K){
            for (g in 1:ndraws){
              yies[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[Economies[i]]][[g]]$Yields[thor, ,h] # All responses to one shock in one horizon
            }
            AllShocksOnePeriodyies[,,h]  <- apply(yies,2, sort) # Ensures that each column is in ascending order

            INFyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[1])
            MEDyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[2])
            SUPyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[3])
          }
        }


        NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Yields$INF <- INFyies
        NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Yields$MED <- MEDyies
        NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Yields$SUP <- SUPyies



      }
    }



  #################################################################################################################
  ###################################### PLOTS ####################################################################
  #################################################################################################################

  ########################################  Factors ###############################################################

  GraphsBinaryFactors <- c(InputsForOutputs[[ModelType]]$IRF$WishGraphs$RiskFactorsBootstrap,
                           InputsForOutputs[[ModelType]]$GIRF$WishGraphs$RiskFactorsBootstrap)

  IdxFacGraphs <- which(GraphsBinaryFactors == 1)



   if (any(GraphsBinaryFactors==1)){

cat(' ** IRFs/GIRFs (Bootstrap) \n' )

    plot_list <- list()

      for (f in 1:C){
        for (d in IdxFacGraphs){
  # Folder Creation
  dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/Model ", Economies[f], "/", LabIRF[d], sep=""))
  dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/Model ", Economies[f], "/", LabIRF[d], "/Factors", sep=""))

          for (b in 1:K){
            for (i in 1:K){
              LL <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Factors$INF[,i,b]
              UU <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Factors$SUP[,i,b]
              MM <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Factors$MED[,i,b]
              PP <- NumOutPE[[LabIRF[d]]][[ModelType]][[Economies[f]]]$Factors[,i,b ] # Point estimate

              ALL <- data.frame(cbind(LL,MM,PP, UU))
              TimeSpan <- 1:InputsForOutputs[[ModelType]][[LabIRF[d]]]$horiz
              ALL$TimeSpan <- TimeSpan

              nmResponse <- dimnames(NumOutPE[[LabIRF[d]]][[ModelType]][[Economies[f]]]$Factors)[[2]] # names of the "response" factor
              nmShock <- dimnames(NumOutPE[[LabIRF[d]]][[ModelType]][[Economies[f]]]$Factors)[[3]] # names of the shock

              p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
                geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                      axis.title.x=element_blank(), axis.title.y=element_blank() ) +
                      ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

              plot_list[[i]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
            suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabIRF[d],"Factors_shock_to_", nmShock[b], ".png"),
                            path= PathsGraphs[[ModelType]][[LabIRF[d]]][[Economies[f]]][["Factors"]]))
          }
        }
      }


  }

  ########################################  Yields ###############################################################

  GraphsBinaryYields <- c(InputsForOutputs[[ModelType]]$IRF$WishGraphs$YieldsBootstrap,
                           InputsForOutputs[[ModelType]]$GIRF$WishGraphs$YieldsBootstrap)

  IdxYieldsGraphs <- which(GraphsBinaryYields == 1)



  if (any(GraphsBinaryYields==1)){

    plot_list <- list()
      for (f in 1:C){
        for (d in IdxYieldsGraphs){

          # Folder Creation
          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/Model ", Economies[f], "/", LabIRF[d], sep=""))
          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/Model ", Economies[f], "/", LabIRF[d], "/Yields", sep=""))

          for (b in 1:K){
            for (i in 1:J){
              LL <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Yields$INF[,i,b]
              UU <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Yields$SUP[,i,b]
              MM <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Yields$MED[,i,b]
              PP <- NumOutPE[[LabIRF[d]]][[ModelType]][[Economies[f]]]$Yields[,i,b ] # Point estimate

              ALL <- data.frame(cbind(LL,MM,PP, UU))
              TimeSpan <- 1:InputsForOutputs[[ModelType]][[LabIRF[d]]]$horiz
              ALL$TimeSpan <- TimeSpan

              nmResponse <- dimnames(NumOutPE[[LabIRF[d]]][[ModelType]][[Economies[f]]]$Yields)[[2]] # names of the "response" factor
              nmShock <- dimnames(NumOutPE[[LabIRF[d]]][[ModelType]][[Economies[f]]]$Yields)[[3]] # names of the shock

              p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
                geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                      axis.title.x=element_blank(), axis.title.y=element_blank() ) +
                ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

              plot_list[[i]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
            suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabIRF[d],"Yields_shock_to_", nmShock[b], ".png"),
                            path=  PathsGraphs[[ModelType]][[LabIRF[d]]][[Economies[f]]][["Yields"]]))
          }
        }
      }

  }


  return(NumOutBounds)

}
##############################################################################################################
#' Creates the confidence bounds and the graphs of FEVDs and GFEVDs after bootstrap ("sep Q" models)
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelBootstrap list containing the complete set of model parameters after bootstrap estimation procedure
#'@param NumOutPE  list of model parameter point estimates
#'@param InputsForOutputs list conataining the desired inputs for the construction of the outputs of interest
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param PathsGraphs path of the folder in which the graphs will be saved
#'
#'@keywords internal


FEVDandGFEVDbs_sepQ <- function(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs){


  ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
  pctg <-   InputsForOutputs[[ModelType]]$Bootstrap$pctg

  pctg_inf <- (100-pctg)/2
  pctg_sup <- 100 - (100-pctg)/2
  quants <- c(pctg_inf, 50, pctg_sup)/100 # Desired quantiles


  # initializarion
  LabFEVD <- c("FEVD","GFEVD")
  OutNames <- names(ModelBootstrap$NumOutDraws)
  C <- length(Economies)
  J <- length(ModelBootstrap$GeneralInputs$mat)


  HorizNumOut <- c(InputsForOutputs[[ModelType]]$IRF$horiz, InputsForOutputs[[ModelType]]$FEVD$horiz,
                   InputsForOutputs[[ModelType]]$GIRF$horiz, InputsForOutputs[[ModelType]]$GFEVD$horiz)

  NumOutBounds <- list()


  for (nn in match(LabFEVD, OutNames) ){
    for (i in 1:C){

      K <- nrow(ModelBootstrap$ParaDraws[[ModelType]][[Economies[i]]][[1]]$ests$K1Z)

      ############################################# Factors ######################################################################

      INFfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # lower bound
      MEDfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # Median
      SUPfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # upper bound

      DimLabelsFac <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[Economies[i]]][[1]]$Factors)
      dimnames(INFfacs) <- DimLabelsFac
      dimnames(MEDfacs) <- DimLabelsFac
      dimnames(SUPfacs) <- DimLabelsFac


      # Allocation
      AllShocksOnePeriodFacs <- array(NA, c(ndraws, K, K))
      Facs <- matrix(NA, nrow = ndraws, ncol = K)

      for (thor in 1:HorizNumOut[[nn]]){
        for (h in 1:K){ # loop through the shocks
          for (g in 1:ndraws){ # # loop through the draws
            Facs[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[Economies[i]]][[g]]$Factors[thor,,h] # All responses to one shock in one horizon
          }

          AllShocksOnePeriodFacs[,,h]  <- apply(Facs,2, sort) # Ensures that each column is in ascending order

          INFfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[1])
          MEDfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[2])
          SUPfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[3])
        }
      }

      NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Factors$INF <- INFfacs
      NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Factors$MED <- MEDfacs
      NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Factors$SUP <- SUPfacs




      ############################################# Yields ######################################################################
      INFyies <- array(NA, c(HorizNumOut[[nn]], K, J)) # lower bound
      MEDyies <- array(NA, c(HorizNumOut[[nn]], K, J)) # Median
      SUPyies <- array(NA, c(HorizNumOut[[nn]], K, J)) # upper bound


      DimLabelsYies <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[Economies[i]]][[1]]$Yields)
      dimnames(INFyies) <- DimLabelsYies
      dimnames(MEDyies) <- DimLabelsYies
      dimnames(SUPyies) <- DimLabelsYies


      #Allocation
      AllShocksOnePeriodyies <- array(NA, c(ndraws, K, J))
      yies <- matrix(NA, nrow = ndraws, ncol = K)

      for (thor in 1:HorizNumOut[[nn]]){
        for (h in 1:J){
          for (g in 1:ndraws){
            yies[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[Economies[i]]][[g]]$Yields[thor, ,h] # All responses to one shock in one horizon
          }
          AllShocksOnePeriodyies[,,h]  <- apply(yies,2, sort) # Ensures that each column is in ascending order

          INFyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[1])
          MEDyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[2])
          SUPyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[3])
        }
      }


      NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Yields$INF <- INFyies
      NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Yields$MED <- MEDyies
      NumOutBounds[[OutNames[nn]]][[ModelType]][[Economies[i]]]$Yields$SUP <- SUPyies



    }
  }


  #################################################################################################################
  ###################################### PLOTS ####################################################################
  #################################################################################################################

  ########################################  Factors ###############################################################

  GraphsBinaryFactors <- c(InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactorsBootstrap,
                          InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactorsBootstrap)

  IdxFactorsGraphs <- which(GraphsBinaryFactors == 1)



  if (any(GraphsBinaryFactors==1)){

    cat(' ** FEVDs/GFEVDs (Bootstrap) \n' )

    plot_list <- list()
    for (f in 1:C){
      for (d in IdxFactorsGraphs){
    # Folder Creation
      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/Model ", Economies[f], "/", LabFEVD[d], sep=""))
      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/Model ", Economies[f], "/", LabFEVD[d], "/Factors", sep=""))

        for (b in 1:K){
          for (i in 1:K){
            LL <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Factors$INF[,i,b]
            UU <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Factors$SUP[,i,b]
            MM <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Factors$MED[,i,b]
            PP <- NumOutPE[[LabFEVD[d]]][[ModelType]][[Economies[f]]]$Factors[,i,b ] # Point estimate # Point estimate

            ALL <- data.frame(cbind(LL,MM,PP, UU))
            TimeSpan <- 1:InputsForOutputs[[ModelType]][[LabFEVD[d]]]$horiz
            ALL$TimeSpan <- TimeSpan

            nmResponse <- dimnames(NumOutPE[[LabFEVD[d]]][[ModelType]][[Economies[f]]]$Factors)[[2]] # names of the "response" factor
            nmShock <- dimnames(NumOutPE[[LabFEVD[d]]][[ModelType]][[Economies[f]]]$Factors)[[3]] # names of the shock

            p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
              geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) + theme_classic() +
              theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                    axis.title.x=element_blank(), axis.title.y=element_blank() ) +
                    ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

            plot_list[[i]] <- p
          }
          subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabFEVD[d],"Factors_", nmShock[b], ".png"),
                          path=  PathsGraphs[[ModelType]][[LabFEVD[d]]][[Economies[f]]][["Factors"]]))
        }
      }
    }


  }

  ########################################  Yields ###############################################################

  GraphsBinaryYields <- c(InputsForOutputs[[ModelType]]$FEVD$WishGraphs$YieldsBootstrap,
                           InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$YieldsBootstrap)

  IdxYieldsGraphs <- which(GraphsBinaryYields == 1)



  if (any(GraphsBinaryYields==1)){

    plot_list <- list()
    for (f in 1:C){
      for (d in IdxYieldsGraphs){

        # Folder Creation
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/Model ", Economies[f], "/", LabFEVD[d], sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/Model ", Economies[f], "/", LabFEVD[d], "/Yields", sep=""))

        for (b in 1:J){
        for (i in 1:K){
            LL <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Yields$INF[,i,b]
            UU <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Yields$SUP[,i,b]
            MM <- NumOutBounds[[d]][[ModelType]][[Economies[f]]]$Yields$MED[,i,b]
            PP <- NumOutPE[[LabFEVD[d]]][[ModelType]][[Economies[f]]]$Yields[,i,b ] # Point estimate

            ALL <- data.frame(cbind(LL,MM,PP, UU))
            TimeSpan <- 1:InputsForOutputs[[ModelType]][[LabFEVD[d]]]$horiz
            ALL$TimeSpan <- TimeSpan

            nmResponse <- dimnames(NumOutPE[[LabFEVD[d]]][[ModelType]][[Economies[f]]]$Yields)[[2]] # names of the "response" factor
            nmShock <- dimnames(NumOutPE[[LabFEVD[d]]][[ModelType]][[Economies[f]]]$Yields)[[3]] # names of the shock

            p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
              geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) + theme_classic() +
              theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                    axis.title.x=element_blank(), axis.title.y=element_blank() ) +
              ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

            plot_list[[i]] <- p
          }
          subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabFEVD[d],"Yields_", nmShock[b], ".png"),
                          path=  PathsGraphs[[ModelType]][[LabFEVD[d]]][[Economies[f]]][["Yields"]]))
        }
      }
    }

  }


  return(NumOutBounds)

}


######################################################################################################
######################################################################################################
####################### OUTPUTS FOR MODELS IN WHICH THE ESTIMATION ###################################
########################       IS DONE FOR ALL COUNTRIES JOINTLY      ################################
######################################################################################################
######################################################################################################
#' Creates the confidence bounds and the graphs of IRFs and GIRFs after bootstrap ("joint Q" models)
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelBootstrap list containing the complete set of model parameters after bootstrap estimation procedure
#'@param NumOutPE  list of model parameter point estimates
#'@param InputsForOutputs list conataining the desired inputs for the construction of the outputs of interest
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'@param PathsGraphs path of the folder in which the graphs will be saved
#'
#'
#'@keywords internal



IRFandGIRFbs_jointQ <- function(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs){


  ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
  pctg <-   InputsForOutputs[[ModelType]]$Bootstrap$pctg

  #Define the percentiles
  pctg_inf <- (100-pctg)/2
  pctg_sup <- 100 - (100-pctg)/2
  quants <- c(pctg_inf, 50, pctg_sup)/100 # Desired quantiles

  # initializarion
  LabIRF <- c("IRF","GIRF")
  OutNames <- names(ModelBootstrap$NumOutDraws)
  C <- length(Economies)
  J <-  length(ModelBootstrap$GeneralInputs$mat)
  CJ <- C*J

  HorizNumOut <- c(InputsForOutputs[[ModelType]]$IRF$horiz, InputsForOutputs[[ModelType]]$FEVD$horiz,
                   InputsForOutputs[[ModelType]]$GIRF$horiz, InputsForOutputs[[ModelType]]$GFEVD$horiz)

  NumOutBounds <- list()


    for (nn in match(LabIRF, OutNames) ){

      K <- nrow(ModelBootstrap$ParaDraws[[ModelType]][[1]]$ests$K1Z)
      ############################################# Factors ######################################################################

      INFfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # lower bound
      MEDfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # Median
      SUPfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # upper bound

      if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
        DimLabelsFac <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Factors$NonOrtho)
      }else{
        DimLabelsFac <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Factors)
      }

      dimnames(INFfacs) <- DimLabelsFac
      dimnames(MEDfacs) <- DimLabelsFac
      dimnames(SUPfacs) <- DimLabelsFac


      # Allocation
      AllShocksOnePeriodFacs <- array(NA, c(ndraws, K, K))
      Facs <- matrix(NA, nrow = ndraws, ncol = K)

      for (thor in 1:HorizNumOut[[nn]]){
        for (h in 1:K){ # loop through the shocks
          for (g in 1:ndraws){ # # loop through the draws
            if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
              Facs[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Factors$NonOrtho[thor,,h]
            } else{
              Facs[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Factors[thor,,h] # All responses to one shock in one horizon
            }
          }

          AllShocksOnePeriodFacs[,,h]  <- apply(Facs,2, sort) # Ensures that each column is in ascending order

          INFfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[1])
          MEDfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[2])
          SUPfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[3])
        }
      }

      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$INF <- INFfacs
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$MED <- MEDfacs
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$SUP <- SUPfacs


      ############################################# Yields ######################################################################
      INFyies <- array(NA, c(HorizNumOut[[nn]], CJ, K)) # lower bound
      MEDyies <- array(NA, c(HorizNumOut[[nn]], CJ, K)) # Median
      SUPyies <- array(NA, c(HorizNumOut[[nn]], CJ, K)) # upper bound

      if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
        DimLabelsYies <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Yields$NonOrtho)
      } else{
        DimLabelsYies <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Yields)
      }

      dimnames(INFyies) <- DimLabelsYies
      dimnames(MEDyies) <- DimLabelsYies
      dimnames(SUPyies) <- DimLabelsYies


      #Allocation
      AllShocksOnePeriodyies <- array(NA, c(ndraws, CJ, K))
      yies <- matrix(NA, nrow = ndraws, ncol = CJ)

      for (thor in 1:HorizNumOut[[nn]]){
        for (h in 1:K){
          for (g in 1:ndraws){
            if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
              yies[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Yields$NonOrtho[thor, ,h]
            }else{
              yies[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Yields[thor, ,h] # All responses to one shock in one horizon
            }



          }
          AllShocksOnePeriodyies[,,h]  <- apply(yies,2, sort) # Ensures that each column is in ascending order

          INFyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[1])
          MEDyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[2])
          SUPyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[3])
        }
      }


      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$INF <- INFyies
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$MED <- MEDyies
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$SUP <- SUPyies

    }


  #################################################################################################################
  ###################################### PLOTS ####################################################################
  #################################################################################################################

  ########################################  Factors ###############################################################

  GraphsBinaryFactors <- c(InputsForOutputs[[ModelType]]$IRF$WishGraphs$RiskFactorsBootstrap,
                          InputsForOutputs[[ModelType]]$GIRF$WishGraphs$RiskFactorsBootstrap)

  IdxFacGraphs <- which(GraphsBinaryFactors == 1)


  if (any(GraphsBinaryFactors==1)){

    cat(' ** IRFs/GIRFs (Bootstrap) \n' )

    plot_list <- list()


      for (d in IdxFacGraphs){
        # Folder Creation
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabIRF[d], sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabIRF[d], "/Factors", sep=""))

        for (b in 1:K){
          for (i in 1:K){
            LL <- NumOutBounds[[d]][[ModelType]]$Factors$INF[,i,b]
            UU <- NumOutBounds[[d]][[ModelType]]$Factors$SUP[,i,b]
            MM <- NumOutBounds[[d]][[ModelType]]$Factors$MED[,i,b]
            if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
              PP <- NumOutPE[[LabIRF[d]]][[ModelType]]$Factors$NonOrtho[,i,b ]
            }else{ PP <- NumOutPE[[LabIRF[d]]][[ModelType]]$Factors[,i,b ]} # Point estimate

            ALL <- data.frame(cbind(LL,MM,PP, UU))
            TimeSpan <- 1:InputsForOutputs[[ModelType]][[LabIRF[d]]]$horiz
            ALL$TimeSpan <- TimeSpan

            nmResponse <- DimLabelsFac[[2]] # names of the "response" factor
            nmShock <- DimLabelsFac[[3]] # names of the shock

            p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
              geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) + theme_classic() +
              theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                    axis.title.x=element_blank(), axis.title.y=element_blank() ) +
              ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

            plot_list[[i]] <- p
          }
          subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabIRF[d],"Factors_shock_to_", nmShock[b], ".png"),
                          path= PathsGraphs[[ModelType]][[LabIRF[d]]][["Factors"]]))
        }
      }

  }

  ########################################  Yields ###############################################################

  GraphsBinaryYields <- c(InputsForOutputs[[ModelType]]$IRF$WishGraphs$YieldsBootstrap,
                          InputsForOutputs[[ModelType]]$GIRF$WishGraphs$YieldsBootstrap)

  IdxYieldsGraphs <- which(GraphsBinaryYields == 1)


  if (any(GraphsBinaryYields==1)){

    plot_list <- list()

      for (d in IdxYieldsGraphs){

        # Folder Creation
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabIRF[d], sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabIRF[d], "/Yields", sep=""))

        for (b in 1:K){
          for (i in 1:CJ){
            LL <- NumOutBounds[[d]][[ModelType]]$Yields$INF[,i,b]
            UU <- NumOutBounds[[d]][[ModelType]]$Yields$SUP[,i,b]
            MM <- NumOutBounds[[d]][[ModelType]]$Yields$MED[,i,b]
            if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
              PP <- NumOutPE[[LabIRF[d]]][[ModelType]]$Yields$NonOrtho[,i,b ]
            }else{ PP <- NumOutPE[[LabIRF[d]]][[ModelType]]$Yields[,i,b ] } # Point estimate

            ALL <- data.frame(cbind(LL,MM,PP, UU))
            TimeSpan <- 1:InputsForOutputs[[ModelType]][[LabIRF[d]]]$horiz
            ALL$TimeSpan <- TimeSpan

            nmResponse <- DimLabelsYies[[2]] # names of the "response" factor
            nmShock <- DimLabelsYies[[3]] # names of the shock

            p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
              geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) +  theme_classic() +
              theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                    axis.title.x=element_blank(), axis.title.y=element_blank() ) +
              ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

            plot_list[[i]] <- p
          }
          subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabIRF[d],"Yields_shock_to_", nmShock[b], ".png"),
                          path= PathsGraphs[[ModelType]][[LabIRF[d]]][["Yields"]]))
        }
        }

  }

  return(NumOutBounds)

}


##############################################################################################################
#' Creates the confidence bounds and the graphs of FEVDs and GFEVDs after bootstrap ("joint Q" models)
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelBootstrap list containing the complete set of model parameters after bootstrap estimation procedure
#'@param NumOutPE  list of model parameter point estimates
#'@param InputsForOutputs list conataining the desired inputs for the construction of the outputs of interest
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'@param PathsGraphs path of the folder in which the graphs will be saved
#'
#'@keywords internal


FEVDandGFEVDbs_jointQ <- function(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs){

  ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
  pctg <-   InputsForOutputs[[ModelType]]$Bootstrap$pctg


  pctg_inf <- (100-pctg)/2
  pctg_sup <- 100 - (100-pctg)/2
  quants <- c(pctg_inf, 50, pctg_sup)/100 # Desired quantiles


  # initializarion
  LabFEVD <- c("FEVD","GFEVD")
  OutNames <- names(ModelBootstrap$NumOutDraws)
  C <- length(Economies)
  J <-  length(ModelBootstrap$GeneralInputs$mat)
  CJ <- C*J


  HorizNumOut <- c(InputsForOutputs[[ModelType]]$IRF$horiz, InputsForOutputs[[ModelType]]$FEVD$horiz,
                   InputsForOutputs[[ModelType]]$GIRF$horiz, InputsForOutputs[[ModelType]]$GFEVD$horiz)

  NumOutBounds <- list()


      for (nn in match(LabFEVD, OutNames) ){


      K <- nrow(ModelBootstrap$ParaDraws[[ModelType]][[1]]$ests$K1Z)

      ############################################# Factors ######################################################################

      INFfacs <- array(NA, c( HorizNumOut[[nn]], K,K)) # lower bound
      MEDfacs <- array(NA, c( HorizNumOut[[nn]], K,K)) # Median
      SUPfacs <- array(NA, c( HorizNumOut[[nn]], K,K)) # upper bound
      if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
        DimLabelsFac <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Factors$NonOrtho)
      }else{
        DimLabelsFac <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Factors)
      }

      dimnames(INFfacs) <- DimLabelsFac
      dimnames(MEDfacs) <- DimLabelsFac
      dimnames(SUPfacs) <- DimLabelsFac


      # Allocation
      AllShocksOnePeriodFacs <- array(NA, c(ndraws, K, K))
      Facs <- matrix(NA, nrow = ndraws, ncol = K)

      for (thor in 1: HorizNumOut[[nn]]){
        for (h in 1:K){ # loop through the shocks
          for (g in 1:ndraws){ # # loop through the draws
            if (any(ModelType ==c("JLL original",  "JLL No DomUnit",  "JLL joint Sigma"))){
              Facs[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Factors$NonOrtho[thor,,h] # All responses to one shock in one horizon
            }else{
              Facs[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Factors[thor,,h]
            }

          }

          AllShocksOnePeriodFacs[,,h]  <- apply(Facs,2, sort) # Ensures that each column is in ascending order

          INFfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[1])
          MEDfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[2])
          SUPfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[3])
        }
      }

      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$INF <- INFfacs
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$MED <- MEDfacs
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$SUP <- SUPfacs


      ############################################# Yields ######################################################################
      INFyies <- array(NA, c(HorizNumOut[[nn]], K, CJ)) # lower bound
      MEDyies <- array(NA, c(HorizNumOut[[nn]], K, CJ)) # Median
      SUPyies <- array(NA, c(HorizNumOut[[nn]], K, CJ)) # upper bound

      if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
        DimLabelsYies <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Yields$NonOrtho)
      } else{
        DimLabelsYies <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Yields)
      }
      dimnames(INFyies) <- DimLabelsYies
      dimnames(MEDyies) <- DimLabelsYies
      dimnames(SUPyies) <- DimLabelsYies

      #Allocation
      AllShocksOnePeriodyies <- array(NA, c(ndraws, K, CJ))
      yies <- matrix(NA, nrow = ndraws, ncol = K)

      for (thor in 1:HorizNumOut[[nn]]){
        for (h in 1:CJ){
          for (g in 1:ndraws){

            if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
              yies[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Yields$NonOrtho[thor, ,h]
            }else{
              yies[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Yields[thor, ,h] # All responses to one shock in one horizon
            }

          }
          AllShocksOnePeriodyies[,,h]  <- apply(yies,2, sort) # Ensures that each column is in ascending order

          INFyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[1])
          MEDyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[2])
          SUPyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[3])
        }
      }


      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$INF <- INFyies
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$MED <- MEDyies
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$SUP <- SUPyies

    }


  #################################################################################################################
  ###################################### PLOTS ####################################################################
  #################################################################################################################

  ########################################  Factors ###############################################################
  GraphsBinaryFactors <- c(InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactorsBootstrap,
                           InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactorsBootstrap)

  IdxFacGraphs <- which(GraphsBinaryFactors == 1)



  if (any(GraphsBinaryFactors==1)){

    cat(' ** FEVDs/GFEVDs (Bootstrap) \n\n' )

    plot_list <- list()


      for (d in IdxFacGraphs){
        # Folder Creation
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabFEVD[d], sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabFEVD[d], "/Factors", sep=""))


        for (b in 1:K){
          for (i in 1:K){
            LL <- NumOutBounds[[d]][[ModelType]]$Factors$INF[,i,b]
            UU <- NumOutBounds[[d]][[ModelType]]$Factors$SUP[,i,b]
            MM <- NumOutBounds[[d]][[ModelType]]$Factors$MED[,i,b]
            if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
              PP <- NumOutPE[[LabFEVD[d]]][[ModelType]]$Factors$NonOrtho[,i,b ]
            }else{  PP <- NumOutPE[[LabFEVD[d]]][[ModelType]]$Factors[,i,b ]} # Point estimate}


            ALL <- data.frame(cbind(LL,MM,PP, UU))
            TimeSpan <- 1:InputsForOutputs[[ModelType]][[LabFEVD[d]]]$horiz
            ALL$TimeSpan <- TimeSpan

            nmResponse <- DimLabelsFac[[2]] # names of the "response" factor
            nmShock <- DimLabelsFac[[3]] # names of the shock

            p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
              geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) +  theme_classic() +
              theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                    axis.title.x=element_blank(), axis.title.y=element_blank() ) +
              ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

            plot_list[[i]] <- p
          }
          subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabFEVD[d],"Factors_", nmShock[b], ".png"),
                          path= PathsGraphs[[ModelType]][[LabFEVD[d]]][["Factors"]]))
        }
      }

  }

  ########################################  Yields ###############################################################

  GraphsBinarYields <- c(InputsForOutputs[[ModelType]]$FEVD$WishGraphs$YieldsBootstrap,
                           InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$YieldsBootstrap)

  IdxYieldGraphs <- which(GraphsBinarYields == 1)


  if (any(GraphsBinarYields==1)){


    plot_list <- list()

      for (d in IdxYieldGraphs){
        # Folder Creation
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabFEVD[d], sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabFEVD[d], "/Yields", sep=""))

        for (b in 1:CJ){
          for (i in 1:K){
            LL <- NumOutBounds[[d]][[ModelType]]$Yields$INF[,i,b]
            UU <- NumOutBounds[[d]][[ModelType]]$Yields$SUP[,i,b]
            MM <- NumOutBounds[[d]][[ModelType]]$Yields$MED[,i,b]
            if (any(ModelType ==c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
              PP <- NumOutPE[[LabFEVD[d]]][[ModelType]]$Yields$NonOrtho[,i,b ]
            }else{  PP <- NumOutPE[[LabFEVD[d]]][[ModelType]]$Yields[,i,b ] } # Point estimate

            ALL <- data.frame(cbind(LL,MM,PP, UU))
            TimeSpan <- 1:InputsForOutputs[[ModelType]][[LabFEVD[d]]]$horiz
            ALL$TimeSpan <- TimeSpan

            nmResponse <- DimLabelsYies[[2]] # names of the "response" factor
            nmShock <- DimLabelsYies[[3]] # names of the shock

            p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
              geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) + theme_classic() +
              theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                    axis.title.x=element_blank(), axis.title.y=element_blank() ) +
              ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

            plot_list[[i]] <- p
          }
          subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabFEVD[d],"Yields_", nmShock[b], ".png"),
                          path= PathsGraphs[[ModelType]][[LabFEVD[d]]][["Yields"]]))
        }
      }


  }

  return(NumOutBounds)

}


##############################################################################################################
########################  OUTPUTS FOR ORTHOGONALIZED FACTORS #################################################
##############################################################################################################
#' Creates the confidence bounds and the graphs of IRFs and GIRFs after bootstrap (JLL-based models)
#'
#'@param ModelType  string-vector containing the label of the model to be estimated
#'@param ModelBootstrap list containing the complete set of model parameters after bootstrap estimation procedure
#'@param NumOutPE  list of model parameter point estimates
#'@param InputsForOutputs list conataining the desired inputs for the construction of the outputs of interest
#'@param Economies  string-vector containing the names of the economies which are part of the economic system
#'@param PathsGraphs path of the folder in which the graphs will be saved
#'
#'
#'@keywords internal


IRFandGIRFbs_jointQ_Ortho <- function(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs){


  ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
  pctg <-   InputsForOutputs[[ModelType]]$Bootstrap$pctg


  #Define the percentiles
  pctg_inf <- (100-pctg)/2
  pctg_sup <- 100 - (100-pctg)/2
  quants <- c(pctg_inf, 50, pctg_sup)/100 # Desired quantiles

  # initializarion
  LabIRF <- c("IRF","GIRF")
  OutNames <- names(ModelBootstrap$NumOutDraws)
  C <- length(Economies)
  J <-  length(ModelBootstrap$GeneralInputs$mat)
  CJ <- C*J


  HorizNumOut <- c(InputsForOutputs[[ModelType]]$IRF$horiz, InputsForOutputs[[ModelType]]$FEVD$horiz,
                   InputsForOutputs[[ModelType]]$GIRF$horiz, InputsForOutputs[[ModelType]]$GFEVD$horiz)

  NumOutBounds <- list()


  for (nn in match(LabIRF, OutNames) ){


      K <- nrow(ModelBootstrap$ParaDraws[[ModelType]][[1]]$ests$K1Z)
      ############################################# Factors ######################################################################

      INFfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # lower bound
      MEDfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # Median
      SUPfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # upper bound

      DimLabelsFac <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Factors$Ortho)
      dimnames(INFfacs) <- DimLabelsFac
      dimnames(MEDfacs) <- DimLabelsFac
      dimnames(SUPfacs) <- DimLabelsFac


      # Allocation
      AllShocksOnePeriodFacs <- array(NA, c(ndraws, K, K))
      Facs <- matrix(NA, nrow = ndraws, ncol = K)

      for (thor in 1:HorizNumOut[[nn]]){
        for (h in 1:K){ # loop through the shocks
          for (g in 1:ndraws){ # # loop through the draws
            Facs[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Factors$Ortho[thor,,h] # All responses to one shock in one horizon
          }

          AllShocksOnePeriodFacs[,,h]  <- apply(Facs,2, sort) # Ensures that each column is in ascending order

          INFfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[1])
          MEDfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[2])
          SUPfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[3])
        }
      }

      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$INF <- INFfacs
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$MED <- MEDfacs
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$SUP <- SUPfacs


      ############################################# Yields ######################################################################
      INFyies <- array(NA, c(HorizNumOut[[nn]], CJ, K)) # lower bound
      MEDyies <- array(NA, c(HorizNumOut[[nn]], CJ, K)) # Median
      SUPyies <- array(NA, c(HorizNumOut[[nn]], CJ, K)) # upper bound


      DimLabelsYies <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Yields$Ortho)
      dimnames(INFyies) <- DimLabelsYies
      dimnames(MEDyies) <- DimLabelsYies
      dimnames(SUPyies) <- DimLabelsYies


      #Allocation
      AllShocksOnePeriodyies <- array(NA, c(ndraws, CJ, K))
      yies <- matrix(NA, nrow = ndraws, ncol = CJ)

      for (thor in 1:HorizNumOut[[nn]]){
        for (h in 1:K){
          for (g in 1:ndraws){
            yies[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Yields$Ortho[thor, ,h] # All responses to one shock in one horizon
          }
          AllShocksOnePeriodyies[,,h]  <- apply(yies,2, sort) # Ensures that each column is in ascending order

          INFyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[1])
          MEDyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[2])
          SUPyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[3])
        }
      }


      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$INF <- INFyies
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$MED <- MEDyies
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$SUP <- SUPyies

    }


  #################################################################################################################
  ###################################### PLOTS ####################################################################
  #################################################################################################################

  ########################################  Factors ###############################################################
  GraphsBinarFactors <- c(InputsForOutputs[[ModelType]]$IRF$WishGraphsOrtho$RiskFactorsBootstrap,
                         InputsForOutputs[[ModelType]]$GIRF$WishGraphsOrtho$RiskFactorsBootstrap)

  IdxFactorsGraphs <- which(GraphsBinarFactors == 1)


  if (any(GraphsBinarFactors==1)){

    cat(' ** IRFs/GIRFs-Ortho (Bootstrap) \n' )

    plot_list <- list()


      for (d in IdxFactorsGraphs){
        # Folder Creation
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabIRF[d], sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabIRF[d], "/Factors", sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabIRF[d], "/Factors/Ortho", sep=""))

        for (b in 1:K){
          for (i in 1:K){
            LL <- NumOutBounds[[d]][[ModelType]]$Factors$INF[,i,b]
            UU <- NumOutBounds[[d]][[ModelType]]$Factors$SUP[,i,b]
            MM <- NumOutBounds[[d]][[ModelType]]$Factors$MED[,i,b]
            PP <- NumOutPE[[LabIRF[d]]][[ModelType]]$Factors$Ortho[,i,b ] # Point estimate

            ALL <- data.frame(cbind(LL,MM,PP, UU))
            TimeSpan <- 1:HorizNumOut[[nn]]
            ALL$TimeSpan <- TimeSpan

            nmResponse <- DimLabelsFac[[2]] # names of the "response" factor
            nmShock <- DimLabelsFac[[3]] # names of the shock

            p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
              geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) + theme_classic() +
              theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                    axis.title.x=element_blank(), axis.title.y=element_blank() ) +
              ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

            plot_list[[i]] <- p
          }
          subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabIRF[d],"Factors_shock_to_", nmShock[b],"ORTHO", ".png"),
                          path= PathsGraphs[[ModelType]][[LabIRF[d]]][["Factors Ortho"]]))
        }
      }

  }

  ########################################  Yields ###############################################################

  GraphsBinarYields <- c(InputsForOutputs[[ModelType]]$IRF$WishGraphsOrtho$YieldsBootstrap,
                          InputsForOutputs[[ModelType]]$GIRF$WishGraphsOrtho$YieldsBootstrap)

  IdxYieldssGraphs <- which(GraphsBinarYields == 1)


  if (any(GraphsBinarYields==1)){


    plot_list <- list()


      for (d in IdxYieldssGraphs){

        # Folder Creation
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabIRF[d], sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabIRF[d], "/Yields", sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabIRF[d], "/Yields/Ortho", sep=""))

        for (b in 1:K){
          for (i in 1:CJ){
            LL <- NumOutBounds[[d]][[ModelType]]$Yields$INF[,i,b]
            UU <- NumOutBounds[[d]][[ModelType]]$Yields$SUP[,i,b]
            MM <- NumOutBounds[[d]][[ModelType]]$Yields$MED[,i,b]
            PP <- NumOutPE[[LabIRF[d]]][[ModelType]]$Yields$Ortho[,i,b ] # Point estimate

            ALL <- data.frame(cbind(LL,MM,PP, UU))
            TimeSpan <- 1:HorizNumOut[[nn]]
            ALL$TimeSpan <- TimeSpan


            nmResponse <- DimLabelsYies[[2]] # names of the "response" factor
            nmShock <- DimLabelsYies[[3]] # names of the shock

            p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
              geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) + theme_classic() +
              theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                    axis.title.x=element_blank(), axis.title.y=element_blank() ) +
              ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

            plot_list[[i]] <- p
          }
          subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabIRF[d],"Yields_shock_to_", nmShock[b], "ORTHO", ".png"),
                          path= PathsGraphs[[ModelType]][[LabIRF[d]]][["Yields Ortho"]]))
        }
      }



  }

  return(NumOutBounds)

}


##############################################################################################################
#' Creates the confidence bounds and the graphs of FEVDs and GFEVDs after bootstrap (JLL-based models)
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelBootstrap list containing the complete set of model parameters after bootstrap estimation procedure
#'@param NumOutPE  list of model parameter point estimates
#'@param InputsForOutputs list conataining the desired inputs for the construction of the outputs of interest
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'@param PathsGraphs path of the folder in which the graphs will be saved
#'
#'
#'@keywords internal


FEVDandGFEVDbs_jointQ_Ortho <- function(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs){

  ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
  pctg <-   InputsForOutputs[[ModelType]]$Bootstrap$pctg

  # function
  pctg_inf <- (100-pctg)/2
  pctg_sup <- 100 - (100-pctg)/2
  quants <- c(pctg_inf, 50, pctg_sup)/100 # Desired quantiles


  # initializarion
  LabFEVD <- c("FEVD","GFEVD")
  OutNames <- names(ModelBootstrap$NumOutDraws)
  C <- length(Economies)
  J <-  length(ModelBootstrap$GeneralInputs$mat)
  CJ <- C*J


  HorizNumOut <- c(InputsForOutputs[[ModelType]]$IRF$horiz, InputsForOutputs[[ModelType]]$FEVD$horiz,
                   InputsForOutputs[[ModelType]]$GIRF$horiz, InputsForOutputs[[ModelType]]$GFEVD$horiz)

  NumOutBounds <- list()


    for (nn in match(LabFEVD, OutNames) ){

      K <- nrow(ModelBootstrap$ParaDraws[[ModelType]][[1]]$ests$K1Z)

      ############################################# Factors ######################################################################

      INFfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # lower bound
      MEDfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # Median
      SUPfacs <- array(NA, c(HorizNumOut[[nn]], K,K)) # upper bound

      DimLabelsFac <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Factors$Ortho)
      dimnames(INFfacs) <- DimLabelsFac
      dimnames(MEDfacs) <- DimLabelsFac
      dimnames(SUPfacs) <- DimLabelsFac


      # Allocation
      AllShocksOnePeriodFacs <- array(NA, c(ndraws, K, K))
      Facs <- matrix(NA, nrow = ndraws, ncol = K)

      for (thor in 1:HorizNumOut[[nn]]){
        for (h in 1:K){ # loop through the shocks
          for (g in 1:ndraws){ # # loop through the draws
            Facs[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Factors$Ortho[thor,,h] # All responses to one shock in one horizon
          }

          AllShocksOnePeriodFacs[,,h]  <- apply(Facs,2, sort) # Ensures that each column is in ascending order

          INFfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[1])
          MEDfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[2])
          SUPfacs[thor, ,h] <- apply(AllShocksOnePeriodFacs[,,h], 2, stats::quantile, probs = quants[3])
        }
      }

      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$INF <- INFfacs
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$MED <- MEDfacs
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Factors$SUP <- SUPfacs


      ############################################# Yields ######################################################################
      INFyies <- array(NA, c(HorizNumOut[[nn]], K, CJ)) # lower bound
      MEDyies <- array(NA, c(HorizNumOut[[nn]], K, CJ)) # Median
      SUPyies <- array(NA, c(HorizNumOut[[nn]], K, CJ)) # upper bound

      DimLabelsYies <- dimnames(ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[1]]$Yields$Ortho)
      dimnames(INFyies) <- DimLabelsYies
      dimnames(MEDyies) <- DimLabelsYies
      dimnames(SUPyies) <- DimLabelsYies

      #Allocation
      AllShocksOnePeriodyies <- array(NA, c(ndraws, K, CJ))
      yies <- matrix(NA, nrow = ndraws, ncol = K)

      for (thor in 1:HorizNumOut[[nn]]){
        for (h in 1:CJ){
          for (g in 1:ndraws){
            yies[g,] <- ModelBootstrap$NumOutDraws[[OutNames[nn]]][[ModelType]][[g]]$Yields$Ortho[thor, ,h] # All responses to one shock in one horizon
          }
          AllShocksOnePeriodyies[,,h]  <- apply(yies,2, sort) # Ensures that each column is in ascending order

          INFyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[1])
          MEDyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[2])
          SUPyies[thor, ,h] <- apply(AllShocksOnePeriodyies[,,h], 2, stats::quantile, probs = quants[3])
        }
      }


      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$INF <- INFyies
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$MED <- MEDyies
      NumOutBounds[[OutNames[nn]]][[ModelType]]$Yields$SUP <- SUPyies

    }


  #################################################################################################################
  ###################################### PLOTS ####################################################################
  #################################################################################################################

  ########################################  Factors ###############################################################

  GraphsBinarFactors <- c(InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$RiskFactorsBootstrap,
                         InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$RiskFactorsBootstrap)

  IdxFactorsGraphs <- which(GraphsBinarFactors == 1)



  if (any(GraphsBinarFactors==1)){

    cat(' ** FEVDs/GFEVDs-Ortho (Bootstrap) \n' )

    plot_list <- list()

      for (d in IdxFactorsGraphs){
        # Folder Creation
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabFEVD[d], sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabFEVD[d], "/Factors", sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabFEVD[d], "/Factors/Ortho", sep=""))


        for (b in 1:K){
          for (i in 1:K){
            LL <- NumOutBounds[[d]][[ModelType]]$Factors$INF[,i,b]
            UU <- NumOutBounds[[d]][[ModelType]]$Factors$SUP[,i,b]
            MM <- NumOutBounds[[d]][[ModelType]]$Factors$MED[,i,b]
            PP <- NumOutPE[[LabFEVD[d]]][[ModelType]]$Factors$Ortho[,i,b ] # Point estimate

            ALL <- data.frame(cbind(LL,MM,PP, UU))
            TimeSpan <- 1:HorizNumOut[[nn]]
            ALL$TimeSpan <- TimeSpan

            nmResponse <- DimLabelsFac[[2]] # names of the "response" factor
            nmShock <- DimLabelsFac[[3]] # names of the shock

            p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
              geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) + theme_classic() +
              theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                    axis.title.x=element_blank(), axis.title.y=element_blank() ) +
              ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

            plot_list[[i]] <- p
          }
          subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabFEVD[d],"Factors_", nmShock[b], "ORTHO", ".png"),
                          path= PathsGraphs[[ModelType]][[LabFEVD[d]]][["Factors Ortho"]]))
        }
      }

  }

  ########################################  Yields ###############################################################

  GraphsBinarYields <- c(InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$YieldsBootstrap,
                         InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$YieldsBootstrap)

  IdxYieldsGraphs <- which(GraphsBinarYields == 1)



    if (any(GraphsBinarYields==1)){

    plot_list <- list()


      for (d in IdxYieldsGraphs){
        # Folder Creation
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabFEVD[d], sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabFEVD[d], "/Yields", sep=""))
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Bootstrap/", LabFEVD[d], "/Yields/Ortho", sep=""))


        for (b in 1:CJ){
          for (i in 1:K){
            LL <- NumOutBounds[[d]][[ModelType]]$Yields$INF[,i,b]
            UU <- NumOutBounds[[d]][[ModelType]]$Yields$SUP[,i,b]
            MM <- NumOutBounds[[d]][[ModelType]]$Yields$MED[,i,b]
            PP <- NumOutPE[[LabFEVD[d]]][[ModelType]]$Yields$Ortho[,i,b ] # Point estimate

            ALL <- data.frame(cbind(LL,MM,PP, UU))
            TimeSpan <- 1:HorizNumOut[[nn]]
            ALL$TimeSpan <- TimeSpan

            nmResponse <- DimLabelsYies[[2]] # names of the "response" factor
            nmShock <- DimLabelsYies[[3]] # names of the shock

            p <- ggplot(data = ALL, aes(x= TimeSpan )) +  geom_line(aes(y = LL), color = 'blue') + geom_line(aes(y = MM), color = 'green') +
              geom_line(aes(y = UU), color = 'red') +  geom_line(aes(y = PP)) + theme_classic() +
              theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                    axis.title.x=element_blank(), axis.title.y=element_blank() )    +
              ggtitle( nmResponse[i]) + geom_hline(yintercept=0)

            plot_list[[i]] <- p
          }
          subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3)
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(LabFEVD[d],"Yields_", nmShock[b], "ORTHO", ".png"),
                          path= PathsGraphs[[ModelType]][[LabFEVD[d]]][["Yields Ortho"]]))
        }
      }


  }

  return(NumOutBounds)

}

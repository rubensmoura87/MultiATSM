#' Generate the graphical outputs for the selected models (Point estimate)
#'
#'@param ModelType A character vector indicating the model type to be estimated.
#'@param ModelPara List of model parameter estimates (See the "Optimization" function)
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, GFEVDs and Term Premia
#'@param InputsForOutputs list containing the desired inputs for the construction of the desired output
#'@param Economies A character vector containing the names of the economies included in the system.
#'@param FactorLabels A list of character vectors with labels for all variables in the model.
#'
#'@keywords internal



GraphicalOutputs <- function(ModelType, ModelPara, NumOut, InputsForOutputs, Economies, FactorLabels){

# Generate the graph paths and the graph folders
  # Create the folder output folder in which the outputs will be stored
  dir.create(paste(tempdir(), "/Outputs", sep=""))
  dir.create(paste(tempdir(), "/Outputs/",  ModelType, sep=""))
  dir.create(paste(tempdir(), "/Outputs/",  ModelType, "/Point Estimate", sep=""))


  PathsGraphs <- FolderCreationPoint(ModelType, Economies)

# 1) Plot the set of risk factors
  RiskFactorsGraphs(ModelType, ModelPara, Economies, FactorLabels)

  cat(" 2.3) Generating the graphs of interest \n")

  # 2) Model fit
  Fitgraphs(ModelType, InputsForOutputs[[ModelType]]$Fit$WishGraphs, ModelPara, NumOut, Economies, PathsGraphs)

  if (any(ModelType == c("JPS original", "JPS global",  "GVAR single"))){
  for (i in 1:length(Economies)){dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], sep=""))}}

  # 3) IRF and GIRF
  if (any(ModelType == c("JLL original",  "JLL No DomUnit",  "JLL joint Sigma"))){
    OutType <- c("IRF", "GIRF","IRF Ortho", "GIRF Ortho") } else{ OutType <- c("IRF", "GIRF")}


  for(j in 1:length(OutType)){
  WG <- Wished_Graphs_IRFandGIRF(InputsForOutputs, OutType[j], ModelType)
  IRFandGIRFgraphs(ModelType, NumOut, WG$RiskGraphs, WG$YieldGraphs, InputsForOutputs[[ModelType]]$IRF$horiz,
                   PathsGraphs, OutputType = OutType[j], Economies)
     }

# 1) Models for which the estimation is done on a country-by-country basis
if (any(ModelType == c("JPS original", "JPS global",  "GVAR single"))){

  # FEVD and GFEVD
  FEVDgraphsSep(ModelType, NumOut, InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactors,
                InputsForOutputs[[ModelType]]$FEVD$WishGraphs$Yields, InputsForOutputs[[ModelType]]$FEVD$horiz,
                PathsGraphs, Economies)

  GFEVDgraphsSep(ModelType, NumOut, InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactors,
                 InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$Yields, InputsForOutputs[[ModelType]]$GFEVD$horiz,
                 PathsGraphs, Economies)

  # Term premia decomposition
  TPDecompGraphSep(ModelType, NumOut, ModelPara, InputsForOutputs[[ModelType]]$RiskPremia$WishGraphs,
                   InputsForOutputs$UnitMatYields, Economies, PathsGraphs)
  }

# 2) Models for which the estimation is done on jointly for all the countries
if ( any(ModelType == c("GVAR multi", "JPS multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma"))){

# FEVD and GFEVD
FEVDgraphsJoint(ModelType, NumOut, InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactors,
                InputsForOutputs[[ModelType]]$FEVD$WishGraphs$Yields, InputsForOutputs[[ModelType]]$FEVD$horiz, PathsGraphs)

GFEVDgraphsJoint(ModelType, NumOut, InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactors,
                InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$Yields, InputsForOutputs[[ModelType]]$GFEVD$horiz, PathsGraphs)


TPDecompGraphJoint(ModelType, NumOut, ModelPara, InputsForOutputs[[ModelType]]$RiskPremia$WishGraphs,
                   InputsForOutputs$UnitMatYields, Economies, PathsGraphs)

# 2.1) JLL-based models
if (any(ModelType == c("JLL original",  "JLL No DomUnit",  "JLL joint Sigma"))){
#  FEVD and GFEVD

FEVDgraphsJLLOrtho(ModelType, NumOut, InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$RiskFactors,
                  InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$Yields, InputsForOutputs[[ModelType]]$FEVD$horiz,
                  PathsGraphs)

GFEVDgraphsJLLOrtho(ModelType, NumOut, InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$RiskFactors,
                    InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$Yields, InputsForOutputs[[ModelType]]$GFEVD$horiz,
                    PathsGraphs)

}

}

  cat(paste("Desired graphs are saved in your temporary directory. Please, check:",tempdir(), "\n\n"))
}


#########################################################################################################
#################################### RISK FACTORS GRAPHS ################################################
#########################################################################################################
#'  Spanned and unspanned factors plot
#'
#'@param ModelType string-vector containing the label of the model to be estimated
#'@param ModelOutputs list of model parameter estimates (see the "Optimization" function)
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'
#'@importFrom ggplot2 ggplot theme_classic scale_x_date element_rect
#'
#'
#'@keywords internal

RiskFactorsGraphs <- function(ModelType, ModelOutputs, Economies, FactorLabels){


  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)

  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
    X <- c()
    for (i in 1:C){
      if (i ==1){
        X <- ModelOutputs[[ModelType]][[Economies[i]]]$inputs$AllFactors
        Factors <- X
      } else{
       IdxNonG <- if (G > 0) -seq_len(G) else TRUE
        X <- ModelOutputs[[ModelType]][[Economies[i]]]$inputs$AllFactors[IdxNonG , ]
        Factors <- rbind(Factors, X)
        }
         }
  }else{ Factors <- ModelOutputs[[ModelType]]$inputs$AllFactors }

  K <- nrow(Factors)
  T <- ncol(Factors)

  TitleGraphs <- c(FactorLabels$Global, FactorLabels$Domestic)
  Folder2save <-  paste(tempdir(), "/Outputs", "/", ModelType, sep="")

  # 0) Preliminary work
  Factors <- t(Factors)
  # Set a dataframe for the factors
  SS <- data.frame(Factors)
  TimeSpan <- 1:T
  SS$TimeSpan <- TimeSpan
  nmFactors <- colnames(SS)

  #Initialize list of plots
  plot_list <- list()
  plot_list_no_legend <- list()


  ######################################### Global Factors ############################################
  for (i in seq_len(G)){
    g <- ggplot(data = SS, aes(x= TimeSpan )) + geom_line(aes_string(y = (nmFactors[i])))
    g <- g + geom_hline(yintercept=0) + ggtitle(TitleGraphs[i])  + theme_classic()
    g <- g + theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                   axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=10))


    plot_list[[i]] <- g
  }
  plot_list_no_legend <- plot_list

  ######################################### Domestic Factors ##########################################
  # Index of the factors in the list of factors
  IDX <- matrix(NA, nrow= N+M, ncol= C)
  for (j in 1:(N+M)){
    for (i in 1:C){
      IDX[j,i] <- (G+j) +(N+M)*(i-1)
    }
  }

  # Plot the graphs
  for (j in 1:(N+M)){
    d <- ggplot(data = SS, aes(x= TimeSpan))
      for (i in 1:C){  d <- d +  geom_line(aes_string(y = nmFactors[IDX[j,i]], color = shQuote(Economies[i])))}

    d <- d + geom_hline(yintercept=0) + ggtitle(TitleGraphs[G+j]) + theme_classic()
    d <- d  +  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
                     axis.title.y=element_blank(),  axis.text.x = element_text(size=10))
    d <- d + labs(colour = "Legend")

    plot_list[[G+j]] <- d # All plots with legends

    plot_list_no_legend[[G+j]] <- d + theme(legend.position = "none")  # Hide legends from all plots
  }



  ##################################### All factors in a subplot ########################################
  # Legend settings:
  LegendSubPlot <- cowplot::get_legend(plot_list[[G+1]] +
                                theme(legend.direction = "horizontal", legend.position="bottom", legend.justification="center",
                                      legend.title = element_text(size=10, face="bold"),
                                      legend.text = element_text( size = 8),
                                      legend.box.background = element_rect(colour = "black")) )

  # Build the graph subplots:
  subplots <- cowplot::plot_grid(plotlist= plot_list_no_legend, ncol=3)

  #Create a new subplot with graphs + legend:
  FactorsGraph <- cowplot::plot_grid(LegendSubPlot, subplots , ncol=1,  rel_heights = c(.1, 1))

  # Save final graph at the overleaf folder:
  suppressMessages(ggplot2::ggsave(FactorsGraph, filename =paste0("RiskFactors", ".png"), path =Folder2save))

}


######################################################################################################
##################################### 1) FIT ########################################################
#####################################################################################################
#' Model fit graphs for ("sep Q" models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param WishFitgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param ModelPara List of model parameter estimates (See the "Optimization" function)
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'@importFrom ggplot2 ggplot theme_classic scale_x_date element_rect
#'
#'
#'@keywords internal


Fitgraphs  <- function(ModelType, WishFitgraphs, ModelPara, NumOut, Economies, PathsGraphs){


  if (WishFitgraphs== 0){cat("No graphs for bond yields fit were generated \n")
  }else{

    cat(' ** Bond yields Fit \n' )

    C <- length(Economies)

    # 1) Models estimated individually
    if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
    for( i in 1:C){
      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], "/Fit", sep=""))


    mat <- (ModelPara[[ModelType]][[Economies[i]]]$inputs$mat)*12
    J <- length(mat)

    YieldData <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y
    ModelFit <- NumOut$Fit[[ModelType]][[Economies[i]]]$`Yield Fit`
    ModelImplied <- NumOut$Fit[[ModelType]][[Economies[i]]]$`Yield Model Implied`

    # Make graphs
    Fit_Subplot(YieldData, ModelFit, ModelImplied, J, mat, Economies[i], ModelType, PathsGraphs)

}
    } else{

      # 2) Models estimated jointly
      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Fit", sep=""))  # Folder creation
      mat <- (ModelPara[[ModelType]]$inputs$mat)*12
      J <- length(mat)

      Idx0 <- 0
      for (i in 1:C){
        Idx <- Idx0
        IdxGrpahs <- (Idx+1):(Idx+J)

        YieldData <- ModelPara[[ModelType]]$inputs$Y[IdxGrpahs, ]
        ModelFit <- NumOut$Fit$`Yield Fit`[IdxGrpahs, ]
        ModelImplied <- NumOut$Fit$`Yield Model Implied`[IdxGrpahs, ]

        # Make graphs
        Fit_Subplot(YieldData, ModelFit, ModelImplied, J, mat, Economies[i], ModelType, PathsGraphs)

        Idx0 <- (Idx+J)
      }
  }
  }
}


######################################################################################################
##################################### 2) IRF and GIRF ################################################
#####################################################################################################
#' IRFs graphs for all models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param IRFhoriz  single numerical vector containing the desired horizon of analysis for the IRFs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'@param OutputType Available options are 'IRF' and 'GIRF'
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'@importFrom ggplot2 ggplot geom_line labs geom_hline theme ggtitle theme_classic
#'
#'@keywords internal


IRFandGIRFgraphs <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, IRFhoriz,
                             PathsGraphs, OutputType, Economies) {

  if (WishPdynamicsgraphs==0 & WishYieldsgraphs ==0){ cat(paste( "No graphs for", OutputType, "were generated \n"))

    if(any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))
       && any(ModelType == c("IRF Ortho","GIRF Ortho"))){
      cat(paste("No graphs for" , OutputType, "were generated (orthogonolized version) \n"))}

    }else {

    if(OutputType == "IRF"){  cat(' ** IRFs \n');
      Lab_Fac <-  "IRFFactors_shock_to_"; Lab_Yield <- "IRFYields_shock_to_"
    } else if( OutputType == "GIRF"){ cat(' ** GIRFs \n');
      Lab_Fac <-  "GIRFFactors_shock_to_"; Lab_Yield <- "GIRFYields_shock_to_"
    } else if(OutputType == "IRF Ortho"){  cat(' ** IRFs-Ortho \n')
      Lab_Fac <-  "IRFFactors_shock_to_"; Lab_Yield <- "IRFYields_shock_to_"
    } else{ cat(' ** GIRFs-Ortho \n')
      Lab_Fac <-  "GIRFFactors_shock_to_" ; Lab_Yield <- "GIRFYields_shock_to_"}


    ################ 1) Estimation done for countries individually ################
    if (any(ModelType == c("JPS original", "JPS global",  "GVAR single"))){
    K <- dim(NumOut$IRF[[ModelType]][[Economies[1]]]$Factors)[2] # Total number of risk factors
    C <- length(Economies)

    for (i in 1:C){
      # a) Folder Creation
      FolderPrep_IRFs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies[i], ModelType)

      # b) Recast IRFs or GIRFs
      IRFset <- BuildIRFlist(NumOut, Economies[i], ModelType, IRFhoriz, K, OutputType)
      nmFactors <- names(IRFset$IRFFactors[[1]])[-1] # Factor names
      nmYields <- names(IRFset$IRFYields[[1]])[-1] # Factor names

      # c) Graph of Factors
      if (WishPdynamicsgraphs == 1){
      PathAdj <- AdjustPathIRFs(OutputType, "Factors", PathsGraphs, Economies[i], ModelType)
      for(h in 1:K){
        p_list <- IRFandGIRFs_Format_Fac(IRFset$IRFFactors[[h]])
        subplots <- cowplot::plot_grid(plotlist= p_list, ncol = 3) # Gather the graphs in sub-plots
        suppressMessages(ggplot2::ggsave(subplots, file=paste0(Lab_Fac, nmFactors[[h]],"_Model_", Economies[i],".png"),
                                           path= PathAdj)) # Save sub-plots
        }
      }

      # d) Graph of yields
      if (WishYieldsgraphs == 1){
      PathAdj <- AdjustPathIRFs(OutputType, "Yields", PathsGraphs, Economies[i], ModelType)
      for(l in 1:K){
      p_list <- IRFandGIRFs_Format_Yields(IRFset$IRFYields[[l]])
      subplots <- cowplot::plot_grid(plotlist= p_list, ncol=3) # Gather the graphs in sub-plots
      suppressMessages(ggplot2::ggsave(subplots, file=paste0(Lab_Yield, nmFactors[[l]],"_Model_", Economies[i], ".png"),
                                         path= PathAdj)) # Save sub-plots
      }
  }
}
  } else {

    ################ 2) Estimation done for countries jointly ###############################
    if (any(OutputType == c("IRF", "GIRF"))){
    # Total number of risk factors
    if(any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){ K <- dim(NumOut$IRF[[ModelType]]$Factors$Ortho)[2]
    }else{K <- dim(NumOut$IRF[[ModelType]]$Factors)[2]}

      # a) Folder Creation
      FolderPrep_IRFs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType)

      # b) Recast IRFs or GIRFs
      IRFset <- BuildIRFlist(NumOut, Economies, ModelType, IRFhoriz, K, OutputType)

      nmFactors <- names(IRFset$IRFFactors[[1]])[-1] # Factor names
      nmYields <- names(IRFset$IRFYields[[1]])[-1] # Factor names

      # c) Graph of Factors
      if (WishPdynamicsgraphs == 1){
      PathAdj <- AdjustPathIRFs(OutputType, "Factors", PathsGraphs, Economies, ModelType)
      for(h in 1:K){
        p_list <- IRFandGIRFs_Format_Fac(IRFset$IRFFactors[[h]])
        subplots <- cowplot::plot_grid(plotlist= p_list, ncol=3) # Gather the graphs in sub-plots
        suppressMessages(ggplot2::ggsave(subplots, file=paste0(Lab_Fac, nmFactors[[h]],".png"),
                                         path= PathAdj)) # Save sub-plots
      }
      }

      # d) Graph of yields
      if (WishYieldsgraphs == 1){
        PathAdj <- AdjustPathIRFs(OutputType, "Yields", PathsGraphs, Economies, ModelType)
        for(l in 1:K){
          p_list <- IRFandGIRFs_Format_Yields(IRFset$IRFYields[[l]])
          subplots <- cowplot::plot_grid(plotlist= p_list, ncol=3) # Gather the graphs in sub-plots
          suppressMessages(ggplot2::ggsave(subplots, file=paste0(Lab_Yield, nmFactors[[l]], ".png"), path= PathAdj)) # Save sub-plots
        }
        }
    }
      ###################### 3) JLL orthoonalized ############################################################
      if(any(OutputType == c("IRF Ortho", "GIRF Ortho"))){

        K <- dim(NumOut$IRF[[ModelType]]$Factors$Ortho)[2]

        # a) Folder Creation
        FolderPrep_IRFs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType)

        # b) Recast IRFs or GIRFs
        IRFset <- BuildIRFlist(NumOut, Economies, ModelType, IRFhoriz, K, OutputType)

        nmFactors <- names(IRFset$IRFFactors[[1]])[-1] # Factor names
        nmYields <- names(IRFset$IRFYields[[1]])[-1] # Factor names

        # a.2) Graph of Factors
        if (WishPdynamicsgraphs == 1){
          PathAdj <- AdjustPathIRFs(OutputType, "Factors", PathsGraphs, Economies, ModelType)
          for(h in 1:K){
            p_list <- IRFandGIRFs_Format_Fac(IRFset$IRFFactors[[h]])
            subplots <- cowplot::plot_grid(plotlist= p_list, ncol=3) # sub-plots
            suppressMessages(ggplot2::ggsave(subplots, file=paste0(Lab_Fac, nmFactors[[h]],"ORTHO",".png"),
                                             path= PathAdj))
          }
        }

        # d) Graph of yields
        if (WishYieldsgraphs == 1){
          PathAdj <- AdjustPathIRFs(OutputType, "Yields", PathsGraphs, Economies, ModelType)
          for(l in 1:K){
            p_list <- IRFandGIRFs_Format_Yields(IRFset$IRFYields[[l]])
            subplots <- cowplot::plot_grid(plotlist= p_list, ncol=3) # Gather the graphs in sub-plots
            suppressMessages(ggplot2::ggsave(subplots, file=paste0(Lab_Yield, nmFactors[[l]], "ORTHO", ".png"),
                                             path= PathAdj)) # Save sub-plots
          }
        }
        }
      }
    }
}

######################################################################################################
######################################################################################################
####################### OUTPUTS FOR MODELS IN WHICH THE ESTIMATION ###################################
########################    IS DONE ON A COUNTRY-BY-COUNTRY BASIS       #############################
######################################################################################################
######################################################################################################

######################################################################################################
##################################### 3) FEVD ########################################################
#####################################################################################################
#' FEVDs graphs for ("sep Q" models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param FEVDhoriz  single numerical vector conataining the desired horizon of analysis for the FEVDs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'@importFrom ggplot2 ggplot theme_classic aes element_text theme labs ggtitle element_blank aes_string geom_bar geom_hline
#'@importFrom zoo index
#'
#'@keywords internal



FEVDgraphsSep <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, FEVDhoriz,
                          PathsGraphs, Economies) {


  if (WishPdynamicsgraphs == 0 & WishYieldsgraphs ==0){cat("No graphs for FEVDs were generated \n")
  } else {

    cat(' ** FEVD \n' )

    C <- length(Economies)

      for (h in 1:C){

        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[h], "/FEVD", sep=""))
  ######################################### Factors #########################################################
      if (WishPdynamicsgraphs == 1){
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[h], "/FEVD/Factors", sep=""))

          FEVDFactors <- list()
          K <- dim(NumOut$FEVD[[ModelType]][[Economies[h]]]$Factors)[2] # Total number of risk factors

          for (i in 1:K){ # Recast the model in the data-frame format
            FEVDFactors[[i]] <- data.frame(NumOut$FEVD[[ModelType]][[Economies[h]]]$Factors[,,i])
            INDX <- cbind(index=rownames(FEVDFactors[[i]]), FEVDFactors[[i]])
            FEVDFactors[[i]] <- suppressMessages(reshape2::melt(INDX))
            FEVDFactors[[i]]$index <- factor(FEVDFactors[[i]]$index, levels= 1:(FEVDhoriz))
          }

          nmFactors <- colnames(NumOut$FEVD[[ModelType]][[Economies[h]]]$Factors) # Factor names
          ## Generate plots:
          for (i in 1:K){
            index <- FEVDFactors[[i]]$index
            value <- FEVDFactors[[i]]$value
            variable <- FEVDFactors[[i]]$variable

            p <- ggplot(FEVDFactors[[i]], aes(x= index, y=value , fill=variable)) +
              geom_bar(stat="identity", width = 0.25) +
              labs(title= paste0("FEVD_",nmFactors[i])) + theme_classic() +
            theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                  plot.title = element_text(size = 10, face = "bold", hjust = 0.5) )
            suppressMessages(ggplot2::ggsave(p, file=paste0("FEVDFactors_", nmFactors[[i]],".png"),
                            path= PathsGraphs[[ModelType]]$FEVD[[Economies[h]]][["Factors"]], width= 11, height= 7.98))
          }

}
          ############################################ Yields #########################################################
          if (WishYieldsgraphs == 1){
            dir.create( paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/Model ", Economies[h], "/FEVD/Yields", sep=""))

            FEVDYields <- list()
            J <- dim(NumOut$FEVD[[ModelType]][[Economies[h]]]$Yields)[3] # Total number of yields of the system

            for (i in 1:J){ # Recast the outputs to the data-frame format
              FEVDYields[[i]] <- data.frame(NumOut$FEVD[[ModelType]][[Economies[h]]]$Yields[,,i])
              INDX <- cbind(index=rownames(FEVDYields[[i]]), FEVDYields[[i]])
              FEVDYields[[i]] <- suppressMessages(reshape2::melt(INDX))
              FEVDYields[[i]]$index <- factor(FEVDYields[[i]]$index, levels= 1:(FEVDhoriz))
            }

            nmYields <- dimnames(NumOut$FEVD[[ModelType]][[Economies[h]]]$Yields)[[3]] # Yield names

            ## Generate plots:
            for (i in 1:J){
              index <- FEVDYields[[i]]$index
              value <- FEVDYields[[i]]$value
              variable <- FEVDYields[[i]]$variable

              p <- ggplot(FEVDYields[[i]], aes(x= index, y=value , fill=variable)) +
                geom_bar(stat="identity", width = 0.25) +
                labs(title= paste0("FEVD_", nmYields[i])) + theme_classic() +
                theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                      plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
              suppressMessages(ggplot2::ggsave(p, file=paste0("FEVDYields_", nmYields[i],".png"),
                              path= PathsGraphs[[ModelType]]$FEVD[[Economies[h]]][["Yields"]], width= 11, height= 7.98))
            }
          }
        }

    }
  }
######################################################################################################
##################################### 5) GFEVD ########################################################
#####################################################################################################
#' GFEVDs graphs for ("sep Q" models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param GFEVDhoriz  single numerical vector conataining the desired horizon of analysis for the GFEVDs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@importFrom ggplot2 ggplot aes element_text theme theme_classic labs ggtitle element_blank aes_string geom_bar geom_hline
#'@importFrom zoo index
#'
#'@keywords internal


GFEVDgraphsSep <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, GFEVDhoriz,
                           PathsGraphs, Economies) {

  if (WishPdynamicsgraphs==0 & WishYieldsgraphs==0){cat("No graphs for GFEVDs were generated \n")
    } else{

      cat(' ** GFEVDs \n' )

      C <- length(Economies)
      for (h in 1:C){

        dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[h], "/GFEVD", sep=""))

          ######################################### Factors #########################################################
          if (WishPdynamicsgraphs == 1){

            dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[h], "/GFEVD/Factors", sep=""))

            GFEVDFactors <- list()
          K <- dim(NumOut$GFEVD[[ModelType]][[Economies[h]]]$Factors)[2] # Total number of risk factors

          for (i in 1:K){ # Recast the outputs to the data-frame format
            GFEVDFactors[[i]] <- data.frame(NumOut$GFEVD[[ModelType]][[Economies[h]]]$Factors[,,i])
            GFEVDFactors[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(GFEVDFactors[[i]]), GFEVDFactors[[i]])))
            GFEVDFactors[[i]]$index <- factor(GFEVDFactors[[i]]$index, levels= 1:(GFEVDhoriz))
          }

          nmFactors <- colnames(NumOut$GFEVD[[ModelType]][[Economies[h]]]$Factors) # Factor names
          ## Generate plots:
          for (i in 1:K){
            index <- GFEVDFactors[[i]]$index
            value <- GFEVDFactors[[i]]$value
            variable <- GFEVDFactors[[i]]$variable

            p <- ggplot(GFEVDFactors[[i]], aes(x= index, y=value , fill=variable)) +
              geom_bar(stat="identity", width = 0.25) +
              labs(title= paste0("GFEVD_",nmFactors[i])) + theme_classic() +
              theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                    plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
            suppressMessages(ggplot2::ggsave(p, file=paste0("GFEVDFactors_", nmFactors[[i]],".png"),
                            path= PathsGraphs[[ModelType]]$GFEVD[[Economies[h]]][["Factors"]], width= 11, height= 7.98))
          }

}
          ############################################ Yields #########################################################
          if (WishYieldsgraphs == 1){

      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[h], "/GFEVD/Yields", sep=""))

            GFEVDYields <- list()
            J <- dim(NumOut$GFEVD[[ModelType]][[Economies[h]]]$Yields)[3] # Total number of yields of the system

            for (i in 1:J){ # Recast the outputs to the data-frame format
              GFEVDYields[[i]] <- data.frame(NumOut$GFEVD[[ModelType]][[Economies[h]]]$Yields[,,i])
              GFEVDYields[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(GFEVDYields[[i]]), GFEVDYields[[i]])))
              GFEVDYields[[i]]$index <- factor(GFEVDYields[[i]]$index, levels= 1:(GFEVDhoriz))
            }

            nmYields <- dimnames(NumOut$GFEVD[[ModelType]][[Economies[h]]]$Yields)[[3]] # Yield names

            ## Generate plots:
            for (i in 1:J){
              index <- GFEVDYields[[i]]$index
              value <- GFEVDYields[[i]]$value
              variable <- GFEVDYields[[i]]$variable

              p <- ggplot(GFEVDYields[[i]], aes(x= index, y=value , fill=variable)) +
                geom_bar(stat="identity", width = 0.25) +
                labs(title= paste0("GFEVD_", nmYields[i])) + theme_classic() +
                theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                      plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
              suppressMessages(ggplot2::ggsave(p, file=paste0("GFEVDYields_", nmYields[i],".png"),
                              path= PathsGraphs[[ModelType]]$GFEVD[[Economies[h]]][["Yields"]], width= 11, height= 7.98))
            }
          }
        }

    }
  }


#####################################################################################################################
######################################## 6) TERM PREMIA DECOMPOSITION ##################################################
#####################################################################################################################
#' Term Premia decomposition graphs for "joint Q" models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, GFEVDs and Risk premia
#'@param ModelPara list of model parameter estimates (See the "Optimization" function)
#'@param WishRPgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param UnitYields (i) "Month": if maturity of yields are expressed in months or
#'                  (ii) "Year": if maturity of yields are expressed in years
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'
#'@keywords internal


TPDecompGraphSep <- function(ModelType, NumOut, ModelPara, WishRPgraphs, UnitYields, Economies, PathsGraphs){



  if (WishRPgraphs ==0){cat("No graphs for term-premia were generated \n")
  } else {

    cat(' ** Term premia \n')


    # Preliminary work
    dt <- ModelPara[[ModelType]][[Economies[1]]]$inputs$dt
    mat <- ModelPara[[ModelType]][[Economies[1]]]$inputs$mat

    J <- length(mat)
    C <- length(Economies)
    T <- ncol(ModelPara[[ModelType]][[Economies[1]]]$inputs$Y)

    TPdecomp <- NumOut$TermPremiaDecomp


    if (UnitYields== "Month"){
      k <- 12
      YLab <- "Months"
    }
    if (UnitYields== "Year"){
      k <- 1
      YLab <- "Years"
    }
    matAdjUnit <- mat*k


    # Graph title, legends and axes
    Dates <- 1:T

    matRPLab <- paste(matAdjUnit, YLab)
    GraphLegend <- c("Data", "Expected Component", "Term Premium")


    # Prepare plots
    subplots_RP <- list()
    YieldData <- list()

    plot_list_RP <- list()
    plot_list_no_legend_RP <- list()

    for (i in 1:C){ # per country

      Yields <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y
      YieldData[[Economies[i]]] <- t(Yields)*100

      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], "/TermPremia", sep=""))  # Folder creation


      p <- c()

      for (h in 1:length(mat)){  # per maturity

        Data <- YieldData[[Economies[i]]][,h]
        ExpCom <- TPdecomp$RiskPremia[["Expected Component"]][[Economies[i]]][,h]
        TP <- TPdecomp$RiskPremia[["Term Premia"]][[Economies[i]]][,h]
        TimeSpan <- Dates

        DataGraph <- data.frame(Data, ExpCom, TP, TimeSpan)

        # Graph: Fit x Data
        p <- ggplot(DataGraph) + geom_line(aes(x=TimeSpan, y=Data, color = GraphLegend[1])) +
          geom_line(aes(x=TimeSpan, y = ExpCom, color = GraphLegend[2])) +
          geom_line(aes(x=TimeSpan, y = TP, color = GraphLegend[3])) +
          labs(color = "Legend") + ggtitle( matRPLab[h]) +  theme_classic() +
          theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
                axis.title.y=element_blank(), axis.text.x = element_text(size=6) ) +
          geom_hline(yintercept=0)


        plot_list_RP[[h]] <- p # All plots with legends

        plot_list_no_legend_RP[[h]] <- p + theme(legend.position = "none")  # Hide legends from all plots

      }

      # Common title for each model
      title <- cowplot::ggdraw() + cowplot::draw_label(Economies[i], fontface='bold', size=14)

      # Generate legend
      # Legend settings:
      LegendSubPlot_RP <- cowplot::get_legend(plot_list_RP[[1]] +
                                                theme(legend.direction = "horizontal", legend.position="bottom", legend.justification="center",
                                                      legend.title = element_text(size=10, face="bold"),
                                                      legend.text = element_text( size = 8),
                                                      legend.box.background = element_rect(colour = "black")) )


      # Build the graph subplots:
      subplots_CS <- cowplot::plot_grid(plotlist= plot_list_no_legend_RP, ncol=3)
      subplots2save <- cowplot::plot_grid(LegendSubPlot_RP, subplots_CS, ncol=1, rel_heights=c(0.2, 1))
      suppressMessages(ggplot2::ggsave(subplots2save, file=paste0("TermPremia_", Economies[i],"_", ModelType, ".png"),
                      path = PathsGraphs[[ModelType]]$TermPremia[[Economies[i]]]))

    }

  }
}
######################################################################################################
######################################################################################################
####################### OUTPUTS FOR MODELS IN WHICH THE ESTIMATION ###################################
########################       IS DONE FOR ALL COUNTRIES JOINTLY   ###################################
########################     (FOR THE JLL MODELS, THE OUTPUTS ARE  ###################################
########################      GENERATED FOR THE NON-ORTHOGONALIZED ###################################
######################              FACTORS ONLY)                ###################################
######################################################################################################
######################################################################################################

######################################################################################################
##################################### 3) FEVD ########################################################
#####################################################################################################
#' FEVDs graphs for ("joint Q" models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param FEVDhoriz  single numerical vector conataining the desired horizon of analysis for the FEVDs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'@importFrom ggplot2 ggplot theme_classic aes element_text theme labs ggtitle element_blank aes_string geom_bar geom_hline
#'@importFrom zoo index
#'
#'@keywords internal


FEVDgraphsJoint <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, FEVDhoriz, PathsGraphs) {



    if (WishPdynamicsgraphs ==0 & WishYieldsgraphs==0){cat("No graphs for FEVDs were generated \n")
      }else {


        cat(' ** FEVDs  \n' )

        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD", sep="")) # Folder Creation

 ######################################### Factors #########################################################
        if (WishPdynamicsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Factors", sep="")) # Folder Creation

          FEVDFactors <- list()


          if ( any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
            K <- dim(NumOut$FEVD[[ModelType]]$Factors$NonOrtho)[2] # Total number of risk factors
            for (i in 1:K){
              FEVDFactors[[i]] <- data.frame(NumOut$FEVD[[ModelType]]$Factors$NonOrtho[,,i])
              FEVDFactors[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDFactors[[i]]), FEVDFactors[[i]])))
              FEVDFactors[[i]]$index <- factor(FEVDFactors[[i]]$index, levels= 1:(FEVDhoriz))
              nmFactors <- colnames(NumOut$FEVD[[ModelType]]$Factors$NonOrtho) # Factor names
            }
          }else{
            K <- dim(NumOut$FEVD[[ModelType]]$Factors)[2] # Total number of risk factors
            for (i in 1:K){

              FEVDFactors[[i]] <- data.frame(NumOut$FEVD[[ModelType]]$Factors[,,i])
              FEVDFactors[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDFactors[[i]]), FEVDFactors[[i]])))
              FEVDFactors[[i]]$index <- factor(FEVDFactors[[i]]$index, levels= 1:(FEVDhoriz))
            }
            nmFactors <- colnames(NumOut$FEVD[[ModelType]]$Factors) # Factor names
          }

          ## Generate plots:
          for (i in 1:K){
            index <- FEVDFactors[[i]]$index
            value <- FEVDFactors[[i]]$value
            variable <- FEVDFactors[[i]]$variable


            p <- ggplot(FEVDFactors[[i]], aes(x= index, y=value , fill=variable)) +
              geom_bar(stat="identity", width = 0.25) +
              labs(title= paste0("FEVD_",nmFactors[i])) + theme_classic() +
              theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                    plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
            suppressMessages(ggplot2::ggsave(p, file=paste0("FEVDFactors_", nmFactors[[i]],".png"),
                            path= PathsGraphs[[ModelType]]$FEVD[["Factors"]], width= 11, height= 7.98))
          }

        }
        ############################################ Yields #########################################################
        if (WishYieldsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Yields", sep="")) # Folder Creation

          FEVDYields <- list()

          if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){

            CJ <- dim(NumOut$FEVD[[ModelType]]$Yields$NonOrtho)[3] # Total number of yields of the system

            for (i in 1:CJ){ # Store data in data-frame format
              FEVDYields[[i]] <- data.frame(NumOut$FEVD[[ModelType]]$Yields$NonOrtho[,,i])
              FEVDYields[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDYields[[i]]), FEVDYields[[i]])))
              FEVDYields[[i]]$index <- factor(FEVDYields[[i]]$index, levels= 1:(FEVDhoriz))
            }

            nmYields <- dimnames(NumOut$FEVD[[ModelType]]$Yields$NonOrtho)[[3]] # Yield names

          }else{
            CJ <- dim(NumOut$FEVD[[ModelType]]$Yields)[3] # Total number of yields of the system

            for (i in 1:CJ){ # Store data in data-frame format
              FEVDYields[[i]] <- data.frame(NumOut$FEVD[[ModelType]]$Yields[,,i])
              FEVDYields[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDYields[[i]]), FEVDYields[[i]])))
              FEVDYields[[i]]$index <- factor(FEVDYields[[i]]$index, levels= 1:(FEVDhoriz))
            }

            nmYields <- dimnames(NumOut$FEVD[[ModelType]]$Yields)[[3]] # Yield names
          }
          ## Generate plots:
          for (i in 1:CJ){
            index <- FEVDYields[[i]]$index
            value <- FEVDYields[[i]]$value
            variable <- FEVDYields[[i]]$variable


            p <- ggplot(FEVDYields[[i]], aes(x= index, y=value , fill=variable)) +
              geom_bar(stat="identity", width = 0.25) +
              labs(title= paste0("FEVD_", nmYields[i])) + theme_classic() +
              theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                    plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
            suppressMessages(ggplot2::ggsave(p, file=paste0("FEVDYields_", nmYields[i],".png"),
                            path= PathsGraphs[[ModelType]]$FEVD[["Yields"]], width= 11, height= 7.98))
          }
        }


  }
}

######################################################################################################
##################################### 5) GFEVD ########################################################
#####################################################################################################
#' GFEVDs graphs for "joint Q" models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param GFEVDhoriz  single numerical vector conataining the desired horizon of analysis for the GFEVDs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'
#'@importFrom ggplot2 ggplot aes element_text theme theme_classic labs ggtitle element_blank aes_string geom_bar geom_hline
#'@importFrom zoo index
#'
#'@keywords internal



GFEVDgraphsJoint <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, GFEVDhoriz, PathsGraphs){


    if (WishPdynamicsgraphs ==0 & WishYieldsgraphs ==0){cat("No graphs for GFEVDs were generated \n")
      } else {


        cat(' ** GFEVD \n' )


        dir.create( paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/GFEVD", sep="")) # Folder Creation
        ######################################### Factors #########################################################
        if (WishPdynamicsgraphs == 1){
          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Factors", sep="")) # Folder Creation

          GFEVDFactors <- list()


          if ( any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
            K <- dim(NumOut$GFEVD[[ModelType]]$Factors$NonOrtho)[2] # Total number of risk factors
            for (i in 1:K){ # Outputs in data-frame
              GFEVDFactors[[i]] <- data.frame(NumOut$GFEVD[[ModelType]]$Factors$NonOrtho[,,i])
              GFEVDFactors[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(GFEVDFactors[[i]]), GFEVDFactors[[i]])))
              GFEVDFactors[[i]]$index <- factor(GFEVDFactors[[i]]$index, levels= 1:(GFEVDhoriz))
            }
            nmFactors <- colnames(NumOut$GFEVD[[ModelType]]$Factors$NonOrtho) # Factor names
          }else{
            K <- dim(NumOut$GFEVD[[ModelType]]$Factors)[2] # Total number of risk factors
            for (i in 1:K){ # Outputs in data-frame
              GFEVDFactors[[i]] <- data.frame(NumOut$GFEVD[[ModelType]]$Factors[,,i])
              GFEVDFactors[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(GFEVDFactors[[i]]), GFEVDFactors[[i]])))
              GFEVDFactors[[i]]$index <- factor(GFEVDFactors[[i]]$index, levels= 1:(GFEVDhoriz))
            }
            nmFactors <- colnames(NumOut$GFEVD[[ModelType]]$Factors) # Factor names
          }

          ## Generate plots:
          for (i in 1:K){

            index <- GFEVDFactors[[i]]$index
            value <- GFEVDFactors[[i]]$value
            variable <- GFEVDFactors[[i]]$variable

            p <- ggplot(GFEVDFactors[[i]], aes(x= index, y=value , fill=variable)) +
              geom_bar(stat="identity", width = 0.25) +
              labs(title= paste0("GFEVD_",nmFactors[i])) + theme_classic() +
              theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                    plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
            suppressMessages(ggplot2::ggsave(p, file=paste0("GFEVDFactors_", nmFactors[[i]],".png"),
                            path= PathsGraphs[[ModelType]]$GFEVD[["Factors"]], width= 11, height= 7.98))
          }

        }
        ############################################ Yields #########################################################
        if (WishYieldsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Yields", sep="")) # Folder Creation

          GFEVDYields <- list()

          if (ModelType == "JLL original"  || ModelType== "JLL No DomUnit" || ModelType =="JLL joint Sigma"){
            CJ <- dim(NumOut$GFEVD[[ModelType]]$Yields$NonOrtho)[3] # Total number of yields of the system
            for (i in 1:CJ){ # Outputs in data-frame
              GFEVDYields[[i]] <- data.frame(NumOut$GFEVD[[ModelType]]$Yields$NonOrtho[,,i])
              GFEVDYields[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(GFEVDYields[[i]]), GFEVDYields[[i]])))
              GFEVDYields[[i]]$index <- factor(GFEVDYields[[i]]$index, levels= 1:(GFEVDhoriz))
            }
            nmYields <- dimnames(NumOut$GFEVD[[ModelType]]$Yields$NonOrtho)[[3]] # Yield names

          }else{
            CJ <- dim(NumOut$GFEVD[[ModelType]]$Yields)[3] # Total number of yields of the system
            for (i in 1:CJ){ # Outputs in data-frame
              GFEVDYields[[i]] <- data.frame(NumOut$GFEVD[[ModelType]]$Yields[,,i])
              GFEVDYields[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(GFEVDYields[[i]]), GFEVDYields[[i]])))
              GFEVDYields[[i]]$index <- factor(GFEVDYields[[i]]$index, levels= 1:(GFEVDhoriz))
            }
            nmYields <- dimnames(NumOut$GFEVD[[ModelType]]$Yields)[[3]] # Yield names
          }

          ## Generate plots:
          for (i in 1:CJ){
            index <- GFEVDYields[[i]]$index
            value <- GFEVDYields[[i]]$value
            variable <- GFEVDYields[[i]]$variable

             p <- ggplot(GFEVDYields[[i]], aes(x= index, y=value , fill=variable)) +
              geom_bar(stat="identity", width = 0.25) +
              labs(title= paste0("GFEVD_", nmYields[i])) + theme_classic() +
               theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                     plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
             suppressMessages(ggplot2::ggsave(p, file=paste0("GFEVDYields_", nmYields[i],".png"),
                            path= PathsGraphs[[ModelType]]$GFEVD[["Yields"]], width= 11, height= 7.98))
          }
        }

    }

}

#####################################################################################################################
######################################## 6) TERM PREMIA DECOMPOSITION ##################################################
#####################################################################################################################
#' Term Premia decomposition graphs for "joint Q" models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, GFEVDs and Risk premia
#'@param ModelPara list of model parameter estimates (See the "Optimization" function)
#'@param WishRPgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param UnitYields (i) "Month": if maturity of yields are expressed in months or
#'                  (ii) "Year": if maturity of yields are expressed in years
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'@keywords internal


TPDecompGraphJoint <- function(ModelType, NumOut, ModelPara, WishRPgraphs, UnitYields, Economies, PathsGraphs){



  if (WishRPgraphs ==0){cat("No graphs for term-premia were generated \n")
  } else {

    cat(' ** Term premia \n' )


    dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/TermPremia", sep=""))  # Folder creation


  # Preliminary work
  dt <- ModelPara[[ModelType]]$inputs$dt
  mat <- ModelPara[[ModelType]]$inputs$mat
  Yields <- ModelPara[[ModelType]]$inputs$Y

  TPdecomp <- NumOut$TermPremiaDecomp

  J <- length(mat)
  C <- length(Economies)
  T <- ncol(Yields)

  if (UnitYields== "Month"){
    k <- 12
    YLab <- "Months"
    }
  if (UnitYields== "Year"){
    k <- 1
    YLab <- "Years"
    }
  matAdjUnit <- mat*k



  YieldData <- list()
  for (i in 1:C){
    IdxRP <- 1:J + J*(i-1)
    YieldData[[Economies[i]]] <- t(Yields[IdxRP, ]*100)
  }

  # Graph title, legends and axes
  Dates <- 1:T

  matRPLab <- paste(matAdjUnit, YLab)
  GraphLegend <- c("Data", "Expected Component", "Term Premium")


    # Prepare plots
    subplots_RP <- list()

    for (i in 1:C){ # per country

      p <- c()
      plot_list_RP <- list()
      plot_list_no_legend_RP <- list()

      for (h in 1:length(matRPLab)){  # per maturity

        Data <- YieldData[[Economies[i]]][,h]
        ExpCom <- TPdecomp$RiskPremia[["Expected Component"]][[Economies[i]]][,h]
        TP <- TPdecomp$RiskPremia[["Term Premia"]][[Economies[i]]][,h]
        TimeSpan <- Dates

        DataGraph <- data.frame(Data, ExpCom, TP, TimeSpan)

        # Graph: Fit x Data
        p <- ggplot(DataGraph) + geom_line(aes(x=TimeSpan, y=Data, color = GraphLegend[1])) +
          geom_line(aes(x=TimeSpan, y = ExpCom, color = GraphLegend[2])) +
          geom_line(aes(x=TimeSpan, y = TP, color = GraphLegend[3])) +
          labs(color = "Legend") + ggtitle( matRPLab[h]) +  theme_classic() +
          theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
                axis.title.y=element_blank(), axis.text.x = element_text(size=6) ) +
          geom_hline(yintercept=0)


        plot_list_RP[[h]] <- p # All plots with legends

        plot_list_no_legend_RP[[h]] <- p + theme(legend.position = "none")  # Hide legends from all plots

      }

      # Common title for each model
      title <- cowplot::ggdraw() + cowplot::draw_label(Economies[i], fontface='bold', size=14)

      # Generate legend
      # Legend settings:
      LegendSubPlot_RP <- cowplot::get_legend(plot_list_RP[[1]] +
                                                theme(legend.direction = "horizontal", legend.position="bottom", legend.justification="center",
                                                      legend.title = element_text(size=10, face="bold"),
                                                      legend.text = element_text( size = 8),
                                                      legend.box.background = element_rect(colour = "black")) )


      # Build the graph subplots:
      subplots_CS <- cowplot::plot_grid(plotlist= plot_list_no_legend_RP, ncol=3)
      subplots2save <- cowplot::plot_grid(LegendSubPlot_RP, subplots_CS, ncol=1, rel_heights=c(0.2, 1))

      suppressMessages(ggplot2::ggsave(subplots2save, file=paste0("TermPremia_", Economies[i],"_", ModelType, ".png"),
                path = PathsGraphs[[ModelType]]$TermPremia))

        }

}
}

##########################################################################################################
##################### GRAPHICAL OUTPUTS FOR THE ORTHOGONALIZED FACTORS OF JLL ############################
##########################################################################################################

######################################################################################################
##################################### 2) FEVD ########################################################
#####################################################################################################
#' FEVDs graphs for orthogonalized risk factors of JLL-based models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param FEVDhoriz  single numerical vector conataining the desired horizon of analysis for the FEVDs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'
#'@importFrom ggplot2 ggplot theme_classic aes element_text theme labs ggtitle element_blank aes_string geom_bar geom_hline
#'@importFrom zoo index
#'
#'@keywords internal

FEVDgraphsJLLOrtho <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, FEVDhoriz, PathsGraphs) {


  if (WishPdynamicsgraphs== 0 & WishYieldsgraphs==0){
    cat("No graphs for FEVDs were generated (orthogonolized version) \n")
  }else{

    cat(' ** FEVDs-Ortho \n')

        ######################################### Factors #########################################################
        if (WishPdynamicsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Factors/Ortho", sep="")) # Folder Creation

          FEVDFactors <- list()
          K <- dim(NumOut$FEVD[[ModelType]]$Factors$Ortho)[2] # Total number of risk factors

          for (i in 1:K){ # Recast the outputs in data-frame format
            FEVDFactors[[i]] <- data.frame(NumOut$FEVD[[ModelType]]$Factors$Ortho[,,i])
            FEVDFactors[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDFactors[[i]]), FEVDFactors[[i]])))
            FEVDFactors[[i]]$index <- factor(FEVDFactors[[i]]$index, levels= 1:(FEVDhoriz))
          }

          nmFactors <- colnames(NumOut$FEVD[[ModelType]]$Factors$Ortho) # Factor names
          ## Generate plots:
          for (i in 1:K){
            index <- FEVDFactors[[i]]$index
            value <- FEVDFactors[[i]]$value
            variable <- FEVDFactors[[i]]$variable


            p <- ggplot(FEVDFactors[[i]], aes(x= index, y=value , fill=variable)) +
              geom_bar(stat="identity", width = 0.25) +
              labs(title= paste0("FEVD_",nmFactors[i])) + theme_classic() +
              theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                    plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
            suppressMessages(ggplot2::ggsave(p, file=paste0("FEVDFactors_", nmFactors[[i]], "ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$FEVD[["Factors Ortho"]], width= 11, height= 7.98))
          }

        }
        ############################################ Yields #########################################################
        if (WishYieldsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Yields/Ortho", sep="")) # Folder Creation

          FEVDYields <- list()
          CJ <- dim(NumOut$FEVD[[ModelType]]$Yields$Ortho)[3] # Total number of yields of the system

          for (i in 1:CJ){ # Recast the outputs in data-frame format
            FEVDYields[[i]] <- data.frame(NumOut$FEVD[[ModelType]]$Yields$Ortho[,,i])
            FEVDYields[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDYields[[i]]), FEVDYields[[i]])))
            FEVDYields[[i]]$index <- factor(FEVDYields[[i]]$index, levels= 1:(FEVDhoriz))
          }

          nmYields <- dimnames(NumOut$FEVD[[ModelType]]$Yields$Ortho)[[3]] # Yield names

          ## Generate plots:
          for (i in 1:CJ){

            index <- FEVDYields[[i]]$index
            value <- FEVDYields[[i]]$value
            variable <- FEVDYields[[i]]$variable

            p <- ggplot(FEVDYields[[i]], aes(x= index, y=value , fill=variable)) +
              geom_bar(stat="identity", width = 0.25) +
              labs(title= paste0("FEVD_", nmYields[i])) + theme_classic() +
              theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                    plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
            suppressMessages(ggplot2::ggsave(p, file=paste0("FEVDYields_", nmYields[i], "ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$FEVD[["Yields Ortho"]], width= 11, height= 7.98))
          }
        }

    }

}
######################################################################################################
##################################### 4) GFEVD ########################################################
#####################################################################################################
#' GFEVDs graphs for orthogonalized risk factors of JLL-based models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param GFEVDhoriz  single numerical vector conataining the desired horizon of analysis for the GFEVDs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'
#'@importFrom ggplot2 ggplot aes element_text theme theme_classic labs ggtitle element_blank aes_string geom_bar geom_hline
#'@importFrom zoo index
#'
#'@keywords internal


GFEVDgraphsJLLOrtho <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, GFEVDhoriz, PathsGraphs) {


  if (WishPdynamicsgraphs== 0 & WishYieldsgraphs==0){
    cat("No graphs for GFEVDs  were generated (orthogonolized version) \n")
  }else{

    cat(' ** GEFVDs-Ortho \n\n' )

  ######################################### Factors #########################################################
  if (WishPdynamicsgraphs == 1){

    dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Factors/Ortho", sep="")) # Folder Creation

        GFEVDFactors <- list()
          K <- dim(NumOut$GFEVD[[ModelType]]$Factors$Ortho)[2] # Total number of risk factors

          for (i in 1:K){ # Recast outputs in data-frame format
            GFEVDFactors[[i]] <- data.frame(NumOut$GFEVD[[ModelType]]$Factors$Ortho[,,i])
            GFEVDFactors[[i]] <-  suppressMessages(reshape2::melt(cbind(index=rownames(GFEVDFactors[[i]]), GFEVDFactors[[i]])))
            GFEVDFactors[[i]]$index <- factor(GFEVDFactors[[i]]$index, levels= 1:(GFEVDhoriz))
          }

          nmFactors <- colnames(NumOut$GFEVD[[ModelType]]$Factors$Ortho) # Factor names
          ## Generate plots:
          for (i in 1:K){
            index <- GFEVDFactors[[i]]$index
            value <- GFEVDFactors[[i]]$value
            variable <- GFEVDFactors[[i]]$variable

            p <- ggplot(GFEVDFactors[[i]], aes(x= index, y=value , fill=variable)) +
              geom_bar(stat="identity", width = 0.25) +
              labs(title= paste0("GFEVD_",nmFactors[i])) + theme_classic() +
              theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8) )
            suppressMessages(ggplot2::ggsave(p, file=paste0("GFEVDFactors_", nmFactors[[i]], "ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$GFEVD[["Factors Ortho"]], width= 11, height= 7.98))
          }

        }
        ############################################ Yields #########################################################
        if (WishYieldsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Yields/Ortho", sep="")) # Folder Creation

          GFEVDYields <- list()
          CJ <- dim(NumOut$GFEVD[[ModelType]]$Yields$Ortho)[3] # Total number of yields of the system

          for (i in 1:CJ){ # Recast outputs in data-frame format
            GFEVDYields[[i]] <- data.frame(NumOut$GFEVD[[ModelType]]$Yields$Ortho[,,i])
            GFEVDYields[[i]] <- suppressMessages(reshape2::melt(cbind(index=rownames(GFEVDYields[[i]]), GFEVDYields[[i]])))
            GFEVDYields[[i]]$index <- factor(GFEVDYields[[i]]$index, levels= 1:(GFEVDhoriz))
          }

          nmYields <- dimnames(NumOut$GFEVD[[ModelType]]$Yields$Ortho)[[3]] # Yield names

          ## Generate plots:
          for (i in 1:CJ){
            index <- GFEVDYields[[i]]$index
            value <- GFEVDYields[[i]]$value
            variable <- GFEVDYields[[i]]$variable

            p <- ggplot(GFEVDYields[[i]], aes(x= index, y=value , fill=variable)) +
              geom_bar(stat="identity", width = 0.25) +
              labs(title= paste0("GFEVD_", nmYields[i])) + theme_classic() +
              theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8) )
            suppressMessages(ggplot2::ggsave(p, file=paste0("GFEVDYields_", nmYields[i], "ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$GFEVD[["Yields Ortho"]], width= 11, height= 7.98))
          }
        }

  }

}

################################################################################################################
#'Build subplot for fitted yields
#'
#'@param YieldData Time series of bond yields
#'@param ModelFit Time series of fitted bond yields
#'@param ModelImplied Time series of model-implied bond yields
#'@param MatLength number of country-specific maturities
#'@param mat vector of maturities
#'@param Economies Economies of the economic system
#'@param ModelType Desired estimated model
#'@param PathsGraphs Path to save the graphs
#'
#'@keywords internal

Fit_Subplot <- function(YieldData, ModelFit, ModelImplied, MatLength, mat, Economies, ModelType, PathsGraphs){

  if (any(ModelType == c("JPS original", "JPS global",  "GVAR single"))){
    AdjPath <- PathsGraphs[[ModelType]]$Fit[[Economies]]
    } else{ AdjPath <- PathsGraphs[[ModelType]]$Fit }


  #Graph titles
  matmonths <- paste(mat, "months - ")
  GraphTitles <- paste(matmonths, Economies)

  GraphLegend <- c("Data", "Model Fit", "Model-Implied")

  T <- ncol(YieldData)
  # Prepare plots
  p <- c()
  plot_list_FIT <- list()
  plot_list_no_legend_FIT <- list()


  for (j in 1:MatLength){
    Data <- YieldData[j,-1]
    Modfit <- ModelFit[j,-1]
    ModImp <- ModelImplied[j,-1]
    TimeSpan <- 1:(T-1)

    YieldCompar <- data.frame(Data, Modfit, ModImp, TimeSpan)


    # Graph: Fit x Data
    p <- ggplot(data = YieldCompar, aes(x= TimeSpan )) +
      geom_line(aes_string(y = Data, color = shQuote(GraphLegend[1])), size = 0.5) +
      geom_line(aes_string(y = Modfit, color = shQuote(GraphLegend[2])), size = 0.5) +
      geom_line(aes_string(y = ModImp, color = shQuote(GraphLegend[3])), linetype = "dashed") +
      labs(color = "Legend") + ggtitle( GraphTitles[j]) +  theme_classic() +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
            axis.title.y=element_blank(), axis.text.x = element_text(size=8) ) +
      geom_hline(yintercept=0)


    plot_list_FIT[[j]] <- p # All plots with legends

    plot_list_no_legend_FIT[[j]] <- p + theme(legend.position = "none")  # Hide legends from all plots

  }

  # Generate legend
  # Legend settings:
  LegendSubPlot_FIT <- cowplot::get_legend(plot_list_FIT[[1]] +
                                             theme(legend.direction = "horizontal", legend.position="bottom", legend.justification="center",
                                                   legend.title = element_text(size=10, face="bold"),
                                                   legend.text = element_text( size = 8),
                                                   legend.box.background = element_rect(colour = "black")) )


  # Build the graph subplots:
  subplots_FIT <- cowplot::plot_grid(plotlist= plot_list_no_legend_FIT, ncol=3)
  subplots2save_FIT <- cowplot::plot_grid(LegendSubPlot_FIT, subplots_FIT , ncol=1,  rel_heights = c(.1, 1))
  suppressMessages(ggplot2::ggsave(subplots2save_FIT, file=paste0("Yields_Fit_", Economies, ".png"),
                                   path = AdjPath))

}
#########################################################################################################
#'Create folders for storing IRFs and GIRFs
#'
#'@param OutputType available options are "IRF", "GIRF", "IRF Ortho" and "GIRF Ortho"
#'@param WishPdynamicsgraphs binary variable specifing whether the user whishes IRFs and/or GIRFs of risk factors
#'@param WishYieldsgraphs binary variable specifing whether the user whishes IRFs and/or GIRFs of bond yields
#'@param Economies Set of economies that are part of the economic system
#'@param ModelType Desired modem type
#'
#'@keywords internal

FolderPrep_IRFs <- function(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType){


  ################# 1) SINGLE-COUNTRY CONTRY MODELS #################
    if (any(ModelType == c("JPS original", "JPS global",  "GVAR single"))){

  # A) General Folder Creation
  if(OutputType == "IRF"){
    dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies, "/IRF", sep=""))
  }else{ dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies, "/GIRF", sep="")) }


  # a.1) Folders for graph of Factors
  if (WishPdynamicsgraphs == 1){
    if(OutputType == "IRF"){
      dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
                       "/IRF/Factors", sep=""))
    }else{
      dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
                       "/GIRF/Factors", sep=""))
    }
  }

  # a.2) Folders for graph of yields
  if (WishYieldsgraphs == 1){

    if(OutputType == "IRF"){
      dir.create( paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
                        "/IRF/Yields", sep=""))
    }else{
      dir.create( paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
                        "/GIRF/Yields", sep=""))
    }
  }
  } else if(any(ModelType == c("JPS multi", "GVAR multi"))){

    ################# 2) JOINT COUNTRY CONTRY MODELS #################
    # A) General Folder Creation
    if(OutputType == "IRF"){
      dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/IRF", sep=""))
    }else{ dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/GIRF", sep=""))}

    # a.1) Folders for graph of Factors
    if (WishPdynamicsgraphs == 1){
    if(OutputType == "IRF"){
      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF/Factors", sep=""))
    }else{ dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GIRF/Factors", sep="")) }
    }

    # a.2) Folders for graph of yields
    if (WishYieldsgraphs == 1){
    if(OutputType == "IRF"){
      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF/Yields", sep=""))
    }else{ dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GIRF/Yields", sep="")) }
    }

    ################# 3) JLL SPECIFIC MODELS #################
    } else{
      # A) General Folder Creation
      if(OutputType == "IRF"){
        dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/IRF", sep=""))
      }else if(OutputType == "GIRF"){ dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/GIRF", sep=""))}

      ########### NON-ORTHO #####################################
      if(any(OutputType == c("IRF", "GIRF"))){
      # a.1) Folders for graph of Factors
      if (WishPdynamicsgraphs == 1){
        if(OutputType == "IRF"){
          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF/Factors", sep=""))
        }else{ dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GIRF/Factors", sep="")) }
      }

      # a.2) Folders for graph of yields
      if (WishYieldsgraphs == 1){
        if(OutputType == "IRF"){
          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF/Yields", sep=""))
        }else{ dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GIRF/Yields", sep="")) }
      }
    } else{

    ########### ORTHO #####################################
      if (WishPdynamicsgraphs == 1){
      if(OutputType == "IRF Ortho"){
    dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF/Factors/Ortho", sep=""))
    }else{ dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GIRF/Factors/Ortho", sep=""))  }
    }

      # a.2) Folders for graph of yields
      if (WishYieldsgraphs == 1){
      if(OutputType == "IRF Ortho"){
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF/Yields/Ortho", sep=""))
      }else{  dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GIRF/Yields/Ortho", sep="")) }
      }
      }
}
}

#########################################################################################################
#'Generate paths to save IRFs/GIRFs graphs
#'
#'@param OutputType available options are "IRF" and "GIRF"
#'@param ElementType available options are "Factors" and "Yields"
#'@param PathsGraphs desired path to save the graphs
#'@param Economies Economies of the economic system
#'@param ModelType Desired estimated model type
#'
#'@keywords internal

AdjustPathIRFs <- function(OutputType, ElementType, PathsGraphs, Economies, ModelType){

  # 1) Models estimated seperatly
  if (any(ModelType == c("JPS original", "JPS global",  "GVAR single"))){
    if(OutputType == "IRF"){
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$IRF[[Economies]]$Factors} # Factors
      else{ PathAdj <- PathsGraphs[[ModelType]]$IRF[[Economies]]$Yields } # Yields

    }else{
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$GIRF[[Economies]]$Factors}
      else{ PathAdj <- PathsGraphs[[ModelType]]$GIRF[[Economies]]$Yields }
    }

    # 2) Models estimated jointly
  }else{
    if(OutputType == "IRF"){ # IRF
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$IRF$Factors} # Factors
      else{ PathAdj <- PathsGraphs[[ModelType]]$IRF$Yields } # Yields

    }else if(OutputType == "GIRF") { # GIRF
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$GIRF$Factors} # Factors
      else{ PathAdj <- PathsGraphs[[ModelType]]$GIRF$Yields } # Yields

    # 3) Exclusively for JLL models
      # IRF ortho
    }else if(OutputType == "IRF Ortho"){
    if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$IRF[["Factors Ortho"]] } # Factors
    else{ PathAdj <- PathsGraphs[[ModelType]]$IRF[["Yields Ortho"]] } # Yields

      # GIRF ortho
    }else {
    if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$GIRF[["Factors Ortho"]] } # Factors
    else{ PathAdj <- PathsGraphs[[ModelType]]$GIRF[["Yields Ortho"]] } # Yields
  }
  }

  return(PathAdj)
}
##########################################################################################################
#' Build the list of IRF and GIRF for both factors and bond yields
#'
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, GFEVDs and Term Premia
#'@param Economies Economies of the economic system
#'@param ModelType Desired model type
#'@param IRFhoriz time-horizon of the IRF and GIRF
#'@param FacDim dimension of the risk factor vector
#'@param OutputType available option are 'IRF' and 'GIRF'
#'
#'@keywords internal

BuildIRFlist <- function(NumOut, Economies, ModelType, IRFhoriz, FacDim, OutputType){

  Horiz <- 0:(IRFhoriz-1)
  IRFFactors <- list()
  IRFYields <- list()

  # 1) Extract IRFs
  if (OutputType == "IRF"){
  # a) Models estimated individually
    if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
    for (g in 1:FacDim){ # Recast the IRFs into a data-frame
      IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]][[Economies]]$Factors[,,g]))
      IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]][[Economies]]$Yields[,,g]))
    }
  # b) Models estimated jointly
} else if (any(ModelType == c("GVAR single", "GVAR multi"))){
    for (g in 1:FacDim){ # Outputs in data-frame format
      IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Factors[,,g]))
      IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Yields[,,g]))
    }
} else{
  # c) JLL-based setups (Non-Orthogonalized version)
  for (g in 1:FacDim){
    IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Factors$NonOrtho[,,g]))
    IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Yields$NonOrtho[,,g]))
  }
}
  # 2) Extract GIRFs
} else if (OutputType == "GIRF"){

    # a) Models estimated individually
    if ( any(ModelType == c("JPS original", "JPS globnal", "GVAR single"))){
      for (g in 1:FacDim){ # Recast the IRFs into a data-frame
        IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]][[Economies]]$Factors[,,g]))
        IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]][[Economies]]$Yields[,,g]))
      }

    # b) Models estimated jointly
    } else if (any(ModelType == c("GVAR single", "GVAR multi"))){
      for (g in 1:FacDim){ # Outputs in data-frame format
        IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Factors[,,g]))
        IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Yields[,,g]))
      }

      # c) JLL-based setups (non-Orthogonalized version)
    } else{
        for (g in 1:FacDim){
          IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Factors$NonOrtho[,,g]))
          IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Yields$NonOrtho[,,g]))
        }
  }


  # JLL-based setups (orthogonalized version)
    }else if(OutputType == "IRF Ortho"){
      for (g in 1:FacDim){ # Recast outputs in a data-frame
      IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Factors$Ortho[,,g]))
      IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Yields$Ortho[,,g]))
    }
    } else{ # GIRF Ortho
      for (g in 1:FacDim){ # Recast outputs in a data-frame
        IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Factors$Ortho[,,g]))
        IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Yields$Ortho[,,g]))
      }
  }


      nmFactors <- names(IRFFactors[[1]]) # Factor names
      nmYields <- names(IRFYields[[1]]) # Yield names

  return(list(IRFFactors = IRFFactors, IRFYields = IRFYields))
    }


##########################################################################################################
#'Gather data for IRFs and GIRFs grahs (version "Factors")
#'
#'@param IRFFac  Data-frame with basic features of a single IRF for risk factors
#'
#'@keywords internal

IRFandGIRFs_Format_Fac <- function(IRFFac){

  nmFactors <- names(IRFFac) # Factor names
  K <- length(nmFactors) - 1
  plot_list <- list()

  for (x in 1:K){
    p <- c()
    p <- ggplot(IRFFac, aes_string(x=nmFactors[1], y= nmFactors[x+1])) + geom_line() +  ggtitle( nmFactors[x+1])
    p <- p  + theme_classic() + theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                  axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                   axis.text.y = element_text(size=4))
    plot_list[[x]] <- p
  }

  return(plot_list)
}

#####################################################################################################
#'Gather data for IRFs and GIRFs grahs (version "Yields")
#'
#'@param IRFYields Data-frame with basic features of a single IRF for yields
#'
#'@keywords internal

IRFandGIRFs_Format_Yields <- function(IRFYields){

  nmYields <- names(IRFYields) # Yields names
  J <- length(nmYields) - 1

  plot_list <- list()

for (x in 1:J){ # Generate graph-by-graph
  p <- c()
  p <- ggplot(IRFYields, aes_string(x=nmYields[1], y= nmYields[x+1])) + geom_line() + labs(title=nmYields[x+1])
  p <- p  + theme_classic() + theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                                    axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                                    axis.text.y = element_text(size=4))
  plot_list[[x]] <- p
}

  return(plot_list)
}

##########################################################################################################
#' Extract list of desired graph features (IRFs anc GIRFs)
#'
#'@param InputsForOutputs List of inputs for outputs
#'@param OutType Output types "IRF", "GIRF" and "IRF Ortho"
#'@param ModelType desired model type
#'
#'@keywords internal

Wished_Graphs_IRFandGIRF <- function(InputsForOutputs, OutType, ModelType){

  # 1) JLL models
  if (any(ModelType == c("JLL original",  "JLL No DomUnit",  "JLL joint Sigma"))){
    if(OutType == "IRF"){
      RiskGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphs$Yields

    }else if(OutType == "GIRF"){
      RiskGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphs$Yields

    } else if (OutType == "IRF Ortho"){
      RiskGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphsOrtho$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphsOrtho$Yields

    } else{ # GIRF ortho
      RiskGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphsOrtho$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphsOrtho$Yields
    }


  }else{
  # 2) All other models
  if(OutType == "IRF"){
    RiskGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphs$RiskFactors
    YieldGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphs$Yields

  }else if(OutType == "GIRF"){
    RiskGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphs$RiskFactors
    YieldGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphs$Yields
  }
  }
  return(list(RiskGraphs = RiskGraphs, YieldGraphs = YieldGraphs))
}

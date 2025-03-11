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
  dirs <- c("Outputs", paste0("Outputs/", ModelType), paste0("Outputs/", ModelType, "/Point Estimate"))
  lapply(file.path(tempdir(), dirs), dir.create, recursive = TRUE, showWarnings = FALSE)

  PathsGraphs <- FolderCreationPoint(ModelType, Economies)

# 1) Plot the set of risk factors
  RiskFactorsGraphs(ModelType, ModelPara, Economies, FactorLabels)

  cat(" 2.3) Generating the graphs of interest \n")

  if (ModelType %in% c("JPS original", "JPS global",  "GVAR single")){
    for (i in 1:length(Economies)){dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], sep=""))}}

  # 2) Model fit
  Fitgraphs(ModelType, InputsForOutputs[[ModelType]]$Fit$WishGraphs, ModelPara, NumOut, Economies, PathsGraphs)

  # 3) IRF and GIRF
  if (any(ModelType == c("JLL original",  "JLL No DomUnit",  "JLL joint Sigma"))){
    OutType <- c("IRF", "GIRF","IRF Ortho", "GIRF Ortho") } else{ OutType <- c("IRF", "GIRF")}

  for(j in 1:length(OutType)){
  WG <- Wished_Graphs_IRFandGIRF(InputsForOutputs, OutType[j], ModelType)
  IRFandGIRFgraphs(ModelType, NumOut, WG$RiskGraphs, WG$YieldGraphs, InputsForOutputs[[ModelType]]$IRF$horiz,
                   PathsGraphs, OutputType = OutType[j], Economies)
     }

  # 4) FEVD and GFEVD
  if (any(ModelType == c("JLL original",  "JLL No DomUnit",  "JLL joint Sigma"))){
    OutType <- c("FEVD", "GFEVD","FEVD Ortho", "GFEVD Ortho") } else{ OutType <- c("FEVD", "GFEVD")}

  for(j in 1:length(OutType)){
    WG <- Wished_Graphs_FEVDandGFEVD(InputsForOutputs, OutType[j], ModelType)
    FEVDandGFEVDgraphs(ModelType, NumOut, WG$RiskGraphs, WG$YieldGraphs, InputsForOutputs[[ModelType]]$FEVD$horiz,
                       PathsGraphs, OutType[j], Economies)
  }

  # 5) Term premia decomposition
  TPDecompGraph(ModelType, NumOut, ModelPara, InputsForOutputs[[ModelType]]$RiskPremia$WishGraphs,
                   InputsForOutputs$UnitMatYields, Economies, PathsGraphs)

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
#' Model fit graphs for all models
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

    cat(' ** Fit of bond yields \n' )

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
#' IRF and GIRF graphs for all models
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
    if(any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
      K <- dim(NumOut$IRF[[ModelType]]$Factors$Ortho)[2]
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
##################################### 3) FEVD and GFEVD ##############################################
#####################################################################################################
#' FEVD and GFEVD graphs for all models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param FEVDhoriz  single numerical vector containing the desired horizon of analysis for the FEVDs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'@param OutputType Available options are 'FEVD' and 'GFEVD'
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'@importFrom ggplot2 ggplot theme_classic aes element_text theme labs ggtitle element_blank aes_string geom_bar geom_hline
#'
#'@keywords internal


FEVDandGFEVDgraphs <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, FEVDhoriz,
                               PathsGraphs, OutputType, Economies){

  if (WishPdynamicsgraphs==0 & WishYieldsgraphs ==0){ cat(paste( "No graphs for", OutputType, "were generated \n"))

    if(any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))
       && any(ModelType == c("FEVD Ortho","GFEVD Ortho"))){
      cat(paste("No graphs for" , OutputType, "were generated (orthogonolized version) \n"))}

  }else {

    if(OutputType == "FEVD"){  cat(' ** FEVDs \n');
      Lab_Fac <-  "FEVDFactors_"; Lab_Yield <- "FEVDYields_"
    } else if( OutputType == "GFEVD"){ cat(' ** GFEVDs \n');
      Lab_Fac <-  "GFEVDFactors_"; Lab_Yield <- "GFEVDYields_"
    } else if(OutputType == "FEVD Ortho"){  cat(' ** FEVDs-Ortho \n')
      Lab_Fac <-  "FEVDFactors_ORTHO_"; Lab_Yield <- "FEVDYields_ORTHO_"
    } else{ cat(' ** GFEVDs-Ortho \n')
      Lab_Fac <-  "GFEVDFactors_ORTHO_" ; Lab_Yield <- "GFEVDYields_ORTHO_"}


    ################ 1) Estimation done for countries individually ################
    if (any(ModelType == c("JPS original", "JPS global",  "GVAR single"))){

      K <- dim(NumOut$FEVD[[ModelType]][[Economies[1]]]$Factors)[[2]] # Total number of risk factors
      J <- dim(NumOut$FEVD[[ModelType]][[Economies[1]]]$Yields)[3]

      C <- length(Economies)

      for (i in 1:C){
        # a) Folder Creation
        FolderPrep_FEVDs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies[i], ModelType)

        # b) Recast FEVDs or GFEVDs
        FEDset <- BuildFEVDlist(NumOut, Economies[i], ModelType, FEVDhoriz, K, J, OutputType)

        # c) Graph of Factors
        if (WishPdynamicsgraphs == 1){
          PathAdj <- AdjustPathFEVDs(OutputType, "Factors", PathsGraphs, Economies[i], ModelType)
          nmFactors <- colnames(NumOut$FEVD[[ModelType]][[Economies[i]]]$Factors) # Factor names

          for(h in 1:K){ FEVDandGFEVDs_Graphs(OutputType, FEDset$FEVDFactors[[h]], nmFactors[h], Lab_Fac,
                                              PathAdj) }
        }

        # d) Graph of Yields
        if (WishYieldsgraphs == 1){
          PathAdj <- AdjustPathFEVDs(OutputType, "Yields", PathsGraphs, Economies[i], ModelType)
          nmYields <- dimnames(NumOut$FEVD[[ModelType]][[Economies[i]]]$Yields)[[3]] # Yield labels

          for(h in 1:J){ FEVDandGFEVDs_Graphs(OutputType, FEDset$FEVDYields[[h]], nmYields[h], Lab_Yield,
                                              PathAdj) }
        }
      }
    } else{

      ################ 2) Estimation done for countries jointly ###############################
      if (any(OutputType == c("FEVD", "GFEVD"))){
        # Total number of risk factors
        if(any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))){
          K <- dim(NumOut$FEVD[[ModelType]]$Factors$NonOrtho)[2]
          CJ <- dim(NumOut$GFEVD[[ModelType]]$Yields$NonOrtho)[3]
          nmFactors <- colnames(NumOut$FEVD[[ModelType]]$Factors$NonOrtho)
          nmYields <- dimnames(NumOut$FEVD[[ModelType]]$Yields$NonOrtho)[[3]]
        }else{
          K <- dim(NumOut$FEVD[[ModelType]]$Factors)[[2]]
          CJ <- dim(NumOut$FEVD[[ModelType]]$Yields)[[3]]
          nmFactors <- colnames(NumOut$FEVD[[ModelType]]$Factors) # Factor names
          nmYields <- dimnames(NumOut$FEVD[[ModelType]]$Yields)[[3]] # Yield labels
        }

        # a) Folder Creation
        FolderPrep_FEVDs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType)

        # b) Recast FEVDs or GFEVDs
        FEDset <- BuildFEVDlist(NumOut, Economies, ModelType, FEVDhoriz, K, CJ, OutputType)

        # c) Graph of Factors
        if (WishPdynamicsgraphs == 1){
          PathAdj <- AdjustPathFEVDs(OutputType, "Factors", PathsGraphs, Economies, ModelType)
          for(h in 1:K){ FEVDandGFEVDs_Graphs(OutputType, FEDset$FEVDFactors[[h]], nmFactors[h], Lab_Fac,
                                              PathAdj) }
        }

        # d) Graph of Yields
        if (WishYieldsgraphs == 1){
          PathAdj <- AdjustPathFEVDs(OutputType, "Yields", PathsGraphs, Economies, ModelType)
          for(h in 1:CJ){ FEVDandGFEVDs_Graphs(OutputType, FEDset$FEVDYields[[h]], nmYields[h], Lab_Yield,
                                               PathAdj) }
        }

      }
      ###################### 3) JLL orthoonalized ############################################################
      if(any(OutputType == c("FEVD Ortho", "GFEVD Ortho"))){

        # Total number of risk factors
        K <- dim(NumOut$FEVD[[ModelType]]$Factors$Ortho)[2]
        CJ <- dim(NumOut$GFEVD[[ModelType]]$Yields$Ortho)[3]
        nmFactors <- colnames(NumOut$FEVD[[ModelType]]$Factors$Ortho)
        nmYields <- dimnames(NumOut$FEVD[[ModelType]]$Yields$Ortho)[[3]]

        # a) Folder Creation
        FolderPrep_FEVDs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType)

        # b) Recast FEVDs or GFEVDs
        FEDset <- BuildFEVDlist(NumOut, Economies, ModelType, FEVDhoriz, K, CJ, OutputType)

        # c) Graph of Factors
        if (WishPdynamicsgraphs == 1){
          PathAdj <- AdjustPathFEVDs(OutputType, "Factors", PathsGraphs, Economies, ModelType)
          for(h in 1:K){ FEVDandGFEVDs_Graphs(OutputType, FEDset$FEVDFactors[[h]], nmFactors[h], Lab_Fac,
                                              PathAdj) }
        }

        # d) Graph of Yields
        if (WishYieldsgraphs == 1){
          PathAdj <- AdjustPathFEVDs(OutputType, "Yields", PathsGraphs, Economies, ModelType)
          for(h in 1:CJ){ FEVDandGFEVDs_Graphs(OutputType, FEDset$FEVDYields[[h]], nmYields[h], Lab_Yield,
                                               PathAdj) }
        }
      }
    }
  }
}

#####################################################################################################################
######################################## 4) TERM PREMIA DECOMPOSITION ##################################################
#####################################################################################################################
#' Term Premia decomposition graphs for all models
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

TPDecompGraph <- function(ModelType, NumOut, ModelPara, WishRPgraphs, UnitYields, Economies, PathsGraphs){


  if (WishRPgraphs ==0){cat("No graphs for term-premia were generated \n")
  } else {

    cat(' ** Term premia \n')

    # Preliminary work
    SepQ_Lab <- c("JPS original", "JPS global",  "GVAR single")

    if(!ModelType %in% SepQ_Lab){
      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/TermPremia", sep=""))  # Folder creation
    }

    if(ModelType %in% SepQ_Lab){
      dt <- ModelPara[[ModelType]][[Economies[1]]]$inputs$dt
      mat <- ModelPara[[ModelType]][[Economies[1]]]$inputs$mat
      T <- ncol(ModelPara[[ModelType]][[Economies[1]]]$inputs$Y)
    }else{
      dt <- ModelPara[[ModelType]]$inputs$dt
      mat <- ModelPara[[ModelType]]$inputs$mat
      T <- ncol(ModelPara[[ModelType]]$inputs$Y)
      Yields <- ModelPara[[ModelType]]$inputs$Y
    }

    J <- length(mat)
    C <- length(Economies)

    TPdecomp <- NumOut$TermPremiaDecomp

    if (UnitYields== "Month"){  k <- 12;  YLab <- "Months"
    }else {   k <- 1;   YLab <- "Years"} # Case  UnitYields== "Year"
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

      if(ModelType %in% SepQ_Lab){
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], "/TermPremia", sep=""))  # Folder creation
        GraphPath <- PathsGraphs[[ModelType]]$TermPremia[[Economies[i]]]

        Yields <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y
        YieldData[[Economies[i]]] <- t(Yields)*100
      }else{
        GraphPath <- PathsGraphs[[ModelType]]$TermPremia

        IdxRP <- 1:J + J*(i-1)
        YieldData[[Economies[i]]] <- t(Yields[IdxRP, ]*100)
      }

      p <- c()

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
                                       path = GraphPath))

    }

  }
}
################################################################################################################
######################### AUXILIARY FUNCTIONS ##################################################################
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
#'@param WishPdynamicsgraphs binary variable specifying whether the user whishes IRFs and/or GIRFs of risk factors
#'@param WishYieldsgraphs binary variable specifying whether the user whishes IRFs and/or GIRFs of bond yields
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

  # 1) Models estimated separately
  if (any(ModelType == c("JPS original", "JPS global",  "GVAR single"))){
    if(OutputType == "IRF"){
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$IRF[[Economies]]$Factors} # Factors
      else{ PathAdj <- PathsGraphs[[ModelType]]$IRF[[Economies]]$Yields } # Yields

    }else if(OutputType == "GIRF"){
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
} else if (any(ModelType == c("JPS multi", "GVAR multi"))){
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
    if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
      for (g in 1:FacDim){ # Recast the IRFs into a data-frame
        IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]][[Economies]]$Factors[,,g]))
        IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]][[Economies]]$Yields[,,g]))
      }

    # b) Models estimated jointly
    } else if (any(ModelType == c("JPS multi", "GVAR multi"))){
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


#################################################################################################
#'Create folders for storing IRFs and GIRFs
#'
#'@param OutputType available options are "IRF", "GIRF", "IRF Ortho" and "GIRF Ortho"
#'@param WishPdynamicsgraphs binary variable specifing whether the user whishes IRFs and/or GIRFs of risk factors
#'@param WishYieldsgraphs binary variable specifing whether the user whishes IRFs and/or GIRFs of bond yields
#'@param Economies Set of economies that are part of the economic system
#'@param ModelType Desired modem type
#'
#'@keywords internal

FolderPrep_FEVDs <- function(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType){

  ################# 1) SINGLE-COUNTRY CONTRY MODELS #################
  if (any(ModelType == c("JPS original", "JPS global",  "GVAR single"))){

    # A) General Folder Creation
    if(OutputType == "FEVD"){
      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies, "/FEVD", sep=""))
    }else{ dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies, "/GFEVD", sep="")) }


    # a.1) Folders for graph of Factors
    if (WishPdynamicsgraphs == 1){
      if(OutputType == "FEVD"){
        dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
                         "/FEVD/Factors", sep=""))
      }else{
        dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
                         "/GFEVD/Factors", sep=""))
      }
    }

    # a.2) Folders for graph of yields
    if (WishYieldsgraphs == 1){

      if(OutputType == "FEVD"){
        dir.create( paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
                          "/FEVD/Yields", sep=""))
      }else{
        dir.create( paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
                          "/GFEVD/Yields", sep=""))
      }
    }
  } else if(any(ModelType == c("JPS multi", "GVAR multi"))){

    ################# 2) JOINT COUNTRY CONTRY MODELS #################
    # A) General Folder Creation
    if(OutputType == "FEVD"){
      dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/FEVD", sep=""))
    }else{ dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/GFEVD", sep=""))}

    # a.1) Folders for graph of Factors
    if (WishPdynamicsgraphs == 1){
      if(OutputType == "FEVD"){
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Factors", sep=""))
      }else{ dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Factors", sep="")) }
    }

    # a.2) Folders for graph of yields
    if (WishYieldsgraphs == 1){
      if(OutputType == "FEVD"){
        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Yields", sep=""))
      }else{ dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Yields", sep="")) }
    }

    ################# 3) JLL SPECIFIC MODELS #################
  } else{
    # A) General Folder Creation
    if(OutputType == "FEVD"){
      dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/FEVD", sep=""))
    }else if(OutputType == "GFEVD"){ dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/GFEVD", sep=""))}

    ########### NON-ORTHO #####################################
    if(any(OutputType == c("FEVD", "GFEVD"))){
      # a.1) Folders for graph of Factors
      if (WishPdynamicsgraphs == 1){
        if(OutputType == "FEVD"){
          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Factors", sep=""))
        }else{ dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Factors", sep="")) }
      }

      # a.2) Folders for graph of yields
      if (WishYieldsgraphs == 1){
        if(OutputType == "FEVD"){
          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Yields", sep=""))
        }else{ dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Yields", sep="")) }
      }
    } else{

      ########### ORTHO #####################################
      if (WishPdynamicsgraphs == 1){
        if(OutputType == "FEVD Ortho"){
          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Factors/Ortho", sep=""))
        }else{ dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Factors/Ortho", sep=""))  }
      }

      # a.2) Folders for graph of yields
      if (WishYieldsgraphs == 1){
        if(OutputType == "FEVD Ortho"){
          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Yields/Ortho", sep=""))
        }else{  dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Yields/Ortho", sep="")) }
      }
    }
  }
}

##########################################################################################################
#' Build the list of IRF and GIRF for both factors and bond yields
#'
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, GFEVDs and Term Premia
#'@param Economies Economies of the economic system
#'@param ModelType Desired model type
#'@param FEVDhoriz time-horizon of the FEVD and GFEVD
#'@param FacDim dimension of the risk factor vector
#'@param YieldsDim dimension of the model set of yields
#'@param OutputType available option are 'FEVD' and 'GFEVD'
#'
#'@keywords internal

BuildFEVDlist <- function(NumOut, Economies, ModelType, FEVDhoriz, FacDim, YieldsDim, OutputType){

  Horiz <- 0:(FEVDhoriz-1)
  FEVDFactors <- list()
  FEVDYields <- list()

  # 1) Extract IRFs
  if (OutputType == "FEVD"){
    # a) Models estimated individually
    if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
      # Factors
      for (g in 1:FacDim){
        FEVDFactors[[g]] <- data.frame(NumOut$FEVD[[ModelType]][[Economies]]$Factors[,,g])
        INDX <- cbind(index=rownames(FEVDFactors[[g]]), FEVDFactors[[g]])
        FEVDFactors[[g]] <- suppressMessages(reshape2::melt(INDX))
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels= 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim){
      FEVDYields[[g]] <- data.frame(NumOut$FEVD[[ModelType]][[Economies]]$Yields[,,g])
        INDX <- cbind(index=rownames(FEVDYields[[g]]), FEVDYields[[g]])
        FEVDYields[[g]] <- suppressMessages(reshape2::melt(INDX))
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels= 1:(FEVDhoriz))
      }

      # b) Models estimated jointly
    } else if (any(ModelType == c("JPS multi", "GVAR multi"))){
      # Factors
      for (g in 1:FacDim){
        FEVDFactors[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Factors[ , , g])
        FEVDFactors[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDFactors[[g]]), FEVDFactors[[g]])))
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels= 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim){
        FEVDYields[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Yields[ , , g])
        FEVDYields[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDYields[[g]]), FEVDYields[[g]])))
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels= 1:(FEVDhoriz))
        }
    } else{
      # c) JLL-based setups (Non-Orthogonalized version)
      # Factors
      for (g in 1:FacDim){
        FEVDFactors[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Factors$NonOrtho[,,g])
        FEVDFactors[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDFactors[[g]]), FEVDFactors[[g]])))
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels= 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim){
        FEVDYields[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Yields$NonOrtho[,,g])
        FEVDYields[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDYields[[g]]), FEVDYields[[g]])))
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels= 1:(FEVDhoriz))
      }
    }

    # 2) Extract GFEVD
  } else if (OutputType == "GFEVD"){

    # a) Models estimated individually
    if ( any(ModelType == c("JPS original", "JPS global", "GVAR single"))){
      for (g in 1:FacDim){ # Recast the IRFs into a data-frame
        FEVDFactors[[g]] <- data.frame(NumOut$GFEVD[[ModelType]][[Economies]]$Factors[,,g])
        INDX <- cbind(index=rownames(FEVDFactors[[g]]), FEVDFactors[[g]])
        FEVDFactors[[g]] <- suppressMessages(reshape2::melt(INDX))
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels= 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim){
        FEVDYields[[g]] <- data.frame(NumOut$GFEVD[[ModelType]][[Economies]]$Yields[,,g])
        INDX <- cbind(index=rownames(FEVDYields[[g]]), FEVDYields[[g]])
        FEVDYields[[g]] <- suppressMessages(reshape2::melt(INDX))
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels= 1:(FEVDhoriz))
      }
      # b) Models estimated jointly
    } else if (any(ModelType == c("JPS multi", "GVAR multi"))){

      # Factors
      for (g in 1:FacDim){ # Outputs in data-frame format
        FEVDFactors[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Factors[ , , g])
        FEVDFactors[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDFactors[[g]]), FEVDFactors[[g]])))
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels= 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim){
        FEVDYields[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Yields[ , , g])
        FEVDYields[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDYields[[g]]), FEVDYields[[g]])))
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels= 1:(FEVDhoriz))
      }

      # c) JLL-based setups (non-Orthogonalized version)
    } else{
      # Factors
      for (g in 1:FacDim){
        FEVDFactors[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Factors$NonOrtho[,,g])
        FEVDFactors[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDFactors[[g]]), FEVDFactors[[g]])))
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels= 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim){
        FEVDYields[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Yields$NonOrtho[,,g])
        FEVDYields[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDYields[[g]]), FEVDYields[[g]])))
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels= 1:(FEVDhoriz))
      }
    }


    # JLL-based setups (orthogonalized version)
  }else if(OutputType == "FEVD Ortho"){
    # Factors
    for (g in 1:FacDim){ # Recast outputs in a data-frame
      FEVDFactors[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Factors$Ortho[,,g])
      FEVDFactors[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDFactors[[g]]), FEVDFactors[[g]])))
      FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels= 1:(FEVDhoriz))
    }

    # Yields
    for (g in 1:YieldsDim){
      FEVDYields[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Yields$Ortho[,,g])
      FEVDYields[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDYields[[g]]), FEVDYields[[g]])))
      FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels= 1:(FEVDhoriz))
    }
  } else{ # GFEVD Ortho
    # Factors
    for (g in 1:FacDim){ # Recast outputs in a data-frame
      FEVDFactors[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Factors$Ortho[,,g])
      FEVDFactors[[g]] <-  suppressMessages(reshape2::melt(cbind(index=rownames(FEVDFactors[[g]]), FEVDFactors[[g]])))
      FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels= 1:(FEVDhoriz))
    }

    # Yields
    for (g in 1:YieldsDim){
      FEVDYields[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Yields$Ortho[,,g])
      FEVDYields[[g]] <- suppressMessages(reshape2::melt(cbind(index=rownames(FEVDYields[[g]]), FEVDYields[[g]])))
      FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels= 1:(FEVDhoriz))
    }
  }


  return(list(FEVDFactors = FEVDFactors, FEVDYields = FEVDYields))
}
###################################################################################################
#'Generate paths to save IRFs/GIRFs graphs
#'
#'@param OutputType available options are "IRF" and "GIRF"
#'@param ElementType available options are "Factors" and "Yields"
#'@param PathsGraphs desired path to save the graphs
#'@param Economies Economies of the economic system
#'@param ModelType Desired estimated model type
#'
#'@keywords internal

AdjustPathFEVDs <- function(OutputType, ElementType, PathsGraphs, Economies, ModelType){

  # 1) Models estimated seperatly
  if (any(ModelType == c("JPS original", "JPS global",  "GVAR single"))){
    if(OutputType == "FEVD"){
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$FEVD[[Economies]]$Factors} # Factors
      else{ PathAdj <- PathsGraphs[[ModelType]]$FEVD[[Economies]]$Yields } # Yields

    }else{
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$GFEVD[[Economies]]$Factors}
      else{ PathAdj <- PathsGraphs[[ModelType]]$GFEVD[[Economies]]$Yields }
    }

    # 2) Models estimated jointly
  }else{
    if(OutputType == "FEVD"){ # FEVD
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$FEVD$Factors} # Factors
      else{ PathAdj <- PathsGraphs[[ModelType]]$FEVD$Yields } # Yields

    }else if(OutputType == "GFEVD") { # GFEVD
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$GFEVD$Factors} # Factors
      else{ PathAdj <- PathsGraphs[[ModelType]]$GFEVD$Yields } # Yields

      # 3) Exclusively for JLL models
      # FEVD ortho
    }else if(OutputType == "FEVD Ortho"){
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$FEVD[["Factors Ortho"]] } # Factors
      else{ PathAdj <- PathsGraphs[[ModelType]]$FEVD[["Yields Ortho"]] } # Yields

      # GFEVD ortho
    }else {
      if (ElementType == "Factors"){ PathAdj <- PathsGraphs[[ModelType]]$GFEVD[["Factors Ortho"]] } # Factors
      else{ PathAdj <- PathsGraphs[[ModelType]]$GFEVD[["Yields Ortho"]] } # Yields
    }
  }

  return(PathAdj)
}
##################################################################################################
#'Generates graphs for FEVDs and GFEVDs
#'
#'@param OutputType available options are "FEVD", "GFEVD", "FEVD Ortho" and "GFEVD Ortho"
#'@param FEVDlist list of FEVD and GFEVD outputs
#'@param nmVarInt name of variable of interest. Options: "Factors" and "Yields"
#'@param Lab_Fac label of the model factors
#'@param PathsGraphs Path to save graphs
#'
#'@keywords internal

FEVDandGFEVDs_Graphs <- function(OutputType, FEVDlist, nmVarInt, Lab_Fac, PathsGraphs){


  index <- FEVDlist$index
  value <- FEVDlist$value
  variable <- FEVDlist$variable

  p <- ggplot(FEVDlist, aes(x= index, y=value , fill=variable)) +
    geom_bar(stat="identity", width = 0.25) +
    labs(title= paste0(OutputType," - " , nmVarInt)) + theme_classic() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5) )
  suppressMessages(ggplot2::ggsave(p, file=paste0(Lab_Fac, nmVarInt,".png"),
                                   path= PathsGraphs, width= 11, height= 7.98))

}

#############################################################################################
#' Extract list of desired graph features (IRFs anc GIRFs)
#'
#'@param InputsForOutputs List of inputs for outputs
#'@param OutType Output types "FEVD", "GFEVD" and "FEVD Ortho"
#'@param ModelType desired model type
#'
#'@keywords internal

Wished_Graphs_FEVDandGFEVD <- function(InputsForOutputs, OutType, ModelType){

  # 1) JLL models
  if (any(ModelType == c("JLL original",  "JLL No DomUnit",  "JLL joint Sigma"))){
    if(OutType == "FEVD"){
      RiskGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphs$Yields

    } else if (OutType == "GFEVD"){
      RiskGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$Yields

    } else if (OutType == "FEVD Ortho"){
      RiskGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$Yields

    } else{ # GIRF ortho
      RiskGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$Yields
    }

  }else{
    # 2) All other models
    if(OutType == "FEVD"){
      RiskGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphs$Yields

    }else if(OutType == "GFEVD"){
      RiskGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$Yields
    }
  }
  return(list(RiskGraphs = RiskGraphs, YieldGraphs = YieldGraphs))
}

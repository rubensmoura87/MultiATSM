#' Generate the graphical outputs for the selected models (Point estimate)


#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param ModelPara List of model parameter estimates (See the "Optimization" function)
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param InputsForOutputs list containing the desired inputs for the construction of the desired output
#'@param Economies string-vector containing the names of the economies which are part of the economic system
#'@param FactorLabels string-list based which contains the labels of all the variables present in the model
#'
#'@keywords internal



GraphicalOutputs <- function(ModelType, ModelPara, NumOut, InputsForOutputs, Economies, FactorLabels){

# Generate the graph paths and the graph folders
  # Create the folder output folder in which the outputs will be stored
  dir.create(paste(tempdir(), "/Outputs", sep=""))
  dir.create(paste(tempdir(), "/Outputs/",  ModelType, sep=""))
  dir.create(paste(tempdir(), "/Outputs/",  ModelType, "/Point Estimate", sep=""))


  PathsGraphs <- FolderCreationPoint(ModelType, Economies)

# 0) Plot the set of risk factors
  RiskFactorsGraphs(ModelType, ModelPara, Economies, FactorLabels)


# 1) Models for which the estimation is done on a country-by-country basis
if ( "JPS" %in% ModelType || "JPS jointP" %in% ModelType ||  "GVAR sepQ" %in% ModelType){

    for (i in 1:length(Economies)){dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], sep=""))}

  # Model fit
  FitgraphsSep(ModelType, InputsForOutputs[[ModelType]]$Fit$WishGraphs, ModelPara, NumOut, Economies, PathsGraphs)

  # IRF and FEVD
  IRFgraphsSep(ModelType, NumOut, InputsForOutputs[[ModelType]]$IRF$WishGraphs$RiskFactors,
               InputsForOutputs[[ModelType]]$IRF$WishGraphs$Yields, InputsForOutputs[[ModelType]]$IRF$horiz,
               PathsGraphs, Economies)

  FEVDgraphsSep(ModelType, NumOut, InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactors,
                InputsForOutputs[[ModelType]]$FEVD$WishGraphs$Yields, InputsForOutputs[[ModelType]]$FEVD$horiz,
                PathsGraphs, Economies)

  # GIRF and GFEVD
  GIRFgraphsSep(ModelType, NumOut, InputsForOutputs[[ModelType]]$GIRF$WishGraphs$RiskFactors,
                InputsForOutputs[[ModelType]]$GIRF$WishGraphs$Yields, InputsForOutputs[[ModelType]]$GIRF$horiz,
                PathsGraphs, Economies)

  GFEVDgraphsSep(ModelType, NumOut, InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactors,
                 InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$Yields, InputsForOutputs[[ModelType]]$GFEVD$horiz,
                 PathsGraphs, Economies)

  # Term premia decomposition
  TPDecompGraphSep(ModelType, NumOut, ModelPara, InputsForOutputs[[ModelType]]$RiskPremia$WishGraphs,
                   InputsForOutputs$UnitMatYields, Economies, PathsGraphs)

  }

# 2) Models for which the estimation is done on jointly for all the countries

if ( "GVAR jointQ" %in% ModelType || "VAR jointQ" %in% ModelType || "JLL original" %in% ModelType ||
      "JLL NoDomUnit" %in% ModelType || "JLL jointSigma" %in% ModelType){

# Model fit
FitgraphsJoint(ModelType, InputsForOutputs[[ModelType]]$Fit$WishGraphs, ModelPara, NumOut, Economies, PathsGraphs)

# IRF and FEVD
IRFgraphsJoint(ModelType, NumOut, InputsForOutputs[[ModelType]]$IRF$WishGraphs$RiskFactors,
               InputsForOutputs[[ModelType]]$IRF$WishGraphs$Yields, InputsForOutputs[[ModelType]]$IRF$horiz, PathsGraphs)

FEVDgraphsJoint(ModelType, NumOut, InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactors,
                InputsForOutputs[[ModelType]]$FEVD$WishGraphs$Yields, InputsForOutputs[[ModelType]]$FEVD$horiz, PathsGraphs)

# GIRF and GFEVD
GIRFgraphsJoint(ModelType, NumOut, InputsForOutputs[[ModelType]]$GIRF$WishGraphs$RiskFactors,
                InputsForOutputs[[ModelType]]$GIRF$WishGraphs$Yields, InputsForOutputs[[ModelType]]$GIRF$horiz, PathsGraphs)


GFEVDgraphsJoint(ModelType, NumOut, InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactors,
                InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$Yields, InputsForOutputs[[ModelType]]$GFEVD$horiz, PathsGraphs)


TPDecompGraphJoint(ModelType, NumOut, ModelPara, InputsForOutputs[[ModelType]]$RiskPremia$WishGraphs,
                   InputsForOutputs$UnitMatYields, Economies, PathsGraphs)

# 2.1) JLL-based models
if ("JLL original" %in% ModelType  || "JLL NoDomUnit" %in% ModelType || "JLL jointSigma" %in% ModelType){
# IRF and FEVD
IRFgraphsJLLOrtho(ModelType, NumOut, InputsForOutputs[[ModelType]]$IRF$WishGraphsOrtho$RiskFactors,
                  InputsForOutputs[[ModelType]]$IRF$WishGraphsOrtho$Yields, InputsForOutputs[[ModelType]]$IRF$horiz,
                  PathsGraphs)


FEVDgraphsJLLOrtho(ModelType, NumOut, InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$RiskFactors,
                  InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$Yields, InputsForOutputs[[ModelType]]$FEVD$horiz,
                  PathsGraphs)

# GIRF and GFEVD
GIRFgraphsJLLOrtho(ModelType, NumOut, InputsForOutputs[[ModelType]]$GIRF$WishGraphsOrtho$RiskFactors,
                  InputsForOutputs[[ModelType]]$GIRF$WishGraphsOrtho$Yields, InputsForOutputs[[ModelType]]$GIRF$horiz,
                  PathsGraphs)

GFEVDgraphsJLLOrtho(ModelType, NumOut, InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$RiskFactors,
                    InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$Yields, InputsForOutputs[[ModelType]]$GFEVD$horiz,
                    PathsGraphs)

}

}

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
#'@importFrom wrapr seqi
#'
#'
#'@keywords internal



RiskFactorsGraphs <- function(ModelType, ModelOutputs, Economies, FactorLabels){


  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)


  if (ModelType == "JPS" || ModelType == "JPS jointP" || ModelType == "GVAR sepQ"  ){
    X <- c()
    for (i in 1:C){
      if (i ==1){
        X <- ModelOutputs[[ModelType]][[Economies[i]]]$inputs$AllFactors
        Factors <- X
      } else{
        X <- ModelOutputs[[ModelType]][[Economies[i]]]$inputs$AllFactors[-seqi(1,G),]
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

  for (i in seqi(1,G)){
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
      for (i in 1:C){  d <- d +  geom_line(aes_string(y = (nmFactors[IDX[j,i]]), color = shQuote(Economies[i])))}

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
  ggplot2::ggsave(FactorsGraph, filename =paste0("RiskFactors", ".png"), path =Folder2save)

}



######################################################################################################
######################################################################################################
####################### OUTPUTS FOR MODELS IN WHICH THE ESTIMATION ###################################
########################    IS DONE ON A COUNTRY-BY-COUNTRY BASIS       #############################
######################################################################################################
######################################################################################################

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


FitgraphsSep  <- function(ModelType, WishFitgraphs, ModelPara, NumOut, Economies, PathsGraphs){


  if (WishFitgraphs== 0){print(paste(ModelType,": No fit graphs were generated"))
  }else{

    print('################################# Generating Fit graphs #################################' )


    C <- length(Economies)

    for( i in 1:C){
      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], "/Fit", sep=""))


    mat <- (ModelPara[[ModelType]][[Economies[i]]]$inputs$mat)*12
    J <- length(mat)

    ZZ <- ModelPara[[ModelType]][[Economies[i]]]$inputs$AllFactors
    YieldData <- ModelPara[[ModelType]][[Economies[i]]]$inputs$Y
    ModelFit <- NumOut$Fit[[ModelType]][[Economies[i]]]$`Yield Fit`
    ModelImplied <- NumOut$Fit[[ModelType]][[Economies[i]]]$`Yield Model Implied`

    T <- ncol(YieldData)

    # Graph title, legends and axes
    matmonths <- paste(mat, "months - ")
    GraphTitles <- paste(matmonths, Economies[i])

    GraphLegend <- c("Data", "Model Fit", "Model-Implied")


    # Prepare plots
    p <- c()
    plot_list_FIT <- list()
    plot_list_no_legend_FIT <- list()


    for (j in 1:J){
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
      ggplot2::ggsave(subplots2save_FIT, file=paste0("Yields_Fit_", Economies[i], ".png"),
                      path = PathsGraphs[[ModelType]]$Fit[[Economies[i]]])

}

}
}

######################################################################################################
##################################### 2) IRF ########################################################
#####################################################################################################
#' IRFs graphs for ("sep Q" models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param IRFhoriz  single numerical vector conataining the desired horizon of analysis for the IRFs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@importFrom ggplot2 ggplot geom_line labs geom_hline theme ggtitle theme_classic
#'
#'@keywords internal



IRFgraphsSep <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, IRFhoriz,
                         PathsGraphs, Economies) {


  if (WishPdynamicsgraphs==0 & WishYieldsgraphs ==0){print(paste(ModelType,": No IRFs graphs were generated"))
  }else {


    print('################################# Generating IRFs graphs #################################' )

  C <- length(Economies)


    for (i in 1:C){
      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], "/IRF", sep="")) # Folder Creation

        K <- dim(NumOut$IRF[[ModelType]][[Economies[i]]]$Factors)[2] # Total number of risk factors
        J <- dim(NumOut$IRF[[ModelType]][[Economies[i]]]$Yields)[2] # Total number of yields

        IRFFactors <- list()
        IRFYields <- list()
        Horiz <- 0:(IRFhoriz-1)

        for (g in 1:K){ # Recast the IRFs into a data-frame
          IRFFactors[[g]] <- data.frame(cbind(NumOut$IRF[[ModelType]][[Economies[i]]]$Factors[,,g],Horiz))
          IRFYields[[g]] <- data.frame(cbind(NumOut$IRF[[ModelType]][[Economies[i]]]$Yields[,,g], Horiz))
        }

        nmFactors <- names(IRFFactors[[1]]) # Factor names
        nmYields <- names(IRFYields[[1]]) # Yield names

        ########################### IRF FACTORS ########################################################
        if (WishPdynamicsgraphs == 1){
          dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i],
                           "/IRF/Factors", sep=""))

          plot_list <- list()
          for(h in 1:K){
            for (x in 1:K){ # Generate graph-by-graph
              p <- ggplot(IRFFactors[[h]], aes_string(x=nmFactors[K+1], y=nmFactors[x])) + geom_line() +
                labs(title=dimnames(IRFFactors)[[2]][[x]]) +   geom_hline(yintercept=0) +
                ggtitle( nmFactors[x]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                axis.text.y = element_text(size=4))

              plot_list[[x]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # Gather the graphs in sub-plots
            ggplot2::ggsave(subplots, file=paste0("IRFFactors_shock_to_", nmFactors[[h]],"_Model_", Economies[i],".png"),
                            path= PathsGraphs[[ModelType]]$IRF[[Economies[i]]][["Factors"]] ) # Save sub-plots
          }
        }
        ############################ IRF Yields ########################################################
        if (WishYieldsgraphs == 1){
          dir.create( paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/Model ",
                            Economies[i], "/IRF/Yields", sep=""))

          plot_list <- list()

          for(l in 1:K){
            for (x in 1:J){ # Generate graph-by-graph
              p <- ggplot(IRFYields[[l]], aes_string(x=nmYields[J+1], y=nmYields[x])) + geom_line() +
                labs(title=dimnames(IRFYields)[[2]][[x]]) +  geom_hline(yintercept=0) +
                ggtitle( nmYields[x]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=8),
                axis.text.y = element_text(size=4))
              plot_list[[x]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # Gather the graphs in sub-plots
            ggplot2::ggsave(subplots, file=paste0("IRFYields_shock_to_", nmFactors[[l]],"_Model_", Economies[i], ".png"),
                            path= PathsGraphs[[ModelType]]$IRF[[Economies[i]]][["Yields"]]) # Save sub-plots
          }
        }
      }


  }
}

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


  if (WishPdynamicsgraphs == 0 & WishYieldsgraphs ==0){print(paste(ModelType,": No FEVDs graphs were generated"))
  } else {

    print('################################# Generating FEVD graphs #################################' )

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
            ggplot2::ggsave(p, file=paste0("FEVDFactors_", nmFactors[[i]],".png"),
                            path= PathsGraphs[[ModelType]]$FEVD[[Economies[h]]][["Factors"]], width= 11, height= 7.98)
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
              ggplot2::ggsave(p, file=paste0("FEVDYields_", nmYields[i],".png"),
                              path= PathsGraphs[[ModelType]]$FEVD[[Economies[h]]][["Yields"]], width= 11, height= 7.98)
            }
          }
        }

    }
  }

######################################################################################################
##################################### 4) GIRF ########################################################
#####################################################################################################
#' GIRFs graphs for ("sep Q" models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param GIRFhoriz  single numerical vector conataining the desired horizon of analysis for the GIRFs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'
#'@importFrom ggplot2 ggplot aes element_text theme theme_classic labs ggtitle element_blank aes_string geom_line geom_hline
#'
#'@keywords internal


GIRFgraphsSep <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, GIRFhoriz,
                          PathsGraphs, Economies) {


  if (WishPdynamicsgraphs == 0 & WishYieldsgraphs == 0){print(paste(ModelType,": No GIRFs graphs were generated"))
    } else {

      print('################################# Generating GIRFs graphs #################################' )

  C <- length(Economies)

    for (i in 1:C){

      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], "/GIRF", sep=""))

        K <- dim(NumOut$GIRF[[ModelType]][[Economies[i]]]$Factors)[2] # Total number of risk factors
        J <- dim(NumOut$GIRF[[ModelType]][[Economies[i]]]$Yields)[2] # Total number of yields

        GIRFFactors <- list()
        GIRFYields <- list()

        Horiz <- 0:(GIRFhoriz-1)
        for (g in 1:K){ # Recast the outputs to the data-frame format
          GIRFFactors[[g]] <- data.frame(cbind(NumOut$GIRF[[ModelType]][[Economies[i]]]$Factors[,,g], Horiz))
          GIRFYields[[g]] <- data.frame(cbind(NumOut$GIRF[[ModelType]][[Economies[i]]]$Yields[,,g], Horiz))
        }

        nmFactors <- names(GIRFFactors[[1]]) # Factor names
        nmYields <- names(GIRFYields[[1]]) # Yield names

        ########################### GIRF FACTORS ########################################################
        if (WishPdynamicsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], "/GIRF/Factors", sep=""))

          plot_list <- list()
          for(h in 1:K){
            for (x in 1:K){ # Generate graph-by-graph
              p <- ggplot(GIRFFactors[[h]], aes_string(x=nmFactors[K+1], y=nmFactors[x])) + geom_line() +
                labs(title=dimnames(GIRFFactors)[[2]][[x]]) +   geom_hline(yintercept=0) +
                theme(plot.title = element_text(size=6))+ ggtitle( nmFactors[x]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                      axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=6),
                      axis.text.y = element_text(size=4))

              plot_list[[x]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # subplots
            ggplot2::ggsave(subplots, file=paste0("GIRFFactors_shock_to_", nmFactors[[h]], "_Model_", Economies[i], ".png"),
                            path= PathsGraphs[[ModelType]]$GIRF[[Economies[i]]][["Factors"]])
          }
        }
        ############################ GIRF Yields ########################################################
        if (WishYieldsgraphs == 1){

          dir.create(paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], "/GIRF/Yields", sep=""))

          plot_list <- list()
          for(l in 1:K){
            for (x in 1:J){ # Generate graph-by-graph
              p <- ggplot(GIRFYields[[l]], aes_string(x=nmYields[J+1], y=nmYields[i])) + geom_line() +
                labs(title=dimnames(GIRFYields)[[2]][[x]]) +  geom_hline(yintercept=0) +
                ggtitle( nmYields[x]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
                      axis.title.y=element_blank(), axis.text.x = element_text(size=6), axis.text.y = element_text(size=4) )
              plot_list[[x]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # subplots
            ggplot2::ggsave(subplots, file=paste0("GIRFYields_shock_to_", nmFactors[[l]],"_Model_", Economies[i], ".png"),
                            path= PathsGraphs[[ModelType]]$GIRF[[Economies[i]]][["Yields"]])
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

  if (WishPdynamicsgraphs==0 & WishYieldsgraphs==0){print(paste(ModelType,": No GFEVDs graphs were generated"))
    } else{

      print('################################# Generating GFEVDs graphs #################################' )

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
            ggplot2::ggsave(p, file=paste0("GFEVDFactors_", nmFactors[[i]],".png"),
                            path= PathsGraphs[[ModelType]]$GFEVD[[Economies[h]]][["Factors"]], width= 11, height= 7.98)
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
              ggplot2::ggsave(p, file=paste0("GFEVDYields_", nmYields[i],".png"),
                              path= PathsGraphs[[ModelType]]$GFEVD[[Economies[h]]][["Yields"]], width= 11, height= 7.98)
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



  if (WishRPgraphs ==0){print(paste(ModelType,": No term-premia graphs were generated"))
  } else {

    print('################################# Generating term premia graphs #################################' )


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
      ggplot2::ggsave(subplots2save, file=paste0("TermPremia_", Economies[i],"_", ModelType, ".png"),
                      path = PathsGraphs[[ModelType]]$TermPremia[[Economies[i]]])

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
##################################### 1) FIT ########################################################
#####################################################################################################
#' Model fit graphs for ("joint Q" models)
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
#'@keywords internal


FitgraphsJoint  <- function(ModelType, WishFitgraphs, ModelPara, NumOut, Economies, PathsGraphs){


  if (WishFitgraphs== 0){print(paste(ModelType,": No fit graphs were generated"))
  }else{

    print('################################# Generating Fit graphs #################################' )


    dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/Fit", sep=""))  # Folder creation

  mat <- (ModelPara[[ModelType]]$inputs$mat)*12
  C <- length(Economies)
  J <- length(mat)
  CJ <- C*J

  ZZ <- ModelPara[[ModelType]]$inputs$AllFactors
  YieldData <- ModelPara[[ModelType]]$inputs$Y
  ModelFit <- NumOut$Fit$`Yield Fit`
  ModelImplied <- NumOut$Fit$`Yield Model Implied`

  T <- ncol(YieldData)

  # Graph title, legends and axes
  matmonths <- paste(mat, "months - ")

  for (i in 1:C){
    temp <- paste(matmonths, Economies[i])
    if ( i==1){ GraphTitles <- temp
    } else{
      GraphTitles <- c(GraphTitles, temp)
    }
  }

  GraphLegend <- c("Data", "Model Fit", "Model-Implied")


  # Prepare plots
  p <- c()
  plot_list_FIT <- list()
  plot_list_no_legend_FIT <- list()


  for (i in 1:CJ){
    Data <- YieldData[i,-1]
    Modfit <- ModelFit[i,-1]
    ModImp <- ModelImplied[i,-1]
    TimeSpan <- 1:(T-1)

    YieldCompar <- data.frame(Data, Modfit, ModImp, TimeSpan)


    # Graph: Fit x Data
    p <- ggplot(data = YieldCompar, aes(x= TimeSpan )) +
      geom_line(aes_string(y = Data, color = shQuote(GraphLegend[1])), size = 0.5) +
      geom_line(aes_string(y = Modfit, color = shQuote(GraphLegend[2])), size = 0.5) +
      geom_line(aes_string(y = ModImp, color = shQuote(GraphLegend[3])), linetype = "dashed") +
      labs(color = "Legend") + ggtitle( GraphTitles[i]) +  theme_classic() +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
            axis.title.y=element_blank(), axis.text.x = element_text(size=6) ) +
      geom_hline(yintercept=0)


    plot_list_FIT[[i]] <- p # All plots with legends

    plot_list_no_legend_FIT[[i]] <- p + theme(legend.position = "none")  # Hide legends from all plots

  }

  # Generate legend
  # Legend settings:
  LegendSubPlot_FIT <- cowplot::get_legend(plot_list_FIT[[1]] +
                                    theme(legend.direction = "horizontal", legend.position="bottom", legend.justification="center",
                                          legend.title = element_text(size=10, face="bold"),
                                          legend.text = element_text( size = 8),
                                          legend.box.background = element_rect(colour = "black")) )


  # Build the graph subplots:

  Idx0 <- 0
  for (j in 1:C){
    Idx <- Idx0
    IdxGrpahs <- (Idx+1):(Idx+J)

    subplots_FIT <- cowplot::plot_grid(plotlist= plot_list_no_legend_FIT[IdxGrpahs], ncol=3)
    subplots2save_FIT <- cowplot::plot_grid(LegendSubPlot_FIT, subplots_FIT , ncol=1,  rel_heights = c(.1, 1))
    ggplot2::ggsave(subplots2save_FIT, file=paste0("Yields_Fit_", Economies[j], ".png"),
                    path = PathsGraphs[[ModelType]]$Fit)

    Idx0 <- (Idx+J)
  }


}
}

######################################################################################################
##################################### 2) IRF ########################################################
#####################################################################################################
#' IRFs graphs for ("joint Q" models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param IRFhoriz  single numerical vector conataining the desired horizon of analysis for the IRFs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'
#'@importFrom ggplot2 ggplot geom_line labs geom_hline theme ggtitle theme_classic
#'
#'@keywords internal



IRFgraphsJoint <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs,IRFhoriz, PathsGraphs){



    if (WishPdynamicsgraphs== 0 & WishYieldsgraphs==0){print(paste(ModelType,": No IRFs graphs were generated"))
      }else{


        print('################################# Generating IRFs graphs #################################' )

      dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF", sep="")) # Folder Creation

      IRFFactors <- list()
      IRFYields <- list()

      Horiz <- 0:(IRFhoriz-1)


        if ( ModelType == "JLL original"  || ModelType== "JLL NoDomUnit" || ModelType =="JLL jointSigma"){
          K <- dim(NumOut$IRF[[ModelType]]$Factors$NonOrtho)[2]  # Total number of risk factors
          CJ <- dim(NumOut$IRF[[ModelType]]$Yields$NonOrtho)[2]  # Total number of yields
          for (i in 1:K){ # Outputs in data-frame format
            IRFFactors[[i]] <- data.frame(cbind(NumOut$IRF[[ModelType]]$Factors$NonOrtho[,,i], Horiz))
            IRFYields[[i]] <- data.frame(cbind(NumOut$IRF[[ModelType]]$Yields$NonOrtho[,,i], Horiz))
          }

        }else{
          K <- dim(NumOut$IRF[[ModelType]]$Factors)[2]  # Total number of risk factors
          CJ <- dim(NumOut$IRF[[ModelType]]$Yields)[2]  # Total number of yields
          for (i in 1:K){ # Outputs in data-frame format
            IRFFactors[[i]] <- data.frame(cbind(NumOut$IRF[[ModelType]]$Factors[,,i], Horiz))
            IRFYields[[i]] <- data.frame(cbind(NumOut$IRF[[ModelType]]$Yields[,,i], Horiz))
          }
        }

        nmFactors <- names(IRFFactors[[1]]) # Factor names
        nmYields <- names(IRFYields[[1]]) # Yield names

        ########################### IRF FACTORS ########################################################
        if (WishPdynamicsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF/Factors", sep="")) # Folder Creation

            plot_list <- list()
          for(h in 1:K){
            for (i in 1:K){ # Generate graph-by-graph
              p <- ggplot(IRFFactors[[h]], aes_string(x=nmFactors[K+1], y=nmFactors[i])) + geom_line() +
                labs(title=dimnames(IRFFactors)[[2]][[i]]) +   geom_hline(yintercept=0) +
                ggtitle( nmFactors[i]) +  theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
                      axis.title.y=element_blank(), axis.text.x = element_text(size=6),
                      axis.text.y = element_text(size=4) )
              plot_list[[i]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # store in subplots
            ggplot2::ggsave(subplots, file=paste0("IRFFactors_shock_to_", nmFactors[[h]],".png"),
                            path= PathsGraphs[[ModelType]]$IRF[["Factors"]])
          }
        }
        ############################ IRF Yields ########################################################
        if (WishYieldsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF/Yields", sep="")) # Folder Creation

          plot_list <- list()
          for(l in 1:K){
            for (i in 1:CJ){ # Generate graph-by-graph
              p <- ggplot(IRFYields[[l]], aes_string(x=nmYields[CJ+1], y=nmYields[i])) + geom_line() +
                labs(title=dimnames(IRFYields)[[2]][[i]]) +  geom_hline(yintercept=0) +
                ggtitle( nmYields[i]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
                      axis.title.y=element_blank(), axis.text.x = element_text(size=6),
                      axis.text.y = element_text(size=4))
              plot_list[[i]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # store in subplots
            ggplot2::ggsave(subplots, file=paste0("IRFYields_shock_to_", nmFactors[[l]],".png"),
                            path= PathsGraphs[[ModelType]]$IRF[["Yields"]])
          }

        }

    }

}

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



    if (WishPdynamicsgraphs ==0 & WishYieldsgraphs==0){print(paste(ModelType,": No FEVDs graphs were generated"))
      }else {


        print('################################# Generating FEVDs graph #################################' )

        dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD", sep="")) # Folder Creation

 ######################################### Factors #########################################################
        if (WishPdynamicsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Factors", sep="")) # Folder Creation

          FEVDFactors <- list()


          if ( ModelType == "JLL original"  || ModelType== "JLL NoDomUnit" || ModelType =="JLL jointSigma"){
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
            ggplot2::ggsave(p, file=paste0("FEVDFactors_", nmFactors[[i]],".png"),
                            path= PathsGraphs[[ModelType]]$FEVD[["Factors"]], width= 11, height= 7.98)
          }

        }
        ############################################ Yields #########################################################
        if (WishYieldsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/FEVD/Yields", sep="")) # Folder Creation

          FEVDYields <- list()

          if (ModelType == "JLL original"  || ModelType == "JLL NoDomUnit" || ModelType =="JLL jointSigma"){

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
            ggplot2::ggsave(p, file=paste0("FEVDYields_", nmYields[i],".png"),
                            path= PathsGraphs[[ModelType]]$FEVD[["Yields"]], width= 11, height= 7.98)
          }
        }


  }
}

######################################################################################################
##################################### 4) GIRF ########################################################
#####################################################################################################
#' GIRFs graphs for ("joint Q" models)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param GIRFhoriz  single numerical vector conataining the desired horizon of analysis for the GIRFs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'
#'@importFrom ggplot2 ggplot aes element_text theme theme_classic labs ggtitle element_blank aes_string geom_line geom_hline
#'@keywords internal





GIRFgraphsJoint <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, GIRFhoriz, PathsGraphs){


    if (WishPdynamicsgraphs ==0 & WishYieldsgraphs==0){print(paste(ModelType,": No GIRFs graphs were generated"))
      } else{

        print('################################# Generating GIRFs graphs #################################' )

        dir.create(paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/GIRF", sep="")) # Folder Creation

      GIRFFactors <- list()
      GIRFYields <- list()

      Horiz <- 0:(GIRFhoriz-1)


        if ( ModelType == "JLL original"  || ModelType == "JLL NoDomUnit" || ModelType =="JLL jointSigma"){
          K <- dim(NumOut$GIRF[[ModelType]]$Factors$NonOrtho)[2]  # Total number of risk factors
          CJ <- dim(NumOut$GIRF[[ModelType]]$Yields$NonOrtho)[2]  # Total number of yields
          for (i in 1:K){
            GIRFFactors[[i]] <- data.frame(cbind(NumOut$GIRF[[ModelType]]$Factors$NonOrtho[,,i], Horiz))
            GIRFYields[[i]] <- data.frame(cbind(NumOut$GIRF[[ModelType]]$Yields$NonOrtho[,,i], Horiz))
          }

        }else{
          K <- dim(NumOut$GIRF[[ModelType]]$Factors)[2]  # Total number of risk factors
          CJ <- dim(NumOut$GIRF[[ModelType]]$Yields)[2]  # Total number of yields
          for (i in 1:K){
            GIRFFactors[[i]] <- data.frame(cbind(NumOut$GIRF[[ModelType]]$Factors[,,i], Horiz))
            GIRFYields[[i]] <- data.frame(cbind(NumOut$GIRF[[ModelType]]$Yields[,,i], Horiz))
          }
        }

        nmFactors <- names(GIRFFactors[[1]]) # Factor names
        nmYields <- names(GIRFYields[[1]]) # Yield names

        ########################### GIRF FACTORS ########################################################
        if (WishPdynamicsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GIRF/Factors", sep="")) # Folder Creation

          plot_list <- list()
          for(h in 1:K){
            for (i in 1:K){ # Generate graph-by-graph
              p <- ggplot(GIRFFactors[[h]], aes_string(x=nmFactors[K+1], y=nmFactors[i])) + geom_line() +
                labs(title=dimnames(GIRFFactors)[[2]][[i]]) +   geom_hline(yintercept=0) +
                ggtitle( nmFactors[i]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
                      axis.title.y=element_blank(), axis.text.x = element_text(size=6),
                      axis.text.y = element_text(size=4))
              plot_list[[i]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # sub-plots
            ggplot2::ggsave(subplots, file=paste0("GIRFFactors_shock_to_", nmFactors[[h]],".png"),
                            path= PathsGraphs[[ModelType]]$GIRF[["Factors"]])
          }
        }
        ############################ IRF Yields ########################################################
        if (WishYieldsgraphs == 1){

          dir.create( paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/GIRF/Yields", sep="")) # Folder Creation

          plot_list <- list()
          for(l in 1:K){
            for (i in 1:CJ){ # Generate graph-by-graph
              p <- ggplot(GIRFYields[[l]], aes_string(x=nmYields[CJ+1], y=nmYields[i])) + geom_line() +
                labs(title=dimnames(GIRFYields)[[2]][[i]]) +  geom_hline(yintercept=0) +
                ggtitle( nmYields[i]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                      axis.title.x=element_blank(), axis.title.y=element_blank(),
                      axis.text.x = element_text(size=6), axis.text.y = element_text(size=4))
              plot_list[[i]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # sub-plots
            ggplot2::ggsave(subplots, file=paste0("GIRFYields_shock_to_", nmFactors[[l]],".png"),
                            path= PathsGraphs[[ModelType]]$GIRF[["Yields"]])
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


    if (WishPdynamicsgraphs ==0 & WishYieldsgraphs ==0){print(paste(ModelType,": No GFEVDs graphs were generated"))
      } else {


        print('################################# Generating GFEVD graphs #################################' )


        dir.create( paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/GFEVD", sep="")) # Folder Creation
        ######################################### Factors #########################################################
        if (WishPdynamicsgraphs == 1){
          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Factors", sep="")) # Folder Creation

          GFEVDFactors <- list()


          if ( ModelType == "JLL original"  || ModelType == "JLL NoDomUnit" || ModelType =="JLL jointSigma"){
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
            ggplot2::ggsave(p, file=paste0("GFEVDFactors_", nmFactors[[i]],".png"),
                            path= PathsGraphs[[ModelType]]$GFEVD[["Factors"]], width= 11, height= 7.98)
          }

        }
        ############################################ Yields #########################################################
        if (WishYieldsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GFEVD/Yields", sep="")) # Folder Creation

          GFEVDYields <- list()

          if (ModelType == "JLL original"  || ModelType== "JLL NoDomUnit" || ModelType =="JLL jointSigma"){
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
            ggplot2::ggsave(p, file=paste0("GFEVDYields_", nmYields[i],".png"),
                            path= PathsGraphs[[ModelType]]$GFEVD[["Yields"]], width= 11, height= 7.98)
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



  if (WishRPgraphs ==0){print(paste(ModelType,": No term-premia graphs were generated"))
  } else {

    print('################################# Generating term premia graphs #################################' )


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

      ggplot2::ggsave(subplots2save, file=paste0("TermPremia_", Economies[i],"_", ModelType, ".png"),
                path = PathsGraphs[[ModelType]]$TermPremia)

        }

}
}

##########################################################################################################
##################### GRAPHICAL OUTPUTS FOR THE ORTHOGONALIZED FACTORS OF JLL ############################
##########################################################################################################


######################################################################################################
##################################### 1) IRF ########################################################
#####################################################################################################
#' IRFs graphs for orthogonalized risk factors of JLL-based models
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param IRFhoriz  single numerical vector conataining the desired horizon of analysis for the IRFs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'
#'@importFrom ggplot2 ggplot geom_line labs geom_hline theme ggtitle theme_classic
#'
#'@keywords internal



IRFgraphsJLLOrtho <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, IRFhoriz, PathsGraphs) {


  Horiz <- 0:(IRFhoriz-1)

  if (WishPdynamicsgraphs== 0 & WishYieldsgraphs==0){
    print(paste(ModelType,": No IRFs graphs were generated (orthogonolized version)"))
  }else{

    print('################################# Generating IRFs-Ortho graphs #################################' )

      K <- dim(NumOut$IRF[[ModelType]]$Factors$Ortho)[2]  # Total number of risk factors
      CJ <- dim(NumOut$IRF[[ModelType]]$Yields$Ortho)[2]  # Total number of yields

        IRFFactors <- list()
        IRFYields <- list()
        for (i in 1:K){ # Recast outputs in a data-frame
          IRFFactors[[i]] <- data.frame(cbind(NumOut$IRF[[ModelType]]$Factors$Ortho[,,i], Horiz))
          IRFYields[[i]] <- data.frame(cbind(NumOut$IRF[[ModelType]]$Yields$Ortho[,,i], Horiz))
        }

        nmFactors <- names(IRFFactors[[1]]) # Factor names
        nmYields <- names(IRFYields[[1]]) # Yield names

        ########################### IRF FACTORS ########################################################
        if (WishPdynamicsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF/Factors/Ortho", sep="")) # Folder Creation

          plot_list <- list()
          for(h in 1:K){
            for (i in 1:K){ # generate graph-by-graph
              p <- ggplot(IRFFactors[[h]], aes_string(x=nmFactors[K+1], y=nmFactors[i])) + geom_line() +
                labs(title=dimnames(IRFFactors)[[2]][[i]]) +   geom_hline(yintercept=0) +
                ggtitle( nmFactors[i]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
                       axis.title.y=element_blank(), axis.text.x = element_text(size=6),
                      axis.text.y = element_text(size=4))
              plot_list[[i]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # sub-plots
            ggplot2::ggsave(subplots, file=paste0("IRFFactors_shock_to_", nmFactors[[h]],"ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$IRF[["Factors Ortho"]])
          }
        }
        ############################ IRF Yields ########################################################
        if (WishYieldsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/IRF/Yields/Ortho", sep="")) # Folder Creation

          plot_list <- list()
          for(l in 1:K){
            for (i in 1:CJ){ # generate graph-by-graph
              p <- ggplot(IRFYields[[l]], aes_string(x=nmYields[CJ+1], y=nmYields[i])) + geom_line() +
                labs(title=dimnames(IRFYields)[[2]][[i]]) +  geom_hline(yintercept=0) +
                ggtitle( nmYields[i]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5), axis.title.x=element_blank(),
                      axis.title.y=element_blank(), axis.text.x = element_text(size=6),
                      axis.text.y = element_text(size=4))
              plot_list[[i]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # sub-plots
            ggplot2::ggsave(subplots, file=paste0("IRFYields_shock_to_", nmFactors[[l]], "ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$IRF[["Yields Ortho"]])
          }

        }


  }
}
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
    print(paste(ModelType,": No FEVDs graphs were generated (orthogonolized version)"))
  }else{

    print('################################# Generating FEVDs-Ortho graphs #################################' )

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
            ggplot2::ggsave(p, file=paste0("FEVDFactors_", nmFactors[[i]], "ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$FEVD[["Factors Ortho"]], width= 11, height= 7.98)
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
            ggplot2::ggsave(p, file=paste0("FEVDYields_", nmYields[i], "ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$FEVD[["Yields Ortho"]], width= 11, height= 7.98)
          }
        }

    }

}

######################################################################################################
##################################### 3) GIRF ########################################################
#####################################################################################################
#' GIRFs graphs for orthogonalized risk factors of JLL-based models

#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, and GFEVDs
#'@param WishPdynamicsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param WishYieldsgraphs binary variable: set 1, if the user wishes graphs to be generated; or set 0, otherwise
#'@param GIRFhoriz  single numerical vector conataining the desired horizon of analysis for the GIRFs
#'@param PathsGraphs Path of the folder in which the graphs will be saved
#'
#'
#'@importFrom ggplot2 ggplot aes element_text theme theme_classic labs ggtitle element_blank aes_string geom_line geom_hline
#'
#'@keywords internal



GIRFgraphsJLLOrtho <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, GIRFhoriz, PathsGraphs){


  Horiz <- 0:(GIRFhoriz-1)

  if (WishPdynamicsgraphs== 0 & WishYieldsgraphs==0){
    print(paste(ModelType,": No GIRFs graphs were generated (orthogonolized version)"))
  }else{

    print('################################# Generating GIRFs-Ortho graphs #################################' )

        K <- dim(NumOut$GIRF[[ModelType]]$Factors$Ortho)[2]  # Total number of risk factors
        CJ <- dim(NumOut$GIRF[[ModelType]]$Yields$Ortho)[2]  # Total number of yields

        GIRFFactors <- list()
        GIRFYields <- list()
        for (i in 1:K){ # Recast the outputs in data-frame format
          GIRFFactors[[i]] <- data.frame(cbind(NumOut$GIRF[[ModelType]]$Factors$Ortho[,,i], Horiz))
          GIRFYields[[i]] <- data.frame(cbind(NumOut$GIRF[[ModelType]]$Yields$Ortho[,,i], Horiz))
        }

        nmFactors <- names(GIRFFactors[[1]]) # Factor names
        nmYields <- names(GIRFYields[[1]]) # Yield names

        ########################### GIRF FACTORS ########################################################
        if (WishPdynamicsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GIRF/Factors/Ortho", sep="")) # Folder Creation

          plot_list <- list()
          for(h in 1:K){
            for (i in 1:K){ # Generate graph-by-graph
              p <- ggplot(GIRFFactors[[h]], aes_string(x=nmFactors[K+1], y=nmFactors[i])) + geom_line() +
                labs(title=dimnames(GIRFFactors)[[2]][[i]]) +   geom_hline(yintercept=0) +
                ggtitle( nmFactors[i]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                      axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=6) )
              plot_list[[i]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # sub-plots
            ggplot2::ggsave(subplots, file=paste0("GIRFFactors_shock_to_", nmFactors[[h]],"ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$GIRF[["Factors Ortho"]])
          }
        }
        ############################ IRF Yields ########################################################
        if (WishYieldsgraphs == 1){

          dir.create( paste(tempdir(), "/Outputs/", ModelType, "/Point Estimate/GIRF/Yields/Ortho", sep="")) # Folder Creation

          plot_list <- list()
          for(l in 1:K){
            for (i in 1:CJ){ # Generate graph-by-graph
              p <- ggplot(GIRFYields[[l]], aes_string(x=nmYields[CJ+1], y=nmYields[i])) + geom_line() +
                labs(title=dimnames(GIRFYields)[[2]][[i]]) +  geom_hline(yintercept=0) +
                ggtitle( nmYields[i]) + theme_classic() +
                theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
                      axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=6),
                      axis.text.y = element_text(size=4))
              plot_list[[i]] <- p
            }
            subplots <- cowplot::plot_grid(plotlist= plot_list, ncol=3) # sub-plots
            ggplot2::ggsave(subplots, file=paste0("GIRFYields_shock_to_", nmFactors[[l]], "ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$GIRF[["Yields Ortho"]])
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
    print(paste(ModelType,": No GFEVDs graphs were generated (orthogonolized version)"))
  }else{

    print('################################# Generating GEFVDs-Ortho graph #################################' )

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
            ggplot2::ggsave(p, file=paste0("GFEVDFactors_", nmFactors[[i]], "ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$GFEVD[["Factors Ortho"]], width= 11, height= 7.98)
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
            ggplot2::ggsave(p, file=paste0("GFEVDYields_", nmYields[i], "ORTHO",".png"),
                            path= PathsGraphs[[ModelType]]$GFEVD[["Yields Ortho"]], width= 11, height= 7.98)
          }
        }

  }

}

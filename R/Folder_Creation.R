#' Creates the folders and the path in which the graphical outputs are stored (ponit estimate version)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal
#'


FolderCreationPoint <- function(ModelType, Economies){

C <- length(Economies)

PathsGraphs <- list()

OutputTypeSet <- c("Fit", "IRF", "FEVD", "GIRF", "GFEVD")
VarTypeSet <- c("Factors", "Yields")


  for (h in 1:length(OutputTypeSet)){
    PathsGraphs[[ModelType]][[OutputTypeSet[h]]] <- list()
  }


if (ModelType == "JPS" || ModelType == "JPS jointP" || ModelType == "GVAR sepQ"){
for (h in 1:length(OutputTypeSet)){
  for (v in 1:length(VarTypeSet)){
  for (i in 1:C){
    if (h==1){
    PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[Economies[i]]]  <- paste("Outputs/", ModelType, "/Point Estimate/Model ",
                                                                           Economies[i], "/", OutputTypeSet[h], sep="")
    } else{
      PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[Economies[i]]][[VarTypeSet[v]]]  <- paste("Outputs/", ModelType, "/Point Estimate/Model ",
                                                                                              Economies[i], "/", OutputTypeSet[h], "/",
                                                                                              VarTypeSet[v], sep="")
       }
    }
}
}
} else{

  for (h in 1:length(OutputTypeSet)){
    for (v in 1:length(VarTypeSet)){
        if (h==1){
          PathsGraphs[[ModelType]][[OutputTypeSet[h]]]  <- paste("Outputs/", ModelType, "/Point Estimate/",
                                                                  OutputTypeSet[h], sep="")
        } else{
          PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[VarTypeSet[v]]]  <- paste("Outputs/", ModelType, "/Point Estimate/",
                                                                                  OutputTypeSet[h], "/", VarTypeSet[v], sep="")

        }


  if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma"){
    if (h !=1){
    OrthoLabel <- paste(VarTypeSet,"Ortho")
    PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[OrthoLabel[v]]] <- paste("Outputs/", ModelType, "/Point Estimate/",
                                                                           OutputTypeSet[h], "/", VarTypeSet[v],
                                                                           "/Ortho", sep="")
    }
    }
    }
  }
}


return(PathsGraphs)
}


#######################################################################################################################
#' Creates the folders and the path in which the graphical outputs are stored (Bootstrap version)
#'
#'@param ModelType a string-vector containing the label of the model to be estimated
#'@param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#'@keywords internal
#'


FolderCreationBoot <- function(ModelType, Economies){

  C <- length(Economies)

  PathsGraphs <- list()

  OutputTypeSet <- c("Fit", "IRF", "FEVD", "GIRF", "GFEVD")
  VarTypeSet <- c("Factors", "Yields")


  for (h in 1:length(OutputTypeSet)){
    PathsGraphs[[ModelType]][[OutputTypeSet[h]]] <- list()
  }


  if (ModelType == "JPS" || ModelType == "JPS jointP" || ModelType == "GVAR sepQ"){
    for (h in 1:length(OutputTypeSet)){
      for (v in 1:length(VarTypeSet)){
        for (i in 1:C){
          if (h==1){
            PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[Economies[i]]]  <- paste("Outputs/", ModelType, "/Bootstrap/Model ",
                                                                                   Economies[i], "/", OutputTypeSet[h], sep="")
          } else{
            PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[Economies[i]]][[VarTypeSet[v]]]  <- paste("Outputs/", ModelType, "/Bootstrap/Model ",
                                                                                                    Economies[i], "/", OutputTypeSet[h], "/",
                                                                                                    VarTypeSet[v], sep="")
          }
        }
      }
    }

  } else{

    for (h in 1:length(OutputTypeSet)){
      for (v in 1:length(VarTypeSet)){
        if (h==1){
          PathsGraphs[[ModelType]][[OutputTypeSet[h]]]  <- paste("Outputs/", ModelType, "/Bootstrap/",
                                                                 OutputTypeSet[h], sep="")
        } else{
          PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[VarTypeSet[v]]]  <- paste("Outputs/", ModelType, "/Bootstrap/",
                                                                                  OutputTypeSet[h], "/", VarTypeSet[v], sep="")

        }


        if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma"){
          if (h !=1){
            OrthoLabel <- paste(VarTypeSet,"Ortho")
            PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[OrthoLabel[v]]] <- paste("Outputs/", ModelType, "/Bootstrap/",
                                                                                   OutputTypeSet[h], "/", VarTypeSet[v],
                                                                                   "/Ortho", sep="")
          }
        }
      }
    }
  }


  return(PathsGraphs)
}


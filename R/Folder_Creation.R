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

OutputTypeSet <- c("Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremia")
VarTypeSet <- c("Factors", "Yields")


  for (h in 1:length(OutputTypeSet)){PathsGraphs[[ModelType]][[OutputTypeSet[h]]] <- list()}


if (ModelType == "JPS" || ModelType == "JPS jointP" || ModelType == "GVAR sepQ"){
for (h in 1:length(OutputTypeSet)){
  for (v in 1:length(VarTypeSet)){
  for (i in 1:C){
    if (h==1){
    PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[Economies[i]]]  <- paste(tempdir(), "/Outputs/", ModelType,
                                                                           "/Point Estimate/Model ", Economies[i], "/",
                                                                           OutputTypeSet[h], sep="")
    }
    if (h>=2 & h <=5){
      PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[Economies[i]]][[VarTypeSet[v]]]  <- paste(tempdir(), "/Outputs/",
                                                                                              ModelType, "/Point Estimate/Model ",
                                                                                              Economies[i], "/", OutputTypeSet[h], "/",
                                                                                              VarTypeSet[v], sep="")
       }
    if (h==6){
      PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[Economies[i]]]  <- paste(tempdir(), "/Outputs/", ModelType,
                                                                             "/Point Estimate/Model ", Economies[i], "/",
                                                                             OutputTypeSet[h], sep="")
          }
      }
}
}
  }else{

  for (h in 1:length(OutputTypeSet)){
    for (v in 1:length(VarTypeSet)){
        if (h==1){
          PathsGraphs[[ModelType]][[OutputTypeSet[h]]]  <- paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/",
                                                                  OutputTypeSet[h], sep="")
          }
      if (h>=2 & h <=5){
          PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[VarTypeSet[v]]]  <- paste(tempdir(),"/Outputs/", ModelType,
                                                                                  "/Point Estimate/", OutputTypeSet[h],
                                                                                  "/", VarTypeSet[v], sep="")

          if (ModelType == "JLL original" || ModelType == "JLL NoDomUnit" || ModelType == "JLL jointSigma"){

              OrthoLabel <- paste(VarTypeSet,"Ortho")
              PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[OrthoLabel[v]]] <- paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/",
                                                                                     OutputTypeSet[h], "/", VarTypeSet[v],
                                                                                     "/Ortho", sep="")
          }
      }

      if (h==6){
        PathsGraphs[[ModelType]][[OutputTypeSet[h]]]  <- paste(tempdir(),"/Outputs/", ModelType, "/Point Estimate/",
                                                               OutputTypeSet[h], sep="")
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

  OutputTypeSet <- c("Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremia")
  VarTypeSet <- c("Factors", "Yields")


  for (h in 1:length(OutputTypeSet)){
    PathsGraphs[[ModelType]][[OutputTypeSet[h]]] <- list()
  }


  if (any(ModelType == c("JPS", "JPS jointP", "GVAR sepQ"))){
    for (h in 1:length(OutputTypeSet)){
      for (v in 1:length(VarTypeSet)){
        for (i in 1:C){
          if (h==1){
            PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[Economies[i]]]  <- paste(tempdir(),"/Outputs/", ModelType, "/Bootstrap/Model ",
                                                                                   Economies[i], "/", OutputTypeSet[h], sep="")
          }
          if (h>=2 & h <=5){
            PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[Economies[i]]][[VarTypeSet[v]]]  <- paste(tempdir(),"/Outputs/", ModelType, "/Bootstrap/Model ",
                                                                                                    Economies[i], "/", OutputTypeSet[h], "/",
                                                                                                    VarTypeSet[v], sep="")
          }
          if (h==6){
            PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[Economies[i]]]  <- paste(tempdir(),"/Outputs/", ModelType, "/Bootstrap/Model ",
                                                                                   Economies[i], "/", OutputTypeSet[h], sep="")
          }
        }
      }
    }

  } else{

    for (h in 1:length(OutputTypeSet)){
      for (v in 1:length(VarTypeSet)){
        if (h==1){
          PathsGraphs[[ModelType]][[OutputTypeSet[h]]]  <- paste(tempdir(),"/Outputs/", ModelType, "/Bootstrap/",
                                                                 OutputTypeSet[h], sep="")
        }
        if (h>=2 & h <=5){
          PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[VarTypeSet[v]]]  <- paste(tempdir(),"/Outputs/", ModelType, "/Bootstrap/",
                                                                                  OutputTypeSet[h], "/", VarTypeSet[v], sep="")



        if (any(ModelType == c("JLL original", "JLL NoDomUnit", "JLL jointSigma"))){
            OrthoLabel <- paste(VarTypeSet,"Ortho")
            PathsGraphs[[ModelType]][[OutputTypeSet[h]]][[OrthoLabel[v]]] <- paste(tempdir(),"/Outputs/", ModelType, "/Bootstrap/",
                                                                                   OutputTypeSet[h], "/", VarTypeSet[v],
                                                                                   "/Ortho", sep="")
          }
        }
        if (h==6){
          PathsGraphs[[ModelType]][[OutputTypeSet[h]]]  <- paste(tempdir(),"/Outputs/", ModelType, "/Bootstrap/",
                                                                 OutputTypeSet[h], sep="")

      }
    }
  }

}
  return(PathsGraphs)
}


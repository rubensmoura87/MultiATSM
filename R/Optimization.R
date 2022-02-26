#' Peform the minimization of mean(f)
#'
#' @param f vector-valued objective function (function)
#' @param tol convergence tolerance (scalar). For ML estimation, a reasonable value is tol <- 1e-4
#' @param varargin list containg starting values and constraints:
#'                        for each input argument K (of f), we need four inputs that look like:
#'    \enumerate{
#'        \item a starting value: K0
#'        \item a variable label ('K0') followed by a ':' followed by a type of constraint. The constraint can be:
#'                  \itemize{
#'                      \item 'bounded': bounded matrix;
#'                      \item 'Jordan' or 'Jordan MultiCountry': a matrix of Jordan type;
#'                      \item 'psd': psd matrix;
#'                      \item 'stationary': largest eigenvalue of the risk-neutral feedback matrix is strictly smaller than 1;
#'                      \item 'diag' or 'BlockDiag': a diagonal or block diagonal matrix.
#'                      \item 'JLLstructure': to impose the zero-restrictions on the variance-voriance matrix along
#'                              the lines of the JLL models
#'                          }
#'        \item a lower bound lb (lb <- NULL -> no lower bound)
#'        \item an upper bound ub (ub <- NULL -> no upper bound)
#'        \item Specification of the optimization settings:
#'        \itemize{
#'                  \item 'iter off': hide the printouts of the numerical optimization routines;
#'                  \item 'fminunc only': only uses fminunc for the optimization;
#'                  \item ''fminsearch only': only uses fminsearch for the optimization.
#'
#'   }
#'   }
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param JLLinputs inputs used in the estimation of the JLL-based models; Default is set to NULL
#' @param GVARinputs inputs used in the estimation of the GVAR-based models; Default is set to NULL
#'
#'@examples
#'#' # See examples in the vignette file of this package (Section 4).
#'
#'@return
#' (i) out: list of second output produced by f (the first output of f must be the objective value to be minimized). \cr
#' (ii) x:  list containing parameter estimates
#'
#' @details
#' If a variable name starts with a '@@', it means that that parameter will be analytically concentrated out in
#' the specification of f. In this case, no starting value is needed for this particular parameter (an empty matrix
#' can be provided as a starting value).
#'
#'
#'
#'@importFrom pracma fminunc
#'
#'@references
#' This function is based on the "LS__opt" function by Le and Singleton (2018). \cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'
#'@export

Optimization <- function(f, tol, varargin, FactorLabels, Economies, ModelType, JLLinputs = NULL, GVARinputs= NULL){
  Jmisc::tic()

  print('#########################################################################################################')
  print( paste('###########################', 'Optimization  (Point Estimate) -', ModelType, '#################################' ))
  print('#########################################################################################################')


  NP <- floor(sum(lengths(varargin))/4)
  if (sum(lengths(varargin))>NP*4){
    EstType <- varargin[NP+1]
  }else{
    EstType <- NULL
  }
  iter <- 'iter'

  if (grepl("iter off", EstType) ){ iter <- "off" }

  if (!exists("tol") || is.null(tol) ){ tol <- 1e-4}
  Max_AG_Iteration <- 1e4
  Previous_Optimal_Obj <-  -1e20

  AuxVal <- getx(con = "concentration", varargin, Economies, FactorLabels, JLLinputs) # Transform the initial guesses of K1XQ and SSZ in the auxiliary parameters
                                            # that will NOT be concentrated out of the llk function
                                            #(NOTE: in the case in which all the parameters are concentrated out
                                            # of the llk, the initial other parameters other than K1XQ and SSP will
                                            # not be changed by the function getx)

  # Value of the likelihood function:
  FFvec <- functional::Curry(f_with_vectorized_parameters, sizex= AuxVal$sizex , f, con = 'concentration', varargin, ModelType,
                 FactorLabels, Economies, JLLinputs, GVARinputs, nargout=1)
  FF <-function(x0){ mean(FFvec(x=x0))   }

  options200 <- neldermead::optimset(MaxFunEvals = 200*numel(AuxVal$x0), Display = iter,
                         MaxIter = 200, GradObj='off', TolFun= 10^-8, TolX= 10^-8)

  options1000 <- options200
  options1000$MaxIter <- 1000

  converged<-(tol>1e5)
  oldF <- FF(x=AuxVal$x0)

  scaling_vector <- NULL; count <- 0

  while (!converged){
    if (!contain('fminsearch only', EstType)){

      if (!contain('no rescaling', EstType)){
        if (isempty(scaling_vector)){
          dFFvec <- df__dx(f=FFvec, x=AuxVal$x0) # first order derivative of the llK function for each point in time for the initial guess of the parameters which are NOT concentrated out of the llk.

          vv <- 1/rowMeans(abs(dFFvec))
          vv[is.infinite(vv)] <- max(vv[!is.infinite(vv)])
          vv[vv==0] <- min(vv[vv>0])
          scaling_vector <- t(t(vv))
        }

        FFtemporary <- function(xtemp,scaling_vector){  FF(x=scaling_vector*xtemp)       }

        FFtemp <- functional::Curry(FFtemporary, scaling_vector = scaling_vector)

        x1 <- fminunc(x0=AuxVal$x0/scaling_vector, FFtemp , gr = NULL, tol = options200$TolFun,
                      maxiter = options200$MaxIter , maxfeval = options200$MaxFunEvals )

        x1 <- x1$par*scaling_vector

      } else{
        x1 <- fminunc(x0=AuxVal$x0, FF , gr = NULL, tol = options200$TolFun, maxiter = options200$MaxIter,
                      maxfeval = options200$MaxFunEvals)
        x1 <- x1$par
      }

      if (FF(x=x1)<FF(x=AuxVal$x0)){ AuxVal$x0 <- x1 }
    }
    if (!contain('fminunc only', EstType)){
      x1<-neldermead::fminsearch(FF, AuxVal$x0, options1000)$optbase$xopt
      if (FF(x=x1)<FF(x=AuxVal$x0)){ AuxVal$x0 <- x1 }
    }

    newF <- FF(x=AuxVal$x0)
    print(newF)
    converged <-   (abs(oldF - newF)<tol) ||(count>Max_AG_Iteration && newF>Previous_Optimal_Obj)
    oldF <- newF

    count <- count+1

  }

  # Build the full auxiliary vector, including concentrated parameters
  ud <- update_para(AuxVal$x0, sizex= AuxVal$sizex, ii= NULL, con= 'concentration',
                    FactorLabels, Economies, JLLinputs, GVARinputs, varargin) # update the parameter set which were NOT concentrated out after the optimization.
  FF <- functional::Curry(f_with_vectorized_parameters, sizex= AuxVal$sizex, f, con = 'concentration', varargin=ud,
              ModelType, FactorLabels, Economies, JLLinputs, GVARinputs, nargout=2)



  tryCatch({
    out <- FF(x=AuxVal$x0)$out # computes the estimates of all parameters after the optimization.

    for (i in 1:NP){ # for loop: remove the @ from the parameters' labels.
      if (contain('@', ud[[i]]$Label) ) {
        namexi <- killa(ud[[i]]$Label)
        ud[[i]]$Label <- namexi

        si <- strfind(namexi, ':')
        if (!isempty(si)){
          namexi <- substr(namexi, start = 1, stop = si-1 )
        }

        ud[[i]]$Value <- out$ests[[i]]
      }
    }
  }, error=function(e){}
  )

  # Produce outputs:
  x0 <- getx(con='', ud, Economies, FactorLabels, JLLinputs)$x0
  sizex <- getx(con='', ud, Economies, FactorLabels, JLLinputs)$sizex
  namex <- getx(con='', ud, Economies, FactorLabels, JLLinputs)$namex

  x.vec_ <- x0 # Vector of auxiliary parameters.

  FF <- functional::Curry(f_with_vectorized_parameters, sizex=sizex, f=f, con='', varargin=ud, ModelType,
              FactorLabels, Economies, JLLinputs, GVARinputs, nargout=2)
  out <- FF(x=x0)$out


  # Stack outputs in a list:
  x <- vector(mode = "list", length = length(namex) +1)
  x[[1]] <- x.vec_

  for (i in 1:dim(sizex)[1]){
    Fi <- functional::Curry(update_para, sizex=sizex, ii=i , con='', FactorLabels, Economies,
                            JLLinputs, GVARinputs, varargin= ud)
    if (!isempty(namex[i]) ){
      x[[i+1]] <- Fi(x=x.vec_)
    }


  }


  # Add labels to the variables
  names(x) <- c('vec_',namex)
  outputs <- list(x,out)

  names(outputs) <- c("FinalEstimates", "Summary")

  Jmisc::toc()
  return(outputs)
}


#############################################################################################################
#' Peform the minimization of mean(f) (adapted for the bootstrap setting)
#'
#' @param f vector-valued objective function (function)
#' @param tol convergence tolerance (scalar). For ML estimation, a reasonable value is tol <- 1e-4
#' @param varargin list containg starting values and constraints:
#'                        for each input argument K (of f), we need four inputs that look like:
#'    \enumerate{
#'        \item a starting value: K0
#'        \item a variable label ('K0') followed by a ':' followed by a type of constraint. The constraint can be:
#'                  \itemize{
#'                      \item 'bounded': bounded matrix;
#'                      \item 'Jordan' or 'Jordan MultiCountry': a matrix of Jordan type;
#'                      \item 'psd': psd matrix;
#'                      \item 'stationary': largest eigenvalue of the risk-neutral feedback matrix is strictly smaller than 1;
#'                      \item 'diag' or 'BlockDiag': a diagonal or block diagonal matrix.
#'                      \item 'JLLstructure': to impose the zero-restrictions on the variance-voriance matrix along
#'                              the lines of the JLL models
#'                          }
#'        \item a lower bound lb (lb <- NULL -> no lower bound)
#'        \item an upper bound ub (ub <- NULL -> no upper bound)
#'        \item Specification of the optimization settings:
#'        \itemize{
#'                  \item 'iter off': hide the printouts of the numerical optimization routines;
#'                  \item 'fminunc only': only uses fminunc for the optimization;
#'                  \item ''fminsearch only': only uses fminsearch for the optimization.
#'
#'   }
#'   }
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param JLLinputs inputs used in the estimation of the JLL-based models; Default is set to NULL
#' @param GVARinputs inputs used in the estimation of the GVAR-based models; Default is set to NULL


#'@return
#' (i) out: list of second output produced by f (the first output of f must be the objective value to be minimized)\cr
#' (ii) x:  list containing parameter estimates
#'
#' @details
#' If a variable name starts with a '@@', it means that that parameter will be analytically concentrated out in
#' the specification of f. In this case, no starting value is needed for this particular parameter.An  empty matrix
#' can be provided as a starting value


#'@importFrom pracma fminunc
#'
#'@references
#' This function is based on the "LS__opt" function by Le and Singleton (2018).\cr
#'  "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models." \cr
#'  (Euro Area Business Cycle Network Training School - Term Structure Modelling).
#'  Available at: https://cepr.org/40029
#'



Optimization_Boot <- function(f, tol, varargin, FactorLabels, Economies, ModelType, JLLinputs = NULL, GVARinputs = NULL){



  NP <- floor(sum(lengths(varargin))/4)
  if (sum(lengths(varargin))>NP*4){
    EstType <- varargin[NP+1]
  }else{
    EstType <- NULL
  }
  iter <- 'iter'

  if (grepl("iter off", EstType) ){ iter <- "off" }

  if (!exists("tol") || is.null(tol) ){ tol <- 1e-4}
  Max_AG_Iteration <- 1e4
  Previous_Optimal_Obj <-  -1e20

  AuxVal <- getx(con = "concentration", varargin, Economies, FactorLabels, JLLinputs) # Transform the initial guesses of K1XQ and SSZ
                                                                                    # in the auxiliary parameters that will
                                                                                    # NOT be concentrated out of the llk function


  # Value of the likelihood function:
  FFvec <- functional::Curry(f_with_vectorized_parameters, sizex= AuxVal$sizex , f, con = 'concentration', varargin,
                 ModelType, FactorLabels, Economies, JLLinputs, GVARinputs, nargout=1)
  FF <-function(x0){ mean(FFvec(x=x0))   }

  options200 <- neldermead::optimset(MaxFunEvals = 200*numel(AuxVal$x0), Display = iter,
                         MaxIter = 1000, GradObj='off', TolFun= 10^-8, TolX= 10^-8)

  options1000 <- options200
  options1000$MaxIter <- 1000

  converged<-(tol>1e5)
  oldF <- FF(x=AuxVal$x0)


  scaling_vector <- NULL; count <- 0

  while (!converged){
    if (!contain('fminsearch only', EstType)){

      if (!contain('no rescaling', EstType)){
        if (isempty(scaling_vector)){
          dFFvec <- df__dx(f=FFvec, x=AuxVal$x0) # first order derivative of the llK function for each point in time for the initial guess of the parameters which are NOT concentrated out of the llk.
          vv <- 1/rowMeans(abs(dFFvec))
          vv[is.infinite(vv)] <- max(vv[!is.infinite(vv)])
          vv[vv==0] <- min(vv[vv>0])
          scaling_vector <- t(t(vv))
        }

        FFtemporary <- function(xtemp,scaling_vector){  FF(x=scaling_vector*xtemp)       }

        FFtemp <- functional::Curry(FFtemporary, scaling_vector = scaling_vector)

        x1 <- fminunc(x0=AuxVal$x0/scaling_vector, FFtemp , gr = NULL, tol = options200$TolFun,
                      maxiter = options200$MaxIter , maxfeval = options200$MaxFunEvals )

        x1 <- x1$par*scaling_vector

      } else{
        x1 <- fminunc(x0=AuxVal$x0, FF , gr = NULL, tol = options200$TolFun, maxiter = options200$MaxIter,
                      maxfeval = options200$MaxFunEvals)
        x1 <- x1$par
      }

      if (FF(x=x1)<FF(x=AuxVal$x0)){ AuxVal$x0 <- x1 }
    }
    if (!contain('fminunc only', EstType)){

      x1<-neldermead::fminsearch(FF, AuxVal$x0, options1000)$optbase$xopt
      if (FF(x=x1)<FF(x=AuxVal$x0)){ AuxVal$x0 <- x1 }
    }

    newF <- FF(x=AuxVal$x0)
    # disp(newF)
    converged <-   (abs(oldF - newF)<tol) ||(count>Max_AG_Iteration && newF>Previous_Optimal_Obj)
    oldF <- newF

    count <- count+1

  }

  # now build the full auxiliary vector, including concentrated parameters
  ud <- update_para(AuxVal$x0, sizex= AuxVal$sizex, ii= NULL, con= 'concentration',
                    FactorLabels, Economies, JLLinputs, GVARinputs, varargin) # update the parameter set which were NOT concentrated out after the optimization.
  FF <- functional::Curry(f_with_vectorized_parameters, sizex= AuxVal$sizex, f, con = 'concentration', varargin=ud,
              ModelType, FactorLabels, Economies, JLLinputs, GVARinputs, nargout=2)


  tryCatch({
    out <- FF(x=AuxVal$x0)$out # computes all the estimates of all parameters after the optimization.

    for (i in 1:NP){ # for loop: remove the @ from the parameters' labels.
      if (contain('@', ud[[i]]$Label) ) {
        namexi <- killa(ud[[i]]$Label)
        ud[[i]]$Label <- namexi

        si <- strfind(namexi, ':')
        if (!isempty(si)){
          namexi <- substr(namexi, start = 1, stop = si-1 )
        }

        ud[[i]]$Value <- out$ests[[i]]
      }
    }
  }, error=function(e){}
  )

  ## Produce outputs:
  x0 <- getx(con='', ud, Economies, FactorLabels, JLLinputs)$x0
  sizex <- getx(con='', ud, Economies, FactorLabels, JLLinputs)$sizex
  namex <- getx(con='', ud, Economies, FactorLabels, JLLinputs)$namex

  x.vec_ <- x0 # Vector of auxiliary parameters.

  FF <- functional::Curry(f_with_vectorized_parameters, sizex=sizex, f=f, con='', varargin=ud,
              ModelType, FactorLabels, Economies, JLLinputs, GVARinputs, nargout=2)
  out <- FF(x=x0)$out


  return(out)

}

##########################################################################################################################
#' Print method for ATSMModelInputs objects
#' @param x An object of class 'ATSMModelInputs'
#' @param ... Additional arguments (not used)
#'
#' @method print ATSMModelInputs
#' @usage \method{print}{ATSMModelInputs}(x, ...)
#'
#' @export

print.ATSMModelInputs <- function(x, ...) {

  info <- attr(x, "ModelInfo")

  cat("General Model Inputs: \n")
  cat("------------------------------------------\n")
  cat("Model Type:", info$ModelType, "\n")
  cat("Economic System:", paste(info$Economies, collapse = ", "), "\n")
  cat("Sample Span:", info$InitialSampleDate, "-", info$FinalSampleDate, "\n")
  cat("Data Frequency:", info$DataFrequency, "\n")
  cat("Common Maturities across Countries:", info$Maturities, "years \n")
  cat("------------------------------------------\n\n")

  cat("Key Structure of the Economic System: \n")
  cat("------------------------------------------\n")
  cat("Total amount of spanned factors in the system:", info$NC, "\n")
  cat("Total amount of global unspanned factors in the system:", info$G, "\n")
  cat("Total amount of country-specific unspanned factors in the system:", info$MC, "\n")
  cat("Total amount of risk factors in the system:", info$NC + info$G + info$MC, "\n")
  cat("------------------------------------------\n")

  invisible(x)
}
##########################################################################################################################
#' Summary method for ATSMModelInputs objects
#' @param object An object of class 'ATSMModelInputs'
#' @param ... Additional arguments (not used)
#'
#' @method summary ATSMModelInputs
#' @usage \method{summary}{ATSMModelInputs}(object, ...)
#' @export

summary.ATSMModelInputs <- function(object, ...) {

  info <- attr(object, "ModelInfo")

  # Function to print the summary statistics
  summary_TS <- function(mat) {
    data.frame(
      Mean = round(rowMeans(mat, na.rm = TRUE), 3),
      Std_Dev = round(apply(mat, 1, stats::sd, na.rm = TRUE), 3),
      Min = round(apply(mat, 1, min, na.rm = TRUE), 3),
      Max = round(apply(mat, 1, max, na.rm = TRUE), 3),
      row.names = rownames(mat)
    )
  }

  cat("Summary Statistics from the Time Series Components:\n")
  cat("------------------------------------------------------\n")

  cat("1) Risk Factors:\n")
  # Single-country estimation
  if(info$ModelType %in% c("JPS original", "JPS global", "GVAR single")){
    for(i in 1:length(info$Economies)){
      cat("\n Model", info$Economies[i] , "\n")
      summary_result <- summary_TS(info$RiskFactors[[i]])
      print(summary_result)
    }

    # Multicountry estimation
  }else{
    summary_result <- summary_TS(info$RiskFactors)
    print(summary_result)
  }


  cat("\n 2) Bond Yields:\n")
  # Single-country estimation
  if(info$ModelType %in% c("JPS original", "JPS global", "GVAR single")){

    J <- length(info$Maturities)

    for(i in 1:length(info$Economies)){
      cat("\n Model", info$Economies[i] , "\n")

      # Extract rows for the current country
      country_rows <- ((i - 1) * J + 1):(i * J)
      country_yields <- info$Yields[country_rows, ] * 100

      summary_result <- summary_TS(country_yields)
      print(summary_result)
    }

    # Multicountry estimation
  }else{
    summary_result <- summary_TS((info$Yields)*100)
    print(summary_result)
  }
  cat("------------------------------------------------------\n")

  invisible(object)
}

##########################################################################################################################
#' Summary method for ATSMModelOutputs objects
#' @param object An object of class 'ATSMModelOutputs'
#' @param ... Additional arguments (not used)
#'
#' @method summary ATSMModelOutputs
#' @usage \method{summary}{ATSMModelOutputs}(object, ...)
#' @export

summary.ATSMModelOutputs <- function(object, ...) {

  info <- attr(object, "ModelOutInfo")

  # Risk-neutral eigenvalues
  Get_Q_Eigein <- function(K1XQ, Economies){
    All_Eigen <- diag(K1XQ)
    C <- length(Economies)
    N <- length(All_Eigen)/C

    Eigen_CS <- matrix(NA, nrow = C, ncol = N)
    for (tt in 1:C) {
      start <- (tt - 1) * N + 1
      end <- tt * N
      Eigen_CS[tt , ] <- All_Eigen[start:end]
    }

    rownames(Eigen_CS) <- Economies
    colnames(Eigen_CS) <-paste0("Eigen_", 1:N)
    return(Eigen_CS)
    }

  # r0 function
  Get_r0_All <- function(r0, Economies){
  r0_All <- matrix(r0, ncol = length(Economies), dimnames = list(NULL, Economies))
  row.names(r0_All) <- "r0"
  return(r0_All)
  }

  # Maximum eigenvalue
 Max_Eigen <- function(K1Z){
  Eigen_All   <- abs(eigen(K1Z)$values)
  LargeEigen <- sum(Eigen_All > 1)

  return(LargeEigen)
  }

  # Pre-allocation
  Economies <- info$Economies

  # 1) Single country models
  if(info$ModelType %in% c("JPS original", "JPS global", "GVAR single")){
    for(i in 1:length(info$Economies)){

      K1XQ <- info$Outs[[Economies[i]]]$ests$K1XQ
      r0 <- info$Outs[[Economies[i]]]$ests$r0
      K1Z <- info$Outs[[Economies[i]]]$ests$K1Z

      cat("Model", info$Economies[i] , "\n")
      cat("General Model Outputs: \n")
      cat("------------------------------------------------------\n")
      cat("1) Maximum Log-likelihood value:", mean(info$Outs[[Economies[i]]]$llk[[1]]), "\n")

      cat("\n 2) Key Features from the Q-measure: \n")
      cat(" a) Long run risk-neutral mean (r0):", round(r0, 5), "\n")
      cat(" b) Risk-neutral eigenvalues:", diag(round(K1XQ, 5)), "\n" )
      cat("\n c) Intercept matrix of loadings (Matrix A): \n")
      print(round(info$Outs[[Economies[i]]]$rot$P$A, 6))
      cat("\n d) Slope Matrix of Loadings (Matrix B) : \n")
      print(round(info$Outs[[Economies[i]]]$rot$P$B, 6))

      cat("\n 3) Key Features from the P-measure: \n")
      cat(" - Number of physical-measure eigenvalues:", Max_Eigen(K1Z), "\n")

      if (Max_Eigen(K1Z) > 1){
        cat(" - Some (Generalized) Impulse Response Functions may exhibit explosive behaviour. \n")
      }
      cat("------------------------------------------------------\n")

    }
    } else {

  # 2) Multicountry models
  K1XQ <- info$Outs$ests$K1XQ
  r0 <- info$Outs$ests$r0
  K1Z <- info$Outs$ests$K1Z

  cat("General Model Outputs: \n")
  cat("------------------------------------------------------\n")
  cat("1) Maximum Log-likelihood value:", mean(info$Outs$llk[[1]]), "\n")

  cat("2) Key Features from the Q-measure: \n")
  cat(" \n a) Long run risk-neutral mean (r0): \n")
  print(round( Get_r0_All(r0, Economies), 5))
  cat(" \n b) Risk-neutral eigenvalues: \n" )
  print(Get_Q_Eigein(K1XQ, Economies))
  cat(" \n c) Intercept matrix of loadings (Matrix A): \n")
  print(round(info$Outs$rot$P$A, 6))
  cat(" \n d) Slope Matrix of Loadings (Matrix B) : \n")
  print(round(info$Outs$rot$P$B, 6))

  cat("\n 3) Key Features from the P-measure: \n")
  cat(" - Number of physical-measure eigenvalues:", Max_Eigen(K1Z), "\n")

  if (Max_Eigen(K1Z) > 1){
    cat(" - Some (Generalized) Impulse Response Functions may exhibit explosive behaviour. \n")
  }
  cat("------------------------------------------------------\n")
  }

  invisible(object)
}
#######################################################################################################
#' Plot method for ATSMModelForecast objects
#'
#' @param x An object of class \code{ATSMModelForecast}
#' @param ... Additional arguments (not used)
#'
#' @method plot ATSMModelForecast
#' @usage \method{plot}{ATSMModelForecast}(x, ...)
#'
#' @export


  plot.ATSMModelForecast <- function(x, ...) {

  info <- attr(x, "ModelForecast")

# Preliminary work
SeQ_Lab <- c("JPS original", "JPS global", "GVAR single")
ModelType <- info$ModelType
Economies <- info$Economies

C <- length(Economies)
H <- info$ForHoriz
J <- if (ModelType %in% SeQ_Lab) {  nrow(x$RMSE[[Economies[1]]]) } else { nrow(x$RMSE[[ModelType]]) / C }

# Set up the plot layout to have two columns
graphics::par(mfrow = c(ceiling(J / 2), 2), mar = c(4, 4, 2, 1) + 0.1)


# Format the y-axis labels to show only two decimal places
format_y_axis <- function(x) {
  sprintf("%.2f", x)
}

# Plot the bar charts
for (i in 1:C){

  if( ModelType %in% SeQ_Lab){
  DataGraph_CS <- x$RMSE[[Economies[i]]]
  } else {
    Idx0 <- 1 + (i-1)*J
    IdxF <- Idx0 + J - 1
    DataGraph_CS <- x$RMSE[[ModelType]][Idx0:IdxF, ]
  }

  MatLab <- row.names(DataGraph_CS)

  for (j in 1:J) {
  DataGraph <- 100*DataGraph_CS[j, ]

  ylim_range <- range(0, max(DataGraph))
  graphics::barplot(DataGraph, xlab = "", ylab = "RMSE (%)", col = "lightblue", names.arg = 1:H, cex.names = 0.7, yaxt = "n", ylim = ylim_range)
  graphics::axis(2, at = pretty(ylim_range), labels = format_y_axis(pretty(ylim_range)))
  graphics::mtext(MatLab[j], side = 1, line = 2, cex = 0.8)
  }
  }
}

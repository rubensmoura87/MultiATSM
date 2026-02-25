# Perform the optimization of the log-likelihood function of the chosen ATSM

Perform the optimization of the log-likelihood function of the chosen
ATSM

## Usage

``` r
Optimization(
  MLEinputs,
  StatQ,
  DataFreq,
  FactorLabels,
  Economies,
  ModelType,
  tol = 1e-04,
  EstType = c("BFGS", "Nelder-Mead"),
  TimeCount = TRUE,
  BS_outputs = FALSE,
  verbose = TRUE
)
```

## Arguments

- MLEinputs:

  list. Contains the inputs for building the log-likelihood function
  (see
  [`InputsForOpt`](https://rubensmoura87.github.io/MultiATSM/reference/InputsForOpt.md)).

- StatQ:

  A logical value indicating whether to impose that the largest
  eigenvalue under Q is strictly smaller than 1. Set TRUE to impose this
  restriction.

- DataFreq:

  character vector specifying the data frequency. Available options:
  "Daily All Days", "Daily Business Days", "Weekly", "Monthly",
  "Quarterly", "Annually".

- FactorLabels:

  list. Labels for all variables present in the model, as returned by
  [`LabFac`](https://rubensmoura87.github.io/MultiATSM/reference/LabFac.md).

- Economies:

  character vector. Names of the `C` economies included in the system.

- ModelType:

  character. Model type to be estimated. Permissible choices: "JPS
  original", "JPS global", "GVAR single", "JPS multi", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma".

- tol:

  numeric. Convergence tolerance. The default is 1e-4.

- EstType:

  Available options are"BFGS" and/or "Nelder-Mead".

- TimeCount:

  Logical. If TRUE, computes the time required for model estimation.
  Default is TRUE.

- BS_outputs:

  Logical. If TRUE, generates a simplified output list in the bootstrap
  setting. Default is FALSE.

- verbose:

  Logical flag controlling function messaging. Default is TRUE.

## Value

An object of class 'ATSMModelOutputs' containing model outputs after the
optimization of the chosen ATSM specification.

## Available Methods

\- \`summary(object)\`

## References

- Candelon, C. and Moura, R. (2024). “A Multicountry Model of the Term
  Structures of Interest Rates with a GVAR.” Journal of Financial
  Econometrics 22 (5): 1558–87.

- Jotikasthira, C; Le, A. and Lundblad, C (2015). “Why Do Term
  Structures in Different Currencies Co-Move?” Journal of Financial
  Economics 115: 58–83.

- Joslin, S,; Priebsch, M. and Singleton, K. (2014). “Risk Premiums in
  Dynamic Term Structure Models with Unspanned Macro Risks.” Journal of
  Finance 69 (3): 1197–1233.

- Joslin, S., Singleton, K. and Zhu, H. (2011). "A new perspective on
  Gaussian dynamic term structure models". The Review of Financial
  Studies.

- Le, A. and Singleton, K. (2018). "A Small Package of Matlab Routines
  for the Estimation of Some Term Structure Models." Euro Area Business
  Cycle Network Training School - Term Structure Modelling.

## Examples

``` r
LoadData("CM_2024")
ModelType <- "JPS original"
Economy <- "Brazil"
t0 <- "01-05-2007" # Initial Sample Date (Format: "dd-mm-yyyy")
tF <- "01-12-2018" # Final Sample Date (Format: "dd-mm-yyyy")
N <- 1
GlobalVar <- "Gl_Eco_Act" # Global Variables
DomVar <- "Eco_Act" # Domestic Variables
DataFreq <- "Monthly"
StatQ <- FALSE

FacLab <- LabFac(N, DomVar, GlobalVar, Economy, ModelType)
ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacro, DomMacro,
  FacLab, Economy, DataFreq,
  CheckInputs = FALSE, verbose = FALSE
)

OptPara <- Optimization(ATSMInputs, StatQ, DataFreq, FacLab, Economy, ModelType, verbose = FALSE)
```

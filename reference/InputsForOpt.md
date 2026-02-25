# Generates inputs necessary to build the likelihood function for the ATSM model

Generates inputs necessary to build the likelihood function for the ATSM
model

## Usage

``` r
InputsForOpt(
  InitialSampleDate,
  FinalSampleDate,
  ModelType,
  Yields,
  GlobalMacro,
  DomMacro,
  FactorLabels,
  Economies,
  DataFrequency,
  GVARlist = NULL,
  JLLlist = NULL,
  WishBRW = FALSE,
  BRWlist = NULL,
  UnitYields = "Month",
  CheckInputs = TRUE,
  BS_Adj = FALSE,
  verbose = TRUE
)
```

## Arguments

- InitialSampleDate:

  Start date of the sample period in the format "dd-mm-yyyy"

- FinalSampleDate:

  End date of the sample period in the format "dd-mm-yyyy"

- ModelType:

  character. Model type to be estimated. Permissible choices: "JPS
  original", "JPS global", "GVAR single", "JPS multi", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma".

- Yields:

  numerical matrix with time series of yields (`J x Td` or `CJ x Td`)

- GlobalMacro:

  numerical matrix with time series of the global risk factors
  (`G x Td`)

- DomMacro:

  numerical matrix with time series of the country-specific risk factors
  for all `C` countries ( `C X Td` or `CM x Td`)

- FactorLabels:

  list. Labels for all variables present in the model, as returned by
  [`LabFac`](https://rubensmoura87.github.io/MultiATSM/reference/LabFac.md).

- Economies:

  character vector. Names of the `C` economies included in the system.

- DataFrequency:

  character. Data frequency. Permissible choices: "Daily All Days",
  "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually".

- GVARlist:

  list. Inputs for GVAR model estimation. See details below.

- JLLlist:

  list. Inputs for JLL model estimation. See details below.

- WishBRW:

  logical. Whether to estimate the physical parameter model with bias
  correction (see
  [`Bias_Correc_VAR`](https://rubensmoura87.github.io/MultiATSM/reference/Bias_Correc_VAR.md)).
  Default is FALSE.

- BRWlist:

  list. Inputs for bias-corrected estimation.

- UnitYields:

  character. Maturity unit of yields. Permissible choices: "Month" or
  "Year". Default is "Month".

- CheckInputs:

  logical. Whether to perform a prior check on the consistency of the
  provided input list. Default is TRUE.

- BS_Adj:

  logical. Whether to adjust the global series for the sepQ models in
  the Bootstrap setting. Default is FALSE.

- verbose:

  logical. Print progress messages. Default is TRUE.

## Value

An object of class 'ATSMModelInputs' containing the necessary inputs for
performing the model optimization.

## Permissible options for GVARlist

- `VARXtype`: "unconstrained" or "constrained"

- `W_type`: "Time-varying" or "Sample Mean"

- `t_First_Wgvar`, `t_Last_Wgvar`: year as character

## Permissible options for JLLlist

- `DomUnit`: name of the dominant economy or `None`

- `WishSigmas`: TRUE (estimate variance-covariance matrices) or FALSE

- `SigmaNonOrtho`: NULL or `K x K` matrix

## Permissible options for BRWlist

- `BiasCorrection`: TRUE (bias-corrected) or FALSE

- `flag_mean`: TRUE (mean) or FALSE (median)

- `gamma`: numeric adjustment parameter

- `N_iter`: number of iterations

- `N_burn`: number of burn-in iterations

- `B`: number of bootstrap samples

- `checkBRW`: TRUE or FALSE

- `B_check`: number of bootstrap samples for closeness check

## General Notation

- `Td` model time series dimension.

- `C` number of countries in the system.

- `G` number of global unspanned factors.

- `M` number of country-specific unspanned factors.

- `K` total number of risk factors.

- `J` number of bond yields per country used in estimation.

## Available Methods

\- \`print(object)\` - \`summary(object)\`

## Examples

``` r
# \donttest{
# Example 1:
data(GlobalMacro)
data(DomMacro)
data(Yields)

ModelType <- "JPS original"
Economies <- "Mexico"
t0 <- "01-05-2007" # Initial Sample Date (Format: "dd-mm-yyyy")
tF <- "01-12-2018" # Final Sample Date (Format: "dd-mm-yyyy")
N <- 3
GlobalVar <- c("Gl_Eco_Act") # Global Variables
DomVar <- c("Eco_Act") # Domestic Variables
FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)

DataFreq <- "Monthly"

ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacro, DomMacro,
  FactorLabels, Economies, DataFreq,
  CheckInputs = FALSE, verbose = FALSE
)

# Example 2:
LoadData("CM_2024")

ModelType <- "GVAR multi"

Economies <- c("China", "Brazil", "Mexico", "Uruguay")
t0 <- "01-05-2007" # InitialSampleDate (Format: "dd-mm-yyyy")
tF <- "01-12-2019" # FinalSampleDate (Format: "dd-mm-yyyy")
N <- 2
GlobalVar <- c("Gl_Eco_Act", "Gl_Inflation") # Global Variables
DomVar <- c("Inflation") # Domestic Variables
FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)

DataFreq <- "Monthly"
GVARlist <- list(
  VARXtype = "unconstrained", W_type = "Sample Mean",
  t_First_Wgvar = "2007", t_Last_Wgvar = "2019", DataConnectedness = TradeFlows
)

ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacro, DomMacro,
  FactorLabels, Economies, DataFreq, GVARlist,
  CheckInputs = FALSE, verbose = FALSE
)

# Example 3:
LoadData("CM_2024")

ModelType <- "JLL original"

Economies <- c("China", "Brazil", "Uruguay")
t0 <- "01-05-2007" # InitialSampleDate (Format: "dd-mm-yyyy")
tF <- "01-12-2019" # FinalSampleDate (Format: "dd-mm-yyyy")
N <- 2
GlobalVar <- c("Gl_Eco_Act", "Gl_Inflation") # Global Variables
DomVar <- c("Eco_Act", "Inflation") # Domestic Variables
FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)

JLLinputs <- list(DomUnit = "China")

DataFrequency <- "Monthly"

ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacro, DomMacro,
  FactorLabels, Economies, DataFreq,
  JLLlist = JLLinputs,
  CheckInputs = FALSE, verbose = FALSE
)
# }
```

# Collects the inputs that are used to construct the numerical and graphical outputs

Collects the inputs that are used to construct the numerical and
graphical outputs

## Usage

``` r
InputsForOutputs(
  ModelType,
  Horiz,
  ListOutputWished,
  OutputLabel,
  WishStationarityQ,
  DataFrequency,
  WishGraphYields = FALSE,
  WishGraphRiskFactors = FALSE,
  WishOrthoJLLgraphs = FALSE,
  WishForwardPremia = FALSE,
  LimFP = NULL,
  WishBootstrap = FALSE,
  ListBoot = NULL,
  WishForecast = FALSE,
  ListForecast = NULL,
  UnitYields = "Month"
)
```

## Arguments

- ModelType:

  character. Model type to be estimated. Permissible choices: "JPS
  original", "JPS global", "GVAR single", "JPS multi", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma".

- Horiz:

  numeric scalar. Desired analysis horizon for the outputs.

- ListOutputWished:

  character vector. Desired graphical outputs. Available options:
  "RiskFactors", "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremia",
  "ForwardPremia".

- OutputLabel:

  character. Name of the output label to be stored.

- WishStationarityQ:

  logical. Whether to impose that the largest eigenvalue under Q is
  strictly smaller than 1. TRUE to impose.

- DataFrequency:

  character. Data frequency. Permissible choices: "Daily All Days",
  "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually".

- WishGraphYields:

  logical. Whether to generate graphs for yields. Default is FALSE.

- WishGraphRiskFactors:

  logical. Whether to generate graphs for risk factors. Default is
  FALSE.

- WishOrthoJLLgraphs:

  logical. Whether to generate orthogonalized JLL-based graphs. Default
  is FALSE.

- WishForwardPremia:

  logical. Whether to generate forward premia graphs. Default is FALSE.

- LimFP:

  numeric vector. Maturities associated with the start and end dates of
  the loan.

- WishBootstrap:

  logical. Whether to perform bootstrap-based estimation. Default is
  FALSE.

- ListBoot:

  list. Contains bootstrap settings: methodBS ("bs", "wild", "block"),
  BlockLength (numeric), ndraws (numeric), pctg (numeric).

- WishForecast:

  logical. Whether to generate forecasts. Default is FALSE.

- ListForecast:

  list. Contains forecast settings: ForHoriz (numeric), t0Sample
  (numeric), t0Forecast (numeric), ForType ("Rolling", "Expanding").

- UnitYields:

  character. Maturity unit of yields. Options: "Month" or "Year".
  Default is "Month".

## Value

List of necessary inputs to generate the graphs and outputs of the
desired model.

## Examples

``` r
ModelType <- "JPS original"
Horiz <- 100
DesiredOutputGraphs <- c("Fit", "GIRF", "GFEVD")
OutputLabel <- "Test"
WishStationarityQ <- TRUE
WishGraphRiskFac <- FALSE
WishGraphYields <- TRUE

InputsList <- InputsForOutputs(
  ModelType, Horiz, DesiredOutputGraphs, OutputLabel,
  WishStationarityQ, WishGraphYields, WishGraphRiskFac
)
```

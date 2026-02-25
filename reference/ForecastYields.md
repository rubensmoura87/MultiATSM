# Generates forecasts of bond yields for all model types

Generates forecasts of bond yields for all model types

## Usage

``` r
ForecastYields(
  ModelType,
  ModelPara,
  InputsForOutputs,
  FactorLabels,
  Economies,
  JLLlist = NULL,
  GVARlist = NULL,
  WishBRW = FALSE,
  BRWlist = NULL,
  Folder2save = NULL,
  verbose = TRUE
)
```

## Arguments

- ModelType:

  character. Model type to be estimated. Permissible choices: "JPS
  original", "JPS global", "GVAR single", "JPS multi", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma".

- ModelPara:

  list. Point estimates of the model parameters. See outputs from
  [`Optimization`](https://rubensmoura87.github.io/MultiATSM/reference/Optimization.md).

- InputsForOutputs:

  list. Inputs for generating IRFs, GIRFs, FEVDs, GFEVDs, and Term
  Premia.

- FactorLabels:

  list. Labels for all variables present in the model, as returned by
  [`LabFac`](https://rubensmoura87.github.io/MultiATSM/reference/LabFac.md).

- Economies:

  character vector. Names of the `C` economies included in the system.

- JLLlist:

  list. Inputs for JLL model estimation (see
  [`JLL`](https://rubensmoura87.github.io/MultiATSM/reference/JLL.md)).
  Default is NULL.

- GVARlist:

  list. Inputs for GVAR model estimation (see
  [`GVAR`](https://rubensmoura87.github.io/MultiATSM/reference/GVAR.md)).
  Default is NULL.

- WishBRW:

  logical. Whether to estimate the physical parameter model with bias
  correction (see
  [`Bias_Correc_VAR`](https://rubensmoura87.github.io/MultiATSM/reference/Bias_Correc_VAR.md)).
  Default is FALSE.

- BRWlist:

  list. Inputs for bias-corrected estimation (see
  [`Bias_Correc_VAR`](https://rubensmoura87.github.io/MultiATSM/reference/Bias_Correc_VAR.md)).

- Folder2save:

  character. Folder path where outputs will be stored. Default saves
  outputs in a temporary directory.

- verbose:

  logical. Print progress messages. Default is TRUE.

## Value

An object of class 'ATSMModelForecast' containing the following
elements:

1.  Out-of-sample forecasts of bond yields per forecast horizon

2.  Out-of-sample forecast errors of bond yields per forecast horizon

3.  Root mean square errors per forecast horizon

## Permissible options - forecast list (`InputsForOutputs` input)

- `ForHoriz`: forecast horizon. Must be a positive integer.

- `t0Sample`: initial sample date. Must be a positive integer smaller
  than the time series dimension of the model (`Td`)

- `t0Forecast`: last sample date for the first forecast. Note that
  `Td > t0Forecast + ForHoriz`.

- `ForType`: `"Rolling"` (rolling window forecast) or `"Expanding"` (for
  expanding window forecast)

## Available Methods

\- \`plot(object)\`

## Examples

``` r
# \donttest{
data("ParaSetEx")
data("InpForOutEx")
# Adjust inputs according to the loaded features
ModelType <- "JPS original"
Economy <- "Brazil"
FacLab <- LabFac(N = 1, DomVar = "Eco_Act", GlobalVar = "Gl_Eco_Act", Economy, ModelType)
 # Adjust Forecasting setting
InpForOutEx[[ModelType]]$Forecasting <- list(
  WishForecast = 1, ForHoriz = 12, t0Sample = 1,
  t0Forecast = 143, ForType = "Expanding"
)

Forecast <- ForecastYields(ModelType, ParaSetEx, InpForOutEx, FacLab, Economy,
  WishBRW = FALSE, verbose = TRUE
)
#> 4) OUT-OF-SAMPLE FORECASTING ANALYSIS
#> Out-of-sample forecast for the information set: 01-07-2006 || 01-05-2018 
#> Out-of-sample forecast for the information set: 01-07-2006 || 01-06-2018 
#> Out-of-sample forecast for the information set: 01-07-2006 || 01-07-2018 
#> Elapsed time: 3.34 seconds
# }
```

# Generates the bootstrap-related outputs

Generates the bootstrap-related outputs

## Usage

``` r
Bootstrap(
  ModelType,
  ModelParaPE,
  NumOutPE,
  Economies,
  InputsForOutputs,
  FactorLabels,
  JLLlist,
  GVARlist,
  WishBC = FALSE,
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

- ModelParaPE:

  list. Point estimates of the model parameters. See outputs from
  [`Optimization`](https://rubensmoura87.github.io/MultiATSM/reference/Optimization.md).

- NumOutPE:

  list. Point estimates from numerical outputs. See outputs from
  [`NumOutputs`](https://rubensmoura87.github.io/MultiATSM/reference/NumOutputs.md).

- Economies:

  character vector. Names of the `C` economies included in the system.

- InputsForOutputs:

  list. Inputs for generating IRFs, GIRFs, FEVDs, GFEVDs, and Term
  Premia.

- FactorLabels:

  list. Labels for all variables present in the model, as returned by
  [`LabFac`](https://rubensmoura87.github.io/MultiATSM/reference/LabFac.md).

- JLLlist:

  list. Inputs for JLL model estimation (see
  [`JLL`](https://rubensmoura87.github.io/MultiATSM/reference/JLL.md)).
  Default is NULL.

- GVARlist:

  list. Inputs for GVAR model estimation (see
  [`GVAR`](https://rubensmoura87.github.io/MultiATSM/reference/GVAR.md)).
  Default is NULL.

- WishBC:

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

An object of class 'ATSMModelBoot' containing:

- List of model parameters for each draw

- List of numerical outputs (IRFs, GIRFs, FEVDs, GFEVDs) for each draw

- Confidence bounds for the chosen level of significance

## Permissible options - Bootstrap list in `InputsForOutputs`

- `methodBS` : `"bs"` (standard bootstrap), `"wild"` (wild bootstrap),
  `"block"` (block bootstrap)

- `BlockLength` : required input for the block bootstrap method. Block
  length must be larger than 0 and smallar than the model time series
  dimension (`Td`).

- `ndraws`: number of draws. Must be a positive integer.

- `pctg` : confidence level. Must be a positive integer. Common choices
  are: 68, 90 and 95.

## Available methods

\- `autoplot(object, NumOutPE, type)`

## Examples

``` r
# \donttest{
data("ParaSetEx")
data("InpForOutEx")
data("NumOutEx")
ModelType <- "JPS original"
Economy <- "Brazil"
FacLab <- LabFac(N = 1, DomVar = "Eco_Act", GlobalVar = "Gl_Eco_Act", Economy, ModelType)

# Adjust Forecasting setting
InpForOutEx[[ModelType]]$Bootstrap <- list(
  WishBootstrap = 1, methodBS = "bs", BlockLength = 4,
  ndraws = 5, pctg = 95
)

Boot <- Bootstrap(ModelType, ParaSetEx, NumOutEx, Economy, InpForOutEx, FacLab,
  JLLlist = NULL,
  GVARlist = NULL, WishBC = FALSE, BRWlist = NULL, Folder2save = NULL, verbose = FALSE
)
#> Warning: cannot create dir '/tmp/RtmpPWG5lH/Outputs/JPS original/Bootstrap', reason 'No such file or directory'
#> Warning: cannot create dir '/tmp/RtmpPWG5lH/Outputs/JPS original/Bootstrap/Model Brazil', reason 'No such file or directory'
# }
```

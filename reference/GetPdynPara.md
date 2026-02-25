# Compute the parameters used in the P-dynamics of the model

Compute the parameters used in the P-dynamics of the model

## Usage

``` r
GetPdynPara(
  RiskFactors,
  FactorLabels,
  Economies,
  ModelType,
  BRWinputs,
  GVARinputs,
  JLLinputs,
  CheckInputs = F,
  verbose
)
```

## Arguments

- RiskFactors:

  time series of risk factors (F x T). Could be stored in a list
  depending on the model

- FactorLabels:

  string-list based which contains the labels of all variables present
  in the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- ModelType:

  string-vector containing the label of the model to be estimated

- BRWinputs:

  list of necessary inputs for performing the bias-corrected estimation
  (see "Bias_Correc_VAR" function)

- GVARinputs:

  list of necessary inputs for the estimation of GVAR-based models (see
  "GVAR" function)

- JLLinputs:

  list of necessary inputs for the estimation of JLL-based models (see
  "JLL" function)

- CheckInputs:

  Logical. Whether to perform a prior check on the consistency of the
  provided input list. Default is FALSE.

- verbose:

  Logical flag controlling function messaging.

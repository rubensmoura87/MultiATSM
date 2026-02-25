# Compute P-dynamics parameters using the bias correction method from BRW (2012)

Compute P-dynamics parameters using the bias correction method from BRW
(2012)

## Usage

``` r
GetPdynPara_BC(
  ModelType,
  BRWinputs,
  RiskFactors,
  Economies,
  FactorLabels,
  GVARinputs,
  JLLinputs,
  verbose
)
```

## Arguments

- ModelType:

  string-vector containing the label of the model to be estimated

- BRWinputs:

  list of necessary inputs for performing the bias-corrected estimation
  (see "Bias_Correc_VAR" function)

- RiskFactors:

  time series of risk factors (F x T). Could be stored in a list
  depending on the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- FactorLabels:

  string-list based which contains the labels of all variables present
  in the model

- GVARinputs:

  list of necessary inputs for the estimation of GVAR-based models (see
  "GVAR" function)

- JLLinputs:

  list of necessary inputs for the estimation of JLL-based models (see
  "JLL" function)

- verbose:

  Logical flag controlling function messaging.

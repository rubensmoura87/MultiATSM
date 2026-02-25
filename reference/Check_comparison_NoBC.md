# check how close the mean or median of the bias-corrected aproach is from the non-corrected approach

check how close the mean or median of the bias-corrected aproach is from
the non-corrected approach

## Usage

``` r
Check_comparison_NoBC(
  K1Z_BC,
  K1Z_NoBC,
  B_check,
  RiskFactors,
  GVARinputs,
  JLLinputs,
  FactorLabels,
  Economies,
  ModelType,
  Use_Mean,
  verbose
)
```

## Arguments

- K1Z_BC:

  feedback matrix after the bias correction procedure

- K1Z_NoBC:

  feedback matrix before the bias correction procedure

- B_check:

  number of bootstrap samples used in the closeness check

- RiskFactors:

  time series of the risk factors (F x T)

- GVARinputs:

  inputs used in the estimation of the GVAR-based models (see "GVAR"
  function). Default is set to NULL

- JLLinputs:

  inputs used in the estimation of the JLL-based models (see "JLL"
  function). Default is set to NULL

- FactorLabels:

  string-list based which contains the labels of all variables present
  in the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- ModelType:

  string-vector containing the label of the model to be estimated

- Use_Mean:

  Choose between the mean or the median across the bootstrap iterations.

- verbose:

  Logical flag controlling function messaging.

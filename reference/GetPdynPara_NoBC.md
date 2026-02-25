# Compute P-dynamics parameters without using the bias correction method from BRW (2012)

Compute P-dynamics parameters without using the bias correction method
from BRW (2012)

## Usage

``` r
GetPdynPara_NoBC(
  ModelType,
  RiskFactors,
  Economies,
  N,
  GVARinputs,
  JLLinputs,
  CheckInpts = F
)
```

## Arguments

- ModelType:

  string-vector containing the label of the model to be estimated

- RiskFactors:

  time series of risk factors (F x T). Could be stored in a list
  depending on the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- N:

  number of country-specific spanned factors

- GVARinputs:

  list of necessary inputs for the estimation of GVAR-based models (see
  "GVAR" function)

- JLLinputs:

  list of necessary inputs for the estimation of JLL-based models (see
  "JLL" function)

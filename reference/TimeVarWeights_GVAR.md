# Compute the star variables with time-varying weights

Compute the star variables with time-varying weights

## Usage

``` r
TimeVarWeights_GVAR(
  RiskFactors,
  Economies,
  RiskFactors_List,
  ListFactors,
  Wgvar,
  FactorLabels
)
```

## Arguments

- RiskFactors:

  A matrix of the complete set of risk factors (F x T).

- Economies:

  A character vector containing the names of the economies included in
  the system.

- RiskFactors_List:

  List of domestic risk factors (both spanned and unspanned)

- ListFactors:

  List of risk factors

- Wgvar:

  List of transition matrices

- FactorLabels:

  A list of character vectors with labels for all variables in the
  model.

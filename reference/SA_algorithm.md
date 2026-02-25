# Stochastic approximation algorithm

Stochastic approximation algorithm

## Usage

``` r
SA_algorithm(
  K1Z_NoBC,
  RiskFactors,
  BRWlist,
  GVARinputs,
  JLLinputs,
  FactorLabels,
  Economies,
  ModelType,
  verbose
)
```

## Arguments

- K1Z_NoBC:

  feedback matrix before bias-correction

- RiskFactors:

  A numeric matrix (T x F) representing the time series of risk factors.

- BRWlist:

  A list containing the necessary inputs for the BRW model estimation

- GVARinputs:

  List. Inputs for GVAR model estimation.

- JLLinputs:

  List. Inputs for JLL model estimation.

- FactorLabels:

  A list of character vectors with labels for all variables in the
  model.

- Economies:

  A character vector containing the names of the economies included in
  the system.

- ModelType:

  A character vector indicating the model type to be estimated.

- verbose:

  verbose Logical flag controlling function messaging.

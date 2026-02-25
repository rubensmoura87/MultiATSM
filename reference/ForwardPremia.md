# Compute the forward premia for all models

Compute the forward premia for all models

## Usage

``` r
ForwardPremia(
  ModelPara,
  avexpFP,
  ModelType,
  FactorLabels,
  InputsForOutputs,
  Economies
)
```

## Arguments

- ModelPara:

  list of model parameter estimates

- avexpFP:

  list containing the country-specific expected component of the forward
  period

- ModelType:

  desired model type

- FactorLabels:

  List of factor labels

- InputsForOutputs:

  list containing the desired horizon of analysis for the model fit,
  IRFs, GIRFs, FEVDs, GFEVDs, and risk premia decomposition

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

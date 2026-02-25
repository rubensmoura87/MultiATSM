# Decomposition of yields into the average of expected future short-term interest rate and risk premia for all models

Decomposition of yields into the average of expected future short-term
interest rate and risk premia for all models

## Usage

``` r
TermPremiaDecomp(
  ModelPara,
  FactorLabels,
  ModelType,
  InputsForOutputs,
  Economies
)
```

## Arguments

- ModelPara:

  list of model parameter estimates (see the "Optimization" function)

- FactorLabels:

  string-list based which contains all the labels of all the variables
  present in the model

- ModelType:

  string-vector containing the label of the model to be estimated

- InputsForOutputs:

  list containing the desired horizon of analysis for the model fit,
  IRFs, GIRFs, FEVDs, GFEVDs, and risk premia decomposition

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

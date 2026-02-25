# Numerical outputs (variance explained, model fit, IRFs, GIRFs, FEVDs, GFEVDs, and risk premia decomposition) for all models

Numerical outputs (variance explained, model fit, IRFs, GIRFs, FEVDs,
GFEVDs, and risk premia decomposition) for all models

## Usage

``` r
OutputConstruction(
  ModelType,
  ModelPara,
  InputsForOutputs,
  FactorLabels,
  Economies,
  verbose = TRUE
)
```

## Arguments

- ModelType:

  string-vector containing the label of the model to be estimated

- ModelPara:

  list of model parameter estimates (See the "Optimization" function)

- InputsForOutputs:

  list containing the desired horizon of analysis for the model fit,
  IRFs, GIRFs, FEVDs, GFEVDs, and risk premia decomposition

- FactorLabels:

  string-list based which contains all the labels of all the variables
  present in the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- verbose:

  Logical flag controlling function messaging. Default is TRUE.

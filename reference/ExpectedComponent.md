# Get the expected component of all models

Get the expected component of all models

## Usage

``` r
ExpectedComponent(
  ModelPara,
  InputsForOutputs,
  ModelType,
  Economies,
  FactorLabels,
  WishFP = FALSE
)
```

## Arguments

- ModelPara:

  list of model parameter estimates

- InputsForOutputs:

  list containing the desired horizon of analysis for the model fit,
  IRFs, GIRFs, FEVDs, GFEVDs, and risk premia decomposition

- ModelType:

  desired model type

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- FactorLabels:

  string-list based which contains all the labels of all the variables
  present in the model

- WishFP:

  If users wants to compute the forward premia. Default is FALSE.

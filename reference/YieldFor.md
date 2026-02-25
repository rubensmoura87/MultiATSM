# Compile the bond yield forecast for any model type

Compile the bond yield forecast for any model type

## Usage

``` r
YieldFor(
  ModelParaList,
  ForHoriz,
  Economies,
  FactorLabels,
  ForLabels,
  ModelType
)
```

## Arguments

- ModelParaList:

  List of model parameter estimates

- ForHoriz:

  Forecast horizon (scalar)

- Economies:

  String-vector containing the names of the economies which are part of
  the economic system

- FactorLabels:

  A string-list based which contains all the labels of all the variables
  present in the model

- ForLabels:

  Forecast labels (string-based vector)

- ModelType:

  A string-vector containing the label of the model to be estimated

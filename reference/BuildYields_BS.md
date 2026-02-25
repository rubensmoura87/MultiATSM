# Build the time-series of bond yields for each bootstrap draw

Build the time-series of bond yields for each bootstrap draw

## Usage

``` r
BuildYields_BS(
  ModelParaPE,
  ModelType,
  RiskFactors_BS,
  BFull,
  BS_Set,
  Economies
)
```

## Arguments

- ModelParaPE:

  List of point estimates of the model parameter

- ModelType:

  String-vector containing the label of the model to be estimated

- RiskFactors_BS:

  Time-series of the risk factors (F x T)

- BFull:

  B matrix of loadings

- BS_Set:

  Set of bootstrap inputs

- Economies:

  String-vector containing the names of the economies which are part of
  the economic system

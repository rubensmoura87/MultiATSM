# Perform out-of-sample forecast of bond yields

Perform out-of-sample forecast of bond yields

## Usage

``` r
OOS_Forecast(
  ForHoriz,
  t_Last,
  ModelParaList,
  FactorLabels,
  Yields_FullSample,
  Economies,
  ModelType
)
```

## Arguments

- ForHoriz:

  Forecast horizon-ahead (scalar)

- t_Last:

  Index of the last set of observations in the information set at a
  given forecasting round

- ModelParaList:

  List of model parameter estimates

- FactorLabels:

  A string-list based which contains all the labels of all the variables
  present in the model

- Yields_FullSample:

  Time-series of bond yields, complete set (J x T or CJ x T)

- Economies:

  String-vector containing the names of the economies which are part of
  the economic system

- ModelType:

  A string-vector containing the label of the model to be estimated

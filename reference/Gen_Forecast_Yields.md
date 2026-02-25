# Compute the bond yield forecast for any model type

Compute the bond yield forecast for any model type

## Usage

``` r
Gen_Forecast_Yields(
  K0Z,
  K1Z,
  A,
  Bfull,
  ZZsubsample,
  C,
  J,
  YieldsLabels,
  ForLabels,
  ForHoriz,
  ModelType
)
```

## Arguments

- K0Z:

  Intercept from the P-dynamics (F x 1)

- K1Z:

  Feedback matrix from the P-dynamics (F x F)

- A:

  Intercept of model-implied yields model (J x 1)

- Bfull:

  Slope of model-implied yields model (J x N or CJ x CN)

- ZZsubsample:

  Sub-sample of risk factors (F x t)

- C:

  Number of countries in the economic cohort (scalar)

- J:

  Number of country-specific bond yields

- YieldsLabels:

  Labels of bond yields

- ForLabels:

  Forecast labels (string-based vector)

- ForHoriz:

  Forecast horizon (scalar)

- ModelType:

  A string-vector containing the label of the model to be estimated

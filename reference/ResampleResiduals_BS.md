# Compute the residuals from the original model

Compute the residuals from the original model

## Usage

``` r
ResampleResiduals_BS(
  residPdynOriginal,
  residYieOriginal,
  InputsForOutputs,
  ModelType,
  nlag = 1
)
```

## Arguments

- residPdynOriginal:

  Time-series of the residuals from the P-dynamics equation (T x F)

- residYieOriginal:

  Time-series of the residuals from the observational equation (T x J or
  T x CJ)

- InputsForOutputs:

  List containing the desired inputs for the construction of the
  numerical outputs.

- ModelType:

  A character vector indicating the model type to be estimated

- nlag:

  Number of lags in the P-dynamics. Default is set to 1.

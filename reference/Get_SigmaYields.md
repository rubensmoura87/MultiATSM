# Compute the variance-covariance matrix of the bond yields

Compute the variance-covariance matrix of the bond yields

## Usage

``` r
Get_SigmaYields(YieldsTS, N, mat, WpcaFull, se, ModelType)
```

## Arguments

- YieldsTS:

  matrix of yields used in estimation (J x T or CJ x T)

- N:

  number of country-specific spanned factors

- mat:

  vector of maturities (in years) of yields used in estimation (J x 1)

- WpcaFull:

  composite matrix of weights the portfolios observed with and without
  errors

- se:

  Variance of the portfolio of yields observed with error (scalar).
  Default is set to NULL

- ModelType:

  string-vector containing the label of the model to be estimated

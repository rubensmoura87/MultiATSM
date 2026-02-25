# Prepare outputs to export after the model optimization

Prepare outputs to export after the model optimization

## Usage

``` r
OptOutputs(
  Y,
  Z,
  mat,
  N,
  dt,
  Wpca,
  K1XQ,
  SSZ,
  LoadAs,
  LoadBs,
  r0,
  se,
  K0Z,
  K1Z,
  Gy.0,
  VarYields,
  y,
  GVARinputs,
  JLLinputs,
  Economies,
  ModelType,
  BS_out = FALSE
)
```

## Arguments

- Y:

  matrix of yields used in estimation (J x T or CJ x T)

- Z:

  complete set of spanned and unspanned factors (F x T)

- mat:

  vector of maturities (in years) of yields used in estimation (J x 1)

- N:

  number of country-specific spanned factors

- dt:

  time interval unit of the model (scalar)

- Wpca:

  matrix of weights of the portfolios observed without errors (N x J or
  CN x J)

- K1XQ:

  risk-neutral feedback matrix (N x N or CN x CN)

- SSZ:

  variance-covariance matrix (F x F)

- LoadAs:

  list containing the A loadings

- LoadBs:

  list containing the B loadings

- r0:

  long-run interest rate (scalar or vector with length C)

- se:

  Variance of the portfolio of yields observed with error (scalar).

- K0Z:

  intercept from the P-dynamics (F x 1)

- K1Z:

  feedback matrix from the P-dynamics (F x F)

- Gy.0:

  matrix of contemporaneous terms from the P-dynamics (F x F)

- VarYields:

  variance-covariance matrix of the bond yields

- y:

  likelihood of each time series (Tx1)

- GVARinputs:

  List of inputs from GVAR models

- JLLinputs:

  List of inputs from JLL models

- Economies:

  string containing the names of the economy to be estimated

- ModelType:

  string-vector containing the label of the model to be estimated

- BS_out:

  Bootstrap output. Default is FALSE.

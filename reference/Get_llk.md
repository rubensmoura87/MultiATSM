# Compute the log-likelihood function

Compute the log-likelihood function

## Usage

``` r
Get_llk(P, Y, Z, N, mat, We, Wpca, K0Z, K1Z, SSZ, LoadBs, LoadAs, ModelType)
```

## Arguments

- P:

  time-series of spanned factors (N x T or CN x T)

- Y:

  time-series of yields (J x T or CJ x T)

- Z:

  time-series of risk factors (F x T)

- N:

  number of country-specific spanned factors

- mat:

  vector of maturities (in years) of yields used in estimation (J x 1)

- We:

  matrix of weights of the portfolios observed with errors ((J-N) x J or
  C(J-N) x CJ)

- Wpca:

  matrix of weights of the portfolios observed without errors (N x J or
  CN x CJ)

- K0Z:

  matrix of intercepts (P-dynamics)

- K1Z:

  feedback matrix (P-dynamics)

- SSZ:

  variance-covariance matrix (P-dynamics)

- LoadBs:

  list containing the B loadings

- LoadAs:

  list containing the A loadings

- ModelType:

  string-vector containing the label of the model to be estimated

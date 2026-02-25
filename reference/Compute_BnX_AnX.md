# Compute the latent loading AnX and BnX

Compute the latent loading AnX and BnX

## Usage

``` r
Compute_BnX_AnX(mat, N, K1XQ, r0, dX, SSX, Economies, ModelType, Lab_SingleQ)
```

## Arguments

- mat:

  vector of maturities (J x 1). Maturities are in multiples of the
  discrete interval used in the model

- N:

  number of country-specific spanned factors

- K1XQ:

  risk neutral feedback matrix (N x N)

- r0:

  the long run risk neutral mean of the short rate (scalar)

- dX:

  state loadings for the one-period rate (1xN). Default is a vector of
  ones

- SSX:

  the covariance matrix of the errors (N x N)

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- ModelType:

  string-vector containing the label of the model to be estimated

- Lab_SingleQ:

  string-vector containing the labels of the models estimated on a
  country-by-country basis

# Build the log-likelihood function of the P-dynamics from the JLL-based models

Build the log-likelihood function of the P-dynamics from the JLL-based
models

## Usage

``` r
llk_JLL_Sigma(VecPara, res, IdxNONzero, K)
```

## Arguments

- VecPara:

  vector that contains all the non-zero entries of the lower-triangular
  part of the Cholesky factorization

- res:

  residuals from the VAR of the JLL model (K x T)

- IdxNONzero:

  vector that contains indexes of the matrix of the non-zero entries of
  the Cholesky factorization

- K:

  dimensions of the variance-covariance matrix (scalar)

## Value

value of the log-likelihood function (scalar)

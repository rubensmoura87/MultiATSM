# Compute the cross-section loadings of yields of a canonical A0_N model

Compute the cross-section loadings of yields of a canonical A0_N model

## Usage

``` r
Get__BnXAnX(mat, K1XQ, ModelType, r0 = NULL, SSX = NULL, Economies)
```

## Arguments

- mat:

  vector of maturities (J x 1).

- K1XQ:

  risk-neutral feedback matrix (N x N)

- ModelType:

  string-vector containing the label of the model to be estimated

- r0:

  the long run risk neutral (scalar)

- SSX:

  covariance matrix of the latent states (N x N)

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

## References

- Dai, Q. and Singleton, K. (2000). "Specification Analysis of Affine
  Term Structure Models". The Journal of Finance.

- Joslin, S., Singleton, K. and Zhu, H. (2011). "A new perspective on
  Gaussian dynamic term structure models". The Review of Financial
  Studies.

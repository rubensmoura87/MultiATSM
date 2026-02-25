# Estimate numerically the Cholesky-factorization from the JLL-based models

Estimate numerically the Cholesky-factorization from the JLL-based
models

## Usage

``` r
EstimationSigma_Ye(SigmaUnres, res, M, G, Economies, DomUnit)
```

## Arguments

- SigmaUnres:

  unrestricted variance-covariance matrix (K x K)

- res:

  residuals from the VAR of the JLL model (K x T)

- M:

  number of domestic unspanned factors per country (scalar)

- G:

  number of global unspanned factors (scalar)

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- DomUnit:

  Name of the economy which is assigned as the dominant unit.  
  If no dominant unit is assigned, then this variable is defined as
  "none"

## Value

Cholesky-factorization after the maximization (K x K)

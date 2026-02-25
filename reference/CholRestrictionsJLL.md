# Impose the zero-restrictions on the Cholesky-factorization from JLL-based models.

Impose the zero-restrictions on the Cholesky-factorization from
JLL-based models.

## Usage

``` r
CholRestrictionsJLL(SigmaUnres, M, G, N, Economies, DomUnit)
```

## Arguments

- SigmaUnres:

  unrestricted variance-covariance matrix (K X K)

- M:

  number of domestic unspanned factors per country (scalar)

- G:

  number of global unspanned factors (scalar)

- N:

  number of country-specific spanned factors (scalar)

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- DomUnit:

  Name of the economy which is assigned as the dominant unit.  
  If no dominant unit is assigned, then this variable is defined as
  "none"

## Value

restricted version the Cholesky factorization matrix from JLL-based
models (K x K)

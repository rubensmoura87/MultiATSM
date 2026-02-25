# Find the indexes of zero-restrictions from the orthogonalized variance-covariance matrix from the JLL-based models

Find the indexes of zero-restrictions from the orthogonalized
variance-covariance matrix from the JLL-based models

## Usage

``` r
IDXZeroRestrictionsJLLVarCovOrtho(M, N, G, Economies, DomUnit)
```

## Arguments

- M:

  number of country-specific unspanned factors (scalar)

- N:

  number of country-specific spanned factors (scalar)

- G:

  number of global unspanned factors (scalar)

- Economies:

  Set of economies that are part of the economic system (string-vector)

- DomUnit:

  Name of the economy which is assigned as the dominant unit.  
  If no dominant unit is assigned, then this variable is defined as
  "None"

## Value

restricted version of the JLL of the Cholesky factorization (F x F)

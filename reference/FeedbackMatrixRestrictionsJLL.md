# Set the zero-restrictions on the feedback matrix of JLL's P-dynamics

Set the zero-restrictions on the feedback matrix of JLL's P-dynamics

## Usage

``` r
FeedbackMatrixRestrictionsJLL(DomUnit, K, G, M, N)
```

## Arguments

- DomUnit:

  Name of the economy which is assigned as the dominant unit.  
  If no dominant unit is assigned, then this variable is defined as
  "none"

- K:

  Total number of risk factors of the economic system (scalar)

- G:

  Number of global unspanned factors (scalar)

- M:

  Number of country-specific unspanned factors (scalar)

- N:

  Number of country-specific spanned factors (scalar)

## Value

matrix containing the restrictions of the feedback matrix (K x K)

# Estimates a standard VAR(1)

Estimates a standard VAR(1)

## Usage

``` r
VAR(RiskFactors, VARtype, Bcon_Mat = NULL)
```

## Arguments

- RiskFactors:

  numeric matrix (`K x Td`). Time series of risk factors.

- VARtype:

  character. Permissible choices: "unconstrained" or "constrained".

- Bcon_Mat:

  matrix (`K x K + 1`). Constraints matrix (includes intercept). Entries
  containing NAs are treated as free parameters. Default is NULL.

## Value

list. Contains:

- intercept (K x 1)

- feedback matrix (K x K)

- variance-covariance matrix (K x K) of a VAR(1)

## General Notation

- `Td`: model time series dimension

- `N`: number of country-specific spanned factors

- `K`: total number of risk factors

## Examples

``` r
data(RiskFacFull)
# Example 1: unconstrained case
VAR_para1 <- VAR(RiskFacFull, VARtype = "unconstrained")

# Example 2: constrained case
K <- nrow(RiskFacFull)
Bcon_Mat <- matrix(0, nrow = K, ncol = K + 1)
Bcon_Mat[, 1:3] <- NaN
VAR_para2 <- VAR(RiskFacFull, VARtype = "constrained", Bcon_Mat)
```

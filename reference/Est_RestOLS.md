# Estimate a restricted OLS model

Estimate a restricted OLS model

## Usage

``` r
Est_RestOLS(LHS, RHS, Rmat)
```

## Arguments

- LHS:

  left hand side variables (M x T).

- RHS:

  right hand side variables (should include the intercept, if desired)
  (N x T).

- Rmat:

  matrix of constraints (M x N). Entries containing NAs are treated as
  free parameters.

## Value

matrix of coefficient (M x N)

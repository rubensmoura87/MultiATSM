# Estimates the P-dynamics from JLL-based models

Estimates the P-dynamics from JLL-based models

## Usage

``` r
JLL(NonOrthoFactors, N, JLLinputs, CheckInputs = FALSE)
```

## Arguments

- NonOrthoFactors:

  numeric matrix (`K x Td`). Time series of risk factors before
  orthogonalization.

- N:

  positive integer. Number of country-specific spanned factors.

- JLLinputs:

  list. Necessary inputs to estimate JLL models:

  1.  `Economies`: character vector. Set of `C` economies in the system.

  2.  `DomUnit`: character. Name of the dominant economy, or "None" if
      not assigned (for "JLL No DomUnit" model).

  3.  `WishSigmas`: logical. TRUE to estimate variance-covariance
      matrices and Cholesky factorizations; FALSE otherwise.

  4.  `SigmaNonOrtho`: NULL or F x F matrix from non-orthogonalized
      dynamics.

  5.  `JLLModelType`: character. Permissible choices: "JLL original",
      "JLL joint Sigma", "JLL No DomUnit".

- CheckInputs:

  logical. Whether to perform a prior consistency check on the inputs
  provided in `JLLinputs`. Default is FALSE.

## Value

List of model parameters from both the orthogonalized and
non-orthogonalized versions of the JLL-based models

## General Notation

- `Td`: model time series dimension

- `C` number of countries in the system.

- `K`: total number of risk factors

## References

Jotiskhatira, P. ; Le, A. and Lundblad, C. (2015). "Why do interest
rates in different currencies co-move?" (Journal of Financial Economics)

## Examples

``` r
# \donttest{
data(RiskFacFull)
RF_TS <- RiskFacFull
N <- 3
JLLinputs <- list(
  Economies = c("China", "Brazil", "Mexico", "Uruguay"), DomUnit = "China",
  WishSigmas = TRUE, SigmaNonOrtho = NULL, JLLModelType = "JLL original"
)
JLLPara <- JLL(RF_TS, N, JLLinputs)
# }
```

# Optimization routine for the entire selected ATSM

Optimization routine for the entire selected ATSM

## Usage

``` r
OptimizationSetup_ATSM(x0, MLvec, EstType, tol = 1e-04, verbose = FALSE)
```

## Arguments

- x0:

  List containing features for estimation of the risk-neutral
  parameters.

- MLvec:

  Log-likelihood function

- EstType:

  Estimation type. Available options: "BFGS" and "Nelder-Mead".

- tol:

  convergence tolerance (scalar). Default value is set as 1e-4.

- verbose:

  Logical flag controlling function messaging.
